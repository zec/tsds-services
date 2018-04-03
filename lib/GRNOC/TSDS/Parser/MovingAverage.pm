package GRNOC::TSDS::Parser::MovingAverage;

use strict;
use warnings;

use Data::Dumper;
use List::Util qw(min max sum);
use List::MoreUtils qw(all any);
use Date::Format qw(time2str);
use Math::FFT qw();

use GRNOC::Log;

# Top-level moving-average function; returns a reference to an array of
# output (timestamp, value) pairs on success, or undef or an error-message
# string if there was a problem.
sub moving_average {
    my $set = shift;
    my $window = shift;

    my @set = grep { defined($_) && defined($_->[0]) } @$set;

    my $result;

    # the moving average of zero or one points is just the input array:
    if (scalar(@set) < 2){
        $result = \@set;
    }
    else {
        # Sort @set by time:
        @set = sort { $a->[0] <=> $b->[0] } @set;

        # It's possible that the points in @set aren't evenly distributed in
        # time... so we take some precautions calculating the interval:
        my $interval = min( map { $set[$_]->[0] - $set[$_-1]->[0] } (1..(scalar(@set)-1)) );

        # Put the data points on a uniform grid of size $interval:
        my $first_time = $set[0]->[0];
        my $last_time  = $set[scalar(@set) - 1]->[0];
        my $npoints    = int(($last_time - $first_time + 1) / $interval);

        my @grid = (undef) x $npoints;
        foreach my $datum (@set){
            $grid[int(($datum->[0] - $first_time) / $interval)] = $datum;
        }

        # Now that we've done that, compute the windowed average:
        if ($window eq 'ASAP') {
            # Use the ASAP window-size auto-tuning algorithm to find the "best"
            # window size, then use that window:
            $result = _asap(\@set);
        }
        else {
            # If we get here, window is actually a number: the
            # manually-specified window size, in seconds
            $result = _mavg(\@grid, $interval, $window);
        }
    }

    return $result;
}

# Compute the windowed-average of a data set given
# 1) A ref to an array, which has (timestamp, datum) pairs such that
#    each pair is at array index int((timestamp-[first timestamp])/$interval),
# 2) The interval of the aformentioned array, which should usually be the
#    collection interval for the underlying data, in seconds
# 3) The size of the window (in seconds)
sub _mavg {
    my $data           = shift;
    my $interval       = shift;
    my $window_seconds = shift;

    $window_seconds = abs($window_seconds);
    my @data = @$data;
    my $first_time = $data[0]->[0];

    my $window = int($window_seconds / $interval);
    my $center_of_window = ($window * $interval) / 2;
    my $offset = $first_time + $center_of_window;
    my @smoothed_data;

    my $npoints = scalar(@data);

    for (my $i = 0; $i < $npoints - $window; $i += 1){
        my @points_in_window = grep { defined($_) && defined($_->[1]) } @data[$i..($i+$window-1)];
        # Whether we should include $vals[$i+$window] depends on just what
        # timestamp that datapoint has...
        my $last = $data[$i+$window];
        if (defined($last) && defined($last->[1]) &&
            ($last->[0] <= $first_time + ($i * $interval) + $window_seconds)){

            push @points_in_window, $last;
        }

        my $total = sum( map { $_->[1] } @points_in_window );
        my $n_defined = scalar(@points_in_window);
        my $avg = ($n_defined > 0) ? $total / $n_defined : undef;
        push @smoothed_data, [int($offset + ($i * $interval)), $avg];
    }

    return \@smoothed_data;
}

sub _asap {
    my $data = shift; # array of (timestamp, value) pairs, sorted by timestamp

    my $step = $data->[1][0] - $data->[0][0];

    # For now, we only support data sets that are completely regularly-sampled
    # (no changes in sample rate, and no gaps in data), as this is much easier
    # to calculate an autocorrelation for (indeed, someone else has already
    # written the code :) )
    foreach my $i (1..(scalar(@$data) - 1)){
        next if (($data->[$i][0] - $data->[$i-1][0]) == $step)
                && defined($data->[$i][1]) && defined($data->[$i-1][1]);

        my $date = time2str('%Y-%m-%dT%H:%M:%SZ', $data->[$i-1][0], 'UTC');
        return "ASAP smoothing currently only supports data without gaps (first gap near $date)";
    }

    # Remove the timestamps, and get just the data:
    my @values = map { $_->[1] } @$data;

    # So, our data is regularly-sampled, with no gaps. Get the autocorrelation:
    my $corr = _calc_acf(\@values);

    # Get the peaks of the autocorrelation - they are our candidate window sizes:
    my $candidates = _find_peaks($corr);

    my $base_stats = _calc_stats(\@values);

    # Initial state for search:
    my $initial_opt = {
        lower_bound => 1,
        roughness   => $base_stats->{'roughness'},
        window      => 1,
        lfi         => -1, # "largest feasible index"
    };

    my $opt = _search_periodic(\@values, $base_stats, $candidates, $corr, $initial_opt);

    # Search bounds for binary search
    my $hi = max(int(scalar(@values) / 10), 1);
    my $lo = $opt->{'lower_bound'};

    if ($opt->{'lfi'} >= 0) {
        $hi = $candidates->[$opt->{'lfi'}+1] if $opt->{'lfi'} <= scalar(@$candidates) - 2;
        $lo = max($lo, $candidates->[$opt->{'lfi'}] + 1);
    }

    my $window_size = _binary_search(\@values, $base_stats, $lo, $hi, $opt);

    my @Y = @{_sma(\@values, $window_size)};

    my $window_offset = $data->[0][0] + (($window_size - 1) * $step / 2);
    @Y = map { [ $window_offset + ($_ * $step), $Y[$_] ] } (0..$#Y);

    return \@Y;
}

# Calculate the statistics of a timeseries
sub _calc_stats {
    my $data = shift; # An array of (defined) values; assumed to have at least two elements
    my $N = scalar(@$data); # Number of elements in @$data

    # Standard statistical measures
    my ($mean, $variance) = _mean_var($data);
    my $stdev = sqrt($variance);

    my $expected_fourths = sum(map { ($_ - $mean) ** 4 } @$data) / ($N - 1);
    my $kurtosis = $expected_fourths / ($variance * $variance);

    # Roughness, as defined in the ASAP paper. First, get the series of
    # first differences:
    my @diffs = map { $data->[$_+1] - $data->[$_] } (0..(scalar(@$data)-2));
    # Then, get its standard deviation:
    my ($dm, $dvar) = _mean_var(\@diffs);
    my $roughness = sqrt($dvar);

    return {
        mean      => $mean,
        variance  => $variance,
        stdev     => $stdev,
        kurtosis  => $kurtosis,
        roughness => $roughness,
    };
}

# Used by _calc_stats
sub _mean_var {
    my $data = shift; # An array of (defined) values; assumed to have at least one element
    my $N = scalar(@$data);

    my $mean = (sum @$data) / $N;
    my $var = 0; # variance
    if ($N > 1) {
        $var = sum(map { my $x = $_ - $mean; $x * $x } @$data) / ($N - 1);
    }

    return ($mean, $var);
}

# Calculate the autocorrelation of the data, with some help from FFTs:
sub _calc_acf {
    # A list of (defined) values, representing regularly-sampled data:
    my $data = shift;

    my @values = @$data;

    # Zero-pad out values, as (1) FFT uses power-of-two lengths only and
    # (2) the FFT-based algorithm actually calculates a *circular*
    # autocorrelation, so we want a large blank area after the actual data
    # to avoid wraparound doing silly things to our results.
    my $len = scalar(@values);
    my $len_b2 = int(log($len)/log(2))-1; # The "-1" and search is due to paranoia about rounding errors
    while ((2**$len_b2) <= $len){
        $len_b2 += 1;
    }
    push @values, (0 x ((2**$len_b2) - $len));

    # Calculate the (circular) autocorrelation:
    my $fft = Math::FFT->new(\@values);
    my @corr = @{$fft->correl($fft)};

    # For our purposes, we only care about the first $len entries, so we truncate:
    $#corr = $#values;

    # ASAP uses a variance-normalized autocorrelation:
    my $variance = $corr[0];
    my @normed = map { $_ / $variance } @corr;

    return \@normed;
}

# Find the peaks (local maxima) in an array of data,
# *excluding the first element*
sub _find_peaks {
    my $arr = shift;

    my $len = scalar(@$arr);
    my @peaks;

    # Handle small arrays and the second element specially
    return [] if $len < 2;
    return [1] if $len == 2;
    push @peaks, 1 if $arr->[1] > $arr->[2];

    my $i = 2;
    while ($i < ($len-1)){
        if ( ($arr->[$i-1] < $arr->[$i]) && ($arr->[$i] > $arr->[$i+1]) ){
            push @peaks, $i;
        }

        # The above handles cases with a sharp peak:
        #         *
        #     * *   *
        #            *
        # But what about plateaus (or false plateaus)?
        #                                 *
        #       * * *                 * *
        #     *       *     vs.     *
        #   *                     *
        # Well:
        if ( ($arr->[$i-1] < $arr->[$i]) && ($arr->[$i] == $arr->[$i+1]) ){
            my $j = $i;
            while ( ($j < ($len-1)) && ($arr->[$j] == $arr->[$j+1]) ){
                $j += 1;
            }

            # Real plateau; push all those points
            if ( ($j == ($len-1)) || ($arr->[$j] > $arr->[$j+1]) ){
                push @peaks, $i..$j;
            }

            $i = $j;
        }

        $i += 1;
    }

    # last element
    push @peaks, ($len-1) if $arr->[$len-1] > $arr->[$len-2];

    return \@peaks;
}

# The searchPeriodic function of the ASAP paper, using the following arguments:
# * the timeseries itself
# * the results of _calc_stats(timeseries)
# * the candidate window sizes
# * the autocorrelation of the timeseries
# * the initial search state
#
# Returns the final search state.
sub _search_periodic {
    my ($X, $X_stats, $candidates, $acf, $initial_opt) = @_;

    my %opt = %$initial_opt;
    my @acf = @$acf;
    my $max_acf = max(@acf[@$candidates]);

    for (my $i = scalar(@$candidates) - 1; $i >= 0; $i -= 1) {
        my $w = $candidates->[$i];

        last if $w < $opt{'lower_bound'};

        next if ((1-$acf->[$w])             / ($w*$w))
              > ((1-$acf->[$opt{'window'}]) / ($opt{'window'}*$opt{'window'}));

        my $Y = _sma($X, $w);

        my %Y = %{_calc_stats(@$Y)};

        if ($Y{'kurtosis'} >= $X_stats->{'kurtosis'}) {
            if ($Y{'roughness'} < $opt{'roughness'}) {
                $opt{'window'} = $w;
                $opt{'roughness'} = $Y{'roughness'};
            }
            $opt{'lower_bound'} = int( max(
                $opt{'lower_bound'},
                $w * sqrt( (1-$max_acf) / (1-$acf->[$w]) )
            ) );
            $opt{'lfi'} = max($opt{'lfi'}, $i);
        }
    }

    return \%opt;
}

# Binary search on window sizes
sub _binary_search {
    my ($X, $X_stats, $low, $high, $initial_opt) = @_;

    my %opt = %$initial_opt;

    while ($low <= $high) {
        my $w = int(($low + $high) / 2); # test window size

        my $Y = _sma($X, $w);
        my %Y = %{_calc_stats(@$Y)};

        # For uncorrelated data, kurtosis and roughness tend to decrease
        # with increasing window size. We want to find the smoothest window
        # that still has at least as much kurtosis as the original data, so:

        if ($Y{'kurtosis'} >= $X_stats->{'kurtosis'}) {
            if ($Y{'roughness'} < $opt{'roughness'}) {
                $opt{'window'} = $w;
                $opt{'roughness'} = $Y{'roughness'};
            }

            # search upper half of range
            $low = $w + 1;
        } else {
            # search lower half of range
            $high = $w - 1;
        }
    }

    # Return the best window we found
    return $opt{'window'};
}


# Windowed moving average, with regularly-sampled, all-defined data
sub _sma {
    my $X = shift; # timeseries
    my $window = shift; # window size; must be a positive integer

    my @X = @$X;
    my @Y;

    foreach my $i (0..(scalar(@X)-$window)) {
        push @Y, (sum @X[$i..($i+$window-1)]) / $window;
    }

    return \@Y;
}

1;
