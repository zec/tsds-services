package GRNOC::TSDS::Parser::MovingAverage;

use strict;
use warnings;

use Data::Dumper;
use List::Util qw(min max sum);
use List::MoreUtils qw(all any);
use Date::Format qw(time2str);
use Math::FFT qw();

use GRNOC::Log;

# Top-level moving-average function
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
            $result = _asap(\@set, \@grid, $interval);
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
    my $grid     = shift; # $data, reformatted for use by _mavg
    my $interval = shift; # Grid interval used by $grid

    my $step = $data->[1][0] - $data->[0][0];

    # For now, we only support data sets that are completely regularly-sampled
    # (no changes in sample rate, and no gaps in data), as this is much easier
    # to calculate an autocorrelation for (indeed, someone else has already
    # written the code :) )
    foreach my $i (1..(scalar(@$data) - 1)){
        next if (($data->[$i][0] - $data->[$i-1][0]) == $step)
                && defined($data->[$i][1]) && defined($data->[$i-1][1]);

        # TODO: implement actual error
        my $date = time2str('%Y-%m-%dT%H:%M:%SZ', $data->[$i-1][0], 'UTC')
        return "ASAP smoothing currently only supports data without gaps (first gap near $date)";
    }

    # Remove the timestamps, and get just the data:
    my @values = map { $_->[1] } @$data;

    # So, our data is regularly-sampled, with no gaps. Get the autocorrelation:
    my $corr = _calc_acf(\@values);

    # Get the peaks of the autocorrelation - they are our candidate window sizes:
    my $candidates = _find_peaks($corr);

    my $base_stats = _calc_stats(\@values);

    return [];
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
    my @diffs = map { $data->[$i+1] - $data->[$i] } (0..(scalar(@$data)-2));
    # Then, get its standard deviation:
    my ($dm, $dvar) = _mean_var(\@diffs);
    my $roughness = sqrt($dvar);

    my %results = (
        mean      => $mean,
        variance  => $variance,
        stdev     => $stdev,
        kurtosis  => $kurtosis,
        roughness => $roughness,
    );

    return \%results;
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

    return ($mean, $variance);
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

    return \@corr;
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

1;
