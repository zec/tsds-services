package GRNOC::TSDS::Parser::MovingAverage;

use strict;
use warnings;

use Data::Dumper;
use List::Util qw(min max sum);
use List::MoreUtils qw(all);

use GRNOC::Log;

# Top-level moving-average function
sub moving_average {
    my $set = shift;
    my $window_seconds = shift;

    my @set = @$set;
    $window_seconds = abs($window_seconds);

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

        # Actually compute the windowed average:
        $result = _mavg(\@grid, $interval, $window_seconds);
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

1;
