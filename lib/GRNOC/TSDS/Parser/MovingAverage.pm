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
        $window_seconds = abs($window_seconds);
        my $window = int($window_seconds / $interval);

        # Put the data points on a uniform grid of size $interval:
        my $first_time = $set[0]->[0];
        my $last_time  = $set[scalar(@set) - 1]->[0];
        my $npoints    = int(($last_time - $first_time + 1) / $interval);

        my @vals = (undef) x $npoints;
        foreach my $datum (@set){
            $vals[int(($datum->[0] - $first_time) / $interval)] = $datum;
        }

        # Actually compute the windowed average:
        my $center_of_window = ($window * $interval) / 2;
        my $offset = $first_time + $center_of_window;
        my @new_set;

        for (my $i = 0; $i < $npoints - $window; $i += 1){
            my @points_in_window = grep { defined($_) } @vals[$i..($i+$window-1)];
            # Whether we should include $vals[$i+$window] depends on just what
            # time that datapoint has...
            my $last = $vals[$i+$window];
            push @points_in_window, $last if defined($last) && ($last->[0] <= $first_time + ($i * $interval) + $window_seconds);

            my $total = sum( map { $_->[1] } @points_in_window );
            my $n_defined = scalar(@points_in_window);
            my $avg = ($n_defined > 0) ? $total / $n_defined : undef;
            push @new_set, [int($offset + ($i * $interval)), $avg];
        }

        $result = \@new_set;
    }

    return $result;
}

1;
