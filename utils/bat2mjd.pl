#!/usr/bin/env perl

use strict;
use warnings;
use Math::BigFloat;
use Math::BigInt;
use Astro::Time;

my $DUTC = 37;

=item B<bat2mjd>

  my $mjd = bat2mjd($bat);
` my $mjd = bat2mjd($bat, $dUT);

 Convert a bat into mjd
    $bat           BAT value
    $mjd           MJD (double)
    $dUT           Offset in seconds between TAI and UTC
=cut

sub bat2mjd($;$) {
 my $bat = Math::BigInt->new(shift);
 my $dUT = shift;
 $dUT = $DUTC if (!defined $dUT);
 return (Math::BigFloat->new($bat)->bstr()/1e6-$dUT)/60/60/24;
}




sub bat2time($;$$) {
 my $bat = Math::BigInt->new(shift);

 my $dUT = shift;
 my $np = shift;

 my $mjd = bat2mjd($bat, $dUT);

 return mjd2time($mjd, $np);
}


my ($time, $sec);
foreach (@ARGV) {

  my $t = bat2mjd($_);
  printf("%.12f\n", t);
  
#    $time = turn2str($_/(60*60*24), 'H', 0);
#    if (defined $time) {
#      printf "%5d -> %s\n", $_, $time;
#    } else {
#      print "Error converting $_\n";
#    }
}
