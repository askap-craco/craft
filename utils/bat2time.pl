#!/usr/bin/env perl

use strict;

use Astro::Time;
use ATNF::MoniCA;


my ($time, $sec);
foreach (@ARGV) {

  my $t = bat2time($_,undef,3);
  print "$t\n";
  
#    $time = turn2str($_/(60*60*24), 'H', 0);
#    if (defined $time) {
#      printf "%5d -> %s\n", $_, $time;
#    } else {
#      print "Error converting $_\n";
#    }
}
