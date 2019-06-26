#!/usr/bin/perl

use strict;
use warnings;
no warnings 'portable';

my $gotbat = 0;
my $gotframeid = 0;
my ($triggerFrameId, $triggerBAT);

while (<>) {

  if (/\s*TRIGGER_BAT\s+(\S+)/) {
    $gotbat = 1;
    $triggerBAT = hex($1);
  } elsif (/\s*TRIGGER_FRAMEID\s+(\d+)/) {
    $gotframeid = 1;
    $triggerFrameId = $1;
  }
  last if ($gotbat and $gotframeid);
}

if (! ($gotbat and $gotframeid)) {
  die "Did not find TRIGGER_BAT and TRIGGER_FRAMEID\n";
}

my $startBAT = int ($triggerBAT - ($triggerFrameId * (27.0/32.0)));
my $BAT0 = int(($startBAT + 999999)/1e6); # Round to nearest uSec
$BAT0 *= 1e6;

my $frame0 = int(($BAT0-$startBAT)*(32.0/27.0));

open(BAT0, '>', '.bat0') || die "Could not open .bat0 for writing\n";
printf(BAT0 "0x%llx %d\n", $BAT0, $frame0);
close(BAT0);
 
