#!/usr/bin/perl

use strict;
use warnings;

my $gotbat = 0;
my $gotframeid = 0;
my $gotnsamp = 0;
my ($triggerFrameId, $nsamps, $triggerBAT);

my $samplesPerWord = 4; # Should set

while (<>) {

  if (/\s*TRIGGER_BAT\s+(\S+)/) {
    $gotbat = 1;
    $triggerBAT = hex($1);
  } elsif (/\s*TRIGGER_FRAMEID\s+(\d+)/) {
    $gotframeid = 1;
    $triggerFrameId = $1;
  } elsif (/\s*NSAMPS_REQUEST\s+(\d+)/) {
    $gotnsamp = 1;
    $nsamps = $1;
  }
  last if ($gotbat and $gotframeid and $gotnsamp);
}

if (! ($gotbat and $gotframeid and $gotnsamp)) {
  die "Did not find TRIGGER_BAT, TRIGGER_FRAMEID and NSAMPS_REQUEST\n";
}

my $startFrameId = $triggerFrameId - $nsamps + $samplesPerWord;

print "startFrameID=$startFrameId\n";

my $BAT0 = $triggerBAT - ($triggerFrameId * (27.0/32.0));
printf("BAT0= 0x%X\n", $BAT0);
      
