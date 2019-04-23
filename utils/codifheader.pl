#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use Carp;
use POSIX;

sub turn2str ($;$);
sub readheader($);
sub mjd2cal($);
sub cal2mjd($$$;$);

my ($invalid, $complex, $seconds, $frame, $version, $nbits, $refepoch, 
    $representation,  $antid, $sampleblocklength, $nchan, $threadid, $groupid, $period, 
    $numsamples, $sync, $newframelength);

my $once = 0;
my $check = 0;
my $dosync = 0;
my $help = 0;
my $skip = 0;
my $framelength = undef;

GetOptions('once'=>\$once, 'check'=>\$check, 'help'=>\$help, 'skip=i'=>\$skip, 'framelength=i'=>\$framelength, 'sync'=>\$dosync);

my ($lastframe, $lastsec);

if ($help || @ARGV==0) {
  print<<EOF;
Usage: vdifheader.pl [options] <vdiffile>

Options:
   -once          Print only the first header 
   -check         Do check frames increase monotonically with no gaps (single thread only?)
   -sync          Check SYNC word is correct and report byte offset if not
   -skip <bytes>  Skip <bytes> bytes at the start of each file    
EOF
}

foreach (@ARGV) {
  open(CODIF, $_) || die "Could not open $_: $!\n";


  print "Reading $_\n\n";

  my $first = 1;
  while (1) {

    if ($skip>0) {
      my $status = sysseek CODIF, $skip, SEEK_SET;
      if (! defined $status) {
	warn "Error trying to skip $skip bytes - skipping file\n";
	last;
      } elsif ($status != $skip) {
	warn "Failed to skip $skip bytes - is this a CODIF file?\n";
	last;
      }
    }

    ($invalid, $complex, $seconds, $frame, $version, $nbits, $newframelength, $refepoch, 
     $representation,  $antid, $sampleblocklength, $nchan, $threadid, $groupid, $period, 
     $numsamples, $sync) = readheader(*CODIF);

    if (!defined $invalid) {
      print "   empty file\n" if ($first);
      close(CODIF);
      last;
    }

    $framelength = $newframelength if (! defined $framelength);

    my $hexsync = sprintf("0x%X", $sync);
    my $timestr = turn2str(fmod($seconds/60/60/24, 1.0));

    if (!$check && !$dosync) {
      print "-------------------\n" if (!$first);
      my $mjd =  cal2mjd(1, ($refepoch%2)*6+1, 2000 + int($refepoch/2));
      $mjd += int($seconds/(60*60*24));
      my ($day, $month, $year, $ut) = mjd2cal($mjd);
      my $date = sprintf("%02d/%02d/%04d", $day, $month, $year);

      # Derived parameters
      my $samplesperframe = $newframelength/($nbits*$nchan*($complex+1))*8;
      my $frameperperiod = $numsamples/$samplesperframe;
      my $frameuSec = $period/$frameperperiod*1e6;
      my $secoffset = $frame*($frameuSec/1e6);
      my $timestr2 = turn2str(fmod(($seconds+$secoffset)/60/60/24, 1.0),7);
      
      print<<EOF;
INVALID:     $invalid
COMPLEX:     $complex

SECONDS:     $seconds     $timestr   $date

FRAME#:      $frame       $timestr2

VERSION:     $version
NBITS:       $nbits
FRAMELENGTH: $newframelength

EPOCH:       $refepoch
REPR:        $representation
ANTID:       $antid

SAMPLEBLOCK: $sampleblocklength
NCHAN:       $nchan

THREADID:    $threadid
GROUPID:     $groupid

PERIOD:      $period

#SAMPLES:    $numsamples

SYNC:        $hexsync

Samples/Frame = $samplesperframe
Frame uSec = $frameuSec

EOF
    } elsif ($check) {
      if ($first) {
	$lastsec = $seconds;
	$lastframe = $frame;
      } else {
	if ($frame-$lastframe!=1) {
	  if ($seconds==$lastsec) {
	    printf("Skipped %d frames at $timestr/$lastframe--$frame\n", $frame-$lastframe-1);
	  } elsif ($seconds-$lastsec==1) {
	    if ($frame!=0) {
	      printf("Skipped > %d frames at $timestr/$lastframe--$frame\n", $frame);
	    }
	  } else {
	    printf("Skipped  %d seconds at $timestr\n", $seconds-$lastsec);
	  }
	}
	$lastframe = $frame;
	$lastsec = $seconds;
      }
    } elsif ($dosync) {
      if ($sync!=0xABADDEED) {
	my $pos = sysseek(CODIF, 0, 1) - 64; # Account for header
	die "Bad sync ($hexsync) at offset $pos\n";
      }
    }
    
    my $status = sysseek(CODIF, $framelength, 1);
    if (!defined $status) {
      close(CODIF);
      last;
    }
    last if ($once);
    $first = 0;
  }

  print "All good!\n" if $dosync;
}

sub readheader($) {
  my $vdif = shift;
  my $buf;
  my $nread = sysread($vdif, $buf, 64);
  if (! defined($nread)) {
    die("Error reading $_: $!");
    return undef;
  } elsif ($nread==0) { # EOF
    return undef;
  } elsif ($nread!=64) {
    die "Error: Only read $nread bytes from header\n";
  }

  my ($seconds,  $frame, $framelength, $antid, $nchan, $groupid, $period, 
      $reserved1, $numbersamples, $sync, $reserved) = unpack 'VVVVVVVVQVV', $buf;

  my $invalid = $seconds>>31;
  my $complex = ($seconds>>30)&0x1;
  $seconds &= 0x3FFFFFFF;

  my $nbits = (($framelength>>24)&0x1F);
  my $version = $framelength >> 29;
  $framelength &= 0xFFFFFF;
  #$framelength++;
  $framelength *= 8;
 

  my $refepoch = ($antid>>26)&0x3F;
  my $representation = ($antid>>22)&0xF;
  $antid = unpack 'A2', pack('n', $antid&0xFFFF);

  my $sampleblocklength = ($nchan>>24)&0x1F;
  $nchan &= 0xFFFF;
  #$nchan++;

  my $threadid = ($groupid>>16)&0xFFFF;
  $groupid &= 0xFFFF;

  $period &= 0xFFFF;

  return ($invalid, $complex, $seconds, $frame, $version, $nbits, $framelength, $refepoch, $representation, 
	  $antid, $sampleblocklength, $nchan, $threadid, $groupid, $period, $numbersamples, $sync);

}

sub turn2str ($;$) {
  my($turn) = @_;
  my $mode = 'H';
  if (($mode ne 'H') && ($mode ne 'D')) {
    carp 'turn2str: $mode must equal \'H\' or \'D\'';
    return undef;
  }
  my $strsep = ':';

  my ($angle, $str, $sign, $wholesec, $secfract, $min);

  if ($mode eq 'H') {
    $angle = $turn * 24;
  } else {
    $angle = $turn * 360;
  }

  if ($angle < 0.0) {
    $sign = -1;
    $angle = -$angle;
  } else {
    $sign = 1;
  }

  my $wholeangle = (int $angle);

  $angle -= $wholeangle;
  $angle *= 3600;

  # Get second fraction
  $wholesec = int $angle;
  $secfract = $angle - $wholesec;

  my $sig = 0;
  $sig = $_[1] if (defined $_[1]);
  $wholesec %= 60;
  $min = ($angle-$wholesec - $secfract)/60.0;
  $secfract = int ($secfract * 10**$sig + 0.5); # Add 0.5 to ensure rounding

  # Check we have not rounded too far
  if ($secfract >= 10**$sig) {
    $secfract -= 10**$sig;
    $wholesec++;
    if ($wholesec >= 60.0) {
      $wholesec -= 60;
      $min++;
      if ($min >= 60.0) {
	$min -= 60;
	$wholeangle++;
      }
    }
  }

  my $angleform = '%02';

  my ($sep1, $sep2, $sep3);
  if ($strsep eq 'HMS') {
    if ($mode eq 'H') {
      $sep1 = 'H';
    } else {
      $sep1 = 'D';
    }
    $sep2 = 'M';
    $sep3 = 'S';
  } elsif ($strsep eq 'hms') {
    if ($mode eq 'H') {
      $sep1 = 'h';
    } else {
      $sep1 = 'd';
    }
    $sep2 = 'm';
    $sep3 = 's';
  } elsif ($strsep eq 'deg') { # What if $mode eq 'H'??
    $sep1 = 'd';
    $sep2 = "'";
    $sep3 = '"';
  } else {
    $sep1 = $sep2 = $strsep;
    $sep3 = '';
  }

  if ($sig > 0) {
    $str = sprintf("${angleform}d$sep1%02d".
		   "$sep2%02d.%0${sig}d$sep3", 
		   $wholeangle, $min, $wholesec, $secfract);
  } else {
    $str = sprintf("${angleform}d$sep1%02d".
		   "$sep2%02d$sep3", 
		   $wholeangle, $min, $wholesec);
  }

  if ($sign == -1) {
    $str = '-'.$str;
  }
  return $str;
}


# Pinched from Astro::Time

sub cal2mjd($$$;$) {
  my($day, $month, $year, $ut) = @_;
  $ut = 0.0 if (!defined $ut);

  my ($m, $y);
  if ($month <= 2) {
    $m = int($month+9);
    $y = int($year-1);
  } else {
    $m = int($month-3);
    $y = int($year);
  }
  my $c = int($y/100);
  $y = $y-$c*100;
  my $x1 = int(146097.0*$c/4.0);
  my $x2 = int(1461.0*$y/4.0);
  my $x3 = int((153.0*$m+2.0)/5.0);
  return($x1+$x2+$x3+$day-678882+$ut);
}

# Return the nearest integer (ie round)
sub nint ($) {
  my ($x) = @_;
  ($x<0.0) ? return(ceil($x-0.5)) : return(floor($x+0.5))
}

sub mjd2cal($) {

  my $mjd = shift;

  my $ut = fmod($mjd,1.0);
  if ($ut<0.0) {
    $ut += 1.0;
    $mjd -= 1;
  }

  use integer;  # Calculations require integer division and modulation
  # Get the integral Julian Day number
  my $jd = nint($mjd + 2400001);

  # Do some rather cryptic calculations

  my $temp1 = 4*($jd+((6*(((4*$jd-17918)/146097)))/4+1)/2-37);
  my $temp2 = 10*((($temp1-237)%1461)/4)+5;

  my $year = $temp1/1461-4712;
  my $month =(($temp2/306+2)%12)+1;
  my $day = ($temp2%306)/10+1;

  return($day, $month, $year, $ut);
}
