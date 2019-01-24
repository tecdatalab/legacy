#!/usr/bin/perl

if ($#ARGV != 0) {
  die "Usage: extract_xml_resolution_summary.pl <resolutions file>\n";
}

$total = 0;
$n0 = 0;
$n5 = 0;
$n10 = 0;
$n15 = 0;
$n20 = 0;

open( DATA, $ARGV[0]);
foreach $line (<DATA>) {
  @parts = split(/\s+/, $line);
  $resolution = $parts[1];
  ++$total;
    
  if ($resolution > 20 ) {
    ++$n20; 
  } elsif ($resolution > 15 ) {
    ++$n15; 
  } elsif ($resolution > 10 ) {
    ++$n10; 
  } elsif ($resolution > 5 ) {
    ++$n5; 
  } elsif ($resolution > 0 ) {
    ++$n0; 
  } 
}
close(DATA);

$n_resolution = $n0 + $n5 + $n10 + $n15 + $n20;

print $total."\t";
print $n_resolution."\t";
print $n0."\t";
print $n5."\t";
print $n10."\t";
print $n15."\t";
print $n20."\t";
