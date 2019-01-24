#!/usr/bin/perl

#jmolCEPose.cgi
use strict;
use CGI;

my $q = new CGI;

my $cepose = '../cepose/zorro';

my $num = $q->param('n');
my $query = $q->param('q');
my $result_id = $q->param('r');
my $rmsd = $q->param('rmsd');
my $coverage = $q->param('coverage');

print $q->header;

`$cepose tmp/ce\_$result_id\_$num\.txt > tmp/cepose\_$result_id\_$num\.pdb`;

print "jmolcepose`;$query`;$result_id`;$num`;$rmsd`;$coverage";

