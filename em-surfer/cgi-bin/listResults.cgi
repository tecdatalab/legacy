#!/usr/bin/perl

#use strict;
use CGI;
my $q = new CGI;

my $emdbid = $q->param('emdbid');
my $mapfile = $q->param('mapfile');
my $mapvolume = $q->param('mapvolume');
my $volumefilter = $q->param('volumefilter');
my $minresolution = $q->param('minresolution');
my $maxresolution = $q->param('maxresolution');
my $representation = $q->param('representation');

my $query, $isupload;

if ($emdbid) {
	$query = $emdbid;
	$isupload = 0;
} else {
	$query = $mapfile;
	$isupload = 1;
}

my $validatedvolumefilter;

if ($volumefilter && $volumefilter eq "on") {
	$validatedvolumefilter = "-filter";
} else {
	$validatedvolumefilter = "-nofilter";
}

my $validatedmin, $validatedmax;

if ($minresolution && $minresolution >= 0) {
	$validatedmin = $minresolution;
} else {
	$validatedmin = -1;
}

if ($maxresolution && $maxresolution >= 0) {
	$validatedmax = $maxresolution;
} else {
	$validatedmax = -1;
}

my $inv_dir = '../db/recommend_contour/inv';
if ($representation eq "recommend") {
	$inv_dir = '../db/recommend_contour/inv';
}
if ($representation eq "recommendonethird") {
	$inv_dir = '../db/1_contour/inv';
}
if ($representation eq "recommendtwothird") {
	$inv_dir = '../db/2_contour/inv';
}
if ($representation eq "recommendonethirdtwothird") {
	$inv_dir = '../db/1_2_contour/inv';
}
if ($representation eq "recommendonestd") {
	$inv_dir = '../db/1_std_merge/inv';
}

my $db_prog = "./query_top_n.pl";
# Change this value to retrieve more than 50 matches
my $matches_retrieved = 50;

print $q->header;
print "<pre>";
print "***************************************************************\n";
print "Rank\tEMDB_ID\tEUC_D\tRESOLUTION\n";
print "***************************************************************\n";

my $dbrun;

if ($isupload) {
	# We assume that tmp files are stored in /tmp and only the basename is passed
	$dbrun = "$db_prog $matches_retrieved /tmp/$query $validatedvolumefilter $inv_dir $validatedmin $validatedmax -upload $mapvolume";
} else {
	$dbrun = "$db_prog $matches_retrieved $query $validatedvolumefilter $inv_dir $validatedmin $validatedmax";
}

my $rank = 1;
for (`$dbrun`)
{
	print "$rank\t$_";
	$rank = $rank + 1;
}
print "***************************************************************\n";
print "</pre>";

