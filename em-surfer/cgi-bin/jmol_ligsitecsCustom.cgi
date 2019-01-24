#!/usr/bin/perl

use strict;
use CGI;
my $q = new CGI;

#perl jmol_ligsitecs.cgi
# ./pocketligres.pl 7tim-A ../chain/7tim-A.pdb
my $ligsitecs = './pocketligres.pl';
my $pdb_dir = './tmp';

my $convexhull = '../convexhull/runhull.sh';

#my $pdb_id = $q->param('pdbid') || '7tim-A';
my $pdb_id = $q->param('pdbid');
my $type   = $q->param('type');

#$pdb_id = "7tim-A";
#$type = "cavity";

print $q->header;

my $pdb_file = $pdb_id;

#$pdb_file =~ s/\-/\_/;
$pdb_file .= '.pdb';

my $stdout;
my ( @cav, @prot, $flat, $tmpline );
my ( @cav_script, @prot_script, $flat_script );
my $count;

my $begin_script =
"select all; cpk off; wireframe off; cartoon off; isosurface delete; select all; spacefill on; color [220, 220, 220];";

for $stdout (`$ligsitecs $pdb_id $pdb_dir/$pdb_file`) {
	$stdout =~ tr/\r\n//d;
	$count++;
	if ( $count == 1 ) {
		$cav[1] = $stdout;
		$stdout =~ s/\s+/\,/g;
		$stdout =~ s/\,$//g;
		$cav_script[1] =
"isosurface cav1 select ($stdout) molecular; select $stdout; spacefill off; color \$cav1 red;";
	}
	elsif ( $count == 2 ) {
		$cav[2] = $stdout;
		$stdout =~ s/\s+/\,/g;
		$stdout =~ s/\,$//g;
		$cav_script[2] =
"isosurface cav2 select ($stdout) molecular; select $stdout; spacefill off; color \$cav2 green;";
	}
	elsif ( $count == 3 ) {
		$cav[3] = $stdout;
		$stdout =~ s/\s+/\,/g;
		$stdout =~ s/\,$//g;
		$cav_script[3] =
"isosurface cav3 select ($stdout) molecular; select $stdout; spacefill off; color \$cav3 blue;";
	}
}


# --- Convexhull to Find Surface Area and Volume ---

open(FILE,">tmp/cavity\_$$.txt");
#print FILE "$cav[1]\n$cav[2]\n$cav[3]\n";
my ( @tmpcavarray, @tmparray );
@tmpcavarray = split (/\s+/, $cav[1]);
foreach my $i (0..$#tmpcavarray){
	@tmparray = split (/:/, $tmpcavarray[$i]);
	print FILE "$tmparray[0] ";
}
print FILE "\n";

@tmpcavarray = split (/\s+/, $cav[2]);
foreach my $i (0..$#tmpcavarray){
	@tmparray = split (/:/, $tmpcavarray[$i]);
	print FILE "$tmparray[0] ";
}
print FILE "\n";

@tmpcavarray = split (/\s+/, $cav[3]);
foreach my $i (0..$#tmpcavarray){
	@tmparray = split (/:/, $tmpcavarray[$i]);
	print FILE "$tmparray[0] ";
}
print FILE "\n";

close(FILE);

my @line;
my %convexhull;

# ../convexhull/runhull.sh ../chain/7tim-A.pdb tmp/cavity_7tim-A.txt
my $count = 0;
for (`$convexhull $pdb_dir/$pdb_file tmp/cavity\_$$.txt`) {
	next if /^0/;
	chomp;
	@line = split /\s+/;
	$tmpline=$line[1];
	$count++;
	$convexhull{$count}{'cav'}{'sa'} = sprintf("%4.3f",$line[1]);
	$convexhull{$count}{'cav'}{'vol'} = sprintf("%4.3f",$line[2]);
}

# --- Script Output ---

print "ligsitecs`;";

my $cav_out = <<EOT;
RED:
1st Largest Pocket Residues: (Surface Area = $convexhull{1}{'cav'}{'sa'}\, Volume = $convexhull{1}{'cav'}{'vol'}):
$cav[1]

GREEN:
2nd Largest Pocket Residues: (Surface Area = $convexhull{2}{'cav'}{'sa'}\, Volume = $convexhull{2}{'cav'}{'vol'}):
$cav[2]

BLUE:
3rd Largest Pocket Residues: (Surface Area = $convexhull{3}{'cav'}{'sa'}\, Volume = $convexhull{3}{'cav'}{'vol'}):
$cav[3]
EOT

print $begin_script, $cav_script[1], $cav_script[2], $cav_script[3], "`;";

# --- ligsitecs Output ---

print $cav_out, "`;";
print $type,    "`;";

