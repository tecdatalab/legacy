#!/usr/bin/perl

use strict;
use CGI;
my $q = new CGI;

my $visgrid = '../VisGrid/flat/o';
my $convexhull = '../convexhull/runhull.sh';
my $pdb_dir = '../chain';

#my $pdb_id = $q->param('pdbid') || '7tim-A';
my $pdb_id = $q->param('pdbid');
my $type = $q->param('type');

print $q->header;

my $pdb_file = $pdb_id;
#$pdb_file =~ s/\-/\_/;
$pdb_file .= '.pdb';

my $stdout;
my (@cav,@prot,$flat, $tmpline);
my (@cav_script,@prot_script,$flat_script);
my $count;

my $begin_script = "select all; cpk off; wireframe off; cartoon off; isosurface delete; select all; spacefill on; color [220, 220, 220];";

for $stdout (`$visgrid $pdb_dir/$pdb_file`) {
	$stdout =~ tr/\r\n//d;
	$count++;
	if ($count == 1) {
		$cav[1] = $stdout;
		$stdout =~ s/\s+/\,/g;
		$stdout =~ s/\,$//g;
		$cav_script[1] = "isosurface cav1 select ($stdout) molecular; select $stdout; spacefill off; color \$cav1 red;";
	}
	elsif ($count == 2) {
		$cav[2] = $stdout;
		$stdout =~ s/\s+/\,/g;
		$stdout =~ s/\,$//g;
		$cav_script[2] = "isosurface cav2 select ($stdout) molecular; select $stdout; spacefill off; color \$cav2 green;";
	}
	elsif ($count == 3) {
		$cav[3] = $stdout;
		$stdout =~ s/\s+/\,/g;
		$stdout =~ s/\,$//g;
		$cav_script[3] = "isosurface cav3 select ($stdout) molecular; select $stdout; spacefill off; color \$cav3 blue;";
	}
	elsif ($count == 5) {
		$prot[1] = $stdout;
		$stdout =~ s/\s+/\,/g;
		$stdout =~ s/\,$//g;
		$prot_script[1] = "isosurface prot1 select ($stdout) molecular; select $stdout; spacefill off; color \$prot1 red;";
	}
	elsif ($count == 6) {
		$prot[2] = $stdout;
		$stdout =~ s/\s+/\,/g;
		$stdout =~ s/\,$//g;
		$prot_script[2] = "isosurface prot2 select ($stdout) molecular; select $stdout; spacefill off; color \$prot2 green;";
	}
	elsif ($count == 7) {
		$prot[3] = $stdout;
		$stdout =~ s/\s+/\,/g;
		$stdout =~ s/\,$//g;
		$prot_script[3] = "isosurface prot3 select ($stdout) molecular; select $stdout; spacefill off; color \$prot3 blue;";
	}
	elsif ($count == 9) {
		$flat = $stdout;
		$stdout =~ s/\s+/\,/g;
		$stdout =~ s/\,$//g;
		$flat_script = "isosurface flat1 select ($stdout) molecular; select $stdout; spacefill off; color \$flat1 yellow;";
	}
}


# --- Convexhull to Find Surface Area and Volume ---

open(FILE,">tmp/cavity\_$$.txt");
print FILE "$cav[1]\n$cav[2]\n$cav[3]\n";
close(FILE);

open(FILE,">tmp/protrusion\_$$.txt");
print FILE "$prot[1]\n$prot[2]\n$prot[3]\n";
close(FILE);

open(FILE,">tmp/flat\_$$.txt");
print FILE "$flat\n";
close(FILE);

#chdir("tmp");

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

$count = 0;
for (`$convexhull $pdb_dir/$pdb_file tmp/protrusion\_$$.txt`) {
	next if /^0/;
	chomp;
	@line = split /\s+/;
	$count++;
	$convexhull{$count}{'prot'}{'sa'} = sprintf("%4.3f",$line[1]);
	$convexhull{$count}{'prot'}{'vol'} = sprintf("%4.3f",$line[2]);
}

for (`$convexhull $pdb_dir/$pdb_file tmp/flat\_$$.txt`) {
	next if /^0/;
	chomp;
	@line = split /\s+/;
	$convexhull{1}{'flat'}{'sa'} = sprintf("%4.3f",$line[1]);
	$convexhull{1}{'flat'}{'vol'} = sprintf("%4.3f",$line[2]);
}


# --- Script Output ---

print "visgrid`;";

my $cav_out = <<EOT;
RED:
1st Largest Cavity Residues: (Surface Area = $convexhull{1}{'cav'}{'sa'}\, Volume = $convexhull{1}{'cav'}{'vol'}):
$cav[1]

GREEN:
2nd Largest Cavity Residues: (Surface Area = $convexhull{2}{'cav'}{'sa'}\, Volume = $convexhull{2}{'cav'}{'vol'}):
$cav[2]

BLUE:
3rd Largest Cavity Residues: (Surface Area = $convexhull{3}{'cav'}{'sa'}\, Volume = $convexhull{3}{'cav'}{'vol'}):
$cav[3]
EOT


my $prot_out = <<EOT;
RED:
1st Largest Protrusion Residues: (Surface Area = $convexhull{1}{'prot'}{'sa'}\, Volume = $convexhull{1}{'prot'}{'vol'}):
$prot[1]

GREEN:
2nd Largest Protrusion Residues: (Surface Area = $convexhull{2}{'prot'}{'sa'}\, Volume = $convexhull{2}{'prot'}{'vol'}):
$prot[2]

BLUE:
3rd Largest Prostrusion Residues: (Surface Area = $convexhull{3}{'prot'}{'sa'}\, Volume = $convexhull{3}{'prot'}{'vol'}):
$prot[3]
EOT


my $flat_out = <<EOT;
YELLOW:
Flat Region Residues: (Surface Area = $convexhull{1}{'flat'}{'sa'}\, Volume = $convexhull{1}{'flat'}{'vol'}):
$flat
EOT

print $begin_script, $cav_script[1], $cav_script[2], $cav_script[3], "`;";
print $begin_script, $prot_script[1], $prot_script[2], $prot_script[3], "`;";
print $begin_script, $flat_script, "`;";


# --- VisGrid Output ---

print $cav_out, "`;";
print $prot_out, "`;";
print $flat_out, "`;";

print $type, "`;";


