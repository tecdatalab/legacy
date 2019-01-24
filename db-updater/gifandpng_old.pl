#!/usr/bin/perl
# Script: gif.pl
# Description: Generates 3d animated .gif and the first frame in .png format
# Programmer: Yi Xiong
# email: xiongyi128@gmail.com
# Date: Mar.04.2013

# number of horizontal rotation frames
$hframes = 20;
# number of vertical rotation frames
$vframes = 20;
# image width
$iwidth = 200;
# image height
$iheight = 150;

$pydir = "./pymol";

$curr_id = $ARGV[0];
$pdbfile = $curr_id . ".pdb";

# temporary filename; contains pymol instructions
$temp_file = "./$curr_id\_pymol_gifandpng_tmp.py";

open( TEMP, ">./$temp_file" );
print TEMP "cmd.load(\"./$pdbfile\")\n";
print TEMP "cmd.center(\"all\")\n";
print TEMP "cmd.orient(\"all\")\n";
print TEMP "cmd.hide(\"all\")\n";
print TEMP "cmd.show(\"surface\")\n";
print TEMP "cmd.color(\"white\")\n";

for ( $c = 0 ; $c < $hframes ; $c++ ) {
	print TEMP "cmd.rotate(\"y\", ", ( 360 / $hframes ), ")\n";
	print TEMP "cmd.ray($iwidth, $iheight)\n";
	print TEMP "cmd.png(\"./$curr_id\_", sprintf( "%03d", $c ), "\")\n";
}
for ( $c = $hframes; $c < ( $hframes + $vframes ) ; $c++ ) {
	print TEMP "cmd.rotate(\"x\", ", ( 360 / $vframes ), ")\n";
	print TEMP "cmd.ray($iwidth, $iheight)\n";
	print TEMP "cmd.png(\"./$curr_id\_", sprintf( "%03d", $c ), "\")\n";
}
close(TEMP);

#Generate .png images
#print "Running Pymol..";
`$pydir/pymol -qc $temp_file`;
#print "png..done\n";
`rm $temp_file`;

#Generate .gif
#print "Generating gif..\n";
for ( $c = 0 ; $c < ( $hframes + $vframes ) ; $c++ ) {
	$filenum = sprintf( "%03d", $c );
	$filename = "./$curr_id\_$filenum";
	@convertrun = ( "convert", "-type", "palette", "$filename.png", "$filename.gif" );
	system(@convertrun) == 0 || die("error: convert\n");
}

#print "\tConvert still..";
@convertstil = (
	"convert",      "./$curr_id\_000.png",
	"-transparent", "black",
	"./$curr_id.png"
);

#print "..done\n";
#print "\tConvert animation..";
@convertanim = ("convert", "-loop", "0", "-delay", "0", "-transparent", "black", "-dispose", "background", "-type", "palette", "./$curr_id*.gif", "./$curr_id.gif");

@convertgif = (
	"convert",      "./$curr_id.gif",
	"-transparent", "black",
	"./$curr_id.gif"
);

#print "..done\n";
system(@convertstil) == 0 || die("error: convert\n");
system(@convertanim) == 0 || die("error: convert\n");
system(@convertgif) == 0 || die("error: convert\n");

for ( $c = 0 ; $c < ( $hframes + $vframes ) ; $c++ ) {
	$filenum  = sprintf( "%03d", $c );
	$filename = "./$curr_id\_$filenum";
	`rm $filename.png`;
	`rm $filename.gif`;
}
