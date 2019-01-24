#!/usr/bin/perl

# perl 3D_Surfer_Automatic_update_3dzd_db.pl
chdir("/bio/kihara-web/www/db");

#######################################
#complex 3DZD database
#######################################

$inv_directory    = "./inv";
$complex_list_file    = "./pdbidlistupdate/complexidlist.txt";
$pdb2cath_file    = "./pdbidlistupdate/complexid.cath.mapping.txt";
$pdb2length_file    = "./pdbidlistupdate/complexid.length.mapping.txt";
$output_filename  = "pdb_latest_complex";
$output_dbname  = "$output_filename.db";

# the following two filenames are used as intermediate steps; they need to be deleted at the end
$basic_output_filename     = "$output_filename.basic";
$addlength_output_filename = "$output_filename.length";

print "Integrating zernike descriptors information...\n";
`./create_pdb_latest_basic.pl $inv_directory $complex_list_file $basic_output_filename`;

print "Adding main complex length...\n";
`./add_length_complex.pl $pdb2length_file $basic_output_filename $addlength_output_filename`;

print "Adding CATH codes...\n";
`./add_cath_complex.pl $pdb2cath_file $addlength_output_filename $output_filename`;

`rm $basic_output_filename`;
`rm $addlength_output_filename`;

`mv $output_filename ./pdbidlistupdate/$output_filename`;
`./invdb ./pdbidlistupdate/$output_filename ./pdbidlistupdate/$output_dbname`;
`cp ./pdbidlistupdate/$output_dbname /bio/kihara-web/www/3d-surfer/db/$output_dbname`;
`cp ./pdbidlistupdate/$output_filename /bio/kihara-web/www/3d-surfer/db/$output_filename`;

#######################################
#complex_cacn 3DZD database
#######################################

$inv_directory    = "./inv_cacn";
$complex_list_file    = "./pdbidlistupdate/complexidlist.txt";
$pdb2cath_file    = "./pdbidlistupdate/complexid.cath.mapping.txt";
$pdb2length_file    = "./pdbidlistupdate/complexid.length.mapping.txt";
$output_filename  = "pdb_latest_complex_cacn";
$output_dbname  = "$output_filename.db";

# the following two filenames are used as intermediate steps; they need to be deleted at the end
$basic_output_filename     = "$output_filename.basic";
$addlength_output_filename = "$output_filename.length";

print "Integrating zernike descriptors information...\n";
`./create_pdb_latest_basic.pl $inv_directory $complex_list_file $basic_output_filename`;

print "Adding main complex length...\n";
`./add_length_complex.pl $pdb2length_file $basic_output_filename $addlength_output_filename`;

print "Adding CATH codes...\n";
`./add_cath_complex.pl $pdb2cath_file $addlength_output_filename $output_filename`;

`rm $basic_output_filename`;
`rm $addlength_output_filename`;

`mv $output_filename ./pdbidlistupdate/$output_filename`;
`./invdb ./pdbidlistupdate/$output_filename ./pdbidlistupdate/$output_dbname`;
`cp ./pdbidlistupdate/$output_dbname /bio/kihara-web/www/3d-surfer/db/$output_dbname`;
`cp ./pdbidlistupdate/$output_filename /bio/kihara-web/www/3d-surfer/db/$output_filename`;

#######################################
#Generate the finishing time
#######################################

($sec,$min,$hour,$day,$month,$year,$wday,$yday,$isdst)=localtime(time());
#$time = substr($year,1,2)."-".($month+1)."-".$day."-".$hour."-".$min."-".$sec;

$year = 2000 + substr($year,1,2);
$month = $month+1;
if($month == 1){
	$month = "Jan";
} elsif($month == 2){
	$month = "Feb";
} elsif($month == 3){
	$month = "Mar";
} elsif($month == 4){
	$month = "Apr";
} elsif($month == 5){
	$month = "May";
} elsif($month == 6){
	$month = "June";
} elsif($month == 7){
	$month = "July";
} elsif($month == 8){
	$month = "Aug";
} elsif($month == 9){
	$month = "Sep";
} elsif($month == 10){
	$month = "Oct";
} elsif($month == 11){
	$month = "Nov";
} elsif($month == 12){
	$month = "Dec";
}

$date_file    = "./pdbidlistupdate/latest.updated.date";
open(OUT, ">$date_file");   
print OUT "$month $day, $year\n";
close(OUT);
`cp ./pdbidlistupdate/latest.updated.date /bio/kihara-web/www/3d-surfer/db/latest.updated.date`;

print "Success\n";
