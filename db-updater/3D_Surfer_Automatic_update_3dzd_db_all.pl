#!/usr/bin/perl

# perl 3D_Surfer_Automatic_update_3dzd_db_all.pl

chdir("/bio/kihara-web/www/db");

#######################################
#all 3DZD database
#######################################

`cat ./pdbidlistupdate/chainidlist.txt ./pdbidlistupdate/complexidlist.txt ./pdbidlistupdate/domainidlist.txt > ./pdbidlistupdate/allidlist.txt`;
`cat ./pdbidlistupdate/chainid.length.mapping.txt ./pdbidlistupdate/complexid.length.mapping.txt ./pdbidlistupdate/domainid.length.mapping.txt > ./pdbidlistupdate/allid.length.mapping.txt`;
`cat ./pdbidlistupdate/chainid.cath.mapping.txt ./pdbidlistupdate/complexid.cath.mapping.txt ./pdbidlistupdate/domainid.cath.mapping.txt > ./pdbidlistupdate/allid.cath.mapping.txt`;


$inv_directory    = "./inv";
$all_list_file    = "./pdbidlistupdate/allidlist.txt";
$pdb2cath_file    = "./pdbidlistupdate/allid.cath.mapping.txt";
$pdb2length_file    = "./pdbidlistupdate/allid.length.mapping.txt";
$output_filename  = "pdb_latest_all";
$output_dbname  = "$output_filename.db";

# the following two filenames are used as intermediate steps; they need to be deleted at the end
$basic_output_filename     = "$output_filename.basic";
$addlength_output_filename = "$output_filename.length";

print "Integrating zernike descriptors information...\n";
`./create_pdb_latest_basic.pl $inv_directory $all_list_file $basic_output_filename`;

print "Adding main all length...\n";
`./add_length_all.pl $pdb2length_file $basic_output_filename $addlength_output_filename`;

print "Adding CATH codes...\n";
`./add_cath_all.pl $pdb2cath_file $addlength_output_filename $output_filename`;

`rm $basic_output_filename`;
`rm $addlength_output_filename`;

`mv $output_filename ./pdbidlistupdate/$output_filename`;
`./invdb ./pdbidlistupdate/$output_filename ./pdbidlistupdate/$output_dbname`;
`cp ./pdbidlistupdate/$output_dbname /bio/kihara-web/www/3d-surfer/db/$output_dbname`;
`cp ./pdbidlistupdate/$output_filename /bio/kihara-web/www/3d-surfer/db/$output_filename`;

#######################################
#all_cacn 3DZD database
#######################################

$inv_directory    = "./inv_cacn";
$all_list_file    = "./pdbidlistupdate/allidlist.txt";
$pdb2cath_file    = "./pdbidlistupdate/allid.cath.mapping.txt";
$pdb2length_file    = "./pdbidlistupdate/allid.length.mapping.txt";
$output_filename  = "pdb_latest_all_cacn";
$output_dbname  = "$output_filename.db";

# the following two filenames are used as intermediate steps; they need to be deleted at the end
$basic_output_filename     = "$output_filename.basic";
$addlength_output_filename = "$output_filename.length";

print "Integrating zernike descriptors information...\n";
`./create_pdb_latest_basic.pl $inv_directory $all_list_file $basic_output_filename`;

print "Adding main all length...\n";
`./add_length_all.pl $pdb2length_file $basic_output_filename $addlength_output_filename`;

print "Adding CATH codes...\n";
`./add_cath_all.pl $pdb2cath_file $addlength_output_filename $output_filename`;

`rm $basic_output_filename`;
`rm $addlength_output_filename`;

`mv $output_filename ./pdbidlistupdate/$output_filename`;
`./invdb ./pdbidlistupdate/$output_filename ./pdbidlistupdate/$output_dbname`;
`cp ./pdbidlistupdate/$output_dbname /bio/kihara-web/www/3d-surfer/db/$output_dbname`;
`cp ./pdbidlistupdate/$output_filename /bio/kihara-web/www/3d-surfer/db/$output_filename`;

#########
print "Success\n";
