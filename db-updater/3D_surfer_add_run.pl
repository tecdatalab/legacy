#!/usr/bin/perl

# perl 3D_Surfer_Automatic_update_without_complex.pl
use LWP::Simple;
chdir("/bio/kihara-web/www/db/");

###############################################
#Calculate 3DZD for chain, domain and complex
###############################################
print "Calculate 3DZD ...\n";

@added_array = `ls ./pdbidlistupdate/added.*idlist.txt`;
chomp(@added_array);

foreach my $m ( 0 .. $#added_array ) {

	#read added list
	$list_name = $added_array[$m];

	open( IN, $list_name);
	@pdbfiles = <IN>;
	close(IN);
	chomp(@pdbfiles);
	$pdbcount=$#pdbfiles + 1;

	$group_type = substr($list_name, 24, -10);
	
	#calculate 3DZD of all atoms
	$zd_dir  = "./inv";
	$pdbpath = "chain";
	$inv_log = "inv.added.".$group_type."idlist.log";

	foreach my $i ( 0 .. $#pdbfiles ) {
		$pdb_id_name = $pdbfiles[$i];
		$filename = $pdbfiles[$i] . ".pdb";
		`cp $pdbpath/$filename .`;
   		# try msroll if msms fails
		`perl zernike.pl $filename msms >> ./pdbidlistupdate/$inv_log`; 
		`mv $filename.inv $zd_dir`;
		`rm $pdb_id_name*`;
		if ( ( $i % 100 ) == 0 ) { print "$i :: ($pdbcount total)\n"; }
	}
	print "..done\n";
	
	#calculate 3DZD of cacn
	$zd_cacn_dir  = "./inv_cacn";
	$pdbpath_cacn = "chain_cacn";
	$inv_cacn_log = "inv_cacn.added.".$group_type."idlist.log";

	foreach my $i ( 0 .. $#pdbfiles ) {
		$pdb_id_name = $pdbfiles[$i];
		$filename = $pdbfiles[$i] . ".pdb";
		`cp $pdbpath_cacn/$filename .`;
   		# try msroll if msms fails
		`perl zernike.pl $filename msms >> ./pdbidlistupdate/$inv_cacn_log`;
		`mv $filename.inv $zd_cacn_dir`;
		`rm $pdb_id_name*`;
		if ( ( $i % 100 ) == 0 ) { print "$i :: ($pdbcount total)\n"; }
	}
	print "..done\n";
}

#######################################
#chain/domain/complex 3DZD database
#######################################
print "Generate 3DZD database ...\n";

@group_array = ('chain','domain','complex');

foreach my $i ( 0 .. $#group_array ) {
	$group_type = $group_array[$i];

	$chain_list_file = "./pdbidlistupdate/".$group_type."idlist.txt";
	$added_idlist = "./pdbidlistupdate/added.".$group_type."idlist.txt";
	$idlist_tmp = "./pdbidlistupdate/".$group_type."idlist.txt.tmp";
	$pdb2length_file = "./pdbidlistupdate/".$group_type."id.length.mapping.txt";
	$added_length_mapping = "./pdbidlistupdate/added.".$group_type."idlist.length.txt";
	$length_mapping_tmp = "./pdbidlistupdate/".$group_type."id.length.mapping.txt.tmp";

	$inv_directory = "./inv";
	$inv_cacn_directory = "./inv_cacn";
	$pdb2cath_file    = "./pdbidlistupdate/".$group_type."id.cath.mapping.txt";
	$output_filename = "pdb_latest_".$group_type;
	$output_dbname  = "$output_filename.db";
	$output_cacn_filename    = $output_filename."_cacn";
	$output_dbname  = "$output_filename.db";
	$output_cacn_dbname = "$output_cacn_filename.db";

	# the following two filenames are used as intermediate steps; they need to be deleted at the end
	$basic_output_filename     = "$output_filename.basic";
	$basic_cacn_output_filename     = "$output_cacn_filename.basic";
	$addlength_output_filename = "$output_filename.length";
	$addlength_cacn_output_filename = "$output_cacn_filename.length";

	print "Integrating zernike descriptors information...\n";
	`./create_pdb_latest_basic.pl $inv_directory $chain_list_file $basic_output_filename`;
	`./create_pdb_latest_basic.pl $inv_cacn_directory $chain_list_file $basic_cacn_output_filename`;

	print "Adding main chain length...\n";
	`./add_length_$group_type.pl $pdb2length_file $basic_output_filename $addlength_output_filename`;
	`./add_length_$group_type.pl $pdb2length_file $basic_cacn_output_filename $addlength_cacn_output_filename`;

	print "Adding CATH codes...\n";
	`./add_cath_$group_type.pl $pdb2cath_file $addlength_output_filename $output_filename`;
	`./add_cath_$group_type.pl $pdb2cath_file $addlength_cacn_output_filename $output_cacn_filename`;

	`rm $basic_output_filename`;
	`rm $basic_cacn_output_filename`;
	`rm $addlength_output_filename`;
	`rm $addlength_cacn_output_filename`;

	`mv $output_filename ./pdbidlistupdate/$output_filename`;
	`mv $output_cacn_filename ./pdbidlistupdate/$output_cacn_filename`;
	`./invdb ./pdbidlistupdate/$output_filename ./pdbidlistupdate/$output_dbname`;
	`./invdb ./pdbidlistupdate/$output_cacn_filename ./pdbidlistupdate/$output_cacn_dbname`;
	`cp ./pdbidlistupdate/$output_dbname ../3d-surfer/db/$output_dbname`;
	`cp ./pdbidlistupdate/$output_cacn_dbname ../3d-surfer/db/$output_cacn_dbname`;
	`cp ./pdbidlistupdate/$output_filename ../3d-surfer/db/$output_filename`;
	`cp ./pdbidlistupdate/$output_cacn_filename ../3d-surfer/db/$output_cacn_filename`;
	
}

#######################################
#all 3DZD database
#######################################
@group_array = ('chain','domain','complex');

$chain_list_file = "./pdbidlistupdate/allidlist.txt";
$idlist_tmp = "./pdbidlistupdate/allidlist.txt.tmp";
$pdb2length_file = "./pdbidlistupdate/allid.length.mapping.txt";
$length_mapping_tmp = "./pdbidlistupdate/allid.length.mapping.txt.tmp";

$inv_directory = "./inv";
$inv_cacn_directory = "./inv_cacn";
$pdb2cath_file    = "./pdbidlistupdate/allid.cath.mapping.txt";
$output_filename = "pdb_latest_all";
$output_dbname  = "$output_filename.db";
$output_cacn_filename    = $output_filename."_cacn";
$output_dbname  = "$output_filename.db";
$output_cacn_dbname = "$output_cacn_filename.db";

# the following two filenames are used as intermediate steps; they need to be deleted at the end
$basic_output_filename     = "$output_filename.basic";
$basic_cacn_output_filename     = "$output_cacn_filename.basic";
$addlength_output_filename = "$output_filename.length";
$addlength_cacn_output_filename = "$output_cacn_filename.length";

print "Integrating zernike descriptors information...\n";
`./create_pdb_latest_basic.pl $inv_directory $chain_list_file $basic_output_filename`;
`./create_pdb_latest_basic.pl $inv_cacn_directory $chain_list_file $basic_cacn_output_filename`;

print "Adding main chain length...\n";
`./add_length_all.pl $pdb2length_file $basic_output_filename $addlength_output_filename`;
`./add_length_all.pl $pdb2length_file $basic_cacn_output_filename $addlength_cacn_output_filename`;

print "Adding CATH codes...\n";
`./add_cath_all.pl $pdb2cath_file $addlength_output_filename $output_filename`;
`./add_cath_all.pl $pdb2cath_file $addlength_cacn_output_filename $output_cacn_filename`;

`rm $basic_output_filename`;
`rm $basic_cacn_output_filename`;
`rm $addlength_output_filename`;
`rm $addlength_cacn_output_filename`;

`mv $output_filename ./pdbidlistupdate/$output_filename`;
`mv $output_cacn_filename ./pdbidlistupdate/$output_cacn_filename`;
`./invdb ./pdbidlistupdate/$output_filename ./pdbidlistupdate/$output_dbname`;
`./invdb ./pdbidlistupdate/$output_cacn_filename ./pdbidlistupdate/$output_cacn_dbname`;
`cp ./pdbidlistupdate/$output_dbname ../3d-surfer/db/$output_dbname`;
`cp ./pdbidlistupdate/$output_cacn_dbname ../3d-surfer/db/$output_cacn_dbname`;
`cp ./pdbidlistupdate/$output_filename ../3d-surfer/db/$output_filename`;
`cp ./pdbidlistupdate/$output_cacn_filename ../3d-surfer/db/$output_cacn_filename`;
print "Success\n";
