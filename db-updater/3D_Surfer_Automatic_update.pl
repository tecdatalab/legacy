#!/usr/bin/perl

# perl 3D_Surfer_Automatic_update.pl
use LWP::Simple;
chdir("/bio/kihara-web/www/db");

#######################################
# Download pdb_entry_type file
#######################################
my $file = "./pdbidlistupdate/pdb_entry_type";
@ftparray =`curl ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt`;
 
open(OUT, ">$file");   
print OUT @ftparray;
close(OUT);

#######################################
# Extract existed.pdb file
#######################################
`ls ./pdb | grep -E -o ^[[:alnum:]]{4,4} > ./pdbidlistupdate/existed.pdb`;

`awk '{print $1}' ./pdbidlistupdate/pdb_entry_type | grep -E -o ^[[:alnum:]]{4,4} > ./pdbidlistupdate/allnew.pdb`;

#################################################
# Generate added.pdb file in first not in second
#################################################
open(IN, "./pdbidlistupdate/allnew.pdb");
@firstarray = <IN>;
close(IN);
chomp(@firstarray);

open(IN, "./pdbidlistupdate/existed.pdb");
@secondarray = <IN>;
close(IN);
chomp(@secondarray);

#compare two arrays and derive the common
@firstarray = sort (@firstarray);
@secondarray = sort (@secondarray);

%secondarrayhash=();
foreach my $i (0..$#secondarray) {

  $secondarrayhash{$secondarray[$i]}="1";
}

my $outfile = "./pdbidlistupdate/added.pdb";  
open(OUT,">$outfile");
foreach my $j (0..$#firstarray) {

  if(!exists $secondarrayhash{$firstarray[$j]}){
    print OUT $firstarray[$j]."\n";
  }
}
close(OUT);

#######################################
# Classify_protandpna_idlist
#######################################
my $pdb_entry_type_list = "./pdbidlistupdate/pdb_entry_type";
open(IN, $pdb_entry_type_list) or die "Can't open $pdb_entry_type_list\n";
@entry_type_array = <IN>;
close(IN);

%entry_hash;
foreach my $element (@entry_type_array) {

  $element =~ s/\A\s+|\s+\z//g;
  my @elementarray = split(/\t/,$element); 
  # carb
  if($elementarray[1] eq "prot" || $elementarray[1] eq "prot-nuc") {

    $entry_hash{$elementarray[0]}="1";
  }
}

my $infile = "./pdbidlistupdate/added.pdb";
open(IN, $infile) or die "Can't open $infile\n";
@entry_type_array = <IN>;
close(IN);
chomp(@entry_type_array);

$outfile = "./pdbidlistupdate/added.protandpna.pdb";
open(OUT,">$outfile");
foreach my $element (@entry_type_array) {

	if (exists $entry_hash{$element}) {
		print OUT $element."\n";
	}
}
close(OUT);

#######################################
# Download pdb files
#######################################
$dir = "pdb";
$i   = 0;
$downloadList = "./pdbidlistupdate/added.protandpna.pdb";
open( DATA, $downloadList);
foreach $main (<DATA>) {
	$main =~ s/\A\s+|\s+\z//g;
	my $file = "$dir/$main.pdb";
	unless ( -e $file ) {
		my $url =
"http://www.rcsb.org/pdb/download/downloadFile.do?%20fileFormat=pdb&compression=NO&structureId=$main";
		getstore( $url, $file );
		$i++;
		if ( ( $i > 0 ) && ( $i % 100 == 0 ) ) {
			print "$i\tpdb download\n";
		}
	}
}
close(DATA);

#######################################
#Parse the pdb to separate chains
#######################################
$idlist = "./pdbidlistupdate/added.protandpna.pdb";
open( IN, $idlist );
@listarray = <IN>;
close(IN);
chomp(@listarray);

foreach my $i (0..$#listarray) {
	
	$pdb_id = $listarray[$i];
	$file = "pdb/$pdb_id.pdb";

	`rm chain/$pdb_id*.pdb`;
	open( OUT, ">chain/$pdb_id.pdb" );
	open( IN,  $file );
	lastwhile:while (<IN>) {

		# Skip if non atom
		next if !/^ATOM/;
		if(substr( $_, 0, 6 ) eq "ENDMDL"){
			last lastwhile;
		}
		
		if (
			substr( $_, 17, 3 ) =~ /[A-Z]{3,3}/ &&
			( substr( $_, 16, 1 ) eq " " || substr( $_, 16, 1 ) eq "A" ))
		{

			print OUT $_;
			$chain_id = substr( $_, 21, 1 );

			#output PDB chain file
			$file = "$pdb_id-$chain_id.pdb";
			
			open( FILE, ">>chain/$file" );
			print FILE $_;
			close(FILE);
		}
	}
	close(OUT);
	close(IN);
}

#######################################
#Parse the pdb to separate cacn chains
#######################################
$idlist = "./pdbidlistupdate/added.protandpna.pdb";
open( IN, $idlist );
@listarray = <IN>;
close(IN);
chomp(@listarray);

foreach my $i (0..$#listarray) {
	
	$pdb_id = $listarray[$i];
	$file = "pdb/$pdb_id.pdb";

	`rm chain_cacn/$pdb_id*.pdb`;
	open( OUT, ">chain_cacn/$pdb_id.pdb" );
	open( IN,  $file );
	lastwhile:while (<IN>) {

		# Skip if non atom
		next if !/^ATOM/;
		if(substr( $_, 0, 6 ) eq "ENDMDL"){
			last lastwhile;
		}
		
		$atomname = substr( $_, 12, 4 );
		$atomname =~ s/\A\s+|\s+\z//g;
		if (
			substr( $_, 17, 3 ) =~ /[A-Z]{3,3}/ &&
			( substr( $_, 16, 1 ) eq " " || substr( $_, 16, 1 ) eq "A" ) && ($atomname eq "N" || $atomname eq "CA" || $atomname eq "C"))
		{

			print OUT $_;
			$chain_id = substr( $_, 21, 1 );

			#output PDB chain file
			$file = "$pdb_id-$chain_id.pdb";
			
			open( FILE, ">>chain_cacn/$file" );
			print FILE $_;
			close(FILE);
		}
	}
	close(OUT);
	close(IN);
}

##############################################################################
#Parse the pdb files in chain_cacn directory, and extract the id list
##############################################################################
$idlist = "./pdbidlistupdate/added.protandpna.pdb";
open( IN, $idlist);
@listarray = <IN>;
close(IN);
chomp(@listarray);

%complexidhash; # these complexes include only single chain
open( OUTCHAINIDLIST, ">./pdbidlistupdate/added.chainidlist.txt" );
foreach my $i (0..$#listarray) {
	
	$pdb_id = $listarray[$i];
	$file = "chain_cacn/$pdb_id.pdb";
	
	my %chainidhash;
	open( IN,  $file );
	while (<IN>) {
        
        $chain_id = substr( $_, 21, 1 );
        $chainidhash{$chain_id}=1;
	}
	close(IN);
	
	$size = scalar(keys %chainidhash);
	if($size == 1){
		$complexidhash{$pdb_id}=1;
		foreach my $key (sort keys %chainidhash) {
			print OUTCHAINIDLIST $pdb_id."-".$key."\n";
		}	
	} elsif ($size > 1) {
		foreach my $key (sort keys %chainidhash) {
			print OUTCHAINIDLIST $pdb_id."-".$key."\n";
		}
	}
}
close(OUTCHAINIDLIST);

open( OUTCOMPLEXIDLIST, ">./pdbidlistupdate/added.complexidlist.txt" );
foreach my $i (0..$#listarray) {
	if (!exists $complexidhash{$listarray[$i]}){
		print OUTCOMPLEXIDLIST $listarray[$i]."\n";
	}
}
close(OUTCOMPLEXIDLIST);

#######################################
#Calculate 3DZD for chain
#######################################
$zd_dir  = "./inv";
$pdbpath = "chain";

open( IN, "./pdbidlistupdate/added.chainidlist.txt");
@pdbfiles = <IN>;
close(IN);
chomp(@pdbfiles);

$pdbcount=$#pdbfiles + 1;
foreach my $i ( 0 .. $#pdbfiles ) {
	$filename = $pdbfiles[$i] . ".pdb";
	`cp $pdbpath/$filename .`;
	`perl zernike.pl $filename msms >> ./pdbidlistupdate/inv.added.chainidlist.log`;
	`mv $filename.inv $zd_dir`;
	`rm $filename`;
	if ( ( $i % 100 ) == 0 ) { print "$i :: ($pdbcount total)\n"; }
}
print "..done\n";

#######################################
#Calculate 3DZD for chain_cacn
#######################################
$zd_dir  = "./inv_cacn";
$pdbpath = "chain_cacn";

open( IN, "./pdbidlistupdate/added.chainidlist.txt");
@pdbfiles = <IN>;
close(IN);
chomp(@pdbfiles);

$pdbcount=$#pdbfiles + 1;
foreach my $i ( 0 .. $#pdbfiles ) {
	$filename = $pdbfiles[$i] . ".pdb";
	`cp $pdbpath/$filename .`;
	`perl zernike.pl $filename msms >> ./pdbidlistupdate/inv_cacn.added.chainidlist.log`;
	`mv $filename.inv $zd_dir`;
	`rm $filename`;
	if ( ( $i % 100 ) == 0 ) { print "$i :: ($pdbcount total)\n"; }
}
print "..done\n";

#######################################
#Calculate 3DZD for complex
#######################################
$zd_dir  = "./inv";
$pdbpath = "chain";

open( IN, "./pdbidlistupdate/added.complexidlist.txt");
@pdbfiles = <IN>;
close(IN);
chomp(@pdbfiles);

$pdbcount=$#pdbfiles + 1;
foreach my $i ( 0 .. $#pdbfiles ) {
	$filename = $pdbfiles[$i] . ".pdb";
	`cp $pdbpath/$filename .`;
	`perl zernike.pl $filename msms >> ./pdbidlistupdate/inv.added.complexidlist.log`;
	`mv $filename.inv $zd_dir`;
	`rm $filename`;
	if ( ( $i % 100 ) == 0 ) { print "$i :: ($pdbcount total)\n"; }
}
print "..done\n";

#######################################
#Calculate 3DZD for complex_cacn
#######################################
$zd_dir  = "./inv_cacn";
$pdbpath = "chain_cacn";

open( IN, "./pdbidlistupdate/added.complexidlist.txt");
@pdbfiles = <IN>;
close(IN);
chomp(@pdbfiles);

$pdbcount=$#pdbfiles + 1;
foreach my $i ( 0 .. $#pdbfiles ) {
	$filename = $pdbfiles[$i] . ".pdb";
	`cp $pdbpath/$filename .`;
	`perl zernike.pl $filename msms >> ./pdbidlistupdate/inv_cacn.added.complexidlist.log`;
	`mv $filename.inv $zd_dir`;
	`rm $filename`;
	if ( ( $i % 100 ) == 0 ) { print "$i :: ($pdbcount total)\n"; }
}
print "..done\n";

#######################################
#Calculate residue length for chain
#######################################
$dir = "./chain_cacn";
$infile = "./pdbidlistupdate/added.chainidlist.txt";
$outfile = "./pdbidlistupdate/added.chainidlist.length.txt";
open( OUT, ">$outfile" );
open( IN,  $infile );
while (<IN>) {

	chomp($_);
	$file = $dir . "/" . $_ . ".pdb";

	if ( -e $file ) {

		$length = `grep ^ATOM.........CA $file | wc -l`;
		chomp($length);
		print OUT $_ . "\t" . $length . "\n";
	}
}
close(IN);
close(OUT);

#######################################
#Calculate residue length for complex
#######################################
$dir = "./chain_cacn";
$infile = "./pdbidlistupdate/added.complexidlist.txt";
$outfile = "./pdbidlistupdate/added.complexidlist.length.txt";
open( OUT, ">$outfile" );
open( IN,  $infile );
while (<IN>) {

	chomp($_);
	$file = $dir . "/" . $_ . ".pdb";

	if ( -e $file ) {

		$length = `grep ^ATOM.........CA $file | wc -l`;
		chomp($length);
		print OUT $_ . "\t" . $length . "\n";
	}
}
close(IN);
close(OUT);

#######################################
#Build chain/complex id CATH mapping 
#######################################


#######################################
#chain 3DZD database
#######################################

`cat ./pdbidlistupdate/chainidlist.txt ./pdbidlistupdate/added.chainidlist.txt > ./pdbidlistupdate/chainidlist.txt.tmp`;
`mv ./pdbidlistupdate/chainidlist.txt.tmp ./pdbidlistupdate/chainidlist.txt`;

`cat ./pdbidlistupdate/chainid.length.mapping.txt ./pdbidlistupdate/added.chainidlist.length.txt > ./pdbidlistupdate/chainid.length.mapping.txt.tmp`;
`mv ./pdbidlistupdate/chainid.length.mapping.txt.tmp ./pdbidlistupdate/chainid.length.mapping.txt`;

$inv_directory    = "./inv";
$chain_list_file    = "./pdbidlistupdate/chainidlist.txt";
$pdb2cath_file    = "./pdbidlistupdate/chainid.cath.mapping.txt";
$pdb2length_file    = "./pdbidlistupdate/chainid.length.mapping.txt";
$output_filename  = "pdb_latest_chain";
$output_dbname  = "$output_filename.db";

# the following two filenames are used as intermediate steps; they need to be deleted at the end
$basic_output_filename     = "$output_filename.basic";
$addlength_output_filename = "$output_filename.length";

print "Integrating zernike descriptors information...\n";
`./create_pdb_latest_basic.pl $inv_directory $chain_list_file $basic_output_filename`;

print "Adding main chain length...\n";
`./add_length_chain.pl $pdb2length_file $basic_output_filename $addlength_output_filename`;

print "Adding CATH codes...\n";
`./add_cath_chain.pl $pdb2cath_file $addlength_output_filename $output_filename`;

`rm $basic_output_filename`;
`rm $addlength_output_filename`;

`mv $output_filename ./pdbidlistupdate/$output_filename`;
`./invdb ./pdbidlistupdate/$output_filename ./pdbidlistupdate/$output_dbname`;
`cp ./pdbidlistupdate/$output_dbname /bio/kihara-web/www/3d-surfer/db/$output_dbname`;
`cp ./pdbidlistupdate/$output_filename /bio/kihara-web/www/3d-surfer/db/$output_filename`;

#######################################
#chain_cacn 3DZD database
#######################################

$inv_directory    = "./inv_cacn";
$chain_list_file    = "./pdbidlistupdate/chainidlist.txt";
$pdb2cath_file    = "./pdbidlistupdate/chainid.cath.mapping.txt";
$pdb2length_file    = "./pdbidlistupdate/chainid.length.mapping.txt";
$output_filename  = "pdb_latest_chain_cacn";
$output_dbname  = "$output_filename.db";

# the following two filenames are used as intermediate steps; they need to be deleted at the end
$basic_output_filename     = "$output_filename.basic";
$addlength_output_filename = "$output_filename.length";

print "Integrating zernike descriptors information...\n";
`./create_pdb_latest_basic.pl $inv_directory $chain_list_file $basic_output_filename`;

print "Adding main chain length...\n";
`./add_length_chain.pl $pdb2length_file $basic_output_filename $addlength_output_filename`;

print "Adding CATH codes...\n";
`./add_cath_chain.pl $pdb2cath_file $addlength_output_filename $output_filename`;

`rm $basic_output_filename`;
`rm $addlength_output_filename`;

`mv $output_filename ./pdbidlistupdate/$output_filename`;
`./invdb ./pdbidlistupdate/$output_filename ./pdbidlistupdate/$output_dbname`;
`cp ./pdbidlistupdate/$output_dbname /bio/kihara-web/www/3d-surfer/db/$output_dbname`;
`cp ./pdbidlistupdate/$output_filename /bio/kihara-web/www/3d-surfer/db/$output_filename`;

#######################################
#complex 3DZD database
#######################################

`cat ./pdbidlistupdate/complexidlist.txt ./pdbidlistupdate/added.complexidlist.txt > ./pdbidlistupdate/complexidlist.txt.tmp`;
`mv ./pdbidlistupdate/complexidlist.txt.tmp ./pdbidlistupdate/complexidlist.txt`;

`cat ./pdbidlistupdate/complexid.length.mapping.txt ./pdbidlistupdate/added.complexidlist.length.txt > ./pdbidlistupdate/complexid.length.mapping.txt.tmp`;
`mv ./pdbidlistupdate/complexid.length.mapping.txt.tmp ./pdbidlistupdate/complexid.length.mapping.txt`;

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
#all 3DZD database
#######################################

`cat ./pdbidlistupdate/allidlist.txt ./pdbidlistupdate/added.chainidlist.txt > ./pdbidlistupdate/allidlist.txt.tmp`;
`mv ./pdbidlistupdate/allidlist.txt.tmp ./pdbidlistupdate/allidlist.txt`;

`cat ./pdbidlistupdate/allidlist.txt ./pdbidlistupdate/added.complexidlist.txt > ./pdbidlistupdate/allidlist.txt.tmp`;
`mv ./pdbidlistupdate/allidlist.txt.tmp ./pdbidlistupdate/allidlist.txt`;

`cat ./pdbidlistupdate/allid.length.mapping.txt ./pdbidlistupdate/added.chainidlist.length.txt > ./pdbidlistupdate/allid.length.mapping.txt.tmp`;
`mv ./pdbidlistupdate/allid.length.mapping.txt.tmp ./pdbidlistupdate/allid.length.mapping.txt`;

`cat ./pdbidlistupdate/allid.length.mapping.txt ./pdbidlistupdate/added.complexidlist.length.txt > ./pdbidlistupdate/allid.length.mapping.txt.tmp`;
`mv ./pdbidlistupdate/allid.length.mapping.txt.tmp ./pdbidlistupdate/allid.length.mapping.txt`;


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

#######################################
#Generate gifs/pngs for chain
#######################################

$pdbpath = "chain";

open( IN, "./pdbidlistupdate/added.chainidlist.txt");
@idarray = <IN>;
close(IN);
chomp(@idarray);

#print "Generating gifs/pngs..\n";
foreach my $i (0..$#idarray) {
	
	$filename = $idarray[$i].".pdb";
	`cp $pdbpath/$filename .`;
	`perl gifandpng.pl $idarray[$i]`;
	`rm $filename`;
	
	`mv $idarray[$i].gif anim`;
	`mv $idarray[$i].png img`;
}

#######################################
#Generate gifs/pngs for complex
#######################################

$pdbpath = "chain";

open( IN, "./pdbidlistupdate/added.complexidlist.txt");
@idarray = <IN>;
close(IN);
chomp(@idarray);

#print "Generating gifs/pngs..\n";
foreach my $i (0..$#idarray) {
	
	$filename = $idarray[$i].".pdb";
	`cp $pdbpath/$filename .`;
	`perl gifandpng.pl $idarray[$i]`;
	`rm $filename`;
	
	`mv $idarray[$i].gif anim`;
	`mv $idarray[$i].png img`;
}

=pod
=cut

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
