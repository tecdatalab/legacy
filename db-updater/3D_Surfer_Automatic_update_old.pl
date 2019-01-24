#!/usr/bin/perl

# perl 3D_Surfer_Automatic_update_without_complex.pl
use LWP::Simple;
chdir("/bio/kihara-web/www/db/");

#####################################################
# Download pdb_entry_type file and modified pdb file
#####################################################
my $file = "./pdbidlistupdate/pdb_entry_type";
@ftparray =`curl ftp://ftp.wwpdb.org/pub/pdb/derived_data/pdb_entry_type.txt`;
 
open(OUT, ">$file");   
print OUT @ftparray;
close(OUT);

my $modified = "./pdbidlistupdate/modified.pdb";
@ftparray =`curl ftp://ftp.wwpdb.org/pub/pdb/data/status/latest/modified.pdb`;
 
open(OUT, ">$modified");   
print OUT @ftparray;
close(OUT);

my $obsolete = "./pdbidlistupdate/obsolete.pdb";
@ftparray =`curl ftp://ftp.wwpdb.org/pub/pdb/data/status/latest/obsolete.pdb`;
 
open(OUT, ">$obsolete");   
print OUT @ftparray;
close(OUT);

######################################
# Download latest CATH classification
######################################
sub download_cath {
  my ($local_file, $ftp_file) = @_;
  my $wget_result = system("wget -O $local_file -- $ftp_file");
  if ($wget_result != 0) {
    print "Error downloading $ftp_file\n";
    `rm -f $local_file`;
  }
}

my $cathchain = "./pdbidlistupdate/cath-chain-list.txt";
my $ftpchain = "ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/cath-classification-data/cath-chain-list-v*.txt";
download_cath($cathchain, $ftpchain);

my $domainseq = "./pdbidlistupdate/cath-domain-boundaries-seqreschopping.txt";
my $ftpseq = "ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/cath-classification-data/cath-domain-boundaries-seqreschopping-v*.txt";
download_cath($domainseq, $ftpseq);

my $domainboundary = "./pdbidlistupdate/cath-domain-boundaries.txt";
my $ftpboundary = "ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/cath-classification-data/cath-domain-boundaries-v*.txt";
download_cath($domainboundary, $ftpboundary);

my $domaindescribe= "./pdbidlistupdate/cath-domain-description-file.txt";
my $ftpdescribe = "ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/cath-classification-data/cath-domain-description-file-v*.txt";
download_cath($domaindescribe, $ftpdescribe);

my $domainlist = "./pdbidlistupdate/cath-domain-list.txt";
my $ftplist= "ftp://orengoftp.biochem.ucl.ac.uk/cath/releases/latest-release/cath-classification-data/cath-domain-list-v*.txt";
download_cath($domainlist, $ftplist);

#######################################
# Extract existed.pdb file
#######################################
`ls ./pdb | grep -E -o ^[[:alnum:]]{4,4} > ./pdbidlistupdate/existed.pdb`;

`awk '{print $1}' ./pdbidlistupdate/pdb_entry_type | grep -E -o ^[[:alnum:]]{4,4} > ./pdbidlistupdate/allnew.pdb`;

#############################################################################
# Generate added.pdb file in first not in second, and remove obsolete files
#############################################################################
print "Generate added.pdb ...\n";

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
my @new_pdb_array;

foreach my $j (0..$#firstarray) {

  if(!exists $secondarrayhash{$firstarray[$j]}){
	push( @new_pdb_array, $firstarray[$j] );
  }
}

$addlist = join( "\n", @new_pdb_array );
print OUT $addlist;
print OUT "\n";
close(OUT);

#add modified.pdb to added.pdb
`cat ./pdbidlistupdate/added.pdb ./pdbidlistupdate/modified.pdb | sort -u > ./pdbidlistupdate/added.pdb.edit`;
`mv ./pdbidlistupdate/added.pdb.edit ./pdbidlistupdate/added.pdb`;

print "Success\n";

#update list files for modified pdb
open(IN, "./pdbidlistupdate/modified.pdb");
@modifyarray = <IN>;
close(IN);
chomp(@modifyarray);

@group_array = ('chain','domain','complex', 'all');

foreach my $i (0..$#modifyarray) {
	
	$pdb_id = $modifyarray[$i];

	#clean up list files
	foreach my $i ( 0 .. $#group_array ) {
		$group_type = $group_array[$i];

		$chain_list_file = "./pdbidlistupdate/".$group_type."idlist.txt";
		$idlist_tmp = "./pdbidlistupdate/".$group_type."idlist.txt.tmp";
		`sed '/$pdb_id/d' $chain_list_file > $idlist_tmp`;
		`mv $idlist_tmp $chain_list_file`;

		$pdb2length_file = "./pdbidlistupdate/".$group_type."id.length.mapping.txt";
		$length_mapping_tmp = "./pdbidlistupdate/".$group_type."id.length.mapping.txt.tmp";
		`sed '/$pdb_id/d' $pdb2length_file > $length_mapping_tmp`;
		`mv $length_mapping_tmp $pdb2length_file`;
	}

	`sed '/$pdb_id/d' ./pdbidlistupdate/complex_chain.txt > ./pdbidlistupdate/complex_chain.txt.tmp`;
	`mv ./pdbidlistupdate/complex_chain.txt.tmp ./pdbidlistupdate/complex_chain.txt`;
}


#remove obsolete files
print "Remove obsolete files ...\n";

open(IN, "./pdbidlistupdate/obsolete.pdb");
@obsoletearray = <IN>;
close(IN);
chomp(@obsoletearray);

foreach my $i (0..$#obsoletearray) {
	
	$pdb_id = $obsoletearray[$i];
	`mv chain/$pdb_id*.pdb obsolete/chain/`; 
	`mv chain_cacn/$pdb_id*.pdb obsolete/chain_cacn/`; 
	`mv inv/$pdb_id*.pdb.inv obsolete/inv/`; 
	`mv inv_cacn/$pdb_id*.pdb.inv obsolete/inv_cacn/`; 
	`mv anim/$pdb_id*.gif obsolete/anim/`; 
	`mv img/$pdb_id*.png obsolete/img/`;

	#clean up list files
	foreach my $i ( 0 .. $#group_array ) {
		$group_type = $group_array[$i];

		$chain_list_file = "./pdbidlistupdate/".$group_type."idlist.txt";
		$idlist_tmp = "./pdbidlistupdate/".$group_type."idlist.txt.tmp";
		`sed '/$pdb_id/d' $chain_list_file > $idlist_tmp`;
		`mv $idlist_tmp $chain_list_file`;

		$pdb2length_file = "./pdbidlistupdate/".$group_type."id.length.mapping.txt";
		$length_mapping_tmp = "./pdbidlistupdate/".$group_type."id.length.mapping.txt.tmp";
		`sed '/$pdb_id/d' $pdb2length_file > $length_mapping_tmp`;
		`mv $length_mapping_tmp $pdb2length_file`;
	}

	`sed '/$pdb_id/d' ./pdbidlistupdate/complex_chain.txt > ./pdbidlistupdate/complex_chain.txt.tmp`;
	`mv ./pdbidlistupdate/complex_chain.txt.tmp ./pdbidlistupdate/complex_chain.txt`;
}

print "Success\n";

#######################################
# Classify_protandpna_idlist
#######################################
print "Generate protandpna list ...\n";

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
my @added_pdb_array;

foreach my $element (@entry_type_array) {

	if (exists $entry_hash{$element}) {
		push( @added_pdb_array, $element );
	}
}
$add_pdb_list = join( "\n", @added_pdb_array );
print OUT $add_pdb_list;
print OUT "\n";
close(OUT);
print "Success\n";

#######################################
# Download pdb files
#######################################
print "Download pdb files ...\n";

$dir = "pdb";
$i   = 0;
$downloadList = "./pdbidlistupdate/added.protandpna.pdb";
open( DATA, $downloadList);
foreach $main (<DATA>) {
	$main =~ s/\A\s+|\s+\z//g;
	my $file = "$dir/$main.pdb";
	my $url = "https://files.rcsb.org/download/$main.pdb";
	getstore( $url, $file );
	$i++;
	if ( ( $i > 0 ) && ( $i % 100 == 0 ) ) {
			print "$i\tpdb download\n";
	}
}
close(DATA);
print "Success\n";

#######################################
#Parse the pdb to separate chains
#######################################
print "Parse the pdb to separate chains ...\n";

$idlist = "./pdbidlistupdate/added.protandpna.pdb";
open( IN, $idlist );
@listarray = <IN>;
close(IN);
chomp(@listarray);

#generate complex list
open( OUTCHAINIDLIST, ">./pdbidlistupdate/complexpdblist.txt" );
my @complexarray;

foreach my $i (0..$#listarray) {
	
	$pdb_id = $listarray[$i];
	$file = "pdb/$pdb_id.pdb";
	my @array;
	my %hash;
	my $buffer;

	`rm chain/$pdb_id*.pdb`;
	`rm inv/$pdb_id*.pdb.inv`;

	if (-e $file) {
		@array = `cat $file`;

		lastwhile: foreach (@array) {

			if(substr( $_, 0, 6 ) eq "ENDMDL"){
				last lastwhile;
			}
		
			# Skip if non atom
			next if !/^ATOM/;
		
			if (
				substr( $_, 17, 3 ) =~ /[A-Z]{3,3}/ &&
				( substr( $_, 16, 1 ) eq " " || substr( $_, 16, 1 ) eq "A" ))
			{

				$chain_id = substr( $_, 21, 1 );

				#output PDB chain file
				$file = "$pdb_id-$chain_id.pdb";

				if (exists $hash{$file}){
            		$hash{$file} .= $_;
				}else{
            		$hash{$file} = $_;
           		}

				$buffer .= $_;
			}
		}

		foreach $key ( keys %hash) {
			open( FILE, ">chain/$key" );
			print FILE $hash{$key};
			close(FILE);
		}
	
		if (scalar(keys %hash) >1) {
			open(OUT, ">chain/$pdb_id.pdb");
			print OUT $buffer;
			close(OUT);
			push ( @complexarray, $pdb_id );
		}
	}	
}
$complexidlist = join( "\n", @complexarray );
print OUTCHAINIDLIST $complexidlist;
print OUTCHAINIDLIST "\n";
close(OUTCHAINIDLIST);
print "Success\n";

#######################################
#Parse the pdb to separate cacn chains
#######################################
print "Parse the pdb to separate cacn chains ...\n";

$idlist = "./pdbidlistupdate/added.protandpna.pdb";
open( IN, $idlist );
@listarray = <IN>;
close(IN);
chomp(@listarray);

foreach my $i (0..$#listarray) {
	
	$pdb_id = $listarray[$i];
	$file = "pdb/$pdb_id.pdb";
	my @array;
	my %hash;

	`rm chain_cacn/$pdb_id*.pdb`;
	`rm inv_cacn/$pdb_id*.pdb.inv`;
	`rm anim/$pdb_id*.gif`;
	`rm img/$pdb_id*.png`;

	if (-e $file) {
		@array = `cat $file`;

		lastwhile: foreach (@array) {

			if(substr( $_, 0, 6 ) eq "ENDMDL"){
				last lastwhile;
			}
		
			# Skip if non atom
			next if !/^ATOM/;
		
			$atomname = substr( $_, 12, 4 );
			$atomname =~ s/\A\s+|\s+\z//g;
			if (
				substr( $_, 17, 3 ) =~ /[A-Z]{3,3}/ &&
				( substr( $_, 16, 1 ) eq " " || substr( $_, 16, 1 ) eq "A" ) && ($atomname eq "N" || $atomname eq "CA" || $atomname eq "C"))
			{

				$chain_id = substr( $_, 21, 1 );

				#output PDB chain file
				$file = "$pdb_id-$chain_id.pdb";
			
				if (exists $hash{$file}){
            		$hash{$file} .= $_;
				}else{
            		$hash{$file} = $_;
            	}
			}
		}
		#close(OUT);
		foreach $key ( keys %hash) {
			open( FILE, ">chain_cacn/$key" );
			print FILE $hash{$key};
			close(FILE);
		}
	}
}
print "Success\n";

#######################################
#Parse the pdb to domains
#######################################
print "Parse the pdb to domains ...\n";

$idlist = "./pdbidlistupdate/added.protandpna.pdb";
open( IN, $idlist );
@listarray = <IN>;
close(IN);
chomp(@listarray);

$dir = "./pdb";
$domaindir = "./chain";

foreach my $i ( 0 .. $#listarray ) {
	$pdbid = $listarray[$i];
	
	#check if pdb exists in domain records
	@inarray = `grep $pdbid ./pdbidlistupdate/cath-domain-boundaries-seqreschopping.txt`;
 
	if (@inarray) {
		foreach my $m ( 0 .. $#inarray ) {
			chomp( $inarray[$m] );
			$chainid   = substr( $inarray[$m], 4, 1 );
			$pdbdomainid = $pdbid."-".$chainid."-".substr( $inarray[$m], 5, 2);
			$pdbFile = "$dir/$pdbid.pdb";
			$chainFile = "$domaindir/$pdbid"."-".$chainid.".pdb";
			$domainpdbFile = "$domaindir/$pdbdomainid.pdb";
	
			@linearray=split (/\t/,$inarray[$m]);
			$segmentnum=$linearray[2];

			if ((-e $chainFile) && (-e $pdbFile)) {
	
				open( OUTPDB, ">$domainpdbFile" );
				open( INPDB, $pdbFile );
				lastwhile:while (<INPDB>) {

					if(substr( $_, 0, 6 ) eq "ENDMDL"){
						last lastwhile;
					}
			
					# Skip if non atom
					next if !/^ATOM/;

					if (
						substr( $_, 21, 1 ) eq $chainid &&
						( substr( $_, 16, 1 ) eq " " || substr( $_, 16, 1 ) eq "A" ) ) {

						$resSeq=int (substr( $_, 22, 4 ));
						for (my $j=3;$j < $#linearray;$j = $j+2){
							if($resSeq >= $linearray[$j] && $resSeq <= $linearray[$j+1]) {
								print OUTPDB $_;
							}	
						}
					}
				}
			}
		}
	}
	close(OUTPDB);
	close(INPDB);
}
print "Success\n";

#######################################
#Parse the pdb to cacn domains
#######################################
print "Parse the pdb to cacn domains ...\n";

$idlist = "./pdbidlistupdate/added.protandpna.pdb";
open( IN, $idlist );
@listarray = <IN>;
close(IN);
chomp(@listarray);

$dir = "./pdb";
$domaindir = "./chain_cacn";

foreach my $i ( 0 .. $#listarray ) {
	$pdbid = $listarray[$i];
	
	#check if pdb exists in domain records
	@inarray = `grep $pdbid ./pdbidlistupdate/cath-domain-boundaries-seqreschopping.txt`;
 
	if (@inarray) {
		foreach my $m ( 0 .. $#inarray ) {
			chomp( $inarray[$m] );
			$pdbid       = substr( $inarray[$m], 0, 4 );
			$chainid     = substr( $inarray[$m], 4, 1 );
			$pdbdomainid = $pdbid . "-" . $chainid . "-" . substr( $inarray[$m], 5, 2 );
			$pdbFile     = "$dir/$pdbid.pdb";
			$chainFile = "$domaindir/$pdbid"."-".$chainid.".pdb";
			$domainpdbFile = "$domaindir/$pdbdomainid.pdb";

			@linearray  = split( /\t/, $inarray[$m] );
			$segmentnum = $linearray[2];

			if ((-e $chainFile) && (-e $pdbFile)) {
				open( OUTPDB, ">$domainpdbFile" );
				open( INPDB,  $pdbFile );
  				lastwhile: while (<INPDB>) {

  					if ( substr( $_, 0, 6 ) eq "ENDMDL" ) {
						last lastwhile;
					}

					# Skip if non atom
					next if !/^ATOM/;

					if ( substr( $_, 21, 1 ) eq $chainid
						&& ( substr( $_, 16, 1 ) eq " " || substr( $_, 16, 1 ) eq "A" ) )
					{

						$str = substr( $_, 12, 4 );
						$str =~ s/\A\s+|\s+\z//g;
						if ( $str eq "N" || $str eq "CA" || $str eq "C" ) {

							$resSeq = int( substr( $_, 22, 4 ) );
							for ( my $j = 3 ; $j < $#linearray ; $j = $j + 2 ) {
								if (   $resSeq >= $linearray[$j]
									&& $resSeq <= $linearray[ $j + 1 ] )
								{
									print OUTPDB $_;
								}
							}
						}
					}
				}
			}
		}
	}
	close(OUTPDB);
	close(INPDB);
}
print "Success\n";

#######################################
#Parse the pdb to complexes
#######################################
print "Parse the pdb to complexes ...\n";

#generate contacting pairs
`python complex_multiprocess.py 2 > ./pdbidlistupdate/complex.log`;
#generate complex component
`./clustercomplex.pl`;

@pdbarray = `ls ./distance/*.distance`;
chomp @pdbarray;
foreach my $element (@pdbarray) { 
	$pdb_id = substr($element, 11, 4);
	`rm chain/$pdb_id.pdb`;
}
`rm ./distance/*.distance`;

#update complex chain list
`cat ./pdbidlistupdate/complex_chain.txt ./pdbidlistupdate/complexcomponentidchainlist.txt | sort -u > ./pdbidlistupdate/complex_chain.txt.n`;
`mv ./pdbidlistupdate/complex_chain.txt.n ./pdbidlistupdate/complex_chain.txt`;

$idchainlist = "./pdbidlistupdate/complexcomponentidchainlist.txt";
open( IN, $idchainlist );
@idchainlistarray = <IN>;
close(IN);
chomp(@idchainlistarray);

foreach my $i (0 .. $#idchainlistarray) {
	
	my @tmparr=split(/\s+/, $idchainlistarray[$i]);
	$pdbid=substr($tmparr[0], 0, 4);
	$complexid = $tmparr[0];
	$file = "./pdb/$pdbid.pdb";

	open( OUT, ">./chain/$complexid.pdb" );
	open( IN,  $file );
	lastwhile:while (<IN>) {
		if(substr( $_, 0, 6 ) eq "ENDMDL"){
			last lastwhile;
		}
		# Skip if non atom
		next if !/^ATOM/;
		
		if (
			substr( $_, 17, 3 ) =~ /[A-Z]{3,3}/ &&
			( substr( $_, 16, 1 ) eq " " || substr( $_, 16, 1 ) eq "A" ) && index($tmparr[1], substr( $_, 21, 1 )) >=0 ) {

			print OUT $_;
		}
	}
	close(OUT);
	close(IN);
}    
print "Success\n";

#######################################
#Parse the pdb to cacn complexes
#######################################
print "Parse the pdb to cacn complexes ...\n";

$idchainlist = "./pdbidlistupdate/complexcomponentidchainlist.txt";
open( IN, $idchainlist );
@idchainlistarray = <IN>;
close(IN);
chomp(@idchainlistarray);

my @complexarray;
foreach my $i (0 .. $#idchainlistarray) {
	
	my @tmparr=split(/\s+/, $idchainlistarray[$i]);
	$pdbid=substr($tmparr[0], 0, 4);
	$complexid = $tmparr[0];
	$file = "./pdb/$pdbid.pdb";
	push( @complexarray, $complexid );

	open( OUT, ">./chain_cacn/$complexid.pdb" );
	open( IN,  $file );
	lastwhile:while (<IN>) {
		if(substr( $_, 0, 6 ) eq "ENDMDL"){
			last lastwhile;
		}
		# Skip if non atom
		next if !/^ATOM/;
		
		$atomname = substr( $_, 12, 4 );
		$atomname =~ s/\A\s+|\s+\z//g;
		
		if (
			substr( $_, 17, 3 ) =~ /[A-Z]{3,3}/ &&
			( substr( $_, 16, 1 ) eq " " || substr( $_, 16, 1 ) eq "A" ) && ($atomname eq "N" || $atomname eq "CA" || $atomname eq "C") && index($tmparr[1], substr( $_, 21, 1 )) >=0 ) {

			print OUT $_;
		}
	}
	close(OUT);
	close(IN);
}
print "Success\n";

################################################################################################
#Move pdb to short_seq if residue number<10 and extract id list for chain, domain and complex
################################################################################################
print "Move pdb to short_seq and extract id list ...\n";

$idlist = "./pdbidlistupdate/added.protandpna.pdb";
open( IN, $idlist);
@listarray = <IN>;
close(IN);
chomp(@listarray);

open( OUTCHAINIDLIST, ">./pdbidlistupdate/added.chainidlist.txt" );
open( CHAINLENGTHLIST, ">./pdbidlistupdate/added.chainidlist.length.txt" );
open( OUTDOMAINIDLIST, ">./pdbidlistupdate/added.domainidlist.txt" );
open( DOMAINLENGTHLIST, ">./pdbidlistupdate/added.domainidlist.length.txt" );
open( OUTCOMPLEXIDLIST, ">./pdbidlistupdate/added.complexidlist.txt" );
open( COMPLEXLENGTHLIST, ">./pdbidlistupdate/added.complexidlist.length.txt" );

my @chainarray;
my @chain_length;
my @domainarray;
my @domain_length;
my @complexarray;
my @complex_length;

foreach my $i (0..$#listarray) {
	
	$pdb_id = $listarray[$i];
	@file_array = `ls chain_cacn/$pdb_id*.pdb`;
	chomp(@file_array);
	
	$size = scalar(@file_array);
	if ($size > 1) {
		foreach my $m (0 .. $#file_array) {

			$filename = $file_array[$m];
   			$length = `grep ^ATOM.........CA $filename | wc -l`;
			chomp($length);

			$pdb_name = substr($filename, 11, -4);
			$pdb_length = $pdb_name. "\t" . $length;
	
    		# move to short_seq if residue number<10
    		if ($length < 10) {
        		`mv chain_cacn/$pdb_name.pdb short_seq/chain_cacn/`;
        		`mv chain/$pdb_name.pdb short_seq/chain/`;
				next;
    		}
			
			#add name and length to corresponding array
			if ((substr( $pdb_name, 6, 1) eq "") && (length($pdb_name) > 4)) {
				push (@chainarray, $pdb_name);
				push (@chain_length, $pdb_length);
			} elsif (substr( $pdb_name, 6, 1) eq "-") {
				push (@domainarray, $pdb_name);
				push (@domain_length, $pdb_length);
			} else {
				push (@complexarray, $pdb_name);
				push (@complex_length, $pdb_length);
			}
		}
	} elsif ($size == 1) {
		$filename = $file_array[0];
   		$length = `grep ^ATOM.........CA $filename | wc -l`;
		chomp($length);
	
		$pdb_name = substr($filename, 11, -4);
		$pdb_length = $pdb_name. "\t" . $length;		
		
    	# move to short_seq if residue number<10
    	if ($length < 10) {
        	`mv chain_cacn/$pdb_name.pdb short_seq/chain_cacn/`;
        	`mv chain/$pdb_name.pdb short_seq/chain/`;
			next;
    	}
			
		#add name and length to corresponding array
		if (substr( $pdb_name, 6, 1) eq "" && (length($pdb_name) > 4)) {
			push (@chainarray, $pdb_name);
			push (@chain_length, $pdb_length);
		} elsif (substr( $pdb_name, 6, 1) eq "-") {
			push (@domainarray, $pdb_name);
			push (@domain_length, $pdb_length);
		} else {
			push (@complexarray, $pdb_name);
			push (@complex_length, $pdb_length);
		}		
	}
}

if (@chainarray) {

	$chainidlist = join( "\n", @chainarray );
	print OUTCHAINIDLIST $chainidlist;
	print OUTCHAINIDLIST "\n";
	close(OUTCHAINIDLIST);

	$chainlengthlist = join( "\n", @chain_length);
	print CHAINLENGTHLIST $chainlengthlist;
	print CHAINLENGTHLIST "\n";
	close(CHAINLENGTHLIST);
}

if (@domainarray) {

	$domainidlist = join( "\n", @domainarray );
	print OUTDOMAINIDLIST $domainidlist;
	print OUTDOMAINIDLIST "\n";
	close(OUTDOMAINIDLIST);

	$domainlengthlist = join( "\n", @domain_length);
	print DOMAINLENGTHLIST $domainlengthlist;
	print DOMAINLENGTHLIST "\n";
	close(DOMAINLENGTHLIST);
}

if (@complexarray) {

	$complexidlist = join( "\n", @complexarray );
	print OUTCOMPLEXIDLIST $complexidlist;
	print OUTCOMPLEXIDLIST "\n";
	close(OUTCOMPLEXIDLIST);

	$complexlengthlist = join( "\n", @complex_length);
	print COMPLEXLENGTHLIST $complexlengthlist;
	print COMPLEXLENGTHLIST "\n";
	close(COMPLEXLENGTHLIST);
}

print "Success\n";

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

	`cat $chain_list_file $added_idlist | sort -u > $idlist_tmp`;
	`mv $idlist_tmp $chain_list_file`;

	`cat $pdb2length_file $added_length_mapping | sort -u > $length_mapping_tmp`;
	`mv $length_mapping_tmp $pdb2length_file`;

	#get CATH mapping
	`./pdb2cath_$group_type.pl`;

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

foreach my $i ( 0 .. $#group_array ) {
	$group_type = $group_array[$i];

	$added_idlist = "./pdbidlistupdate/added.".$group_type."idlist.txt";
	$added_length_mapping = "./pdbidlistupdate/added.".$group_type."idlist.length.txt";

	`cat $chain_list_file $added_idlist | sort -u > $idlist_tmp`;
	`mv $idlist_tmp $chain_list_file`;

	`cat $pdb2length_file $added_length_mapping | sort -u > $length_mapping_tmp`;
	`mv $length_mapping_tmp $pdb2length_file`;
}

#get CATH mapping
`cat ./pdbidlistupdate/chainid.cath.mapping.txt ./pdbidlistupdate/complexid.cath.mapping.txt ./pdbidlistupdate/domainid.cath.mapping.txt > ./pdbidlistupdate/allid.cath.mapping.txt`;

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

###############################################
#Generate gifs/pngs for chain/domain/complex
###############################################
print "Generate gif and png ...\n";

`python image_multiprocess.py 2 > ./pdbidlistupdate/image.log`; #output structure id and pid
print "Success\n";

`rm *.gif`;
`rm *.png`;

#######################################
#Generate the finishing time
#######################################
print "Generate finish time ...\n";

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
`cp ./pdbidlistupdate/latest.updated.date ../3d-surfer/db/latest.updated.date`;

print "Success\n";
