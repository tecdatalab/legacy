#!/usr/bin/perl

# perl 3D_Surfer_Automatic_update_split.pl
use LWP::Simple;
chdir("/bio/kihara-web/www/db");

#######################################
#Parse the pdb to separate chains
#######################################
$idlist = "./added.protandpna.pdb";
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
		
		if(substr( $_, 0, 6 ) eq "ENDMDL"){
			last lastwhile;
		}

		# Skip if non atom
		next if !/^ATOM/;
		
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
$idlist = "./added.protandpna.pdb";
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
$idlist = "./added.protandpna.pdb";
open( IN, $idlist);
@listarray = <IN>;
close(IN);
chomp(@listarray);

my %singlecomplexidhash; # these complexes with single chain
my %complexidhash; # these complexes include at least two chain
open( OUTCHAINIDLIST, ">./added.chainidlist.txt" );
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
		
		$singlecomplexidhash{$pdb_id}=1;
		foreach my $key (sort keys %chainidhash) {
			print OUTCHAINIDLIST $pdb_id."-".$key."\n";
		}	
	} elsif ($size > 1) {
		
		$complexidhash{$pdb_id}=1;
		foreach my $key (sort keys %chainidhash) {
			print OUTCHAINIDLIST $pdb_id."-".$key."\n";
		}
	}
}
close(OUTCHAINIDLIST);

open( OUTCOMPLEXIDLIST, ">./added.complexidlist.txt" );
foreach my $key (sort keys %complexidhash) {
	
	print OUTCOMPLEXIDLIST $key."\n";
}
close(OUTCOMPLEXIDLIST);

##########################################################
#Delete complexes of single chains in chain or chain_cacn 
##########################################################
foreach my $key (sort keys %singlecomplexidhash) {
	
	$file = "chain/$key.pdb";
	$file_cacn = "chain_cacn/$key.pdb";
	`rm $file`;
	`rm $file_cacn`;
}

print "Success\n";
