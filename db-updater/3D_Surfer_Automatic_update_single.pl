#!/usr/bin/perl

# perl 3D_Surfer_Automatic_update_single.pl
use LWP::Simple;
chdir("/bio/kihara-web/www/db");

=pod
#######################################
# Protandpna_idlist
#######################################
$outfile = "./pdbidlistupdate/added.protandpna.pdb";
open(OUT,">$outfile");
print OUT "1htq\n";
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
$idlist = "./pdbidlistupdate/added.protandpna.pdb";
open( IN, $idlist);
@listarray = <IN>;
close(IN);
chomp(@listarray);

my %singlecomplexidhash; # these complexes with single chain
my %complexidhash; # these complexes include at least two chain
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

open( OUTCOMPLEXIDLIST, ">./pdbidlistupdate/added.complexidlist.txt" );
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

=cut

########################################################
# Calculate distance of any pair of chains in a complex
########################################################
$dir = "chain/";

my $infile = "./pdbidlistupdate/added.complexidlist.txt";
open( IN, $infile );
my @pdbidarray = <IN>;
close(IN);
chomp(@pdbidarray);

foreach my $m ( 0 .. $#pdbidarray ) {
	
	my $pdbid=$pdbidarray[$m];
	my %chainidhash = ();
	open( FILE, "$dir$pdbid\.pdb" );
	@pdbdataarray = <FILE>;
	close(FILE);

	foreach my $n ( 0 .. $#pdbdataarray ) {
		$chainidhash{ substr( $pdbdataarray[$n], 21, 1 ) } = 1;
	}

	my @chainidarray = keys %chainidhash;
	my $lenhash      = keys %chainidhash;
	if ( $lenhash >= 2 ) {
		
		$outfile = ">distance/$pdbid\.distance";
		open( OUT, $outfile );
	  endloop: foreach my $n ( 0 .. $#chainidarray ) {

			$firstchainid = $chainidarray[$n];
			$pinit        = $n + 1;
			foreach my $p ( $pinit .. $#chainidarray ) {
				$secondchainid = $chainidarray[$p];

				#firstresiduearray
				my %ridhash = ();
				$x = 0;
				foreach my $i ( 0 .. $#pdbdataarray ) {
					if (   substr( $pdbdataarray[$i], 0, 4 ) eq "ATOM"
						&& substr( $pdbdataarray[$i], 21, 1 ) eq $firstchainid )
					{
						$rid = substr( $pdbdataarray[$i], 22, 5 );
						$rid =~ s/\A\s+|\s+\z//g;
						$ridhash{$rid} = $x;
						++$x;
					}
				}

				my @firstresiduearray = ();
				@firstresiduearray =
				  sort { $ridhash{$a} <=> $ridhash{$b} } keys %ridhash;

				#secondresiduearray
				my %ridhash = ();
				$x = 0;
				foreach my $i ( 0 .. $#pdbdataarray ) {
					if (   substr( $pdbdataarray[$i], 0, 4 ) eq "ATOM"
						&& substr( $pdbdataarray[$i], 21, 1 ) eq
						$secondchainid )
					{
						$rid = substr( $pdbdataarray[$i], 22, 5 );
						$rid =~ s/\A\s+|\s+\z//g;
						$ridhash{$rid} = $x;
						++$x;
					}
				}

				my @secondresiduearray = ();
				@secondresiduearray =
				  sort { $ridhash{$a} <=> $ridhash{$b} } keys %ridhash;

				##############################################

				my @firstchainatomarray  = ();
				my @secondchainatomarray = ();

				foreach my $i ( 0 .. $#pdbdataarray ) {
					if (   substr( $pdbdataarray[$i], 0, 4 ) eq "ATOM"
						&& substr( $pdbdataarray[$i], 21, 1 ) eq $firstchainid )
					{
						push @firstchainatomarray, $pdbdataarray[$i];
					}
				}

				foreach my $i ( 0 .. $#pdbdataarray ) {
					if (   substr( $pdbdataarray[$i], 0, 4 ) eq "ATOM"
						&& substr( $pdbdataarray[$i], 21, 1 ) eq
						$secondchainid )
					{
						push @secondchainatomarray, $pdbdataarray[$i];
					}
				}

				##################################################
				#first hash
				##################################################
				my %firstcontactresiduehash;
				foreach my $i ( 0 .. $#firstresiduearray ) {

					my @residueatomarray = ();
					for my $j ( 0 .. $#firstchainatomarray ) {

						$rid = substr( $firstchainatomarray[$j], 22, 5 );
						$rid =~ s/\A\s+|\s+\z//g;
						if ( $rid eq $firstresiduearray[$i] ) {
							push @residueatomarray, $firstchainatomarray[$j];
						}
					}

					# the number of atoms included in the considering residue
					for my $j ( 0 .. $#residueatomarray ) {

						$xp = substr( $residueatomarray[$j], 30, 8 );
						$yp = substr( $residueatomarray[$j], 38, 8 );
						$zp = substr( $residueatomarray[$j], 46, 8 );

						for my $d ( 0 .. $#secondchainatomarray ) {

							$xd = substr( $secondchainatomarray[$d], 30, 8 );
							$yd = substr( $secondchainatomarray[$d], 38, 8 );
							$zd = substr( $secondchainatomarray[$d], 46, 8 );

							$distance =
							  ( $xp - $xd )**2 +
							  ( $yp - $yd )**2 +
							  ( $zp - $zd )**2;

							if ( $distance < 4.5 * 4.5 ) {

								$firstcontactresiduehash{ $firstresiduearray[$i]
								  } = 1;
							}
						}    #end for d
					}    #end for j every residue
				}    #end for i in lenfirstrwa

				$firstcontactnum = keys %firstcontactresiduehash;

				##################################################
				#second hash
				##################################################
				my %secondcontactresiduehash;
				foreach my $i ( 0 .. $#secondresiduearray ) {

					my @residueatomarray = ();

					for my $j ( 0 .. $#secondchainatomarray ) {

						$rid = substr( $secondchainatomarray[$j], 22, 5 );
						$rid =~ s/\A\s+|\s+\z//g;
						if ( $rid eq $secondresiduearray[$i] ) {
							push @residueatomarray, $secondchainatomarray[$j];
						}
					}

					# the number of atoms included in the considering residue
					for my $j ( 0 .. $#residueatomarray ) {

						$xp = substr( $residueatomarray[$j], 30, 8 );
						$yp = substr( $residueatomarray[$j], 38, 8 );
						$zp = substr( $residueatomarray[$j], 46, 8 );

						for my $d ( 0 .. $#firstchainatomarray ) {

							$xd = substr( $firstchainatomarray[$d], 30, 8 );
							$yd = substr( $firstchainatomarray[$d], 38, 8 );
							$zd = substr( $firstchainatomarray[$d], 46, 8 );

							$distance =
							  ( $xp - $xd )**2 +
							  ( $yp - $yd )**2 +
							  ( $zp - $zd )**2;

							if ( $distance < 4.5 * 4.5 ) {

								$secondcontactresiduehash{ $secondresiduearray
									  [$i] } = 1;
							}
						}    #end for d
					}    #end for j every residue
				}    #end for i in lenfirstrwa

				$secondcontactnum = keys %secondcontactresiduehash;
				print OUT $firstchainid . "\t"
				  . $firstcontactnum . "\t"
				  . $secondchainid . "\t"
				  . $secondcontactnum . "\n";
			}
		}
		close(OUT);
	}
}

########################################################
# Cluster chains and generate our own defined complexes
########################################################

my $infile = "./pdbidlistupdate/added.complexidlist.txt";
open( IN, $infile );
my @pdbidarray = <IN>;
close(IN);
chomp(@pdbidarray);

open( OUT, ">./pdbidlistupdate/added.complexidlist.cluster.txt" );
foreach my $pid ( 0 .. $#pdbidarray ) {

	my $pdbid = $pdbidarray[$pid];
	my @idpairarray;

	# adjust the threshold of numbers of contacting residues
	open( IN, "distance/$pdbid.distance" );
	my @linearray = <IN>;
	close(IN);

	my %uniprotidpairhash = ();

	foreach my $i ( 0 .. $#linearray ) {

		chomp( $linearray[$i] );
		my @tmparray = split( /\s+/, $linearray[$i] );

		if ( $tmparray[1] > 4 && $tmparray[3] > 4 ) {

			if ( $tmparray[0] lt $tmparray[2] ) {

				push @idpairarray, $tmparray[0] . $tmparray[2];
			}
			else {

				push @idpairarray, $tmparray[2] . $tmparray[0];
			}
		}
	}

	@idpairarray = sort { $a cmp $b } @idpairarray;

	# to deal with the pairs and generated clusters
	my $chainidliststr;
	foreach my $i ( 0 .. $#idpairarray ) {

		if ( index( $chainidliststr, substr( $idpairarray[$i], 0, 1 ) ) == -1 )
		{
			$chainidliststr .= substr( $idpairarray[$i], 0, 1 );
		}
		if ( index( $chainidliststr, substr( $idpairarray[$i], 1, 1 ) ) == -1 )
		{
			$chainidliststr .= substr( $idpairarray[$i], 1, 1 );
		}
	}

	#print OUT "$pdbid\t";

	#print OUT $chainidliststr. "\n";

	# the first time, cluster array separated by blank in the chainidlist string
	my @oldclusterarray = split( //, $chainidliststr );

  # the rest if time, cluster array separated by comma in the chainidlist string
	my $oldclusterstr = join( ',', @oldclusterarray );
	my $newclusterstr = $oldclusterstr;

  endfor: foreach my $i ( 0 .. $#idpairarray ) {

		my $connect = $idpairarray[$i];

		$oldclusterstr = $newclusterstr;
		@oldclusterarray = split( /,/, $oldclusterstr );
		my $firstelement;
		my $secondelement;

		# $j; the $j th cluster
	  endfor1: foreach my $j ( 0 .. $#oldclusterarray ) {
			if ( index( $oldclusterarray[$j], substr( $connect, 0, 1 ) ) != -1 )
			{
				$firstelement = $j;
				last endfor1;
			}
		}

	  endfor2: foreach my $j ( 0 .. $#oldclusterarray ) {
			if ( index( $oldclusterarray[$j], substr( $connect, 1, 1 ) ) != -1 )
			{
				$secondelement = $j;
				last endfor2;
			}
		}

	   # It means the two elements of a given connect are in the orignal cluster
		if ( $firstelement == $secondelement ) {

			$newclusterstr = $oldclusterstr;

		}
		else {

			$newclusterstr = "";
			$newclusterstr .= $oldclusterarray[$firstelement];

			my @tmparr = split( //, $oldclusterarray[$secondelement] );

			foreach my $j ( 0 .. $#tmparr ) {

				if ( index( $newclusterstr, $tmparr[$j] ) == -1 ) {
					$newclusterstr .= $tmparr[$j];
				}
			}
			foreach my $j ( 0 .. $#oldclusterarray ) {
				if ( ( $j != $firstelement ) && ( $j != $secondelement ) ) {
					$newclusterstr .= "," . $oldclusterarray[$j];
				}
			}
		}

		#print OUT $newclusterstr. "\n";

		my @idarr1 = split( //, $chainidliststr );
		my @idarr2 = split( //, $newclusterstr );
		my @idarr1 = sort { $a cmp $b } @idarr1;
		my @idarr2 = sort { $a cmp $b } @idarr2;
		my $idstr1 = join( '', @idarr1 );
		my $idstr2 = join( '', @idarr2 );

		if ( $idstr1 eq $idstr2 ) {

			#print "$pdbid\t$i\t$#idpairarray\tEqual\n";
			#print OUT $pdbid . "\n";
			last endfor;
		}
	}

	#print OUT $newclusterstr. "\n";
	my @tmparr = split( /,/, $newclusterstr );
	foreach my $i ( 0 .. $#tmparr ) {

		my @idarr = split( //, $tmparr[$i] );
		my @idarr = sort { $a cmp $b } @idarr;
		$tmparr[$i] = join( '', @idarr );
	}
	@tmparr = sort { $a cmp $b } @tmparr;
	$newclusterstr = join( ',', @tmparr );

	#print OUT $newclusterstr. "\n";
	#print OUT "$pdbid\t";
	foreach my $i ( 0 .. $#tmparr ) {

		my $tmpi = $i + 1;
		if ( $tmpi < 10 ) {
			if ($#tmparr == 0) {
				print OUT $pdbid . "\t" . $tmparr[$i] . "\n";
			} else {
				print OUT $pdbid . "-C0" . $tmpi . "\t" . $tmparr[$i] . "\n";
			}
		}
		else {
			print OUT $pdbid . "-C" . $tmpi . "\t" . $tmparr[$i] . "\n";
		}
	}
}
close(OUT);

=pod

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

=cut
