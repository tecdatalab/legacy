#! /usr/bin/perl

# Script: perl clustercomplex.pl
# Description:
# Given a set of N items to be clustered, and an N*N distance (or similarity) matrix, the basic process of hierarchical clustering (defined by S.C. Johnson in 1967) is this:
# Step 1: Start by assigning each item to a cluster.
# Step 2: Find the closest (most similar) pair of clusters and merge them into a single cluster, so that now one cluster less each time.
# Step 3: Compute distances (similarities) between the new cluster and each of the old clusters.
# Step 4: Repeat steps 2 and 3 until all items are clustered into a single cluster of size N, or the loop  of items ends.
# Author: Yi Xiong
# Email: xiongyi128@gmail.com
# Date: Apr.19.2013

my $infile = "./pdbidlistupdate/complexpdblist.txt";
open( IN, $infile );
my @pdbidarray = <IN>;
close(IN);
chomp(@pdbidarray);

open( OUT, ">./pdbidlistupdate/complexcomponentidchainlist.txt" );
foreach my $pid ( 0 .. $#pdbidarray ) {

	my $pdbid = $pdbidarray[$pid];
	my @idpairarray;

	# adjust the threshold of numbers of contacting residues
	open( IN, "./distance/$pdbid.distance" );
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

# 1ado
