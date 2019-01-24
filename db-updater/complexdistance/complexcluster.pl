#! /usr/bin/perl

# Script: perl complexcluster.pl 0 0
# Description: atom distance, output the contacting residues
# Author: Yi Xiong
# Email: xiongyi128@gmail.com
# Date: Apr.8.2013

$dir = "/bio/kihara-web/www/db/chain/";

# ids in idlist.txt
my @idarray = ();
open( IN, "falseidlist9" );
@idarray = <IN>;
close(IN);
chomp(@idarray);

open( OUT, ">falseidlist9.cluster" );
# $ARGV[0] .. $ARGV[1]     0 .. $#idarray
foreach my $m ( $ARGV[0] .. $ARGV[1] ) {

    print $m."\t".$idarray[$m]."\n";
	my %chainidhash     = ();
	my %chainidpairhash = ();
	open( FILE, "$dir$idarray[$m]\.pdb" );
	@pdbdataarray = <FILE>;
	close(FILE);

	foreach my $n ( 0 .. $#pdbdataarray ) {
		$chainidhash{ substr( $pdbdataarray[$n], 21, 1 ) } = 1;
	}

	my @chainidarray = sort { $a cmp $b } keys %chainidhash;
	my $lenhash = keys %chainidhash;

	#if ( $lenhash < 2 ) {
	#	print substr( $idarray[$m], 0, 4 ) . "\n";
	#}
	#else {
	if ( $lenhash >= 2 ) {

		my $chainlist;
		foreach my $n ( 0 .. $#chainidarray ) {
			$chainlist .= $chainidarray[$n];
		}
		#print $chainlist. "\n";

		my @idpairarray;

		# two loops
		foreach my $n ( 0 .. $#chainidarray ) {

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
				if ( ( $firstcontactnum >= 5 ) && ( $secondcontactnum >= 5 ) ) {

					print $firstchainid. "\t" . $secondchainid . "\n";
					my $tmpidpair = $firstchainid . $secondchainid;
					push @idpairarray, $tmpidpair;
				}
			}
		}    # end for two loops

	  endout: foreach my $i ( 0 .. $#idpairarray ) {

			$idcluster = $idpairarray[$i];
			foreach my $j ( 0 .. $#idpairarray ) {

				if ( $j != $i ) {

					$clusterelement = $idpairarray[$j];
					if (
						index( $idcluster, substr( $clusterelement, 0, 1 ) ) !=
						"-1" )
					{
						$idcluster .= substr( $clusterelement, 1, 1 );
					}
					elsif (
						index( $idcluster, substr( $clusterelement, 1, 1 ) ) !=
						"-1" )
					{
						$idcluster .= substr( $clusterelement, 0, 1 );
					}

					my @idarr1 = split( //, $idcluster );
					my @idarr2 = split( //, $chainlist );
					my @idarr1 = sort { $a cmp $b } @idarr1;
					my @idarr2 = sort { $a cmp $b } @idarr2;
					my $idstr1 = join( '', @idarr1 );
					my $idstr2 = join( '', @idarr2 );

					if ( $idstr1 eq $idstr2 ) {

						print OUT "$idarray[$m]\n";
						#print $idarray[$m];
						last endout;
					}
				}
			}
		}
	} # end for $lenhash >= 2
	#print "\n";
}
close(OUT);