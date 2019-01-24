#! /usr/bin/perl

# Script: perl idlist.pl 0 5000 > falseidlist1
# Description: atom distance, output the contacting residues
# Author: Yi Xiong
# Email: xiongyi128@gmail.com
# Date: Apr.8.2013

$dir = "/bio/kihara-web/www/db/chain/";

=pod
open(OUT, ">idlist.txt");
opendir(DIR, $dir);
@dir=readdir DIR;
foreach $file (@dir){
	if (length ($file) == 8 ) {
		print OUT $file."\n";
	}	
}
closedir(DIR);
close(OUT);
=cut

# 50389 ids in idlist.txt
my @idarray = ();
open( IN, "idlist.txt" );
@idarray = <IN>;
close(IN);
chomp(@idarray);

#print $#idarray."\n";

foreach my $m ( $ARGV[0] .. $ARGV[1] ) {

	my %chainidhash = ();
	open( FILE, "$dir$idarray[$m]" );
	@pdbdataarray = <FILE>;
	close(FILE);

	foreach my $n ( 0 .. $#pdbdataarray ) {
		$chainidhash{ substr( $pdbdataarray[$n], 21, 1 ) } = 1;
	}

	my @chainidarray = keys %chainidhash;
	my $lenhash      = keys %chainidhash;
	if ( $lenhash < 2 ) {
		print substr( $idarray[$m], 0, 4 ) . "\n";
	}
	else {

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
				if ( $firstcontactnum < 5 ) {
					print substr( $idarray[$m], 0, 4 ) . "\n";
					last endloop;
				}

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
				if ( $secondcontactnum < 5 ) {
					print substr( $idarray[$m], 0, 4 ) . "\n";
					last endloop;
				}
			}
		}
	}
}
