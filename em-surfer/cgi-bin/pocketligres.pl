#!/usr/bin/perl

#./pocketligres.pl 7tim-A ../chain/7tim-A.pdb

$query_pdb_id   = $ARGV[0];
$query_pdb_file = $ARGV[1];

#$query_pdb_id = "2wiw-B";
$pdbFile = "tmppdblcs/$query_pdb_id.pdb";

`rm $pdbFile`;
`cp $query_pdb_file $pdbFile`;

my @loadedarray = ();
open( TMPIN, $pdbFile );
@loadedarray = <TMPIN>;
close(TMPIN);

$ca_atom_counter = 0;

# TODO: check if the loaded pdb file is valid
for my $i ( 0 .. $#loadedarray ) {

	# check format and output the ATOM and HETATM component
	if ( $loadedarray[$i] =~ /^ATOM.........CA/ ) {
		$ca_atom_counter++;
	}
}

if ( $ca_atom_counter > 2000 )    # this is our current limit  2000
{
	`rm $pdbFile`;

	print "The PDB file is too big for the server to handle $ca_atom_counter!";
	exit;
}

if ( !$ca_atom_counter ) {
	`rm $pdbFile`;

	print
	  "The file that you uploaded does not appear to have atoms in the PDB.";
	exit;
}

# Generate the PDB file with ATOM and simulated HETATM
open( TMP, ">$pdbFile" );
for my $i ( 0 .. $#loadedarray ) {

	# check format and output the ATOM component
	if ( $loadedarray[$i] =~ /^ATOM/
		&& substr( $loadedarray[$i], 17, 3 ) =~ /[A-Z]{3,3}/ )
	{
		print TMP "$loadedarray[$i]";
	}
}
close(TMP);

`cp $pdbFile .`;

`./ligsitecs -i $query_pdb_id.pdb`;
$pocketsitefile = $query_pdb_id . "_pocket.pdb";

open( IN, $pocketsitefile );
@pocketsitearray = <IN>;
close(IN);

$pocketclusterfile = $query_pdb_id . "_cluster_pocket.pdb";
open( IN, $pocketclusterfile );
@pocketclusterarray = <IN>;
close(IN);

open( OUT, ">>$pdbFile" );
$sid = 80000;
for my $i ( 0 .. $#pocketsitearray ) {
	my @tmparray = split( /\s+/, $pocketsitearray[$i] );
	$pocketsitearray[$i] = $tmparray[5];

	$mark = 0;
	for my $j ( 0 .. $#pocketclusterarray ) {
		if ( $pocketclusterarray[$j] eq "BEGIN$tmparray[5]\n" ) {
			$mark = 1;
			next;
		}
		if ( $pocketclusterarray[$j] eq "END$tmparray[5]\n" ) {
			$mark = 0;
			next;
		}
		if ( $mark == 1 ) {

			#print OUT $pocketclusterarray[$j];
			@xyzarray = split( /\s+/, $pocketclusterarray[$j] );
			++$sid;
			if ( $i == 0 ) {
				printf OUT
"HETATM%5d  NA  HEM Z9997    %8.3f%8.3f%8.3f  1.00 19.51           N\n",
				  $sid, $xyzarray[0], $xyzarray[1], $xyzarray[2];
			}
			elsif ( $i == 1 ) {
				printf OUT
"HETATM%5d  NA  HEM Z9998    %8.3f%8.3f%8.3f  1.00 19.51           N\n",
				  $sid, $xyzarray[0], $xyzarray[1], $xyzarray[2];
			}
			elsif ( $i == 2 ) {
				printf OUT
"HETATM%5d  NA  HEM Z9999    %8.3f%8.3f%8.3f  1.00 19.51           N\n",
				  $sid, $xyzarray[0], $xyzarray[1], $xyzarray[2];
			}
		}
	}
}
close(OUT);

###
open( IN, $pdbFile );
@pdbarray = <IN>;
close(IN);

@atomarray = ();
@hetarray0 = ();
@hetarray1 = ();
@hetarray2 = ();

$pocketIDs[0] = $query_pdb_id . "-0";
$pocketIDs[1] = $query_pdb_id . "-1";
$pocketIDs[2] = $query_pdb_id . "-2";

#$pocket_res{$pocketIDs[0]} = "$reslist";

foreach my $element (@pdbarray) {
	if ( $element =~ /^ATOM/ ) {
		push @atomarray, $element;
	}
	elsif ( $element =~ /^HETATM...............Z9997/ ) {
		push @hetarray0, $element;
	}
	elsif ( $element =~ /^HETATM...............Z9998/ ) {
		push @hetarray1, $element;
	}
	elsif ( $element =~ /^HETATM...............Z9999/ ) {
		push @hetarray2, $element;
	}
}

#pocket 0
%pocketreshash = ();
foreach my $atomelement (@atomarray) {

	$x1 = substr( $atomelement, 30, 8 );
	$y1 = substr( $atomelement, 38, 8 );
	$z1 = substr( $atomelement, 46, 8 );

  lasthet: foreach my $hetelement (@hetarray0) {
		$x2 = substr( $hetelement, 30, 8 );
		$y2 = substr( $hetelement, 38, 8 );
		$z2 = substr( $hetelement, 46, 8 );
		$dist = ( $x1 - $x2 )**2 + ( $y1 - $y2 )**2 + ( $z1 - $z2 )**2;
		$dist = $dist**( 1 / 2 );
		if ( $dist < 4.5 ) {
			$resSeq = int( substr( $atomelement, 22, 4 ) );
			$chainID             = substr( $atomelement, 21, 1 );
			$tmp                 = $resSeq . ":" . $chainID;
			$pocketreshash{$tmp} = 1;
			last lasthet;
		}
	}
}

$reslist = "";
$len     = scalar( keys %pocketreshash );
$last    = 0;
foreach $key ( sort keys %pocketreshash ) {
	$reslist .= $key;
	++$last;
	$codeDir = "../code";
	if ( $last < $len ) {
		#$reslist .= ",";
		$reslist .= " ";
	}
}
$pocket_res{ $pocketIDs[0] } = $reslist;

#pocket 1
%pocketreshash = ();
foreach my $atomelement (@atomarray) {

	$x1 = substr( $atomelement, 30, 8 );
	$y1 = substr( $atomelement, 38, 8 );
	$z1 = substr( $atomelement, 46, 8 );

  lasthet: foreach my $hetelement (@hetarray1) {
		$x2 = substr( $hetelement, 30, 8 );
		$y2 = substr( $hetelement, 38, 8 );
		$z2 = substr( $hetelement, 46, 8 );
		$dist = ( $x1 - $x2 )**2 + ( $y1 - $y2 )**2 + ( $z1 - $z2 )**2;
		$dist = $dist**( 1 / 2 );
		if ( $dist < 4.5 ) {
			$resSeq = int( substr( $atomelement, 22, 4 ) );
			$chainID             = substr( $atomelement, 21, 1 );
			$tmp                 = $resSeq . ":" . $chainID;
			$pocketreshash{$tmp} = 1;
			last lasthet;
		}
	}
}

$reslist = "";
$len     = scalar( keys %pocketreshash );
$last    = 0;
foreach $key ( sort keys %pocketreshash ) {
	$reslist .= $key;
	++$last;
	if ( $last < $len ) {
		#$reslist .= ",";
		$reslist .= " ";
	}
}
$pocket_res{ $pocketIDs[1] } = $reslist;

#pocket 2
%pocketreshash = ();
foreach my $atomelement (@atomarray) {

	$x1 = substr( $atomelement, 30, 8 );
	$y1 = substr( $atomelement, 38, 8 );
	$z1 = substr( $atomelement, 46, 8 );

  lasthet: foreach my $hetelement (@hetarray2) {
		$x2 = substr( $hetelement, 30, 8 );
		$y2 = substr( $hetelement, 38, 8 );
		$z2 = substr( $hetelement, 46, 8 );
		$dist = ( $x1 - $x2 )**2 + ( $y1 - $y2 )**2 + ( $z1 - $z2 )**2;
		$dist = $dist**( 1 / 2 );
		if ( $dist < 4.5 ) {
			$resSeq = int( substr( $atomelement, 22, 4 ) );
			$chainID             = substr( $atomelement, 21, 1 );
			$tmp                 = $resSeq . ":" . $chainID;
			$pocketreshash{$tmp} = 1;
			last lasthet;
		}
	}
}

$reslist = "";
$len     = scalar( keys %pocketreshash );
$last    = 0;
foreach $key ( sort keys %pocketreshash ) {
	$reslist .= $key;
	++$last;
	if ( $last < $len ) {
		#$reslist .= ",";
		$reslist .= " ";
	}
}
$pocket_res{ $pocketIDs[2] } = $reslist;

#output list
for my $i ( 0 .. $#pocketIDs ) {
	print $pocket_res{ $pocketIDs[$i] } . "\n";

=pod
	my @pocketResArray = split( /,/, $pocket_res{ $pocketIDs[$i] } );
	foreach my $j (0..$#pocketResArray) {
		print "$pocketResArray[$j]\t";
	}
	print "\n";
=cut

}

`rm $query_pdb_id*`;
