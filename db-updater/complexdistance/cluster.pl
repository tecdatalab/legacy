#! /usr/bin/perl

# Script: perl cluster.pl
# Description: atom distance, output the contacting residues
# Author: Yi Xiong
# Email: xiongyi128@gmail.com
# Date: Apr.11.2013

my @idpairarray;

#push @idpairarray, "AD";
#push @idpairarray, "BC";
#push @idpairarray, "BD";
#@idpairarray=("AD", "BC", "BD");
#@idpairarray = ("AB", "AC", "AD", "AE", "DE", "DN", "EF", "FG", "FH", "FI", "FJ", "IJ", "IS", "KL", "KM", "KN", "KO", "NO", "OP", "PQ", "PR", "PS", "PT", "ST");
@idpairarray = ("AB", "BC");

my $chainidliststr;
foreach my $i ( 0 .. $#idpairarray ) {

	if ( index( $chainidliststr, substr( $idpairarray[$i], 0, 1 ) ) == -1 ) {
		$chainidliststr .= substr( $idpairarray[$i], 0, 1 );
	}
	if ( index( $chainidliststr, substr( $idpairarray[$i], 1, 1 ) ) == -1 ) {
		$chainidliststr .= substr( $idpairarray[$i], 1, 1 );
	}
}

print $chainidliststr. "\n";

my @oldclusterarray = split( //, $chainidliststr );
my $oldclusterstr = join( ',', @oldclusterarray );
my $newclusterstr = $oldclusterstr;

endfor:foreach my $i ( 0 .. $#idpairarray ) {

	my $connect = $idpairarray[$i];
	$oldclusterstr = $newclusterstr;
	@oldclusterarray = split( /,/, $oldclusterstr );
	my $firstelement;
	my $secondelement;
  endfor1: foreach my $j ( 0 .. $#oldclusterarray ) {
		if ( index( $oldclusterarray[$j], substr( $connect, 0, 1 ) ) != -1 ) {
			$firstelement = $j;
			last endfor1;
		}
	}

  endfor2: foreach my $j ( 0 .. $#oldclusterarray ) {
		if ( index( $oldclusterarray[$j], substr( $connect, 1, 1 ) ) != -1 ) {
			$secondelement = $j;
			last endfor2;
		}
	}

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
	print $newclusterstr. "\n";

	my @idarr1 = split( //, $chainidliststr );
	my @idarr2 = split( //, $newclusterstr );
	my @idarr1 = sort { $a cmp $b } @idarr1;
	my @idarr2 = sort { $a cmp $b } @idarr2;
	my $idstr1 = join( '', @idarr1 );
	my $idstr2 = join( '', @idarr2 );

	if ( $idstr1 eq $idstr2 ) {

		print "$i\t$#idpairarray\tEqual\n";
		last endfor;
	}
}

=pod
endout: foreach my $i ( 0 .. $#idpairarray ) {

	$idcluster = $idpairarray[$i];
	foreach my $j ( 0 .. $#idpairarray ) {

		if ( $j != $i ) {
			
			$clusterelement = $idpairarray[$j];
			if ( index( $idcluster, substr( $clusterelement, 0, 1 ) ) != "-1" )
			{
				$idcluster .= substr( $clusterelement, 1, 1 );
			}
			elsif (
				index( $idcluster, substr( $clusterelement, 1, 1 ) ) != "-1" )
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

				print "equal the complete chain list $chainlist\n";
				last endout;
			}
		}
	}
}
=cut
