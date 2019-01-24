#!/usr/bin/perl

# perl calculate_from_uploaded_file.pl

$dir = "tmp/";
$mapfile = $dir."emd_1884.map";
$situsfile = $dir."emd_1884.situs";
$recommendcontourlevel=3.16;
open(IN, $situsfile);
@linearray=<IN>;
close(IN);

$maximum = 0;
for my $i (2 .. $#linearray) {

        $line = $linearray[$i];
	$line =~ s/\A\s+|\s+\z//g;
        my @tmparray = split (/\s+/, $line);
        foreach my $element (@tmparray) {
		if ($element > $maximum) {
			$maximum = $element;
                }
        }
}

#$onethirdcontourlevel=$recommendcontourlevel+($maximum-$recommendcontourlevel)/3;
#$dist = sprintf "%.3f", $dist;
#$twothirdcontourlevel=$recommendcontourlevel+(2*($maximum-$recommendcontourlevel)/3);
#print $onethirdcontourlevel."\n";


$twothirdcontourlevel=$recommendcontourlevel+(2*($maximum-$recommendcontourlevel)/3);
print $twothirdcontourlevel."\n";

$result = `./em_volume 3.16 $mapfile`;
@tmparray=split(/\n/, $result);
@tmplinearray=split(/\s+/, $tmparray[0]);
$volume=$tmplinearray[1];
print $volume."\n";
