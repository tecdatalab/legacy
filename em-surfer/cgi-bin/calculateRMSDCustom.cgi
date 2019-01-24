#!/usr/bin/perl

use CGI;

$q = new CGI;

#my $query = $q->param('q');
# 0 = query pdb
# 1-25 = 0 or comparison pdb
my @params = $q->param();
my $a1 = $q->param(@params[0]);
my $a2;
my $i;
#my $ret = "$a1";
my $ret = "rmsd";
my $ce = "";
my $pdb_chains_dir = "../pdb_chains/";
my $compared_filename;
my ($compared_pdb_id, $compared_chain_id);

print $q->header;


for ($i = 1; $i <= 25; $i++)
{
	$a2 = $q->param(@params[$i]);
#	$a1 = $a1 . $q->param(@params[$i]);
#$line = `cd /dragon3/www/3dzernikeDev/ce_distr; /dragon3/www/3dzernikeDev/ce_distr/CE - /bio/corgi5-scratch2/David/images/pymol/pdb_bench/1cd1-A.pdb - /bio/corgi5-scratch2/David/images/pymol/pdb_bench/1gzp-A.pdb - /dragon3/www/3dzernikeDev/ce_distr/scratch`;
#	if ($a2 != 0)
#	{
#		$a1 = $a1 . $a2;
#	}
	if ($a2 != 0)
	{
		($compared_pdb_id,$compared_chain_id) = split(/-/, $a2);
                $compared_filename = lc($compared_pdb_id)."-".$compared_chain_id.".pdb";
                if(!(-e "$pdb_chains_dir$compared_filename")) # then the pdb id may be uppercase
                {
                        $compared_filename = uc($compared_pdb_id)."-".$compared_chain_id.".pdb";
                }
		
		$line = `cd ../ce_distr; ./CE - ../cgi-bin/tmp/$a1.pdb - $pdb_chains_dir$compared_filename - ./scratch > ../cgi-bin/tmp/ce\_$$\_$i.txt; cat ../cgi-bin/tmp/ce\_$$\_$i.txt`;
		$ce = $ce . "`;" . $line;

		($rmsd) = $line =~ /Rmsd = (\d+\.\d+)A/;
		$rmsd = 'N/A' if !$rmsd;
		$ret = $ret . "`;" . $rmsd;

	}
	else
	{
		$ce = $ce . "`;0";
		$ret = $ret . "`;0";
	}
}

#print "RMSD: $rmsd\n";

#print "$a1";
print "$ret$ce`;$a1`;$$";
#print "$ret$ce";
#print "$ret";
