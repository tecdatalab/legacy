#!/usr/bin/perl

use strict;
use DBI();
use DBI qw(:sql_types);

# Connect to the database.
my $dbh = DBI->connect("DBI:mysql:database=KiharaDB;host=coco1",
			"autoupdateuser", "autoupdateuser",{'RaiseError' => 1});

my $sth = $dbh->prepare("call `KiharaDB`.`SP_NewUpdate`('New Update Perl')");

#$rv = $sth->execute();

my @result = $dbh->selectrow_array($sth);

print "The result ";
print @result;
print "\n";


$sth->finish();
$dbh->disconnect();
