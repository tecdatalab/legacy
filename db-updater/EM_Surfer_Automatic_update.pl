#!/usr/bin/perl

# perl EM_Surfer_Automatic_update.pl
use LWP::Simple;
chdir("/bio/kihara-web/www/db/emdb");

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

`./main_emdb_update.sh /bio/kihara-web/www/em-surfer/db 2>&1 1>main.log`;

$date_file    = "./em-surfer.latest.updated.date";
open(OUT, ">$date_file");   
print OUT "$month $day, $year\n";
close(OUT);
`cp ./em-surfer.latest.updated.date /bio/kihara-web/www/em-surfer/db/latest.updated.date`;
