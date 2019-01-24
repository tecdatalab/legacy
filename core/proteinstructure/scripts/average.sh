awk 'BEGIN {isfirst=1;}
{
	if($1 == $1+0) #i.e. if it is a number
	{
		sum+=$1;++n;
		if(isfirst == 1) {max = min = $1; isfirst=0}
		else {
			if($1 < min) {min = $1};
			if($1 > max) {max = $1};
		}
	}
}
END {print "Average "sum/n" Max "max" Min "min}'
