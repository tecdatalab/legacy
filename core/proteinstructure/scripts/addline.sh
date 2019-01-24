#!/bin/bash
# add line numbers to beginning of line

while [ "$*" != "" ]
do

	FILENAME=$1
	TMPFILE="$1.tmp"
	count=0

	cat $FILENAME | while read LINE
	do
		if [ $count -eq 0 ]
		then
			echo "$count	$LINE" > $TMPFILE
		else
			echo "$count	$LINE" >> $TMPFILE
		fi
		let count++
	done

	mv $TMPFILE $FILENAME

	shift
done

