#!/bin/sh
if [ $# -lt 7 ]; then
    echo "Usage: $0 <emdb XMLs directory> <idlist output> <names output> <fullnames output> <contours output> <resolutions output> <resolution summary output>"
    echo "Example: $0 xml/ idlist.txt name.txt fullname.txt contour_levels.txt resolutions.txt resolution_summary.txt"
    exit 0
fi

ls $1/*.xml | awk '{print "basename "$1}' | sh | sed -e 's/emd-//' | sed -e 's/\.xml//' > $2
grep -e '<name>' $1/*.xml | sed -e 's/^.*emd-//' | sed -e 's/\.xml\://' | sed -e 's/<\/*name>//g' | awk '{id=$1;$1="";print id"\t"substr($0,2,40)}' > $3
grep -e '<name>' $1/*.xml | sed -e 's/^.*emd-//' | sed -e 's/\.xml\://' | sed -e 's/<\/*name>//g' | awk '{id=$1;$1="";print id"\t"substr($0,2)}' > $4
grep contourLevel $1/*.xml | sed -e 's/^.*emd-//' | sed -e 's/\.xml\://' | sed -e 's/<\/contourLevel>//g' | sed -e 's/<contourLevel.*>//g' > $5
grep resolutionByAuthor $1/*.xml | sed -e 's/^.*emd-//' | sed -e 's/\.xml\://' | sed -e 's/<\/*resolutionByAuthor>//g' | awk '{print $1"\t"$2}' > $6
./extract_xml_resolution_summary.pl $6 > $7
