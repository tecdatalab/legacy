#!/bin/sh

if [ $?LD_LIBRARY_PATH ]; then
    export LD_LIBRARY_PATH=/bio/kihara-web/www/3d-surfer/vishhull/qhull-2003.1/src/.libs:$LD_LIBRARY_PATH
else
    export LD_LIBRARY_PATH=/bio/kihara-web/www/3d-surfer/vishhull/qhull-2003.1/src/.libs:$LD_LIBRARY_PATH
fi

# get points input is the pdb file followed by visgrid residue number arranged
# on each line
/bio/kihara-web/www/3d-surfer/convexhull/getpoints -p $1 -v $2

# do the hull stuff
for file in results_*.txt
do
/bio/kihara-web/www/3d-surfer/vishhull/qhull-2003.1/src/qconvex FS FA < $file
done

rm results_*.txt

