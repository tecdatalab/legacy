#!/bin/sh
# get points input is the pdb file followed by visgrid residue number arranged
# on each line
./getpoints -p $1 -v $2

# do the hull stuff
for file in results_*.txt
do
./qconvex FS  FA < $file
done

rm results_*.txt

