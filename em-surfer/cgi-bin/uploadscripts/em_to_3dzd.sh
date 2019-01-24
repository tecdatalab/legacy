#!/bin/sh
if [ $# -ne 3 ]; then
    echo "Usage: $0 <em map file> <contour value> <output file>"
    exit 0
fi

MAPFILE=$1
CONTOUR=$2
OUTPUT=$3
BASE=`basename $MAPFILE`
RANDPREFIX=`mktemp`
PREFIX="$RANDPREFIX${BASE%.*}"

#TMP FILES
SITUS=$PREFIX.situs
GRID=$PREFIX.grid
DEFAULTINV=$PREFIX.inv

echo "====EM Map to 3DZD: $1 $2 $3 $SITUS $GRID $DEFAULTINV"
echo "Running map2map..."
echo 2 | map2map $MAPFILE $SITUS

echo "Running situs2grid.pl..."
./situs2grid.pl $SITUS $GRID $CONTOUR +inf

echo "Creating 3DZD inv file..."
./zernike $GRID 20

echo "Generating output file..."
cp $DEFAULTINV $OUTPUT

echo "Removing tmp files..."
rm $SITUS $GRID $DEFAULTINV
