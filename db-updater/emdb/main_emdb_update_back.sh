#!/bin/bash
if [ $# -lt 1 ]; then
  echo "Usage: $0 <production db directory>"
  exit 0
fi

PRODUCTIONDIR=$1

# PREREQUISITES: within extract_xml_information.sh we call extract_xml_resolution_summary.pl
# calculate_volume_and_3dzd.pl uses em_volume and em2zernike programs.
# Make sure that either the relative or absolute paths to these are correct.

# Make an FTP connection to EBI's EMDB and get the latest update timestamps.
echo "Step 1: Get latest timestamps from EBI..."
date
./emdb_get_timestamps.pl new_emdb_timestamps.txt >& new_emdb_timestamps.log

# Only modified or new entries should be downloaded (compared to the latest update).
# Make the date comparison for each entry and download appropriately
echo "Step 2: Download updated files based on timestamps..."
date
./emdb_download_updated_files.pl $PRODUCTIONDIR/emdb_timestamps.txt new_emdb_timestamps.txt wget.log >& download_events.log

# Generate an updated set of data files just from the new timestamps. We will have to generate
# it again for the complete database later, after we have recomputed 3DZD's
echo "Step 3: Generate data files from XML information..."
date
./extract_xml_information.sh xml/ idlist.txt name.txt fullname.txt contour_levels.txt resolutions.txt resolution_summary.txt std_dev.txt

# Generate all 3DZD representations required
echo "Step 4: Generate 3DZDs for new maps..."
date
./calculate_volume_and_3dzd.pl contour_levels.txt std_dev.txt map/ volume.txt ./ >& volume_3dzd.log

#### Start updating production
echo "Step 5: Copy new files to production directories..."
date
cp xml/*.xml $PRODUCTIONDIR/xml
cp img/*.gif $PRODUCTIONDIR/img
cp recommend_contour/inv/* $PRODUCTIONDIR/recommend_contour/inv
cp 1_contour/inv/* $PRODUCTIONDIR/1_contour/inv
cp 2_contour/inv/* $PRODUCTIONDIR/2_contour/inv
cp 1_2_contour/inv/* $PRODUCTIONDIR/1_2_contour/inv
cp 1_std_merge/inv/* $PRODUCTIONDIR/1_std_merge/inv

# Compute the new data files for the production side
echo "Step 6: Regenerate data files from XMLs, for the complete production database..."
date
./extract_xml_information.sh $PRODUCTIONDIR/xml/ idlist.txt name.txt fullname.txt contour_levels.txt resolutions.txt resolution_summary.txt std_dev.txt
# This only computes the volumes (not 3DZD's)
echo "Step 7: Merge volumes with the ones in production..."
date
sort -k 1,1 -g $PRODUCTIONDIR/volume.txt > old_volumes.txt
sort -k 1,1 -g volume.txt > new_volumes.txt
join -a 1 -a 2 old_volumes.txt new_volumes.txt | awk '{if(NF==2) {print $1"\t"$2;} else {print $1"\t"$3;}}' > volume.txt
echo "Step 8: Copy all data files to production..."
date
cp idlist.txt name.txt fullname.txt resolutions.txt resolution_summary.txt volume.txt std_dev.txt $PRODUCTIONDIR

# Update the timestamps file for the next update
echo "Step 9: Update the timestamps file in the production directory..."
date
mv new_emdb_timestamps.txt $PRODUCTIONDIR/emdb_timestamps.txt
