#This script is to adjust the corners of the landcov data (class 1-12). Dataset from EarthEnv.

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov/Florida_clips

for n in $(seq 1 12 ) ; do 

#gdal_edit.py -a_ullr -180   84  180  -90 consensus_full_class_${n}.tif

#Why is the pixel size getting messed up????!
#I don't think these datasets have no data. They just have zero. So those steps can be skipped
#gdal_edit.py -a_nodata -999 consensus_full_class_${n}.tif
#pksetmask -i consensus_full_class_${n}.tif -m consensus_full_class_${n}.tif -o consensus_full_class_${n}.tif --msknodata 255 -nodata -999                                                   

#gdal_translate  -projwin -85 31.5 -79.8 24.0  consensus_full_class_${n}.tif  $OUTDIR/consensus_full_class_${n}_FloridaClip.tif

gdal_translate -of AAIGrid $OUTDIR/consensus_full_class_${n}_FloridaClip.tif $OUTDIR/consensus_full_class_${n}_FloridaClip.asc


done
