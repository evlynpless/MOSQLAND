cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/ABSHUM

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/ABSHUM/Florida_clips

gdal_edit.py -a_ullr -180   90  180  -60 ABS50S.tif

gdal_translate  -projwin -85 31.5 -79.8 24.0  ABS50S.tif  $OUTDIR/ABS50S_FloridaClip.tif

#gdal_translate -of AAIGrid $OUTDIR/ABS50S_FloridaClip.tif $OUTDIR/ABS50S_FloridaClip.asc
