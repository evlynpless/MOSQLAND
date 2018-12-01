cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/SW

#gdal_translate -a_ullr -180 85 180 -60 -co COMPRESS=DEFLATE -co ZLEVEL=9 occurrence_250m.tif occurrence_250m_tmp.tif

#gdalwarp -tr 0.008333333333333 0.008333333333333 occurrence_250m.tif occurrence_1km.tif

#gdal_translate -a_nodata 255 occurrence_1km.tif occurrence_1km.tif

#skip for now
#pksetmask -i occurrence_250m_tmp.tif -m occurrence_250m_tmp.tif -o occurrence_250m_tmp.tif --msknodata -9999 -nodata 255

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/SW/Florida_Clips

gdal_translate  -ot Float32  -projwin -85 31.5 -79.8 24.0 occurrence_1km.tif  $OUTDIR/occurrence_1km_FloridaClip.tif

gdal_calc.py -A $OUTDIR/occurrence_1km_FloridaClip.tif  --outfile=$OUTDIR/occurrence_1km_FloridaClip_positive.tif --NoDataValue=255  --calc="A+0.01"


gdal_translate -of AAIGrid occurrence_1km_FloridaClip_positive.tif occurrence_1km_FloridaClip\
_positive.asc
