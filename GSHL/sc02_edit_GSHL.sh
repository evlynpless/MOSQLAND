cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/GSHL

#gdal_translate -a_ullr -180 85 180 -60 -co COMPRESS=DEFLATE -co ZLEVEL=9 GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84.tif GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84_tmp.tif

#mv GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84_tmp.tif GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS8.tif

#pksetmask -i GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84.tif -m GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84.tif -o GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84.tif --msknodata -1 -nodata -999

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/GSHL/Florida_clips

gdal_translate  -projwin -85 31.5 -79.8 24.0 GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84.tif  $OUTDIR/GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84_FloridaClip.tif

gdal_translate -of AAIGrid GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84.tif GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84.asc

#Also creating a version with no zeros for circuitscape

gdal_calc.py -A  GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84_FloridaClip.tif  --outfile=GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84_FloridaClip_positive.tif --NoDataValue=-999  --calc="A+0.01"

gdal_translate -of AAIGrid GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84_FloridaClip_positive.tif GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84_FloridaClip_positive.asc
