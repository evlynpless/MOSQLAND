#This script is to fix the resolution and other features of the accessibility to cities data (Weiss et al) and the friction data (Malaria project)

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/access

gdalwarp -r bilinear -tr 0.00833333333333333 0.00833333333333333  -te -180 -60 180 90  accessibility_to_cities_2015_v1.0.tif accessibility_to_cities_2015_v1.0_res.tif - -co COMPRESS=DEFLATE -co ZLEVEL=9

pksetmask -i accessibility_to_cities_2015_v1.0_res.tif  -m accessibility_to_cities_2015_v1.0_res.tif  -o accessibility_to_cities_2015_v1.0_res.tif --msknodata -2147483648 -nodata -9999


OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/access/Florida_clips

gdal_translate  -projwin -85 31.5 -79.8 24.0 accessibility_to_cities_2015_v1.0.tif  $OUTDIR/accessibility_to_cities_2015_v1.0_FloridaClip.tif

gdal_translate  -projwin -126 40 -57 15 accessibility_to_cities_2015_v1.0_res.tif  /project/fas/powell/esp38/dataproces/MOSQLAND/consland/access/NAm_clip/accessibility_to_cities_2015_v1.0_NAmClip.tif

gdal_translate  -projwin -113.1 35.8 -100.5 28.7  accessibility_to_cities_2015_v1.0_res.tif    SW_clip/accessibility_to_cities_2015_v1.0_res_SWclip.tif 
#gdal_translate -of AAIGrid $OUTDIR/accessibility_to_cities_2015_v1.0_FloridaClip.tif $OUTDIR/accessibility_to_cities_2015_v1.0_FloridaClip.asc
