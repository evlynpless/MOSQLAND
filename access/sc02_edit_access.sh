#This script is to fix the corners and other features of the accessibility to cities data (Weiss et al) and the friction data (Malaria project)

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/access

gdal_translate -a_ullr -180 84 180 -90 -co COMPRESS=DEFLATE -co ZLEVEL=9 accessibility_to_cities_2015_v1.0.tif 

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/access/Florida_clips

gdal_translate  -projwin -85 31.5 -79.8 24.0 accessibility_to_cities_2015_v1.0.tif  $OUTDIR/accessibility_to_cities_2015_v1.0_FloridaClip.tif

gdal_translate -of AAIGrid $OUTDIR/accessibility_to_cities_2015_v1.0_FloridaClip.tif $OUTDIR/accessibility_to_cities_2015_v1.0_FloridaClip.asc
