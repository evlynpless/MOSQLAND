#This script is to fix the corners and other features of the friction data (Malaria project)

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/friction                                        

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/friction/Florida_Clips

gdal_translate -a_ullr -180 85 +180 -60 -co COMPRESS=DEFLATE -co ZLEVEL=9 friction_surface_2015_v1.0.tif                                            

pksetmask -i friction_surface_2015_v1.0.tif -m friction_surface_2015_v1.0.tif -o friction_surface_2015_v1.0.tif --msknodata -9999 -nodata -999

gdal_translate  -projwin -85 31.5 -79.8 24.0 friction_surface_2015_v1.0.tif  $OUTDIR/friction_surface_2015_v1.0_FloridaClip.tif

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/friction/Florida_Clips

gdal_translate -of AAIGrid friction_surface_2015_v1.0_FloridaClip.tif friction_surface_2015_v1.0_Florida_Clip.asc

gdal_translate  -projwin -113.1 35.8 -100.5 28.7 friction_surface_2015_v1.0.tif SW_clip/friction_surface_2015_v1.0_SWclip.tif
