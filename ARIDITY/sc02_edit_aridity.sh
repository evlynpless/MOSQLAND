INDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/ARIDITY

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/ARIDITY/Florida_clips

gdal_translate  -projwin -85 31.5 -79.8 24.0  $INDIR/AI_annual.tif  $OUTDIR/AI_annual_FloridaClip.tif

pksetmask -i AI_annual_FloridaClip.tif -m AI_annual_FloridaClip.tif -o AI_annual_FloridaClip_tmp.tif --msknodata -1 -nodata -9999   

gdal_translate -of AAIGrid $OUTDIR/AI_annual_FloridaClip.tif $OUTDIR/AI_annual_FloridaClip.asc

gdal_translate  -projwin -113.1 35.8 -100.5 28.7 AI_annual.tif SW_clip/AI_annual_SWclip.tif
