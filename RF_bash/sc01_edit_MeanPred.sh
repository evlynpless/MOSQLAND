cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/LinFST_best_preds
pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m Avg_surface.tif  -i Avg_surface.tif    --msknodata -1 -p "<"  -nodata -1  -o Avg_surface_msk.tif 
gdal_translate -a_nodata -1 Avg_surface_msk.tif Avg_surface_final.tif
