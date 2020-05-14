/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/ResidFST_best_preds

pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m ResidData_Fold1_BestCor2_Pred_it6.tif  -i ResidData_Fold1_BestCor2_Pred_it6.tif    --msknodata -1 -p "<"  -nodata -1  -o ResidData_Fold1_BestCor2_Pred_it6_msk.tif
gdal_translate -a_nodata -1 ResidData_Fold1_BestCor2_Pred_it6_msk.tif ResidData_Fold1_BestCor2_Pred_it6_final.tif

pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m ResidData_Fold2_BestCor2_Pred_it4.tif  -i ResidData_Fold2_BestCor2_Pred_it4.tif    --msknodata -1 -p "<"  -nodata -1  -o ResidData_Fold2_BestCor2_Pred_it4_msk.tif
gdal_translate -a_nodata -1 ResidData_Fold2_BestCor2_Pred_it4_msk.tif ResidData_Fold2_BestCor2_Pred_it4_final.tif

pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m ResidData_Fold3_BestCor2_Pred_it6.tif  -i ResidData_Fold3_BestCor2_Pred_it6.tif    --msknodata -1 -p "<"  -nodata -1  -o ResidData_Fold3_BestCor2_Pred_it6_msk.tif
gdal_translate -a_nodata -1 ResidData_Fold3_BestCor2_Pred_it6_msk.tif ResidData_Fold3_BestCor2_Pred_it6_final.tif

pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m ResidData_Fold4_BestCor2_Pred_it1.tif  -i ResidData_Fold4_BestCor2_Pred_it1.tif    --msknodata -1 -p "<"  -nodata -1  -o ResidData_Fold4_BestCor2_Pred_it1_msk.tif
gdal_translate -a_nodata -1 ResidData_Fold4_BestCor2_Pred_it1_msk.tif ResidData_Fold4_BestCor2_Pred_it1_final.tif

pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m ResidData_Fold5_BestCor2_Pred_it3.tif  -i ResidData_Fold5_BestCor2_Pred_it3.tif    --msknodata -1 -p "<"  -nodata -1  -o ResidData_Fold5_BestCor2_Pred_it3_msk.tif
gdal_translate -a_nodata -1 ResidData_Fold5_BestCor2_Pred_it3_msk.tif  ResidData_Fold5_BestCor2_Pred_it3_final.tif

pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m ResidData_Fold6_BestCor2_Pred_it5.tif  -i ResidData_Fold6_BestCor2_Pred_it5.tif    --msknodata -1 -p "<"  -nodata -1  -o ResidData_Fold6_BestCor2_Pred_it5_msk.tif
gdal_translate -a_nodata -1 ResidData_Fold6_BestCor2_Pred_it5_msk.tif ResidData_Fold6_BestCor2_Pred_it5_final.tif

pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m ResidData_Fold7_BestCor2_Pred_it3.tif  -i ResidData_Fold7_BestCor2_Pred_it3.tif    --msknodata -1 -p "<"  -nodata -1  -o ResidData_Fold7_BestCor2_Pred_it3_msk.tif
gdal_translate -a_nodata -1 ResidData_Fold7_BestCor2_Pred_it3_msk.tif ResidData_Fold7_BestCor2_Pred_it3_final.tif

pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m ResidData_Fold8_BestCor2_Pred_it4.tif  -i ResidData_Fold8_BestCor2_Pred_it4.tif    --msknodata -1 -p "<"  -nodata -1  -o ResidData_Fold8_BestCor2_Pred_it4_msk.tif
gdal_translate -a_nodata -1 ResidData_Fold8_BestCor2_Pred_it4_msk.tif ResidData_Fold8_BestCor2_Pred_it4_final.tif

pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m ResidData_Fold9_BestCor2_Pred_it6.tif  -i ResidData_Fold9_BestCor2_Pred_it6.tif    --msknodata -1 -p "<"  -nodata -1  -o ResidData_Fold9_BestCor2_Pred_it6_msk.tif
gdal_translate -a_nodata -1 ResidData_Fold9_BestCor2_Pred_it6_msk.tif ResidData_Fold9_BestCor2_Pred_it6_final.tif

pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m ResidData_Fold10_BestCor2_Pred_it6.tif  -i ResidData_Fold10_BestCor2_Pred_it6.tif    --msknodata -1 -p "<"  -nodata -1  -o ResidData_Fold10_BestCor2_Pred_it6_msk.tif
gdal_translate -a_nodata -1 ResidData_Fold10_BestCor2_Pred_it6_msk.tif ResidData_Fold10_BestCor2_Pred_it6_final.tif

gdalbuildvrt -separate input.vrt ResidData_*_final.tif 
pkstatprofile -co COMPRESS=DEFLATE -co ZLEVEL=9 -co INTERLEAVE=BAND -nodata -1 -f mean  -f stdev -of GTiff -i  input.vrt  -o output.tif 
gdal_translate -b 1  -co COMPRESS=DEFLATE -co ZLEVEL=9 -co INTERLEAVE=BAND output.tif ResidData_mean.tif 
gdal_translate -b 2  -co COMPRESS=DEFLATE -co ZLEVEL=9 -co INTERLEAVE=BAND output.tif ResidData_stdev.tif 
