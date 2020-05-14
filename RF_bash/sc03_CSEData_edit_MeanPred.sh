cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc15/CSE_best_preds

pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m CSEData_Fold1_BestCor2_Pred_it1.tif  -i CSEData_Fold1_BestCor2_Pred_it1.tif   --msknodata -1 -p "<"  -nodata -1  -o CSEData_Fold1_BestCor2_Pred_it1_msk.tif
gdal_translate -a_nodata -1 CSEData_Fold1_BestCor2_Pred_it1_msk.tif CSEData_Fold1_BestCor2_Pred_it1_final.tif

pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m CSEData_Fold2_BestCor2_Pred_it2.tif  -i CSEData_Fold2_BestCor2_Pred_it2.tif    --msknodata -1 -p "<"  -nodata -1  -o CSEData_Fold2_BestCor2_Pred_it2_msk.tif
gdal_translate -a_nodata -1 CSEData_Fold2_BestCor2_Pred_it2_msk.tif CSEData_Fold2_BestCor2_Pred_it2_final.tif

pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m CSEData_Fold3_BestCor2_Pred_it6.tif   -i CSEData_Fold3_BestCor2_Pred_it6.tif    --msknodata -1 -p "<"  -nodata -1  -o CSEData_Fold3_BestCor2_Pred_it6_msk.tif
gdal_translate -a_nodata -1 CSEData_Fold3_BestCor2_Pred_it6_msk.tif CSEData_Fold3_BestCor2_Pred_it6_final.tif

pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m CSEData_Fold4_BestCor2_Pred_it2.tif  -i CSEData_Fold4_BestCor2_Pred_it2.tif    --msknodata -1 -p "<"  -nodata -1  -o CSEData_Fold4_BestCor2_Pred_it2_msk.tif
gdal_translate -a_nodata -1 CSEData_Fold4_BestCor2_Pred_it2_msk.tif CSEData_Fold4_BestCor2_Pred_it2_final.tif

pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m CSEData_Fold5_BestCor2_Pred_it1.tif  -i CSEData_Fold5_BestCor2_Pred_it1.tif    --msknodata -1 -p "<"  -nodata -1  -o CSEData_Fold5_BestCor2_Pred_it1_msk.tif
gdal_translate -a_nodata -1 CSEData_Fold5_BestCor2_Pred_it1_msk.tif  CSEData_Fold5_BestCor2_Pred_it1_final.tif

pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m CSEData_Fold6_BestCor2_Pred_it0.tif  -i CSEData_Fold6_BestCor2_Pred_it0.tif    --msknodata -1 -p "<"  -nodata -1  -o CSEData_Fold6_BestCor2_Pred_it0_msk.tif
gdal_translate -a_nodata -1 CSEData_Fold6_BestCor2_Pred_it0_msk.tif CSEData_Fold6_BestCor2_Pred_it0_final.tif

pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m CSEData_Fold7_BestCor2_Pred_it1.tif  -i CSEData_Fold7_BestCor2_Pred_it1.tif    --msknodata -1 -p "<"  -nodata -1  -o CSEData_Fold7_BestCor2_Pred_it1_msk.tif
gdal_translate -a_nodata -1 CSEData_Fold7_BestCor2_Pred_it1_msk.tif CSEData_Fold7_BestCor2_Pred_it1_final.tif

pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m CSEData_Fold8_BestCor2_Pred_it1.tif  -i CSEData_Fold8_BestCor2_Pred_it1.tif    --msknodata -1 -p "<"  -nodata -1  -o CSEData_Fold8_BestCor2_Pred_it1_msk.tif
gdal_translate -a_nodata -1 CSEData_Fold8_BestCor2_Pred_it1_msk.tif CSEData_Fold8_BestCor2_Pred_it1_final.tif

pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m CSEData_Fold9_BestCor2_Pred_it4.tif  -i CSEData_Fold9_BestCor2_Pred_it4.tif    --msknodata -1 -p "<"  -nodata -1  -o CSEData_Fold9_BestCor2_Pred_it4_msk.tif
gdal_translate -a_nodata -1 CSEData_Fold9_BestCor2_Pred_it4_msk.tif CSEData_Fold9_BestCor2_Pred_it4_final.tif

pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m CSEData_Fold10_BestCor2_Pred_it2.tif  -i CSEData_Fold10_BestCor2_Pred_it2.tif    --msknodata -1 -p "<"  -nodata -1  -o CSEData_Fold10_BestCor2_Pred_it2_msk.tif
gdal_translate -a_nodata -1 CSEData_Fold10_BestCor2_Pred_it2_msk.tif CSEData_Fold10_BestCor2_Pred_it2_final.tif

gdalbuildvrt -separate input.vrt CSEData_*_final.tif
pkstatprofile -co COMPRESS=DEFLATE -co ZLEVEL=9 -co INTERLEAVE=BAND -nodata -1 -f mean  -f stdev -of GTiff -i  input.vrt  -o output.tif
gdal_translate -b 1  -co COMPRESS=DEFLATE -co ZLEVEL=9 -co INTERLEAVE=BAND output.tif CSEData_mean.tif
gdal_translate -b 2  -co COMPRESS=DEFLATE -co ZLEVEL=9 -co INTERLEAVE=BAND output.tif CSEData_stdev.tif

