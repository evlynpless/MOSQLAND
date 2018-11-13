#!/bin/bash
#SBATCH -p scavenge
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evlyn.pless@yale.edu
#SBATCH -o  /gpfs/scratch60/fas/powell/esp38/stdout/sc09_GPP_mean_stdev.sh.%J.out
#SBATCH -e  /gpfs/scratch60/fas/powell/esp38/stderr/sc09_GPP_mean_stdev.sh.%J.err

##   sbatch  ~/scripts/MOSQLAND/sc09_GPP_mean_stdev.sh  

#Corner coordinates have been checked already
#This code is to create a monthly averages for our GPP data across the years 2000-2016. 

INDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/GPP/mnth/monthly

 

for month in 01 02 03 04 05 06 07 08 09 10 11 12 ; do 
#Create a multiband vrt 
gdalbuildvrt -overwrite -separate ${INDIR}_proces/GPP.VPM.M${month}.v20.CMG.vrt  $INDIR/GPP.VPM.*.M${month}.v20.CMG.tif 

echo  Calculate mean and standard deviation
pkstatprofile -co  COMPRESS=LZW -nodata -9999 -f mean -f stdev  -i ${INDIR}_proces/GPP.VPM.M${month}.v20.CMG.vrt  -o ${INDIR}_proces/GPP.VPM.mean_stdev.M${month}.v20.CMG.tif

echo split tif uniq band 

gdal_translate -a_srs EPSG:4326 -nodata -9999 -b 1 -co COMPRESS=LZW ${INDIR}_proces/GPP.VPM.mean_stdev.M${month}.v20.CMG.tif  ${INDIR}_mean/GPP.VPM.mean.M${month}.v20.CMG.tif

gdal_translate -a_srs EPSG:4326 -nodata -9999 -b 2 -co COMPRESS=LZW ${INDIR}_proces/GPP.VPM.mean_stdev.M${month}.v20.CMG.tif  ${INDIR}_stdev/GPP.VPM.stdev.M${month}.v20.CMG.tif

done 


#for month in 01 02 03 04 05 06 07 08 09 10 11 12 ; do 

#gdal_edit.py -a_srs EPSG:4326 -a_nodata -9999 ${INDIR}_mean/GPP.VPM.mean.M${month}.v20.CMG.tif

#gdal_edit.py -a_nodata -9999 ${INDIR}_mean/GPP.VPM.mean.M${month}.v20.CMG.tif

#pksetmask   -co COMPRESS=DEFLATE -co ZLEVEL=9  -m ${INDIR}_mean/GPP.VPM.mean.M${month}.v20.CMG.tif -msknodata -9999 -p '='  -nodata -9999 -i ${INDIR}_stdev/GPP.VPM.stdev.M${month}.v20.CMG.tif  -o ${INDIR}_stdev/GPP.VPM.stdev.M${month}.v20.CMG_nodata.tif

#mv ${INDIR}_stdev/GPP.VPM.stdev.M${month}.v20.CMG_nodata.tif ${INDIR}_stdev/GPP.VPM.stdev.M${month}.v20.CMG.tif

#done 
