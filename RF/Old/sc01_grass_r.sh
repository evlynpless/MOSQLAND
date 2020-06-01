#!/bin/bash
#SBATCH -p day
#SBATCH -J sc01_grass_r.sh
#SBATCH -n 1 -c 8 -N 1
#SBATCH -t 1:00:00 
#SBATCH -o /gpfs/scratch60/fas/powell/esp38/stdout/sc01_grass_r.sh.%J.out 
#SBATCH -e /gpfs/scratch60/fas/powell/esp38/stderr/sc01_grass_r.sh.%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evlyn.pless@yale.edu 
#SBATCH --mem=40G

####### for point in  $(seq 1 38 )  ; do sbatch --export=point=$point   /home/esp38/scripts/MOSQLAND/RF/sc01_grass_r.sh ; done
####### for point in  $(seq 1 1  )  ; do sbatch --export=point=$point   /home/esp38/scripts/MOSQLAND/RF/sc01_grass_r.sh ; done
######  bash /home/esp38/scripts/MOSQLAND/RF/sc01_grass_r.sh

module load GEOS/3.6.2-foss-2018a-Python-3.6.4
module load GDAL/3.1.0-foss-2018a-Python-3.6.4
module load GRASS/7.8.0-foss-2018a-Python-3.6.4

module load R/3.5.3-foss-2018a-X11-20180131

export IN_TXT=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTesting
export IN_MSQ=/project/fas/powell/esp38/dataproces/MOSQLAND
export RAM=/dev/shm

#### giuseppe output location
export OUT_TXT=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTesting
export point=$point
###  export point=1

###  spliting in training and testing. Select only one point (all the pairwise from that point)  for the testing

#Evie: I can  update FST_list_NAmRF3.csv so it has the new Houston data

echo "index,V1,V2,locality1,locality2,lat1,long1,lat2,long2,FST_lin,CSE,Resd" > $OUT_TXT/FST_list_NAmRF3_Test$point.csv
awk -v point=$point  -F ","  '{ if ($2==point  || $3==point  ) print   }' $IN_TXT/FST_list_NAmRF4_Houston_linux.csv  >> $OUT_TXT/FST_list_NAmRF3_Test$point.csv
awk -v point=$point  -F ","  '{ if ($2!=point  && $3!=point  ) print   }' $IN_TXT/FST_list_NAmRF4_Houston_linux.csv  > $OUT_TXT/FST_list_NAmRF3_Trai$point.csv

#### create the start-end points for each line

awk -F "," '{ if(NR>1) printf ("%f %f\n%f %f\n%s\n" , $(NF-5), $(NF-6), $(NF-3), $(NF-4), "NaN NaN" ) }' $OUT_TXT/FST_list_NAmRF3_Test$point.csv > $OUT_TXT/FST_line_NAmRF3_Test$point.txt
awk -F "," '{ if(NR>1) printf ("%f %f\n%f %f\n%s\n" , $(NF-5), $(NF-6), $(NF-3), $(NF-4), "NaN NaN" ) }' $OUT_TXT/FST_list_NAmRF3_Trai$point.csv > $OUT_TXT/FST_line_NAmRF3_Trai$point.txt 

### enter in grass and get the mean under straight line
## grass78  -f -text --tmp-location  -c $IN_MSQ/consland/ARIDITY/NAm_clip/AI_annual_NAmClip2_Int16.tif    <<'EOF'

## copy files to the ram 

rm -fr $RAM/grassdb$point
mkdir  $RAM/grassdb$point

ls $IN_MSQ/consland/ARIDITY/NAm_clip/AI_annual_NAmClip2_Int16.tif $IN_MSQ/consland/access/NAm_clip/accessibility_to_cities_2015_v1.0_NAmClip2_Int16.tif $IN_MSQ/consland/chelsa/bio12/NAm_clip/bio12_mean_NAmClip2_Int16.tif $IN_MSQ/consland/chelsa/bio1/NAm_clip/bio1_NAmClip2_Int16.tif $IN_MSQ/consland/GSHL/NAm_clip/GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84_NAmClip2.tif $IN_MSQ/consland/friction/NAm_clip/friction_surface_2015_v1.0_NAmClip2.tif $IN_MSQ/consland/chelsa/bio6/NAm_clip/bio6_mean_NAmClip2_Int16.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_1_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_2_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_3_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_4_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_5_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_6_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_7_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_8_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_9_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_10_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_11_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_12_NAmClip2.tif $IN_MSQ/consland/MERIT/slope/NAm_clip/slope_1KMmedian_MERIT_NAmClip2.tif $IN_MSQ/consland/MERIT/altitude/NAm_clip/altitude_1KMmedian_MERIT_NAmClip2_Int16.tif $IN_MSQ/consland/PET/NAm_clip/pet_mean_NAmClip2.tif $IN_MSQ/consland/chelsa/bio2/NAm_clip/bio2_mean_NAmClip2_Int16.tif $IN_MSQ/consland/chelsa/bio5/NAm_clip/bio5_mean_NAmClip2_Int16.tif $IN_MSQ/consland/chelsa/bio7/NAm_clip/bio7_mean_NAmClip2_Int16.tif $IN_MSQ/consland/chelsa/bio13/NAm_clip/bio13_mean_NAmClip2_Int16.tif $IN_MSQ/consland/chelsa/bio14/NAm_clip/bio14_mean_NAmClip2.tif $IN_MSQ/consland/GPP/mnth/monthly_mean/NAm_clip/GPP_mean_NAmClip2_Int16.tif $IN_MSQ/consland/kernel/KernelRas_100m_fnl.tif | xargs -n 1 -P 4 bash -c $' cp $1 $RAM/grassdb$point  ' _ 



grass78 -f -text -c $RAM/grassdb$point/AI_annual_NAmClip2_Int16.tif  $RAM/grassdb$point/loc$point   <<'EOF'

r.external  input=$RAM/grassdb$point/AI_annual_NAmClip2_Int16.tif  output=arid    --overwrite
r.external  input=$RAM/grassdb$point/accessibility_to_cities_2015_v1.0_NAmClip2_Int16.tif  output=access    --overwrite
r.external  input=$RAM/grassdb$point/bio12_mean_NAmClip2_Int16.tif  output=prec    --overwrite
r.external  input=$RAM/grassdb$point/bio1_NAmClip2_Int16.tif  output=meantemp    --overwrite
r.external  input=$RAM/grassdb$point/GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84_NAmClip2.tif  output=humandensity    --overwrite
r.external  input=$RAM/grassdb$point/friction_surface_2015_v1.0_NAmClip2.tif  output=friction    --overwrite
r.external  input=$RAM/grassdb$point/bio6_mean_NAmClip2_Int16.tif  output=mintemp    --overwrite
r.external  input=$RAM/grassdb$point/consensus_full_class_1_NAmClip2.tif  output=Needleleaf    --overwrite
r.external  input=$RAM/grassdb$point/consensus_full_class_2_NAmClip2.tif  output=EvBroadleaf    --overwrite
r.external  input=$RAM/grassdb$point/consensus_full_class_3_NAmClip2.tif  output=DecBroadleaf    --overwrite
r.external  input=$RAM/grassdb$point/consensus_full_class_4_NAmClip2.tif  output=MiscTrees    --overwrite
r.external  input=$RAM/grassdb$point/consensus_full_class_5_NAmClip2.tif  output=Shrubs    --overwrite
r.external  input=$RAM/grassdb$point/consensus_full_class_6_NAmClip2.tif  output=Herb    --overwrite
r.external  input=$RAM/grassdb$point/consensus_full_class_7_NAmClip2.tif  output=Crop    --overwrite
r.external  input=$RAM/grassdb$point/consensus_full_class_8_NAmClip2.tif  output=Flood    --overwrite
r.external  input=$RAM/grassdb$point/consensus_full_class_9_NAmClip2.tif  output=Urban    --overwrite
r.external  input=$RAM/grassdb$point/consensus_full_class_10_NAmClip2.tif  output=Snow    --overwrite
r.external  input=$RAM/grassdb$point/consensus_full_class_11_NAmClip2.tif  output=Barren    --overwrite
r.external  input=$RAM/grassdb$point/consensus_full_class_12_NAmClip2.tif  output=Water    --overwrite
r.external  input=$RAM/grassdb$point/slope_1KMmedian_MERIT_NAmClip2.tif  output=Slope    --overwrite
r.external  input=$RAM/grassdb$point/altitude_1KMmedian_MERIT_NAmClip2_Int16.tif  output=Altitude    --overwrite
r.external  input=$RAM/grassdb$point/pet_mean_NAmClip2.tif  output=PET    --overwrite
r.external  input=$RAM/grassdb$point/bio2_mean_NAmClip2_Int16.tif  output=DailyTempRange    --overwrite
r.external  input=$RAM/grassdb$point/bio5_mean_NAmClip2_Int16.tif  output=maxtemp    --overwrite
r.external  input=$RAM/grassdb$point/bio7_mean_NAmClip2_Int16.tif  output=AnnualTempRange    --overwrite
r.external  input=$RAM/grassdb$point/bio13_mean_NAmClip2_Int16.tif  output=precwet    --overwrite
r.external  input=$RAM/grassdb$point/bio14_mean_NAmClip2.tif  output=precdry    --overwrite
r.external  input=$RAM/grassdb$point/GPP_mean_NAmClip2_Int16.tif  output=GPP    --overwrite
r.external  input=$RAM/grassdb$point/KernelRas_100m_fnl.tif  output=kernel100    --overwrite

## add all the predictors 

v.in.lines input=$OUT_TXT/FST_line_NAmRF3_Test$point.txt output=FST_line_NAmRF3_Test$point  separator=" " --overwrite
v.in.lines input=$OUT_TXT/FST_line_NAmRF3_Trai$point.txt output=FST_line_NAmRF3_Trai$point  separator=" " --overwrite

echo  extract mean for FST_line_NAmRF3_Trai$point

v.db.addtable map=FST_line_NAmRF3_Trai$point   

echo "cat,arid,access,prec,meantemp,humandensity,friction,mintemp,Needleleaf,EvBroadleaf,DecBroadleaf,MiscTrees,Shrubs,Herb,Crop,Flood,Urban,Snow,Barren,Water,Slope,Altitude,PET,DailyTempRange,maxtemp,AnnualTempRange,precwet,precdry,GPP,kernel100" > $OUT_TXT/FST_list_NAmRF3_PredictTrai$point.csv 

awk -F "," '{ if (NR>1)  print $1   }'   $OUT_TXT/FST_list_NAmRF3_Trai$point.csv | xargs -n 1 -P 8 bash -c $' 
LINE=$1 
v.to.rast  input=FST_line_NAmRF3_Trai1 where="cat == $LINE " output=raster$LINE use="cat" --o   2>/dev/null 

echo $LINE","$(for rast in arid access prec meantemp humandensity friction mintemp Needleleaf EvBroadleaf DecBroadleaf MiscTrees Shrubs Herb Crop Flood Urban Snow Barren Water Slope Altitude PET DailyTempRange maxtemp AnnualTempRange precwet precdry GPP ; do  r.univar -t map=$rast    zones=raster$LINE separator=comma 2>/dev/null | awk  -F , \' { if (NR==2)  printf  ("%s," , $8 ) } \'  ; done )$(  r.univar -t map=kernel100   zones=raster$LINE  separator=comma 2>/dev/null  | awk  -F ,   \' { if (NR==2)  printf  ("%s\\n",$8 ) } \' )

' _   >> $OUT_TXT/FST_list_NAmRF3_PredictTrai$point.csv 


echo  extract mean for FST_line_NAmRF3_Test$point

v.db.addtable map=FST_line_NAmRF3_Test$point   

echo "cat,arid,access,prec,meantemp,humandensity,friction,mintemp,Needleleaf,EvBroadleaf,DecBroadleaf,MiscTrees,Shrubs,Herb,Crop,Flood,Urban,Snow,Barren,Water,Slope,Altitude,PET,DailyTempRange,maxtemp,AnnualTempRange,precwet,precdry,GPP,kernel100" > $OUT_TXT/FST_list_NAmRF3_PredictTest$point.csv 

awk -F "," '{ if (NR>1)  print $1   }'   $OUT_TXT/FST_list_NAmRF3_Test$point.csv | xargs -n 1 -P 8 bash -c $' 
LINE=$1 
v.to.rast  input=FST_line_NAmRF3_Test1 where="cat == $LINE " output=raster$LINE use="cat" --o   2>/dev/null 

echo $LINE","$(for rast in arid access prec meantemp humandensity friction mintemp Needleleaf EvBroadleaf DecBroadleaf MiscTrees Shrubs Herb Crop Flood Urban Snow Barren Water Slope Altitude PET DailyTempRange maxtemp AnnualTempRange precwet precdry GPP ; do  r.univar -t map=$rast    zones=raster$LINE separator=comma 2>/dev/null | awk  -F , \' { if (NR==2)  printf  ("%s," , $8 ) } \'  ; done )$(  r.univar -t map=kernel100   zones=raster$LINE  separator=comma 2>/dev/null  | awk  -F ,   \' { if (NR==2)  printf  ("%s\\n",$8 ) } \' )

' _   >> $OUT_TXT/FST_list_NAmRF3_PredictTest$point.csv 


EOF

rm -fr $RAM/grassdb$point


R --vanilla --no-readline   -q  << 'EOF'


library("randomForest")
library("rgdal")
library("raster")
library("sp")
library("gdistance")

# library("tidyverse")

point <- Sys.getenv(c('point'))
point=1

Env.table.train <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTesting/FST_list_NAmRF3_PredictTrai" , point , ".csv"), sep=",", header=T)
Gen.table.train <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTesting/FST_list_NAmRF3_Trai", point         , ".csv"), sep=",", header=T)

Env.table.test <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTesting/FST_list_NAmRF3_PredictTest" , point , ".csv"), sep=",", header=T)
Gen.table.test <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTesting/FST_list_NAmRF3_Test", point         , ".csv"), sep=",", header=T)

###  Rename columns 

names(Env.table.train) <- c("index","arid","access","prec","mean.temp","human.density","friction","min.temp","Needleleaf","EvBroadleaf","DecBroadleaf","MiscTrees","Shrubs","Herb","Crop","Flood","Urban","Snow","Barren","Water","Slope","Altitude","PET","DailyTempRange","max.temp","AnnualTempRange","prec.wet","prec.dry","GPP","kernel100") 

names(Env.table.test) <-  c("index","arid","access","prec","mean.temp","human.density","friction","min.temp","Needleleaf","EvBroadleaf","DecBroadleaf","MiscTrees","Shrubs","Herb","Crop","Flood","Urban","Snow","Barren","Water","Slope","Altitude","PET","DailyTempRange","max.temp","AnnualTempRange","prec.wet","prec.dry","GPP","kernel100") 

### Add genetic distance column to the new data frame

Env.table.test$CSE   <- Gen.table.test$CSE  
Env.table.train$CSE  <- Gen.table.train$CSE  

set.seed(NULL)

load(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12_rasterstack_image_withKernel.RData")

#Tuning parameters for random forest
tune_x <- Env.table.train[,names(env)]
tune_y <- Env.table.train[,c("CSE")]
bestmtry <- tuneRF(tune_x, tune_y, stepFactor=1.5, improve=1e-5, ntree=500)
mtry_opt <- bestmtry[,"mtry"][which.min(bestmtry[,"OOBError"])]

#Run random forest 
#Predictor variables are mean along straight lines through the rasters
#Response variable is genetic distance

Straight_RF = randomForest(CSE ~ arid + access  +   prec  +   mean.temp  +   human.density  +   friction + min.temp + Needleleaf + EvBroadleaf + DecBroadleaf + MiscTrees + Shrubs + Herb + Crop + Flood + Urban + Snow + Barren + Water + Slope + Altitude + PET + DailyTempRange + max.temp + AnnualTempRange + prec.wet + prec.dry + GPP + kernel100, importance=TRUE, mtry = mtry_opt, na.action=na.omit, data=Env.table.test)

pred.cond <- 1 / predict(env, Straight_RF) 
writeRaster(pred.cond,'/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTesting/pred.cond.tif',options=c(COMPRESS=DEFLATE,"ZLEVEL=9" ))


trNAm1 <- transition(pred.cond, transitionFunction=mean, directions=8) #make transitional matrix

writeRaster(trNAm1,'/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTesting/trNAm1.tif',options=c(COMPRESS=DEFLATE,"ZLEVEL=9" ))


EOF
 
exit 
