#!/bin/bash
#SBATCH -p day
#SBATCH -J sc01_grass_r_jobarray.sh
#SBATCH -n 1 -c 16 -N 1
#SBATCH -t 24:00:00 
#SBATCH -o /gpfs/scratch60/fas/powell/esp38/stdout/sc01_grass_r_jobarray.sh.%A_%a.out
#SBATCH -e /gpfs/scratch60/fas/powell/esp38/stderr/sc01_grass_r_jobarray.sh.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evlyn.pless@yale.edu
#SBATCH --array=1
#SBATCH --mem=80G

####### sbatch  /home/fas/powell/esp38/scripts/MOSQLAND/RF/sc01_grass_r_jobarray.sh
###### 

module load GEOS/3.6.2-foss-2018a-Python-3.6.4
module load GDAL/3.1.0-foss-2018a-Python-3.6.4
module load GSL/2.3-GCCcore-6.4.0
module load Boost/1.66.0-foss-2018a
module load PKTOOLS/2.6.7.6-foss-2018a-Python-3.6.4
module load Armadillo/8.400.0-foss-2018a-Python-3.6.4
module load GRASS/7.8.0-foss-2018a-Python-3.6.4

module load R/3.5.3-foss-2018a-X11-20180131

export IN_TXT=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3
export IN_MSQ=/project/fas/powell/esp38/dataproces/MOSQLAND
export RAM=/dev/shm
export CPU=$SLURM_CPUS_ON_NODE

#### evie output location
export OUT_TXT=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc
export point=$SLURM_ARRAY_TASK_ID
###  export point=1

###  spliting in training and testing. Select only one point (all the pairwise from that point)  for the testing

#Evie: I can  update FST_list_NAmRF3.csv so it has the new Houston data

echo "index,V1,V2,locality1,locality2,lat1,long1,lat2,long2,FST_lin,CSE,Resd" > ${OUT_TXT}_${point}/FST_list_NAmRF3_Test$point.csv
awk -v point=$point  -F ","  '{ if ($2==point || $3==point ) print }' $IN_TXT/FST_list_NAmRF3_linux.csv  >> ${OUT_TXT}_${point}/FST_list_NAmRF3_Test$point.csv
awk -v point=$point  -F ","  '{ if ($2!=point && $3!=point ) print }' $IN_TXT/FST_list_NAmRF3_linux.csv  > ${OUT_TXT}_${point}/FST_list_NAmRF3_Trai$point.csv

#### create the start-end points for each line

awk -F "," '{ if(NR>1) printf ("%f %f\n%f %f\n%s\n" , $(NF-5), $(NF-6), $(NF-3), $(NF-4), "NaN NaN" ) }' ${OUT_TXT}_${point}/FST_list_NAmRF3_Test$point.csv > ${OUT_TXT}_${point}/FST_line_NAmRF3_Test$point.txt
awk -F "," '{ if(NR>1) printf ("%f %f\n%f %f\n%s\n" , $(NF-5), $(NF-6), $(NF-3), $(NF-4), "NaN NaN" ) }' ${OUT_TXT}_${point}/FST_list_NAmRF3_Trai$point.csv > ${OUT_TXT}_${point}/FST_line_NAmRF3_Trai$point.txt 

### enter in grass and get the mean under straight line

## copy files to the ram to speed-up read and write 

rm -fr $RAM/grassdb$point
mkdir  $RAM/grassdb$point

ls $IN_MSQ/consland/ARIDITY/NAm_clip/AI_annual_NAmClip2_Int16.tif $IN_MSQ/consland/access/NAm_clip/accessibility_to_cities_2015_v1.0_NAmClip2_Int16.tif $IN_MSQ/consland/chelsa/bio12/NAm_clip/bio12_mean_NAmClip2_Int16.tif $IN_MSQ/consland/chelsa/bio1/NAm_clip/bio1_NAmClip2_Int16.tif $IN_MSQ/consland/GSHL/NAm_clip/GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84_NAmClip2.tif $IN_MSQ/consland/friction/NAm_clip/friction_surface_2015_v1.0_NAmClip2.tif $IN_MSQ/consland/chelsa/bio6/NAm_clip/bio6_mean_NAmClip2_Int16.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_1_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_2_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_3_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_4_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_5_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_6_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_7_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_8_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_9_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_10_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_11_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_12_NAmClip2.tif $IN_MSQ/consland/MERIT/slope/NAm_clip/slope_1KMmedian_MERIT_NAmClip2.tif $IN_MSQ/consland/MERIT/altitude/NAm_clip/altitude_1KMmedian_MERIT_NAmClip2_Int16.tif $IN_MSQ/consland/PET/NAm_clip/pet_mean_NAmClip2.tif $IN_MSQ/consland/chelsa/bio2/NAm_clip/bio2_mean_NAmClip2_Int16.tif $IN_MSQ/consland/chelsa/bio5/NAm_clip/bio5_mean_NAmClip2_Int16.tif $IN_MSQ/consland/chelsa/bio7/NAm_clip/bio7_mean_NAmClip2_Int16.tif $IN_MSQ/consland/chelsa/bio13/NAm_clip/bio13_mean_NAmClip2_Int16.tif $IN_MSQ/consland/chelsa/bio14/NAm_clip/bio14_mean_NAmClip2.tif $IN_MSQ/consland/GPP/mnth/monthly_mean/NAm_clip/GPP_mean_NAmClip2_Int16.tif $IN_MSQ/consland/kernel/KernelRas_100m_fnl.tif | xargs -n 1 -P $CPU bash -c $' cp $1 $RAM/grassdb$point  ' _ 


# rm -r $OUT_TXT/grassdb$point/loc$point
# grass78 -f -text -c $RAM/grassdb$point/AI_annual_NAmClip2_Int16.tif  $OUT_TXT/grassdb$point/loc$point   <<'EOF'

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

v.in.lines input=${OUT_TXT}_${point}/FST_line_NAmRF3_Test$point.txt output=FST_line_NAmRF3_Test$point  separator=" " --overwrite
v.in.lines input=${OUT_TXT}_${point}/FST_line_NAmRF3_Trai$point.txt output=FST_line_NAmRF3_Trai$point  separator=" " --overwrite

echo addtable for Test and Trai 

v.db.addtable map=FST_line_NAmRF3_Trai$point   
v.db.addtable map=FST_line_NAmRF3_Test$point   

echo  extract mean for straight line in FST_line_NAmRF3_Test$point and FST_line_NAmRF3_Trai$point

for VAR in Trai Test ; do 

export VAR=$VAR
echo "index,arid,access,prec,meantemp,humandensity,friction,mintemp,Needleleaf,EvBroadleaf,DecBroadleaf,MiscTrees,Shrubs,Herb,Crop,Flood,Urban,Snow,Barren,Water,Slope,Altitude,PET,DailyTempRange,maxtemp,AnnualTempRange,precwet,precdry,GPP,kernel100" > ${OUT_TXT}_${point}/FST_list_NAmRF3_Predict${VAR}$point.csv 

awk -F "," '{ if (NR>1)  print NR-1 , $1   }'   ${OUT_TXT}_${point}/FST_list_NAmRF3_${VAR}$point.csv | xargs -n 2 -P $CPU  bash -c $' 
CAT=$1
INDEX=$2 
v.to.rast  input=FST_line_NAmRF3_${VAR}$point where="cat == $CAT " output=raster$CAT  use="cat" --o   2>/dev/null 

echo $INDEX","$(for rast in arid access prec meantemp humandensity friction mintemp Needleleaf EvBroadleaf DecBroadleaf MiscTrees Shrubs Herb Crop Flood Urban Snow Barren Water Slope Altitude PET DailyTempRange maxtemp AnnualTempRange precwet precdry GPP ; do  r.univar -t map=$rast    zones=raster$CAT separator=comma 2>/dev/null | awk  -F , \' { if (NR==2)  printf  ("%s," , $8 ) } \'  ; done )$(r.univar -t map=kernel100   zones=raster$CAT  separator=comma 2>/dev/null   | awk  -F ,   \' { if (NR==2)  printf  ("%s\\n",$8 ) } \' )
g.remove -f  type=raster name=raster$CAT  2>/dev/null
' _   >> ${OUT_TXT}_${point}/FST_list_NAmRF3_Predict${VAR}$point.csv 

done

EOF



##### extract kernel values at point level 

paste -d ","  <(awk -F , '{  if(NR>1) print $2 }' ${OUT_TXT}_${point}/FST_list_NAmRF3_Test$point.csv | uniq )   <(gdallocationinfo -geoloc -wgs84  -valonly   $IN_MSQ/consland/kernel/KernelRas_100m_fnl.tif  $(awk -F , '{  if(NR>1) print  $(NF-5), $(NF-6) }'   ${OUT_TXT}_${point}/FST_list_NAmRF3_Test$point.csv | uniq )) >  ${OUT_TXT}_${point}/FST_list_NAmRF3_KernelTest$point.csv 

awk -F , '{ if(NR>1) print  $(NF-5), $(NF-6) }'   ${OUT_TXT}_${point}/FST_list_NAmRF3_Trai$point.csv | uniq > ${OUT_TXT}_${point}/FST_list_NAmRF3_LatLongTrai$point.txt 
paste -d ","  <(awk -F , '{ if(NR>1) print $2 }'  ${OUT_TXT}_${point}/FST_list_NAmRF3_Trai$point.csv | uniq) <(gdallocationinfo -geoloc -wgs84  -valonly   $IN_MSQ/consland/kernel/KernelRas_100m_fnl.tif  < ${OUT_TXT}_${point}/FST_list_NAmRF3_LatLongTrai$point.txt )   >  ${OUT_TXT}_${point}/FST_list_NAmRF3_KernelTrai$point.csv
rm ${OUT_TXT}_${point}/FST_list_NAmRF3_LatLongTrai$point.txt  

echo start to process in R 

R --vanilla --no-readline   -q  <<'EOF'

library("randomForest")
library("rgdal")
library("raster")
library("gdistance")

point <- Sys.getenv(c('point'))

Env.table.train <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_PredictTrai" , point , ".csv"), sep=",", header=T)
## select only index and CSE
Gen.table.train <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_Trai", point         , ".csv"), sep=",", header=T)[,c( 1, 11)]

Env.table.test <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfscr_", point, "/FST_list_NAmRF3_PredictTest" , point , ".csv"), sep=",", header=T)
Gen.table.test <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_Test", point         , ".csv"), sep=",", header=T)[,c( 1, 11)]

###  Rename columns 

names(Env.table.train) <- c("index","arid","access","prec","mean.temp","human.density","friction","min.temp","Needleleaf","EvBroadleaf","DecBroadleaf","MiscTrees",
"Shrubs","Herb","Crop","Flood","Urban","Snow","Barren","Water","Slope","Altitude","PET","DailyTempRange","max.temp","AnnualTempRange","prec.wet","prec.dry","GPP","kernel100") 
names(Env.table.test) <-  c("index","arid","access","prec","mean.temp","human.density","friction","min.temp","Needleleaf","EvBroadleaf","DecBroadleaf","MiscTrees",
"Shrubs","Herb","Crop","Flood","Urban","Snow","Barren","Water","Slope","Altitude","PET","DailyTempRange","max.temp","AnnualTempRange","prec.wet","prec.dry","GPP","kernel100") 

### Add genetic distance column to the new data frame

Env.table.test  = merge (Env.table.test  ,  Gen.table.test  , by = "index" )
Env.table.train = merge (Env.table.train ,  Gen.table.train , by = "index" )

head(Env.table.train)  ### doble check if the merging is done correctly
head(Env.table.test)   ### doble check if the merging is done correctly

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

Straight_RF = randomForest(CSE ~ arid + access  +   prec  +   mean.temp  +   human.density  +   friction + min.temp + Needleleaf + EvBroadleaf + DecBroadleaf + MiscTrees + 
Shrubs + Herb + Crop + Flood + Urban + Snow + Barren + Water + Slope + Altitude + PET + DailyTempRange + max.temp + AnnualTempRange + prec.wet + prec.dry + GPP + kernel100, 
importance=TRUE, mtry = mtry_opt, na.action=na.omit, data=Env.table.train)

Straight_RF

Cor1 = cor(Straight_RF$predicted, Env.table.train$CSE)
Cor2 = cor(predict(Straight_RF, Env.table.test), Env.table.test$CSE)

paste ("Cor1" , Cor1 )
paste ("Cor2" , Cor2 )

pred.cond <-  1 /  predict(env, Straight_RF)  

#come back to this
writeRaster(pred.cond, paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/prediction.tif"), options=c("COMPRESS=DEFLATE","ZLEVEL=9") , format="GTiff", overwrite=TRUE  )

# writeRaster(trNAm1,"/gpfs/loomis/project/sbsc/ga254/dataproces/MOSQLAND/TrainingTesting/trNAm1.tif",options=c("COMPRESS=DEFLATE","ZLEVEL=9") , format="GTiff", overwrite=TRUE  )

EOF

pksetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -m ${OUT_TXT}_${point}/prediction.tif -msknodata -1 -p "<" -nodata -1 -i ${OUT_TXT}_${point}/prediction.tif -o ${OUT_TXT}_${point}/prediction_msk0.tif
rm ${OUT_TXT}_${point}/prediction.tif 
cp ${OUT_TXT}_${point}/prediction_msk0.tif ${OUT_TXT}_${point}/prediction_msk.tif


### Test and Train  start and stop coordinates 
awk -F "," '{ if(NR>1) print $1 , $(NF-5),  $(NF-6) ,  $(NF-3),  $(NF-4) }' ${OUT_TXT}_${point}/FST_list_NAmRF3_Test$point.csv | uniq > ${OUT_TXT}_${point}/FST_line_NAmRF3_StartStopTest$point.txt
awk -F "," '{ if(NR>1) print $1 , $(NF-5),  $(NF-6) ,  $(NF-3),  $(NF-4) }' ${OUT_TXT}_${point}/FST_list_NAmRF3_Trai$point.csv | uniq > ${OUT_TXT}_${point}/FST_line_NAmRF3_StartStopTrai$point.txt

for ITER in $(seq 1 10 ) ; do  
export ITER=$ITER

echo Start Iteration $ITER  GRASS 

grass78 -f -text $RAM/grassdb$point/loc$point/PERMANENT   <<'EOF'
r.external  input=${OUT_TXT}_${point}/prediction_msk.tif    output=prediction   --overwrite

for VAR in Trai Test ; do 
export VAR=$VAR
echo "index,arid,access,prec,meantemp,humandensity,friction,mintemp,Needleleaf,EvBroadleaf,DecBroadleaf,MiscTrees,Shrubs,Herb,Crop,Flood,Urban,Snow,Barren,Water,Slope,Altitude,PET,DailyTempRange,maxtemp,AnnualTempRange,precwet,precdry,GPP,kernel100" > ${OUT_TXT}_${point}/FST_list_NAmRF3_Iter${ITER}LeastPath${VAR}$point.csv

cat ${OUT_TXT}_${point}/FST_line_NAmRF3_StartStop${VAR}$point.txt   | xargs -n 5 -P $CPU  bash -c $'
INDEX=$1 

r.cost -n -k  input=prediction  output=cost$INDEX  outdir=dir$INDEX    start_coordinates=$2,$3  memory=200  --o --q
g.remove -f  type=raster name=cost$INDEX  2>/dev/null
r.path input=dir$INDEX  raster_path=least_path$INDEX start_coordinates=$4,$5 --o --q 
g.remove -f  type=raster name=dir$INDEX  2>/dev/null

echo $INDEX","$(for rast in arid access prec meantemp humandensity friction mintemp Needleleaf EvBroadleaf DecBroadleaf MiscTrees Shrubs Herb Crop Flood Urban Snow Barren Water Slope Altitude PET DailyTempRange maxtemp AnnualTempRange precwet precdry GPP ; do  r.univar -t map=$rast    zones=least_path$INDEX separator=comma 2>/dev/null | awk  -F , \' { if (NR==2)  printf  ("%s,", $8) } \' ; done )$(r.univar -t map=kernel100 zones=least_path$INDEX separator=comma 2>/dev/null | awk  -F , \' { if (NR==2)  printf  ("%s\\n",$8 ) } \') 
g.remove -f  type=raster name=least_path$INDEX  2>/dev/null

' _  >> ${OUT_TXT}_${point}/FST_list_NAmRF3_Iter${ITER}LeastPath${VAR}$point.csv
done

EOF

echo Start Iteration $ITER   R 

R --vanilla --no-readline   -q  <<'EOF'

library("randomForest")
library("rgdal")
library("raster")
library("gdistance")

point <- Sys.getenv(c('point'))
ITER  <- Sys.getenv(c('ITER'))

Env.table.train <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_Iter",ITER,"LeastPathTrai",point,".csv"), sep=",", header=T)
Gen.table.train <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_Trai",point,".csv"), sep=",", header=T)[,c( 1, 11)]

Env.table.test <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_Iter",ITER,"LeastPathTest",point,".csv"), sep=",", header=T) 
Gen.table.test <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_Test",point, ".csv"), sep=",", header=T)[,c( 1, 11)]

### Add genetic distance column to the new data frame

###  Rename columns 

names(Env.table.train) <- c("index","arid","access","prec","mean.temp","human.density","friction","min.temp","Needleleaf","EvBroadleaf","DecBroadleaf","MiscTrees",
"Shrubs","Herb","Crop","Flood","Urban","Snow","Barren","Water","Slope","Altitude","PET","DailyTempRange","max.temp","AnnualTempRange","prec.wet","prec.dry","GPP","kernel100") 
names(Env.table.test) <-  c("index","arid","access","prec","mean.temp","human.density","friction","min.temp","Needleleaf","EvBroadleaf","DecBroadleaf","MiscTrees",
"Shrubs","Herb","Crop","Flood","Urban","Snow","Barren","Water","Slope","Altitude","PET","DailyTempRange","max.temp","AnnualTempRange","prec.wet","prec.dry","GPP","kernel100") 

### Add genetic distance column to the new data frame

Env.table.test =  merge (Env.table.test , Gen.table.test  , by = "index" )
Env.table.train = merge (Env.table.train, Gen.table.train , by = "index" )

head(Env.table.train)  ### doble check if the merging is done correctly
head(Env.table.test)  ### doble check if the merging is done correctly

set.seed(NULL)

load(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12_rasterstack_image_withKernel.RData")

#Tuning parameters for random forest
tune_x <- Env.table.train[,names(env)]
tune_y <- Env.table.train[,c("CSE")]
bestmtry <- tuneRF(tune_x, tune_y, stepFactor=1.5, improve=1e-5, ntree=500)
mtry_opt <- bestmtry[,"mtry"][which.min(bestmtry[,"OOBError"])]

#Run random forest 
#Predictor variables are mean along the least cost path lines through the rasters
#Response variable is genetic distance

LeastPath_RF = randomForest(CSE ~ arid + access  +   prec  +   mean.temp  +   human.density  +   friction + min.temp + Needleleaf + EvBroadleaf + DecBroadleaf + MiscTrees + 
Shrubs + Herb + Crop + Flood + Urban + Snow + Barren + Water + Slope + Altitude + PET + DailyTempRange + max.temp + AnnualTempRange + prec.wet + prec.dry + GPP + kernel100, 
importance=TRUE, mtry = mtry_opt, na.action=na.omit, data=Env.table.train)

LeastPath_RF

Cor1 = cor(LeastPath_RF$predicted, Env.table.train$CSE)                         
Cor2 = cor(predict(LeastPath_RF, Env.table.test), Env.table.test$CSE)           
                                                                               
paste ( "ITER" , ITER , "Cor1" , Cor1 )                                                                           
paste ( "ITER" , ITER , "Cor2" , Cor2 )                              

pred.cond <- 1 / predict(env, LeastPath_RF) 

#come back to this
writeRaster(pred.cond, paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/prediction.tif"), options=c("COMPRESS=DEFLATE","ZLEVEL=9") , format="GTiff", overwrite=TRUE  )
EOF

pksetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -m ${OUT_TXT}_${point}/prediction.tif -msknodata -1 -p "<" -nodata -1 -i ${OUT_TXT}_${point}/prediction.tif -o ${OUT_TXT}_${point}/prediction_msk$ITER.tif
rm ${OUT_TXT}_${point}/prediction.tif 
cp ${OUT_TXT}_${point}/prediction_msk$ITER.tif ${OUT_TXT}_${point}/prediction_msk.tif

done   # close the iteration loop  

rm ${OUT_TXT}_${point}/prediction_msk.tif 
rm -rf ${OUT_TXT}_${point}/grassdb$point/loc$point 
cp -r $RAM/grassdb$point/loc$point  ${OUT_TXT}_${point}/grassdb$point/
rm -rf $RAM/grassdb$point/loc$point
