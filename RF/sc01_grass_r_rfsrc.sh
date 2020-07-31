#!/bin/bash
#SBATCH -p day
#SBATCH -J sc01_grass_r_rfsrc.sh
#SBATCH -n 1 -c 16 -N 1
#SBATCH -t 24:00:00 
#SBATCH -o /gpfs/scratch60/fas/powell/esp38/stdout/sc01_grass_r_rfsrc.sh.%A_%a.out
#SBATCH -e /gpfs/scratch60/fas/powell/esp38/stderr/sc01_grass_r_rfsrc.sh.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evlyn.pless@yale.edu
#SBATCH --array=1
#SBATCH --mem=80G

####### sbatch  /home/fas/powell/esp38/scripts/MOSQLAND/RF/sc01_grass_r_rfsrc.sh
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

echo "index,V1,V2,locality1,locality2,lat1,long1,lat2,long2,FST_lin,CSE,Resd" > ${OUT_TXT}_${point}/FST_list_NAmRF3_Test$point.csv
awk -v point=$point  -F ","  '{ if ($2==point || $3==point ) print }' $IN_TXT/FST_list_NAmRF4.csv  >> ${OUT_TXT}_${point}/FST_list_NAmRF3_Test$point.csv
awk -v point=$point  -F ","  '{ if ($2!=point && $3!=point ) print }' $IN_TXT/FST_list_NAmRF4.csv  > ${OUT_TXT}_${point}/FST_list_NAmRF3_Trai$point.csv

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

awk -F , '{ if(NR>1) print  $(NF-5), $(NF-6) } END { print $(NF-3), $(NF-4)  }  '   ${OUT_TXT}_${point}/FST_list_NAmRF3_Trai$point.csv | uniq > ${OUT_TXT}_${point}/FST_list_NAmRF3_LatLongTrai$point.txt 

#paste -d ","  <(awk -F , '{ if(NR>1) print $2 }'  ${OUT_TXT}_${point}/FST_list_NAmRF3_Trai$point.csv | uniq) <(gdallocationinfo -geoloc -wgs84  -valonly   $IN_MSQ/consland/kernel/KernelRas_100m_fnl.tif  < ${OUT_TXT}_${point}/FST_list_NAmRF3_LatLongTrai$point.txt )   >  ${OUT_TXT}_${point}/FST_list_NAmRF3_KernelTrai$point.csv

if [  $point -le 37    ] ; then  

paste -d "," <(awk -F , '{ if(NR>1) print $2 } END {print 38 } ' ${OUT_TXT}_${point}/FST_list_NAmRF3_Trai$point.csv | uniq) <(gdallocationinfo -geoloc -wgs84 -valonly $IN_MSQ/consland/kernel/KernelRas_100m_fnl.tif < ${OUT_TXT}_${point}/FST_list_NAmRF3_LatLongTrai$point.txt ) > ${OUT_TXT}_${point}/FST_list_NAmRF3_KernelTrai$point.csv

else

paste -d "," <(awk -F , '{ if(NR>1) print $2 } END {print 37 } ' ${OUT_TXT}_${point}/FST_list_NAmRF3_Trai$point.csv | uniq) <(gdallocationinfo -geoloc -wgs84 -valonly $IN_MSQ/consland/kernel/KernelRas_100m_fnl.tif < ${OUT_TXT}_${point}/FST_list_NAmRF3_LatLongTrai$point.txt ) > ${OUT_TXT}_${point}/FST_list_NAmRF3_KernelTrai$point.csv

paste  -d "," <(awk -F , '{ if(NR>1) print $2 } END {print 38 } '  ${OUT_TXT}_1/FST_list_NAmRF3_Trai1.csv | uniq  )  <(awk -F , '{ if(NR>1) print $2 } END {print 37 } '  ${OUT_TXT}_38/FST_list_NAmRF3_Trai38.csv | uniq   )


rm ${OUT_TXT}_${point}/FST_list_NAmRF3_LatLongTrai$point.txt  


echo start to process in R 

R --vanilla --no-readline   -q  <<'EOF'

library("randomForest")
library("rgdal")
library("raster")
library("gdistance")
library("randomForestSRC")

point <- Sys.getenv(c('point'))

Env.table.train <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_PredictTrai" , point , ".csv"), sep=",", header=T)
## select only index and CSE
Gen.table.train <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_Trai", point         , ".csv"), sep=",", header=T)[,c( 1, 11)]

Env.table.test <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_PredictTest" , point , ".csv"), sep=",", header=T)
Gen.table.test <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_Test", point         , ".csv"), sep=",", header=T)[,c( 1, 11)]

###  Rename columns 

names(Env.table.train) <- c("index","arid","access","prec","mean.temp","human.density","friction","min.temp","Needleleaf","EvBroadleaf","DecBroadleaf","MiscTrees",
"Shrubs","Herb","Crop","Flood","Urban","Snow","Barren","Water","Slope","Altitude","PET","DailyTempRange","max.temp","AnnualTempRange","prec.wet","prec.dry","GPP","kernel100") 
names(Env.table.test) <-  c("index","arid","access","prec","mean.temp","human.density","friction","min.temp","Needleleaf","EvBroadleaf","DecBroadleaf","MiscTrees",
"Shrubs","Herb","Crop","Flood","Urban","Snow","Barren","Water","Slope","Altitude","PET","DailyTempRange","max.temp","AnnualTempRange","prec.wet","prec.dry","GPP","kernel100") 

### Add genetic distance column to the new data frame

Env.table.test  = merge (Env.table.test  ,  Gen.table.test  , by = "index" )
Env.table.train = merge (Env.table.train ,  Gen.table.train , by = "index" )

#head(Env.table.train)  ### doble check if the merging is done correctly
#head(Env.table.test)   ### doble check if the merging is done correctly

set.seed(NULL)

load(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc16_rasterstack.RData")


#Run random forest 
#Predictor variables are mean along straight lines through the rasters
#Response variable is genetic distance

Kernel.train <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_KernelTrai" , point , ".csv"), sep=",", header=F) 
names(Kernel.train) <- c('V1', 'kernel') 
#if (point < 38) {
#Kernel.train[37,1] = 38
#} else { Kernel.train[37,1] = 37
#}
#Kernel.train[37,1] = 37 
nrow(Kernel.train)
tail(Kernel.train)

Pairs.train <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_Trai" , point , ".csv"), sep=",", header=T)
tail(Pairs.train)

Merge.train <- merge (Kernel.train , Pairs.train  , by = "V1" )
Merge.train.sort <- Merge.train[order(Merge.train$index),]
tail(Merge.train.sort)
nrow(Merge.train)

Kernel.Vector <- as.vector(Merge.train.sort[,"kernel"])
tail(Kernel.Vector)
str(Kernel.Vector)

Kernel.train2 <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_KernelTrai" , point , ".csv"), sep=",", header=F)
names(Kernel.train2) <- c('V2', 'kernel')
#if (point < 38) {
#Kernel.train2[37,1] = 38
#} else { Kernel.train2[37,1] = 37
#}
#Kernel.train2[37,1] = 37 
tail(Kernel.train2)


Pairs.train2 <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_Trai" , point , ".csv"), sep=",", header=T)
head(Pairs.train2)

Merge.train2 <- merge (Kernel.train2 , Pairs.train2  , by = "V2" )
Merge.train2.sort <- Merge.train2[order(Merge.train2$index),]
tail(Merge.train2.sort)
nrow(Merge.train2)

Kernel.Vector2 <- as.vector(Merge.train2.sort[,"kernel"])
tail(Kernel.Vector2)

Kernel.Vector.Final <- 1/(pmin(Kernel.Vector, Kernel.Vector2))
head(Kernel.Vector.Final)
tail(Kernel.Vector.Final)
length(Kernel.Vector.Final)

Straight_RF_tune = tune(CSE ~ arid + access  +   prec  +   mean.temp  +   human.density  +   friction + min.temp + Needleleaf + EvBroadleaf + DecBroadleaf + MiscTrees + 
Shrubs + Herb + Crop + Flood + Urban + Snow + Barren + Water + Slope + Altitude + PET + DailyTempRange + max.temp + AnnualTempRange + prec.wet + prec.dry + GPP + kernel100, 
importance=TRUE, na.action=c("na.omit"), case.wt=Kernel.Vector.Final, data=Env.table.train)

Straight_RF_tune$optimal[["mtry"]]
Straight_RF_tune$optimal[["nodesize"]]                  

Straight_RF = rfsrc(CSE ~ arid + access  +   prec  +   mean.temp  +   human.density  +   friction + min.temp + Needleleaf + EvBroadleaf + DecBroadleaf + MiscTrees + 
Shrubs + Herb + Crop + Flood + Urban + Snow + Barren + Water + Slope + Altitude + PET + DailyTempRange + max.temp + AnnualTempRange + prec.wet + prec.dry + GPP + kernel100, 
importance=TRUE, na.action=c("na.omit"), case.wt=Kernel.Vector.Final, mtry = Straight_RF_tune$optimal[["mtry"]], nodesize =  Straight_RF_tune$optimal[["nodesize"]], data=Env.table.train)

Straight_RF

pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/ErrVIMP_iter0.pdf"), 7, 7)
plot(Straight_RF, m.target = NULL, plots.one.page = TRUE, sorted = TRUE, verbose = TRUE)
dev.off()

Rtrain = cor(Straight_RF$predicted, Env.table.train$CSE)
paste (" ITER 0 RtrainVar " , Rtrain ) 

#Turns out this is same as Rtrain
#Rtrain2 = cor((predict.rfsrc(Straight_RF, Env.table.train))$predicted, Env.table.train$CSE)
#paste ("ITER 0 Rtrain2Var " , Rtrain2 )

Rtest = cor((predict.rfsrc(Straight_RF, Env.table.test))$predicted, Env.table.test$CSE)
paste (" ITER 0 RtestVar " , Rtest ) 

RMSEtrain = sqrt(mean((Straight_RF$predicted - Env.table.train$CSE)^2))
paste (" ITER 0 RMSEtrainVar " , RMSEtrain)

RMSEtest = sqrt(mean((predict.rfsrc(Straight_RF, Env.table.test)$predicted - Env.table.test$CSE)^2))
paste (" ITER 0 RMSEtestVar " , RMSEtest) 

MAEtrain =  mean(abs(predict.rfsrc(Straight_RF, Env.table.train)$predicted - Env.table.train$CSE))
paste (" ITER 0 MAEtrainVar " , MAEtrain) 

MAEtest = mean(abs(predict.rfsrc(Straight_RF, Env.table.test)$predicted - Env.table.test$CSE))
paste (" ITER 0  MAEtestVar " , MAEtest) 

pred = predict.rfsrc(Straight_RF, value.raster, na.action = c("na.impute"))

predict.rast=raster(vals=as.vector(pred$predicted),  nrows= 1500 , ncols=4140 , xmn=-113.5, xmx=-79, ymn=24, ymx=36.5)

predict.rast.mask <- mask(predict.rast, arid)

pred.cond <- 1/predict.rast

pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/TrainingScatter_iter0.pdf"), 5, 5)
plot(Env.table.train$CSE, Straight_RF$predicted,  xlab ="Observed CSE (training)", ylab="Predicted CSE")
legend("bottomright", legend=c(paste0("Pearson correlation = ", round(Rtrain,3))), cex=0.7)
dev.off()

#pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/Straight_TestingScatter2_pt", point, ".pdf"), 5, 5)
#plot(Env.table.train$CSE, predict.rfsrc(Straight_RF, Env.table.train)$predicted,  xlab ="Observed CSE (training)", ylab="Predicted CSE")
#legend("bottomright", legend=c(paste0("Pearson correlation = ", round(Rtrain2,3))), cex=0.7)
#dev.off()

pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/TestingScatter_iter0.pdf"), 5, 5)
plot(Env.table.test$CSE, predict.rfsrc(Straight_RF, Env.table.test)$predicted,  xlab ="Observed CSE (training)", ylab="Predicted CSE")
legend("bottomright", legend=c(paste0("Pearson correlation = ", round(Rtest,3))), cex=0.7)
dev.off() 

writeRaster(pred.cond, paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/prediction.tif"), options=c("COMPRESS=DEFLATE","ZLEVEL=9") , format="GTiff", overwrite=TRUE  )

EOF

#Giuseppe suggestion
#pksetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -m ${OUT_TXT}_${point}/prediction.tif -msknodata -1 -p "<" -nodata -1  -m  path/arid.tif  -msknodata -9999   -nodata -1  -i ${OUT_TXT}_${point}/prediction.tif -o ${OUT_TXT}_${point}/prediction_msk0.tif

pksetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -m ${OUT_TXT}_${point}/prediction.tif -msknodata -1 -p "<" -nodata -1 -i ${OUT_TXT}_${point}/prediction.tif -o ${OUT_TXT}_${point}/prediction_msk0.tif

rm ${OUT_TXT}_${point}/prediction.tif 
cp ${OUT_TXT}_${point}/prediction_msk0.tif ${OUT_TXT}_${point}/prediction_msk.tif


### Test and Train  start and stop coordinates 
awk -F "," '{ if(NR>1) print $1 , $(NF-5),  $(NF-6) ,  $(NF-3),  $(NF-4) }' ${OUT_TXT}_${point}/FST_list_NAmRF3_Test$point.csv | uniq > ${OUT_TXT}_${point}/FST_line_NAmRF3_StartStopTest$point.txt
awk -F "," '{ if(NR>1) print $1 , $(NF-5),  $(NF-6) ,  $(NF-3),  $(NF-4) }' ${OUT_TXT}_${point}/FST_list_NAmRF3_Trai$point.csv | uniq > ${OUT_TXT}_${point}/FST_line_NAmRF3_StartStopTrai$point.txt

for ITER in $(seq 1 2 ) ; do  
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
library("randomForestSRC")

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

#head(Env.table.train)  ### doble check if the merging is done correctly
#head(Env.table.test)  ### doble check if the merging is done correctly

set.seed(NULL)

load(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc16_rasterstack.RData")

#Add tuning here?

#Run random forest 
#Predictor variables are mean along the least cost path lines through the rasters
#Response variable is genetic distance

Kernel.train <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_KernelTrai" , point , ".csv"), sep=",", header=F)
names(Kernel.train) <- c('V1', 'kernel')
#if (point < 38) {
 # Kernel.train[37,1] = 38
#} else { Kernel.train[37,1] = 37
#}                      
#Kernel.train[37,1] = 37                                                                                                                                           
tail(Kernel.train)                                                                                                                                                                                                      
Pairs.train <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_Trai" , point , ".csv"), sep=",", header=T)
tail(Pairs.train)

Merge.train <- merge (Kernel.train , Pairs.train  , by = "V1" )
Merge.train.sort <- Merge.train[order(Merge.train$index),]
tail(Merge.train.sort)

Kernel.Vector <- as.vector(Merge.train.sort[,"kernel"])
tail(Kernel.Vector)

Kernel.train2 <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_KernelTrai" , point , ".csv"), sep=",", header=F)
names(Kernel.train2) <- c('V2', 'kernel')
#if (point < 38) {
 # Kernel.train2[37,1] = 38
#} else { Kernel.train2[37,1] = 37
#}
#Kernel.train2[37,1] = 37
tail(Kernel.train2)

Pairs.train2 <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_Trai" , point , ".csv"), sep=",", header=T)
head(Pairs.train2)

Merge.train2 <- merge (Kernel.train2 , Pairs.train2  , by = "V2" )
Merge.train2.sort <- Merge.train2[order(Merge.train2$index),]
tail(Merge.train2.sort)

Kernel.Vector2 <- as.vector(Merge.train2.sort[,"kernel"])

Kernel.Vector.Final <- 1/(pmin(Kernel.Vector + Kernel.Vector2))

LeastPath_RF_tune = tune(CSE ~ arid + access  +   prec  +   mean.temp  +   human.density  +   friction + min.temp + Needleleaf + EvBroadleaf + DecBroadleaf + MiscTrees + 
Shrubs + Herb + Crop + Flood + Urban + Snow + Barren + Water + Slope + Altitude + PET + DailyTempRange + max.temp + AnnualTempRange + prec.wet + prec.dry + GPP + kernel100, 
importance=TRUE, na.action=c("na.omit"), data=Env.table.train)

LeastPath_RF_tune$optimal[["mtry"]]
LeastPath_RF_tune$optimal[["nodesize"]]

LeastPath_RF = rfsrc(CSE ~ arid + access  +   prec  +   mean.temp  +   human.density  +   friction + min.temp + Needleleaf + EvBroadleaf + DecBroadleaf + MiscTrees + 
Shrubs + Herb + Crop + Flood + Urban + Snow + Barren + Water + Slope + Altitude + PET + DailyTempRange + max.temp + AnnualTempRange + prec.wet + prec.dry + GPP + kernel100, 
importance=TRUE, na.action=c("na.omit"), case.wt=Kernel.Vector.Final, mtry = LeastPath_RF_tune$optimal[["mtry"]], nodesize = LeastPath_RF_tune$optimal[["nodesize"]], data=Env.table.train)

LeastPath_RF

pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/ErrVIMP_iter", ITER, ".pdf"), 7, 7)
plot(LeastPath_RF, m.target = NULL, plots.one.page = TRUE, sorted = TRUE, verbose = TRUE)
dev.off()  

Rtrain  = cor(LeastPath_RF$predicted, Env.table.train$CSE)                         
paste ( "ITER" , ITER , "RtrainVar" , Rtrain )

#Rtrain2 = cor((predict.rfsrc(LeastPath_RF, Env.table.train))$predicted, Env.table.train$CSE)
#paste ( "ITER" , ITER , "RtrainVar2" , Rtrain2 )  

Rtest  = cor((predict.rfsrc(LeastPath_RF, Env.table.test))$predicted, Env.table.test$CSE)
paste ( "ITER" , ITER , "RtestVar" , Rtest )                              

RMSEtrain = sqrt(mean((LeastPath_RF$predicted - Env.table.train$CSE)^2))
paste ( "ITER" , ITER ,  "RMSEtrainVar" , RMSEtrain )

RMSEtest = sqrt(mean((predict.rfsrc(LeastPath_RF, Env.table.test)$predicted - Env.table.test$CSE)^2))
paste ( "ITER" , ITER , "RMSEtestVar" , RMSEtest)

MAEtrain =  mean(abs(predict.rfsrc(LeastPath_RF, Env.table.train)$predicted - Env.table.train$CSE))
paste ("ITER" , ITER , "MAEtrainVar" , MAEtrain)

MAEtest = mean(abs(predict.rfsrc(LeastPath_RF, Env.table.test)$predicted - Env.table.test$CSE))
paste ("ITER" , ITER , "MAEtestVar" , MAEtest) 

pred = predict.rfsrc(LeastPath_RF, value.raster, na.action = c("na.impute"))

predict.rast=raster(vals=as.vector(pred$predicted),  nrows= 1500 , ncols=4140 , xmn=-113.5, xmx=-79, ymn=24, ymx=36.5)                                                                                          

predict.rast.mask <- mask(predict.rast, arid)     

pred.cond <- 1/predict.rast

#pred.cond <- 1/raster(vals=as.vector(pred$predicted),  nrows= 1500 , ncols=4140 , xmn=-113.5, xmx=-79, ymn=24, ymx=36.5) 

pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/TrainingScatter_iter", ITER, ".pdf"), 5, 5)
plot(Env.table.train$CSE, LeastPath_RF$predicted,  xlab ="Observed CSE (training)", ylab="Predicted CSE")
legend("bottomright", legend=c(paste0("Pearson correlation = ", round(Rtrain,3))), cex=0.7)
dev.off()

#pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/TrainingScatter2_iter", ITER, ".pdf"), 5, 5)
#plot(Env.table.train$CSE, predict.rfsrc(LeastPath_RF, Env.table.train)$predicted,  xlab ="Observed CSE (training)", ylab="Predicted CSE")
#legend("bottomright", legend=c(paste0("Pearson correlation = ", round(Rtrain2,3))), cex=0.7)
#dev.off()

pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/TestingScatter_iter", ITER, ".pdf"), 5, 5)
plot(Env.table.test$CSE, predict.rfsrc(LeastPath_RF, Env.table.test)$predicted,  xlab ="Observed CSE (training)", ylab="Predicted CSE")
legend("bottomright", legend=c(paste0("Pearson correlation = ", round(Rtest,3))), cex=0.7)
dev.off() 

writeRaster(pred.cond, paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/prediction.tif"), options=c("COMPRESS=DEFLATE","ZLEVEL=9") , format="GTiff", overwrite=TRUE  )
EOF

#pksetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -m ${OUT_TXT}_${point}/prediction.tif -msknodata -1 -p "<" -nodata -1  -m  path/arid.tif  -msknodata -9999   -nodata -1  -i ${OUT_TXT}_${point}/prediction.tif -o ${OUT_TXT}_${point}/prediction_msk0.tif

pksetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -m ${OUT_TXT}_${point}/prediction.tif -msknodata -1 -p "<" -nodata -1 -i ${OUT_TXT}_${point}/prediction.tif -o ${OUT_TXT}_${point}/prediction_msk$ITER.tif

rm ${OUT_TXT}_${point}/prediction.tif 
cp ${OUT_TXT}_${point}/prediction_msk$ITER.tif ${OUT_TXT}_${point}/prediction_msk.tif

done   # close the iteration loop  

rm ${OUT_TXT}_${point}/prediction_msk.tif 
rm -rf ${OUT_TXT}_${point}/grassdb$point/loc$point 
cp -r $RAM/grassdb$point/loc$point  ${OUT_TXT}_${point}/grassdb$point/
rm -rf $RAM/grassdb$point/loc$point
