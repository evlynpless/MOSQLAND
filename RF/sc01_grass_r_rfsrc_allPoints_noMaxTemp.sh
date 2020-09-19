#!/bin/bash
#SBATCH -p day
#SBATCH -J sc01_grass_r_rfsrc_allPoints_noMaxTemp.sh
#SBATCH -n 1 -c 16 -N 1
#SBATCH -t 24:00:00 
#SBATCH -o /gpfs/scratch60/fas/powell/esp38/stdout/sc01_grass_r_rfsrc_allPoints_noMaxTemp.sh.%A_%a.out
#SBATCH -e /gpfs/scratch60/fas/powell/esp38/stderr/sc01_grass_r_rfsrc_allPoints_noMaxTemp.sh.%A_%a.err
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
export OUT_TXT=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_all
export point=$SLURM_ARRAY_TASK_ID
###  export point=1

###  spliting in training and testing. Select only one point (all the pairwise from that point)  for the testing

#echo "index,V1,V2,locality1,locality2,lat1,long1,lat2,long2,FST_lin,CSE,Resd" > ${OUT_TXT}/FST_list_NAmRF3_all.csv
#awk -v point=$point  -F ","  '{ if ($2==point || $3==point ) print }' $IN_TXT/FST_list_NAmRF4.csv  >> ${OUT_TXT}_${point}/FST_list_NAmRF3_Test$point.csv
#awk -v point=$point  -F ","  '{ if ($2!=point && $3!=point ) print }' $IN_TXT/FST_list_NAmRF4.csv  > ${OUT_TXT}_${point}/FST_list_NAmRF3_Trai$point.csv

awk -v point=$point  -F ","  '{ if ( 1 > 0 ) print }' $IN_TXT/FST_list_NAmRF4.csv  > ${OUT_TXT}/FST_list_NAmRF3_all.csv 

#### create the start-end points for each line


awk -F "," '{ if(NR>1) printf ("%f %f\n%f %f\n%s\n" , $(NF-5), $(NF-6), $(NF-3), $(NF-4), "NaN NaN" ) }' ${OUT_TXT}/FST_list_NAmRF3_all.csv > ${OUT_TXT}/FST_line_NAmRF3_all.txt

## awk -F "," '{ if(NR>1) printf ("%f %f\n%f %f\n%s\n" , $(NF-5), $(NF-6), $(NF-3), $(NF-4), "NaN NaN" ) }' ${OUT_TXT}_${point}/FST_list_NAmRF3_Trai$point.csv > ${OUT_TXT}_${point}/FST_line_NAmRF3_Trai$point.txt 

### enter in grass and get the mean under straight line

## copy files to the ram to speed-up read and write 

rm -fr $RAM/grassdAll
mkdir  $RAM/grassdbAll

ls $IN_MSQ/consland/ARIDITY/NAm_clip/AI_annual_NAmClip2_Int16.tif $IN_MSQ/consland/access/NAm_clip/accessibility_to_cities_2015_v1.0_NAmClip2_Int16.tif $IN_MSQ/consland/chelsa/bio12/NAm_clip/bio12_mean_NAmClip2_Int16.tif $IN_MSQ/consland/chelsa/bio1/NAm_clip/bio1_NAmClip2_Int16.tif $IN_MSQ/consland/GSHL/NAm_clip/GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84_NAmClip2.tif $IN_MSQ/consland/friction/NAm_clip/friction_surface_2015_v1.0_NAmClip2.tif $IN_MSQ/consland/chelsa/bio6/NAm_clip/bio6_mean_NAmClip2_Int16.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_1_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_2_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_3_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_4_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_5_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_6_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_7_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_8_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_9_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_10_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_11_NAmClip2.tif $IN_MSQ/consland/landcov/NAm_clip/consensus_full_class_12_NAmClip2.tif $IN_MSQ/consland/MERIT/slope/NAm_clip/slope_1KMmedian_MERIT_NAmClip2.tif $IN_MSQ/consland/MERIT/altitude/NAm_clip/altitude_1KMmedian_MERIT_NAmClip2_Int16.tif $IN_MSQ/consland/PET/NAm_clip/pet_mean_NAmClip2.tif $IN_MSQ/consland/chelsa/bio2/NAm_clip/bio2_mean_NAmClip2_Int16.tif $IN_MSQ/consland/chelsa/bio5/NAm_clip/bio5_mean_NAmClip2_Int16.tif $IN_MSQ/consland/chelsa/bio7/NAm_clip/bio7_mean_NAmClip2_Int16.tif $IN_MSQ/consland/chelsa/bio13/NAm_clip/bio13_mean_NAmClip2_Int16.tif $IN_MSQ/consland/chelsa/bio14/NAm_clip/bio14_mean_NAmClip2.tif $IN_MSQ/consland/GPP/mnth/monthly_mean/NAm_clip/GPP_mean_NAmClip2_Int16.tif $IN_MSQ/consland/kernel/KernelRas_100m_fnl.tif | xargs -n 1 -P $CPU bash -c $' cp $1 $RAM/grassdbAll  ' _ 


# rm -r $OUT_TXT/grassdb$point/loc$point
# grass78 -f -text -c $RAM/grassdb$point/AI_annual_NAmClip2_Int16.tif  $OUT_TXT/grassdb$point/loc$point   <<'EOF'

grass78 -f -text -c $RAM/grassdbAll/AI_annual_NAmClip2_Int16.tif  $RAM/grassdbAll/locAll   <<'EOF'

r.external  input=$RAM/grassdbAll/AI_annual_NAmClip2_Int16.tif  output=arid    --overwrite 
r.external  input=$RAM/grassdbAll/accessibility_to_cities_2015_v1.0_NAmClip2_Int16.tif  output=access    --overwrite
r.external  input=$RAM/grassdbAll/bio12_mean_NAmClip2_Int16.tif  output=prec    --overwrite 
r.external  input=$RAM/grassdbAll/bio1_NAmClip2_Int16.tif  output=meantemp    --overwrite
r.external  input=$RAM/grassdbAll/GHS_BUILT_LDS2014_GLOBE_R2016A_54009_1k_v1_0_WGS84_NAmClip2.tif  output=humandensity    --overwrite
r.external  input=$RAM/grassdbAll/friction_surface_2015_v1.0_NAmClip2.tif  output=friction    --overwrite 
r.external  input=$RAM/grassdbAll/bio6_mean_NAmClip2_Int16.tif  output=mintemp    --overwrite
r.external  input=$RAM/grassdbAll/consensus_full_class_1_NAmClip2.tif  output=Needleleaf    --overwrite
r.external  input=$RAM/grassdbAll/consensus_full_class_2_NAmClip2.tif  output=EvBroadleaf    --overwrite 
r.external  input=$RAM/grassdbAll/consensus_full_class_3_NAmClip2.tif  output=DecBroadleaf    --overwrite
r.external  input=$RAM/grassdbAll/consensus_full_class_4_NAmClip2.tif  output=MiscTrees    --overwrite
r.external  input=$RAM/grassdbAll/consensus_full_class_5_NAmClip2.tif  output=Shrubs    --overwrite 
r.external  input=$RAM/grassdbAll/consensus_full_class_6_NAmClip2.tif  output=Herb    --overwrite
r.external  input=$RAM/grassdbAll/consensus_full_class_7_NAmClip2.tif  output=Crop    --overwrite
r.external  input=$RAM/grassdbAll/consensus_full_class_8_NAmClip2.tif  output=Flood    --overwrite  
r.external  input=$RAM/grassdbAll/consensus_full_class_9_NAmClip2.tif  output=Urban    --overwrite
r.external  input=$RAM/grassdbAll/consensus_full_class_10_NAmClip2.tif  output=Snow    --overwrite
r.external  input=$RAM/grassdbAll/consensus_full_class_11_NAmClip2.tif  output=Barren    --overwrite 
r.external  input=$RAM/grassdbAll/consensus_full_class_12_NAmClip2.tif  output=Water    --overwrite
r.external  input=$RAM/grassdbAll/slope_1KMmedian_MERIT_NAmClip2.tif  output=Slope    --overwrite
r.external  input=$RAM/grassdbAll/altitude_1KMmedian_MERIT_NAmClip2_Int16.tif  output=Altitude    --overwrite
r.external  input=$RAM/grassdbAll/pet_mean_NAmClip2.tif  output=PET    --overwrite 
r.external  input=$RAM/grassdbAll/bio2_mean_NAmClip2_Int16.tif  output=DailyTempRange    --overwrite
r.external  input=$RAM/grassdbAll/bio5_mean_NAmClip2_Int16.tif  output=maxtemp    --overwrite
r.external  input=$RAM/grassdbAll/bio7_mean_NAmClip2_Int16.tif  output=AnnualTempRange    --overwrite 
r.external  input=$RAM/grassdbAll/bio13_mean_NAmClip2_Int16.tif  output=precwet    --overwrite
r.external  input=$RAM/grassdbAll/bio14_mean_NAmClip2.tif  output=precdry    --overwrite
r.external  input=$RAM/grassdbAll/GPP_mean_NAmClip2_Int16.tif  output=GPP    --overwrite
r.external  input=$RAM/grassdbAll/KernelRas_100m_fnl.tif  output=kernel100    --overwrite

## add all the predictors 

echo predictors added

v.in.lines input=${OUT_TXT}/FST_line_NAmRF3_all.txt output=FST_line_NAmRF3_all  separator=" " --overwrite
## v.in.lines input=${OUT_TXT}_${point}/FST_line_NAmRF3_Trai$point.txt output=FST_line_NAmRF3_Trai$point  separator=" " --overwrite

echo addtable

v.db.addtable map=FST_line_NAmRF3_all   
## v.db.addtable map=FST_line_NAmRF3_Test$point   

echo  extract mean for straight line in FST_line_NAmRF3_all.txt

echo "index,arid,access,prec,meantemp,humandensity,friction,mintemp,Needleleaf,EvBroadleaf,DecBroadleaf,MiscTrees,Shrubs,Herb,Crop,Flood,Urban,Snow,Barren,Water,Slope,Altitude,PET,DailyTempRange,AnnualTempRange,precwet,precdry,GPP,kernel100" > ${OUT_TXT}/FST_list_NAmRF3_PredictAll.csv 

awk -F "," '{ if (NR>1)  print NR-1 , $1   }'   ${OUT_TXT}/FST_list_NAmRF3_all.csv | xargs -n 2 -P $CPU  bash -c $' 
CAT=$1
INDEX=$2 
v.to.rast  input=FST_line_NAmRF3_all where="cat == $CAT " output=raster$CAT  use="cat" --o   2>/dev/null 

echo $INDEX","$(for rast in arid access prec meantemp humandensity friction mintemp Needleleaf EvBroadleaf DecBroadleaf MiscTrees Shrubs Herb Crop Flood Urban Snow Barren Water Slope Altitude PET DailyTempRange AnnualTempRange precwet precdry GPP ; do  r.univar -t map=$rast    zones=raster$CAT separator=comma 2>/dev/null | awk  -F , \' { if (NR==2)  printf  ("%s," , $8 ) } \'  ; done )$(r.univar -t map=kernel100   zones=raster$CAT  separator=comma 2>/dev/null   | awk  -F ,   \' { if (NR==2)  printf  ("%s\\n",$8 ) } \' )
g.remove -f  type=raster name=raster$CAT  2>/dev/null
' _   >> ${OUT_TXT}/FST_list_NAmRF3_PredictAll.csv 

done

EOF



##### extract kernel values at point level 

#paste -d ","  <(awk -F , '{  if(NR>1) print $2 }' ${OUT_TXT}/FST_list_NAmRF3_all.csv | uniq )   <(gdallocationinfo -geoloc -wgs84  -valonly   $IN_MSQ/consland/kernel/KernelRas_100m_fnl.tif  $(awk -F , '{  if(NR>1) print  $(NF-5), $(NF-6) }'   ${OUT_TXT}/FST_list_NAmRF3_all.csv | uniq )) >  ${OUT_TXT}/FST_list_NAmRF3_KernelAll.csv

#awk -F , '{ if(NR>1) print  $(NF-5), $(NF-6) } END { print $(NF-3), $(NF-4)  }  '   ${OUT_TXT}/FST_list_NAmRF3_all.csv | uniq > ${OUT_TXT}/FST_list_NAmRF3_LatLonAll.txt 

#paste -d ","  <(awk -F , '{ if(NR>1) print $2 }'  ${OUT_TXT}/FST_list_NAmRF3_all.csv | uniq) <(gdallocationinfo -geoloc -wgs84  -valonly   $IN_MSQ/consland/kernel/KernelRas_100m_fnl.tif  < ${OUT_TXT}_${point}/FST_list_NAmRF3_LatLongAll.txt )   >  ${OUT_TXT}/FST_list_NAmRF3_KernelAll.csv


#rm ${OUT_TXT}/FST_list_NAmRF3_LatLonAll.txt  


echo start to process in R 

R --vanilla --no-readline   -q  <<'EOF'

library("randomForest")
library("rgdal")
library("raster")
library("gdistance")
library("randomForestSRC")

point <- Sys.getenv(c('point'))

Env.table <- read.table("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_all/FST_list_NAmRF3_PredictAll.csv", sep=",", header=T)
## select only index and CSE
Gen.table <- read.table("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_all/FST_list_NAmRF3_all.csv", sep=",", header=T)[,c( 1, 11)]

#Env.table.test <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_PredictTest" , point , ".csv"), sep=",", header=T)
#Gen.table.test <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_Test", point         , ".csv"), sep=",", header=T)[,c( 1, 11)]

###  Rename columns 

names(Env.table) <- c("index","arid","access","prec","mean.temp","human.density","friction","min.temp","Needleleaf","EvBroadleaf","DecBroadleaf","MiscTrees","Shrubs","Herb","Crop","Flood","Urban","Snow","Barren","Water","Slope","Altitude","PET","DailyTempRange","AnnualTempRange","prec.wet","prec.dry","GPP","kernel100") 
#names(Env.table.test) <-  c("index","arid","access","prec","mean.temp","human.density","friction","min.temp","Needleleaf","EvBroadleaf","DecBroadleaf","MiscTrees",
#"Shrubs","Herb","Crop","Flood","Urban","Snow","Barren","Water","Slope","Altitude","PET","DailyTempRange","max.temp","AnnualTempRange","prec.wet","prec.dry","GPP","kernel100") 

### Add genetic distance column to the new data frame

Env.table  = merge (Env.table  ,  Gen.table  , by = "index" )
#Env.table.train = merge (Env.table.train ,  Gen.table.train , by = "index" )

head(Env.table)  ### doble check if the merging is done correctly
#head(Env.table.test)   ### doble check if the merging is done correctly

set.seed(NULL)

load(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc17_rasterstack_noMaxTemp.RData")


#Run random forest 
#Predictor variables are mean along straight lines through the rasters
#Response variable is genetic distance

Kernel <- read.table("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/KernelAll.csv",  sep=",", header=F)
names(Kernel) <- c('V1', 'kernel') 
nrow(Kernel)
tail(Kernel)

Pairs <- read.table("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_all/FST_list_NAmRF3_all.csv", sep=",", header=T)
tail(Pairs)

Merge <- merge (Kernel , Pairs  , by = "V1" )
head(Merge)
Merge.sort <- Merge[order(Merge$index),]
tail(Merge.sort)

Kernel.Vector <- as.vector(Merge.sort[,"kernel"])
tail(Kernel.Vector)
str(Kernel.Vector)

Kernel2 <- read.table("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/KernelAll.csv",  sep=",", header=F)
names(Kernel2) <- c('V2', 'kernel')
tail(Kernel2)

Pairs2 <- read.table("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_all/FST_list_NAmRF3_all.csv", sep=",", header=T)
head(Pairs2)

Merge2 <- merge (Kernel2 , Pairs2  , by = "V2" )
Merge2.sort <- Merge2[order(Merge2$index),]
tail(Merge2.sort)
nrow(Merge2)

Kernel.Vector2 <- as.vector(Merge2.sort[,"kernel"])
tail(Kernel.Vector2)

Kernel.Vector.Final <- 1/(pmin(Kernel.Vector, Kernel.Vector2))
head(Kernel.Vector.Final)
tail(Kernel.Vector.Final)
length(Kernel.Vector.Final)

Straight_RF_tune = tune(CSE ~ arid + access  +   prec  +   mean.temp  +   human.density  +   friction + min.temp + Needleleaf + EvBroadleaf + DecBroadleaf + MiscTrees + 
Shrubs + Herb + Crop + Flood + Urban + Snow + Barren + Water + Slope + Altitude + PET + DailyTempRange + AnnualTempRange + prec.wet + prec.dry + GPP + kernel100, 
importance=TRUE, na.action=c("na.omit"), case.wt=Kernel.Vector.Final,  data=Env.table)

Straight_RF_tune$optimal[["mtry"]]
Straight_RF_tune$optimal[["nodesize"]]                  

Straight_RF = rfsrc(CSE ~ arid + access  +   prec  +   mean.temp  +   human.density  +   friction + min.temp + Needleleaf + EvBroadleaf + DecBroadleaf + MiscTrees + 
Shrubs + Herb + Crop + Flood + Urban + Snow + Barren + Water + Slope + Altitude + PET + DailyTempRange +  AnnualTempRange + prec.wet + prec.dry + GPP + kernel100, 
importance=TRUE, na.action=c("na.omit"),  case.wt=Kernel.Vector.Final,  mtry = Straight_RF_tune$optimal[["mtry"]], nodesize =  Straight_RF_tune$optimal[["nodesize"]], data=Env.table)

Straight_RF

pdf("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_all/ErrVIMP_iter0.pdf", 7, 7)
plot(Straight_RF, m.target = NULL, plots.one.page = TRUE, sorted = TRUE, verbose = TRUE)
dev.off()

R = cor(Straight_RF$predicted.oob, Env.table$CSE)
paste ("ITER 0 RVar" , R )  

RMSE = sqrt(mean((Straight_RF$predicted.oob - Env.table$CSE)^2))
paste ("ITER 0 RMSEVar" , RMSE) 

MAE =  mean(abs(Straight_RF$predicted.oob - Env.table$CSE))
paste ("ITER 0 MAEVar" , MAE) 

pred = predict.rfsrc(Straight_RF, value.raster, na.action = c("na.impute"))

predict.rast=raster(vals=as.vector(pred$predicted),  nrows= 1500 , ncols=4140 , xmn=-113.5, xmx=-79, ymn=24, ymx=36.5)

predict.rast.mask <- mask(predict.rast, arid)

pred.cond <- 1/predict.rast.mask

pdf("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_all/Scatter_iter0.pdf", 5, 5)
plot(Env.table$CSE, Straight_RF$predicted.oob, xlab ="Observed CSE", ylab="Predicted CSE")
abline(a=0, b=1)
abline(lm(Straight_RF$predicted.oob ~ Env.table$FST_lin), col="red")
legend("bottomright", legend=c(paste0("Pearson correlation = ", round(R,3))), cex=0.7)
dev.off()

writeRaster(pred.cond, "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_all/prediction.tif", options=c("COMPRESS=DEFLATE","ZLEVEL=9") , format="GTiff", overwrite=TRUE  )

EOF

#Giuseppe suggestion
#pksetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -m ${OUT_TXT}/prediction.tif -msknodata -1 -p "<" -nodata -1  -m  path/arid.tif  -msknodata -9999   -nodata -1  -i ${OUT_TXT}/prediction.tif -o ${OUT_TXT}/prediction_msk0.tif

pksetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -m ${OUT_TXT}/prediction.tif -msknodata -1 -p "<" -nodata -1 -i ${OUT_TXT}/prediction.tif -o ${OUT_TXT}/prediction_msk0.tif

rm ${OUT_TXT}/prediction.tif 
cp ${OUT_TXT}/prediction_msk0.tif ${OUT_TXT}/prediction_msk.tif


### Test and Train  start and stop coordinates 
awk -F "," '{ if(NR>1) print $1 , $(NF-5),  $(NF-6) ,  $(NF-3),  $(NF-4) }' ${OUT_TXT}/FST_list_NAmRF3_all.csv | uniq > ${OUT_TXT}/FST_line_NAmRF3_StartStopAll.txt
#awk -F "," '{ if(NR>1) print $1 , $(NF-5),  $(NF-6) ,  $(NF-3),  $(NF-4) }' ${OUT_TXT}_${point}/FST_list_NAmRF3_Trai$point.csv | uniq > ${OUT_TXT}_${point}/FST_line_NAmRF3_StartStopTrai$point.txt

for ITER in $(seq 1 10 ) ; do  
export ITER=$ITER

echo Start Iteration $ITER  GRASS 

grass78 -f -text $RAM/grassdbAll/locAll/PERMANENT   <<'EOF'
r.external  input=${OUT_TXT}/prediction_msk.tif    output=prediction   --overwrite

echo "index,arid,access,prec,meantemp,humandensity,friction,mintemp,Needleleaf,EvBroadleaf,DecBroadleaf,MiscTrees,Shrubs,Herb,Crop,Flood,Urban,Snow,Barren,Water,Slope,Altitude,PET,DailyTempRange,maxtemp,AnnualTempRange,precwet,precdry,GPP,kernel100" > ${OUT_TXT}/FST_list_NAmRF3_Iter${ITER}LeastPathAll.csv

cat ${OUT_TXT}/FST_line_NAmRF3_StartStopAll.txt   | xargs -n 5 -P $CPU  bash -c $'
INDEX=$1 

r.cost -n -k  input=prediction  output=cost$INDEX  outdir=dir$INDEX    start_coordinates=$2,$3  memory=200  --o --q
g.remove -f  type=raster name=cost$INDEX  2>/dev/null
r.path input=dir$INDEX  raster_path=least_path$INDEX start_coordinates=$4,$5 --o --q 
g.remove -f  type=raster name=dir$INDEX  2>/dev/null

echo $INDEX","$(for rast in arid access prec meantemp humandensity friction mintemp Needleleaf EvBroadleaf DecBroadleaf MiscTrees Shrubs Herb Crop Flood Urban Snow Barren Water Slope Altitude PET DailyTempRange AnnualTempRange precwet precdry GPP ; do  r.univar -t map=$rast    zones=least_path$INDEX separator=comma 2>/dev/null | awk  -F , \' { if (NR==2)  printf  ("%s,", $8) } \' ; done )$(r.univar -t map=kernel100 zones=least_path$INDEX separator=comma 2>/dev/null | awk  -F , \' { if (NR==2)  printf  ("%s\\n",$8 ) } \') 
g.remove -f  type=raster name=least_path$INDEX  2>/dev/null

' _  >> ${OUT_TXT}/FST_list_NAmRF3_Iter${ITER}LeastPathAll.csv
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

Env.table <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_all/FST_list_NAmRF3_Iter",ITER,"LeastPathAll.csv"), sep=",", header=T)
head(Env.table)
Gen.table <- read.table("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_all/FST_list_NAmRF3_all.csv", sep=",", header=T)[,c( 1, 11)]
head(Gen.table)

#Env.table.test <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_Iter",ITER,"LeastPathTest",point,".csv"), sep=",", header=T) 
#Gen.table.test <- read.table(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_", point, "/FST_list_NAmRF3_Test",point, ".csv"), sep=",", header=T)[,c( 1, 11)]

### Add genetic distance column to the new data frame

###  Rename columns 

names(Env.table) <- c("index","arid","access","prec","mean.temp","human.density","friction","min.temp","Needleleaf","EvBroadleaf","DecBroadleaf","MiscTrees", "Shrubs","Herb","Crop","Flood","Urban","Snow","Barren","Water","Slope","Altitude","PET","DailyTempRange","AnnualTempRange","prec.wet","prec.dry","GPP","kernel100") 
#names(Env.table.test) <-  c("index","arid","access","prec","mean.temp","human.density","friction","min.temp","Needleleaf","EvBroadleaf","DecBroadleaf","MiscTrees",
#"Shrubs","Herb","Crop","Flood","Urban","Snow","Barren","Water","Slope","Altitude","PET","DailyTempRange","max.temp","AnnualTempRange","prec.wet","prec.dry","GPP","kernel100") 

### Add genetic distance column to the new data frame

Env.table =  merge (Env.table , Gen.table  , by = "index" )
#Env.table.train = merge (Env.table.train, Gen.table.train , by = "index" )

#head(Env.table.train)  ### doble check if the merging is done correctly
#head(Env.table.test)  ### doble check if the merging is done correctly

set.seed(NULL)

load(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc17_rasterstack_noMaxTemp.RData")

#Run random forest 
#Predictor variables are mean along the least cost path lines through the rasters
#Response variable is genetic distance

#Add kernel stuff here

Kernel <- read.table("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/KernelAll.csv",  sep=",", header=F)
names(Kernel) <- c('V1', 'kernel')
nrow(Kernel)
tail(Kernel)

Pairs <- read.table("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_all/FST_list_NAmRF3_all.csv", sep=",", header=T)
tail(Pairs)

Merge <- merge (Kernel , Pairs  , by = "V1" )
head(Merge)                                                                                                                                                         

Merge.sort <- Merge[order(Merge$index),]
tail(Merge.sort)

Kernel.Vector <- as.vector(Merge.sort[,"kernel"])
tail(Kernel.Vector)
str(Kernel.Vector)

Kernel2 <- read.table("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/KernelAll.csv",  sep=",", header=F)
names(Kernel2) <- c('V2', 'kernel')
tail(Kernel2)

Pairs2 <- read.table("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_all/FST_list_NAmRF3_all.csv", sep=",", header=T)
head(Pairs2)

Merge2 <- merge (Kernel2 , Pairs2  , by = "V2" )
Merge2.sort <- Merge2[order(Merge2$index),]
tail(Merge2.sort)
nrow(Merge2)

Kernel.Vector2 <- as.vector(Merge2.sort[,"kernel"])
tail(Kernel.Vector2)

Kernel.Vector.Final <- 1/(pmin(Kernel.Vector, Kernel.Vector2))
head(Kernel.Vector.Final)
tail(Kernel.Vector.Final)
length(Kernel.Vector.Final)


LeastPath_RF_tune = tune(CSE ~ arid + access  +   prec  +   mean.temp  +   human.density  +   friction + min.temp + Needleleaf + EvBroadleaf + DecBroadleaf + MiscTrees + 
Shrubs + Herb + Crop + Flood + Urban + Snow + Barren + Water + Slope + Altitude + PET + DailyTempRange + AnnualTempRange + prec.wet + prec.dry + GPP + kernel100, 
importance=TRUE, na.action=c("na.omit"), case.wt=Kernel.Vector.Final, data=Env.table)

LeastPath_RF_tune$optimal[["mtry"]]
LeastPath_RF_tune$optimal[["nodesize"]]

LeastPath_RF = rfsrc(CSE ~ arid + access  +   prec  +   mean.temp  +   human.density  +   friction + min.temp + Needleleaf + EvBroadleaf + DecBroadleaf + MiscTrees + 
Shrubs + Herb + Crop + Flood + Urban + Snow + Barren + Water + Slope + Altitude + PET + DailyTempRange + AnnualTempRange + prec.wet + prec.dry + GPP + kernel100, 
importance=TRUE, na.action=c("na.omit"), mtry = LeastPath_RF_tune$optimal[["mtry"]], nodesize = LeastPath_RF_tune$optimal[["nodesize"]], case.wt=Kernel.Vector.Final, data=Env.table)

LeastPath_RF

pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_all/ErrVIMP_iter", ITER, ".pdf"), 7, 7)
plot(LeastPath_RF, m.target = NULL, plots.one.page = TRUE, sorted = TRUE, verbose = TRUE)
dev.off()  

R  = cor(LeastPath_RF$predicted.oob, Env.table$CSE)                         
paste ( "ITER" , ITER , "RVar" , R )                              

RMSE = sqrt(mean((LeastPath_RF$predicted.oob - Env.table$CSE)^2))
paste ( "ITER" , ITER ,  "RMSEVar" , RMSE )

MAE =  mean(abs(LeastPath_RF$predicted.oob - Env.table$CSE))
paste ("ITER" , ITER , "MAEVar" , MAE)

pred = predict.rfsrc(LeastPath_RF, value.raster, na.action = c("na.impute"))

predict.rast=raster(vals=as.vector(pred$predicted),  nrows= 1500 , ncols=4140 , xmn=-113.5, xmx=-79, ymn=24, ymx=36.5)                                                                                          

predict.rast.mask <- mask(predict.rast, arid)     

pred.cond <- 1/predict.rast.mask

#pred.cond <- 1/raster(vals=as.vector(pred$predicted),  nrows= 1500 , ncols=4140 , xmn=-113.5, xmx=-79, ymn=24, ymx=36.5) 

pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_all/Scatter_iter", ITER, ".pdf"), 5, 5)
plot(Env.table$CSE, LeastPath_RF$predicted.oob,  xlab ="Observed CSE", ylab="Predicted CSE")
legend("bottomright", legend=c(paste0("Pearson correlation = ", round(R,3))), cex=0.7)
dev.off()

pdf(paste0("/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_all/Scatter2_iter", ITER, ".pdf"), 5, 5)
plot(Env.table$CSE, LeastPath_RF$predicted.oob, xlab ="Observed CSE", ylab="Predicted CSE")
abline(a=0, b=1)
abline(lm(LeastPath_RF$predicted.oob ~ Env.table$CSE), col="red")
legend("bottomright", legend=c(paste0("Pearson correlation = ", round(R,3))), cex=0.7)
dev.off()  

writeRaster(pred.cond, "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/TrainingTestingRfsrc_all/prediction.tif", options=c("COMPRESS=DEFLATE","ZLEVEL=9") , format="GTiff", overwrite=TRUE  )
EOF

#pksetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -m ${OUT_TXT}_${point}/prediction.tif -msknodata -1 -p "<" -nodata -1  -m  path/arid.tif  -msknodata -9999   -nodata -1  -i ${OUT_TXT}_${point}/prediction.tif -o ${OUT_TXT}_${point}/prediction_msk0.tif

pksetmask -co COMPRESS=DEFLATE -co ZLEVEL=9 -m ${OUT_TXT}/prediction.tif -msknodata -1 -p "<" -nodata -1 -i ${OUT_TXT}/prediction.tif -o ${OUT_TXT}/prediction_msk$ITER.tif

rm ${OUT_TXT}/prediction.tif 
cp ${OUT_TXT}/prediction_msk$ITER.tif ${OUT_TXT}/prediction_msk.tif

done   # close the iteration loop  

rm ${OUT_TXT}/prediction_msk.tif 
rm -rf ${OUT_TXT}/grassdbAll/locAll 
cp -r $RAM/grassdb/locAll  ${OUT_TXT}/grassdbAll/ #how to edit these last lines?
rm -rf $RAM/grassdb/locAll
