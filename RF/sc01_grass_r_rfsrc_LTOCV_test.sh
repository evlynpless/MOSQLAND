#!/bin/bash
#SBATCH -p day
#SBATCH -J sc01_grass_r_rfsrc_LTOCV.sh
#SBATCH -n 1 -c 16 -N 1
#SBATCH -t 24:00:00 
#SBATCH -o /gpfs/scratch60/fas/powell/esp38/stdout/sc01_grass_r_rfsrc_LTOCV.sh.%A_%a.out
#SBATCH -e /gpfs/scratch60/fas/powell/esp38/stderr/sc01_grass_r_rfsrc_LTOCV.sh.%A_%a.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evlyn.pless@yale.edu
#SBATCH --array=1-6:2
#SBATCH --mem=8G

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

###  spliting in training and testing. Select two points (and all the pairwise from those points)  for the testing

echo "index,V1,V2,locality1,locality2,lat1,long1,lat2,long2,FST_lin,CSE,Resd" > ${OUT_TXT}_${point}/FST_list_NAmRF3_Test$point.csv
awk -v point=$point  -F ","  '{ if ($2==point || $3==point ) print }' $IN_TXT/FST_list_NAmRF4.csv  >> ${OUT_TXT}_${point}/FST_list_NAmRF3_Test$point.csv
awk -v point=$point  -F ","  '{ if ($2==(point+1) || $3==(point+1) ) print }' $IN_TXT/FST_list_NAmRF4.csv  >> ${OUT_TXT}_${point}/FST_list_NAmRF3_Test$point.csv
awk -v point=$point  -F ","  '{ if ($2!=point && $2!=(point+1) && $3!=point && $3!=(point+1) ) print }' $IN_TXT/FST_list_NAmRF4.csv  > ${OUT_TXT}_${point}/FST_list_NAmRF3_Trai$point.csv

