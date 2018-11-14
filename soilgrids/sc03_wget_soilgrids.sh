#!/bin/bash
#SBATCH -p day
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 24:00:00
#SBATCH -o /home/fas/powell/esp38/scratch60/sc01_wget.sh.%J.out #check this?
#SBATCH -e /home/fas/powell/esp38/scratch60/sc01_wget.sh.%J.err #check this?
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email
#SBATCH --job-name=sc03_wget_soilgrids.sh


# sbatch /home/fas/powell/esp38/scripts/MOSQLAND/soilgrids/sc03_wget_soilgrids.sh

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/soilgrids/SWI

#does this part need updating?
for seq in $( seq 1 7) ; do 
    wget ftp://ftp.soilgrids.org/data/recent/AWCtS_M_sl${seq}_250m.tif 
    gdal_edit.py -a_ullr -180 84 180 -56 AWCtS_M_sl${seq}_250m.tif 
done 


# depth   0 15 30 60 100 200 

cd /project/fas/sbsc/ga254/grace0.grace.hpc.yale.internal/dataproces/SoilGrids/TEXMHT
wget ftp://ftp.soilgrids.org/data/recent/TEXMHT_M_sl1_250m.tif
gdal_edit.py -a_ullr -180 84 180 -56 TEXMHT_M_sl1_250m.tif


#Modified based on Giuseppe's script: sc01_wget_soilgrids.sh