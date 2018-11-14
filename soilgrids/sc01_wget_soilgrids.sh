#!/bin/bash
#SBATCH -p day
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 24:00:00
#SBATCH -o /gpfs/scratch60/fas/sbsc/ga254/grace0/stdout/sc01_wget.sh.%J.out
#SBATCH -e /gpfs/scratch60/fas/sbsc/ga254/grace0/stderr/sc01_wget.sh.%J.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=email
#SBATCH --job-name=sc01_wget.sh


# sbatch /gpfs/home/fas/sbsc/ga254/scripts/SoilGrids/sc01_wget.sh

cd /project/fas/sbsc/ga254/grace0.grace.hpc.yale.internal/dataproces/SoilGrids/AWC

for seq in $( seq 1 7) ; do 
    wget ftp://ftp.soilgrids.org/data/recent/AWCtS_M_sl${seq}_250m.tif 
    gdal_edit.py -a_ullr -180 84 180 -56 AWCtS_M_sl${seq}_250m.tif 
done 

# depth   0 15 30 60 100 200 

cd /project/fas/sbsc/ga254/grace0.grace.hpc.yale.internal/dataproces/SoilGrids/TEXMHT
wget ftp://ftp.soilgrids.org/data/recent/TEXMHT_M_sl1_250m.tif
gdal_edit.py -a_ullr -180 84 180 -56 TEXMHT_M_sl1_250m.tif


