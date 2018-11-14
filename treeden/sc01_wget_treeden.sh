#!/bin/bash
#SBATCH -p scavenge
#SBATCH -N 1
#SBATCH -c 1
#SBATCH -t 24:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=evlyn.pless@yale.edu


##   sbatch  ~/scripts/MOSQLAND/sc08_wget_treeden.sh  

#This script is to a Global Tree Density map from Crowther et al 2015

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/treeden
wget http://elischolar.library.yale.edu/context/yale_fes_data/article/1000/type/native/viewcontent
7za x viewcontent