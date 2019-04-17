#!/bin/bash                                         
#SBATCH -p day                                           
#SBATCH -n 1 -c 1  -N 1                                  
#SBATCH --mem-per-cpu=100000
#SBATCH -t 24:00:00                                      
#SBATCH -o /gpfs/scratch60/fas/powell/esp38/stdout/sc07_iterativeRF_NAm1.sh.%J.out  
#SBATCH -e /gpfs/scratch60/fas/powell/esp38/stderr/sc07_iterativeRF_NAm1.sh.%J.err  
#SBATCH --mail-type=ALL                                  
#SBATCH --mail-user=evlyn.pless@yale.edu                                
#SBATCH --job-name=sc07_iterativeRF_NAm1.sh                     


# sbatch   /home/fas/powell/esp38/scripts/MOSQLAND/RF/sc07_iterativeRF_NAm1.sh  


module load Apps/R/3.3.2-generic

module load Rpkgs/RGDAL/1.2-5

R --vanilla -no-readline -q  -f  /home/fas/powell/esp38/scripts/MOSQLAND/RF/sc08_iterativeRF_NAm1.R  


