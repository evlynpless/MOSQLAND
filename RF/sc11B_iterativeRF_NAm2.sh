#!/bin/bash                                         
#SBATCH -p bigmem                                           
#SBATCH --mem=500g
#SBATCH -n 1 -c 8 -N 1
#SBATCH -t 24:00:00                                      
#SBATCH -o /gpfs/scratch60/fas/powell/esp38/stdout/sc11B_iterativeRF_NAm2.sh.%J.out  
#SBATCH -e /gpfs/scratch60/fas/powell/esp38/stderr/sc11B_iterativeRF_NAm2.sh.%J.err  
#SBATCH --mail-type=ALL                                  
#SBATCH --mail-user=evlyn.pless@yale.edu                                
#SBATCH --job-name=sc11B_iterativeRF_NAm2.sh                     


# sbatch   /home/fas/powell/esp38/scripts/MOSQLAND/RF/sc11B_iterativeRF_NAm2.sh  

ulimit -c 0

module load Apps/R/3.3.2-generic

module load Rpkgs/RGDAL/1.2-5

module load Rpkgs/DOPARALLEL/1.0.3
# --slave      use if you only want to see output
R --vanilla --no-readline -q  -f  /home/fas/powell/esp38/scripts/MOSQLAND/RF/sc11B_iterativeRF_NAm2.R  


