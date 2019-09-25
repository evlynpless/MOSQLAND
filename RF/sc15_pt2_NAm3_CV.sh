#!/bin/bash                                         
#SBATCH -p bigmem                                           
#SBATCH --mem=500g
#SBATCH -n 1 -c 8 -N 1
#SBATCH -t 24:00:00                                      
#SBATCH -o /gpfs/scratch60/fas/powell/esp38/stdout/sc15_pt2_NAm3_CV.sh.%J.out  
#SBATCH -e /gpfs/scratch60/fas/powell/esp38/stderr/sc15_pt2_NAm3_CV.sh.%J.err  
#SBATCH --mail-type=ALL                                  
#SBATCH --mail-user=evlyn.pless@yale.edu                                
#SBATCH --job-name=sc15_pt2_NAm3_CV.sh                     

####  for foldnum in 1 2 3 4 5 6 7 8 9 10; do sbatch --export=foldnum=$foldnum /home/fas/powell/esp38/scripts/MOSQLAND/RF/sc15_pt2_NAm3_CV.sh  ; done    

#for testing
####  for foldnum in 1; do sbatch --export=foldnum=$foldnum /home/fas/powell/esp38/scripts/MOSQLAND/RF/sc15_pt2_NAm3_CV.sh  ; done 

ulimit -c 0

module load Apps/R/3.3.2-generic

module load Rpkgs/RGDAL/1.2-5

module load Rpkgs/DOPARALLEL/1.0.3
# --slave      use if you only want to see output

export foldnum=$foldnum

R --vanilla --no-readline -q  -f /home/fas/powell/esp38/scripts/MOSQLAND/RF/sc15_pt2_NAm3_CV.R


