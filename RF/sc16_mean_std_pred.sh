#!/bin/bash                                         
#SBATCH -p day                                           
#SBATCH -n 1 -c 1  -N 1                                  
#SBATCH --mem-per-cpu=50000
#SBATCH -t 10:00:00                                      
#SBATCH -o /gpfs/scratch60/fas/powell/esp38/stdout/sc16_mean_std_pred.sh.%J.out  
#SBATCH -e /gpfs/scratch60/fas/powell/esp38/stderr/sc16_mean_std_pred.sh.%J.err  
#SBATCH --mail-type=ALL                                  
#SBATCH --mail-user=evlyn.pless@yale.edu                                
#SBATCH --job-name=sc16_mean_std_pred.sh                     


# sbatch   /home/fas/powell/esp38/scripts/MOSQLAND/RF/sc16_mean_std_pred.sh  


module load Apps/R/3.3.2-generic

module load Rpkgs/RGDAL/1.2-5

R --vanilla -no-readline -q  -f  /home/fas/powell/esp38/scripts/MOSQLAND/RF/sc16_mean_std_pred.R  
