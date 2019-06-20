#!/bin/bash                                         
#SBATCH -p day                                           
#SBATCH -n 1 -c 1  -N 1                                  
#SBATCH --mem-per-cpu=10000
#SBATCH -t 24:00:00                                      
#SBATCH -o /gpfs/scratch60/fas/powell/esp38/stdout/test_plot.sh.%J.out  
#SBATCH -e /gpfs/scratch60/fas/powell/esp38/stderr/test_plot.sh.%J.err  
#SBATCH --mail-type=ALL                                  
#SBATCH --mail-user=evlyn.pless@yale.edu                                
#SBATCH --job-name=test_plot.sh                     


# sbatch   /home/fas/powell/esp38/scripts/MOSQLAND/RF/test_plot.sh  

module load Apps/R/3.3.2-generic

module load Rpkgs/RGDAL/1.2-5

R --vanilla -no-readline -q  -f  /home/fas/powell/esp38/scripts/MOSQLAND/RF/test_plot.R  
