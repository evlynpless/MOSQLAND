#!/bin/bash                                         
#SBATCH -p day                                           
#SBATCH -n 1 -c 1  -N 1                                                                                 
#SBATCH --mem-per-cpu=50000                                                                            
#SBATCH -t 10:00:00  
#SBATCH -o /gpfs/scratch60/fas/powell/esp38/stdout/sc15_NAm3_CV.sh.%J.out  
#SBATCH -e /gpfs/scratch60/fas/powell/esp38/stderr/sc15_NAm3_CV.sh.%J.err  
#SBATCH --mail-type=ALL                                  
#SBATCH --mail-user=evlyn.pless@yale.edu                                
#SBATCH --job-name=sc15_NAm3_CV.sh                     

# sbatch /home/fas/powell/esp38/scripts/MOSQLAND/RF/sc15_NAm3_CV.sh ; done   

ulimit -c 0

module load Apps/R/3.3.2-generic

module load Rpkgs/RGDAL/1.2-5

module load Rpkgs/DOPARALLEL/1.0.3
# --slave      use if you only want to see output

R --vanilla --no-readline -q  -f /home/fas/powell/esp38/scripts/MOSQLAND/RF/sc15_NAm3_CV.R


