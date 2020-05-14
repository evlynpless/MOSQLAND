#!/bin/bash                                        
#SBATCH -p day                                     
#SBATCH --mem=100g                                
#SBATCH -t 24:00:00                                
#SBATCH -o /gpfs/scratch60/fas/powell/esp38/stdout/sc12_prepare_rasterstack.sh.%J.out                   
#SBATCH -e /gpfs/scratch60/fas/powell/esp38/stderr/sc12_prepare_rasterstack.sh.%J.err                  
#SBATCH --mail-type=ALL                            
#SBATCH --mail-user=evlyn.pless@yale.edu           
#SBATCH --job-name=sc12_prepare_rasterstack.sh                                                                                        

# sbatch   /home/fas/powell/esp38/scripts/MOSQLAND/RF/sc12_prepare_rasterstack.sh                                                     
                 
  

ulimit -c 0
module load Apps/R/3.3.2-generic
module load Rpkgs/RGDAL/1.2-5
module load Rpkgs/DOPARALLEL/1.0.3
R --vanilla -no-readline -q  -f  /home/fas/powell/esp38/scripts/MOSQLAND/RF/sc12_prepare_rasterstack.R
