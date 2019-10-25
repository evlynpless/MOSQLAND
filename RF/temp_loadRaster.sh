#!/bin/bash                                                                                                    
#SBATCH -p day                                                                                                 
#SBATCH -n 1 -c 1  -N 1                                                                                       
#SBATCH -t 10:00:00                                                                                            
#SBATCH -o /gpfs/scratch60/fas/powell/esp38/stdout/temp_loadRaster.sh.%J.out                                   
#SBATCH -e /gpfs/scratch60/fas/powell/esp38/stdout/temp_loadRaster.sh.%J.err                                   
#SBATCH --mail-type=ALL                                                                                        
#SBATCH --mail-user=email                                                                                      
#SBATCH --job-name=temp_loadRaster.sh                                                                          

# sbatch   /home/fas/powell/esp38/scripts/MOSQLAND/RF/temp_loadRaster.sh                                       
                    

module load Apps/R/3.3.2-generic

module load Rpkgs/RGDAL/1.2-5

R --vanilla --no-readline -q  -f /home/fas/powell/esp38/scripts/MOSQLAND/RF/temp_loadRaster.R
