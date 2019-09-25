#!/bin/bash                 
#SBATCH -p day                                                      
#SBATCH -n 1 -c 1 -N 1
#SBATCH -t 1:00:00
#SBATCH -o /gpfs/scratch60/fas/powell/esp38/stdout/test_runnum.sh.%J.out
#SBATCH -e /gpfs/scratch60/fas/powell/esp38/stderr/test_runnum.sh.%J.err
#SBATCH --mail-type=ALL 
#SBATCH --job-name=test_runnum.sh 


####  for runnum in 1:5; do sbatch --export=runnum=$runnum /home/fas/powell/esp38/scripts/MOSQLAND/RF/test_runnum.sh  ; done

ulimit -c 0

module load Apps/R/3.3.2-generic

export runnum=$runnum
R --vanilla --no-readline -q  -f  /home/fas/powell/esp38/scripts/MOSQLAND/RF/test_runnum.R 
