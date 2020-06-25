#!/bin/bash                                        
#SBATCH -p day                                     
#SBATCH --mem=100g                                
#SBATCH -t 24:00:00                                
#SBATCH -o /gpfs/scratch60/fas/powell/esp38/stdout/sc16_prepare_numeric_rasterstack.sh.%J.out                   
#SBATCH -e /gpfs/scratch60/fas/powell/esp38/stderr/sc16_prepare_numeric_rasterstack.sh.%J.err                  
#SBATCH --mail-type=ALL                            
#SBATCH --mail-user=evlyn.pless@yale.edu           
#SBATCH --job-name=sc16_prepare_numeric_rasterstack.sh                                                                                        

# sbatch   /home/fas/powell/esp38/scripts/MOSQLAND/RF/sc16_prepare_numeric_rasterstack.sh

module load R/3.5.3-foss-2018a-X11-20180131


R --vanilla --no-readline   -q  <<'EOF'

load(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc12_rasterstack_image_withKernel.RData")

library("raster")

value.raster = as.data.frame(getValues(arid))
names(value.raster) = "arid"
value.raster$access = getValues(access)
value.raster$prec = getValues(prec)
value.raster$mean.temp = getValues(mean.temp)
value.raster$human.density = getValues(human.density)
value.raster$friction = getValues(friction)
value.raster$min.temp = getValues(min.temp)
value.raster$Needleleaf = getValues(Needleleaf)
value.raster$EvBroadleaf = getValues(EvBroadleaf)
value.raster$DecBroadleaf = getValues(DecBroadleaf)
value.raster$MiscTrees = getValues(MiscTrees)
value.raster$Shrubs = getValues(Shrubs)
value.raster$Herb = getValues(Herb)
value.raster$Crop = getValues(Crop)
value.raster$Flood = getValues(Flood)
value.raster$Urban = getValues(Urban)
value.raster$Snow = getValues(Snow)
value.raster$Barren = getValues(Barren)
value.raster$Water = getValues(Water)
value.raster$Slope = getValues(Slope)
value.raster$Altitude = getValues(Altitude)
value.raster$PET = getValues(PET)
value.raster$DailyTempRange = getValues(DailyTempRange)
value.raster$max.temp = getValues(max.temp)
value.raster$AnnualTempRange = getValues(AnnualTempRange)
value.raster$prec.wet = getValues(prec.wet)
value.raster$prec.dry = getValues(prec.dry)
value.raster$GPP = getValues(GPP)
value.raster$kernel100 = getValues(kernel100)

save.image(file = "/project/fas/powell/esp38/dataproces/MOSQLAND/consland/RF/NAm_RF_3/sc16_rasterstack.RData")

EOF

