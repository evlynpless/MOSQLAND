
# This script is for downloading consensus landcover data from www.earthenv.org/landcover. Each pixel shows percentage for each landcover. 

#(Full version 1.0 with DISCover) using wget link. 1-km resolution. Consensus dataset integrating GlobCover, MODIS, GLC2000, and DISCover. Unsigned 8-bit values, ranging from 0-100, representing consensus prevalence in percentage. All data layers have a spatial extent from 90ºN - 56ºS and from 180ºW - 180ºE, and have a spatial resolution of 30 arc-second per pixel (~1 km per pixel at the equator). What is the projection? "All maps are in Behrmann projection" (but this might just refer to the figures in the paper?).

#1. Evergreen/deciduous needleleaf trees
#2. Evergreen broadleaf trees
#3. Deciduous broadleaf trees
#4. Mixed/other trees
#5. Shrubs
#6. Herbaceous vegetation
#7. Cultivated and managed vegetation
#8. Regularly flooded vegetation
#9. Urban/built-up
#10.Snow/ice
#11.Barren
#12. Open Water

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/landcov
for n in $(seq 1 12 ) ; do 
wget http://data.earthenv.org/consensus_landcover/with_DISCover/consensus_full_class_$n.tif
done 
