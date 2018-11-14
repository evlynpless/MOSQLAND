
# This script is for downloading consensus texture data from www.earthenv.org/texture. I'm using version 1.0, 1-km resolution to match the landcov data. 5-km and 25-km are also available. All data layers are in WGS84 projection and have a spatial extent from 85ºN - 60ºS and from 180ºW - 180ºE. The pixel values of the data layers should be mulitplied by 0.0001 to obtain the actual values of the metrics. Metrics based on EVI (enhances Vegetative Index).

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/texture
http://data.earthenv.org/habitat_heterogeneity/1km/cv_01_05_1km_uint16.tif
http://data.earthenv.org/habitat_heterogeneity/1km/evenness_01_05_1km_uint16.tif
http://data.earthenv.org/habitat_heterogeneity/1km/range_01_05_1km_uint16.tif
http://data.earthenv.org/habitat_heterogeneity/1km/shannon_01_05_1km_uint16.tif
http://data.earthenv.org/habitat_heterogeneity/1km/std_01_05_1km_uint16.tif
http://data.earthenv.org/habitat_heterogeneity/1km/Contrast_01_05_1km_uint32.tif
http://data.earthenv.org/habitat_heterogeneity/1km/Correlation_01_05_1km_int16.tif
http://data.earthenv.org/habitat_heterogeneity/1km/Dissimilarity_01_05_1km_uint32.tif
http://data.earthenv.org/habitat_heterogeneity/1km/Entropy_01_05_1km_uint16.tif
http://data.earthenv.org/habitat_heterogeneity/1km/Homogeneity_01_05_1km_uint16.tif
http://data.earthenv.org/habitat_heterogeneity/1km/Maximum_01_05_1km_uint16.tif
http://data.earthenv.org/habitat_heterogeneity/1km/Uniformity_01_05_1km_uint16.tif
http://data.earthenv.org/habitat_heterogeneity/1km/Variance_01_05_1km_uint32.tif
