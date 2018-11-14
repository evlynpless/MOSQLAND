#This script is to download monthly minimum temperature data from Chelsa. Geographic Coordinate System [EPSG 4326], WGS84.

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/tmin

for n in $(seq 1 12 ) ; do
wget https://www.wsl.ch/lud/chelsa/data/climatologies/temp/integer/tmin/CHELSA_tmin10_${n}_land.7z
7za x CHELSA_tmin10_${n}_land.7z
done

