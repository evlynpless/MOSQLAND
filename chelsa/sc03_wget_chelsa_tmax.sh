#This script is to download monthly maximum temperature data from Chelsa. Geographic Coordinate System [EPSG 4326], WGS84.

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/tmax

for n in $(seq 1 12 ) ; do
wget https://www.wsl.ch/lud/chelsa/data/climatologies/temp/integer/tmax/CHELSA_tmax10_${n}_land.7z
7za x CHELSA_tmax10_${n}_land.7z
done

