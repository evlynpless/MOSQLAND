#This script is to download monthly mean precipitation data from Chelsa. Geographic Coordinate System [EPSG 4326], WGS84.

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/chelsa/prec

for n in $(seq 1 12 ) ; do
wget https://www.wsl.ch/lud/chelsa/data/climatologies/prec/CHELSA_prec_${n}_land.7z
7za x CHELSA_prec_${n}_land.7z
done

