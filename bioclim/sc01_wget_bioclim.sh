#This script is to download BioClim data from Chelsa. Geographic Coordinate System [EPSG 4326], WGS84.                                                                                                
cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/bioclim

#Bio1 = Annual Mean Temperature
wget https://www.wsl.ch/lud/chelsa/data/bioclim/integer/CHELSA_bio10_01_land.7z
7za x CHELSA_bio10_01_land.7z

#Bio2 = Mean Diurnal Range (Mean of monthly (max temp - min temp))
wget https://www.wsl.ch/lud/chelsa/data/bioclim/integer/CHELSA_bio10_02_land.7z
7za x CHELSA_bio10_02_land.7z

#Bio3 = Isothermality (BIO2/BIO7) (* 100)
wget https://www.wsl.ch/lud/chelsa/data/bioclim/integer/CHELSA_bio10_03_land.7z
7za x CHELSA_bio10_03_land.7z

#Bio4 = Temperature Seasonality (standard deviation *100)
wget https://www.wsl.ch/lud/chelsa/data/bioclim/integer/CHELSA_bio10_04_land.7z
7za x CHELSA_bio10_04_land.7z

#Bio5 = Max Temperature of Warmest Month
wget https://www.wsl.ch/lud/chelsa/data/bioclim/integer/CHELSA_bio10_05_land.7z
7za x CHELSA_bio10_05_land.7z

#Bio6 = Min Temperature of Coldest Month
wget https://www.wsl.ch/lud/chelsa/data/bioclim/integer/CHELSA_bio10_06_land.7z
7za x CHELSA_bio10_06_land.7z

#Bio7 = Temperature Annual Range
wget https://www.wsl.ch/lud/chelsa/data/bioclim/integer/CHELSA_bio10_07_land.7z
7za x CHELSA_bio10_07_land.7z

#Bio8 = Mean Temperature of Wettest Quarter
wget https://www.wsl.ch/lud/chelsa/data/bioclim/integer/CHELSA_bio10_08_land.7z
7za x CHELSA_bio10_08_land.7z

#Bio9 = Mean Temperature of Driest Quarter
wget https://www.wsl.ch/lud/chelsa/data/bioclim/integer/CHELSA_bio10_09_land.7z
7za x CHELSA_bio10_09_land.7z

#Bio10 = Mean Temperature of Warmest Quarter
wget https://www.wsl.ch/lud/chelsa/data/bioclim/integer/CHELSA_bio10_10_land.7z
7za x CHELSA_bio10_10_land.7z

#Bio11 = Mean Temperature of Coldest Quarter
wget https://www.wsl.ch/lud/chelsa/data/bioclim/integer/CHELSA_bio10_11_land.7z
7za x CHELSA_bio10_11_land.7z

#Bio12 = Annual Precipitation
wget https://www.wsl.ch/lud/chelsa/data/bioclim/integer/CHELSA_bio10_12_land.7z
7za x CHELSA_bio10_12_land.7z

#Bio13 = Precipitation of Wettest Month
wget https://www.wsl.ch/lud/chelsa/data/bioclim/integer/CHELSA_bio10_13_land.7z
7za x CHELSA_bio10_13_land.7z

#Bio14 = Precipitation of Driest Month
wget https://www.wsl.ch/lud/chelsa/data/bioclim/integer/CHELSA_bio10_14_land.7z
7za x CHELSA_bio10_14_land.7z

#Bio15 = Precipitation Seasonality (Coefficient of Variation)
wget https://www.wsl.ch/lud/chelsa/data/bioclim/integer/CHELSA_bio10_15_land.7z
7za x CHELSA_bio10_15_land.7z

#Bio16 = Precipitation of Wettest Quarter
wget https://www.wsl.ch/lud/chelsa/data/bioclim/integer/CHELSA_bio10_16_land.7z
7za x CHELSA_bio10_16_land.7z

#Bio17 = Precipitation of Driest Quarter
wget https://www.wsl.ch/lud/chelsa/data/bioclim/integer/CHELSA_bio10_17_land.7z
7za x CHELSA_bio10_17_land.7z

#Bio18 = Precipitation of Warmest Quarter
wget https://www.wsl.ch/lud/chelsa/data/bioclim/integer/CHELSA_bio10_18_land.7z
7za x CHELSA_bio10_18_land.7z

#Bio19 = Precipitation of Coldest Quarter
wget https://www.wsl.ch/lud/chelsa/data/bioclim/integer/CHELSA_bio10_19_land.7z
7za x CHELSA_bio10_19_land.7z
