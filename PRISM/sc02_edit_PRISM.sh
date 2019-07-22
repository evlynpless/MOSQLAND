
OUTDIR = /project/fas/powell/esp38/dataproces/MOSQLAND/consland/PRISM

#Editing dewpoint

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/PRISM/tdmean

gdal_translate -of GTiff PRISM_tdmean_30yr_normal_800mM2_annual_asc.asc  PRISM_tdmean_30yr_normal_800mM2_annual_asc.tif

gdalwarp -tr 0.008333333333333 0.008333333333333 PRISM_tdmean_30yr_normal_800mM2_annual_asc.tif PRISM_tdmean_30yr_normal_800mM2_annual_asc2.tif

gdalwarp PRISM_tdmean_30yr_normal_800mM2_annual_asc2.tif PRISM_tdmean_30yr_normal_800mM2_annual_asc3.tif -t_srs "+proj=longlat +ellps=WGS84"


#Calulating vapor pressure (in kPa) from Tetens conversion
#vapor pressure = 0.6108*((17.27*tdmean)/(tdmean+237.3)

gdal_calc.py -A PRISM_tdmean_30yr_normal_800mM2_annual_asc3.tif  --outfile=vapor_pressure.tif --NoDataValue=-9999  --calc="0.6108*((17.27*A)/(A+237.3))"

#Editing temperature

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/PRISM/tmean

gdal_translate -of GTiff PRISM_tmean_30yr_normal_800mM2_annual_asc.asc PRISM_tmean_30yr_normal_800mM2_annual_asc.tif

gdalwarp -tr 0.008333333333333 0.008333333333333 PRISM_tmean_30yr_normal_800mM2_annual_asc.tif PRISM_tmean_30yr_normal_800mM2_annual_asc2.tif

gdalwarp  PRISM_tmean_30yr_normal_800mM2_annual_asc2.tif   PRISM_tmean_30yr_normal_800mM2_annual_asc3.tif -t_srs "+proj=longlat +ellps=WGS84"


#Calculating Absolute humidity

gdal_calc.py -A /project/fas/powell/esp38/dataproces/MOSQLAND/consland/PRISM/tdmean/vapor_pressure.tif -B /project/fas/powell/esp38/dataproces/MOSQLAND/consland/PRISM/tmean/PRISM_tmean_30yr_normal_800mM2_annual_asc3.tif --outfile=$OUTDIR/PRISM_abshum.tif --NoDataValue=-9999  --calc="(2165*A)/(B+273.16)"
