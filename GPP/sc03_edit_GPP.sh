cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/GPP/mnth/monthly_mean

OUTDIR=/project/fas/powell/esp38/dataproces/MOSQLAND/consland/GPP/mnth/monthly_mean/Florida_clips

for n in 01 02 03 04 05 06 07 08 09 10 11 12 ; do                                                         

gdal_edit.py -a_ullr -180   84  180  -90 GPP.VPM.mean.M${n}.v20.CMG.tif                                                                                                                                                                                                              
gdalwarp -tr 0.008333333333333 -0.008333333333333 GPP.VPM.mean.M${n}.v20.CMG.tif GPP.VPM.mean.M${n}.v20.CMG_res.tif

done

#Make an annual average
gdalbuildvrt -separate -overwrite  GPP.vrt    GPP.VPM.mean.M01.v20.CMG_res.tif     GPP.VPM.mean.M02.v20.CMG_res.tif     GPP.VPM.mean.M03.v20.CMG_res.tif     GPP.VPM.mean.M04.v20.CMG_res.tif     GPP.VPM.mean.M05.v20.CMG_res.tif      GPP.VPM.mean.M06.v20.CMG_res.tif     GPP.VPM.mean.M07.v20.CMG_res.tif     GPP.VPM.mean.M08.v20.CMG_res.tif    GPP.VPM.mean.M09.v20.CMG_res.tif    GPP.VPM.mean.M10.v20.CMG_res.tif    GPP.VPM.mean.M11.v20.CMG_res.tif    GPP.VPM.mean.M12.v20.CMG_res.tif   

pkstatprofile  -f mean  -i  GPP.vrt  -o  GPP_mean.tif

#Make a Florida average
gdal_translate  -projwin -85 31.5 -79.8 24.0  GPP_mean.tif   $OUTDIR/GPP_mean_Florida_clip.tif                             

gdal_translate -of AAIGrid $OUTDIR/GPP_mean_Florida_clip.tif  $OUTDIR/GPP_mean_Florida_clip.asc 
