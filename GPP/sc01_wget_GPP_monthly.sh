#This script is to download monthly GPP (0.05 degree) from the Zhang et al 2017 publication.
#It also changes the no data value from infinite to -9999



cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/GPP/mnth
wget https://ndownloader.figshare.com/files/8945542
7za x 8945542

cd /project/fas/powell/esp38/dataproces/MOSQLAND/consland/GPP/mnth/monthly
for year in 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016; do
	for month in 01 02 03 04 05 06 07 08 09 10 11 12 ; do 

	pksetmask -ot Float32  -co COMPRESS=DEFLATE -co ZLEVEL=9    -m   GPP.VPM.${year}.M${month}.v20.CMG.tif  -msknodata -100000000 -nodata -9999 -p "<"    -i   GPP.VPM.${year}.M${month}.v20.CMG.tif  -o  GPP.VPM.${year}.M${month}.v20.CMG_test.tif 
	mv  GPP.VPM.${year}.M${month}.v20.CMG_test.tif  GPP.VPM.${year}.M${month}.v20.CMG.tif 
	done
done
