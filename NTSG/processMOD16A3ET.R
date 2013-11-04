# title         : processMOD16A3ET.R
# purpose       : processing of the MODIS based Evapotranspiration product;
# reference     : [http://worldgrids.org/doku.php?id=wiki:layers#modis_products]
# producer      : Prepared by T. Hengl
# addresss      : In Wageningen, NL.
# inputs        : data available for download via FTP from [http://www.ntsg.umt.edu/project/mod16];
# outputs       : geotiff images projected in the "+proj=longlat +datum=WGS84" system;
# remarks 1     : First download and install FWtools [http://fwtools.maptools.org]; 

library(rgdal)
library(RSAGA)
myenv <- rsaga.env(path="C:/saga_vc")
library(XML)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"

## list all geotifs:
tif.lst <- list.files(pattern=glob2rx("MOD16A3_ET_*.tif$"))
tif.lst
## convert to SAGA GIS Format:
for(k in 2:length(tif.lst)){
  rsaga.geoprocessor(lib="io_gdal", 0, param=list(FILES=tif.lst[k], GRIDS=set.file.extension(tif.lst[k], ".sgrd"), TRANSFORM=FALSE, INTERPOL=0), env=myenv)
}

sdat.lst <- list.files(pattern=glob2rx("MOD16A3_ET_*.sgrd$"))
## derive mean value for all years:
rsaga.geoprocessor(lib="geostatistics_grid", module=4, param=list(GRIDS=paste(sdat.lst, collapse=";"), MEAN="ETMNTS3a.sgrd", MIN="tmp.sgrd", MAX="tmp.sgrd", STDDEV="ETSNTS3a.sgrd"), show.output.on.console = FALSE, env=myenv)
rsaga.geoprocessor(lib="io_gdal", module=2, param=list(GRIDS="ETMNTS3a.sgrd", FILE="tmp.tif"), env=myenv)
system(paste(gdalwarp, "tmp.tif ETMNTS3a.tif -t_srs \"+proj=longlat +datum=WGS84\" -r near -ot \"UInt16\" -dstnodata \"0\" -r near -te -180 -90 180 90 -tr", 1/120, 1/120))
GDALinfo("ETMNTS3a.tif")

## Resample to other resolutions using the majority value!
# 2.5 km:
rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT="ETMNTS3a.sgrd", USER_GRID="ETMNTS2a.sgrd", KEEP_TYPE=FALSE, SCALE_UP_METHOD=9, USER_SIZE=1/40, USER_XMIN=-180+1/40*0.5, USER_XMAX=180-1/40*0.5, USER_YMIN=-90+1/40*0.5, USER_YMAX=90-1/40*0.5), env=myenv) 
system(paste(gdal_translate, "ETMNTS2a.sdat ETMNTS2a.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"UInt16\" -a_nodata \"0\""))
## 5 km:
rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT="ETMNTS2a.sgrd", USER_GRID="ETMNTS1a.sgrd", KEEP_TYPE=FALSE, SCALE_UP_METHOD=9, USER_SIZE=1/20, USER_XMIN=-180+1/20*0.5, USER_XMAX=180-1/20*0.5, USER_YMIN=-90+1/20*0.5, USER_YMAX=90-1/20*0.5), env=myenv)
system(paste(gdal_translate, "ETMNTS1a.sdat ETMNTS1a.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"UInt16\" -a_nodata \"0\""))
## 20 km:
rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT="ETMNTS1a.sgrd", USER_GRID="ETMNTS0a.sgrd", KEEP_TYPE=FALSE, SCALE_UP_METHOD=9, USER_SIZE=1/5, USER_XMIN=-180+1/5*0.5, USER_XMAX=180-1/5*0.5, USER_YMIN=-90+1/5*0.5, USER_YMAX=90-1/5*0.5), env=myenv)
system(paste(gdal_translate, "ETMNTS0a.sdat ETMNTS0a.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"UInt16\" -a_nodata \"0\""))

## compress:  
for(i in 0:3){
    outname = paste("ETMNTS", i, 'a', sep="")
    system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), set.file.extension(outname, ".tif")))
    system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(outname, ".tif.gz"))
}


## end of script;