# title         : getISG.R
# purpose       : download and resampling of global surface geology;
# reference     : [http://worldgrids.org]
# producer      : Prepared by E. Ribeiro and T. Hengl
# address       : In Wageningen, NL.
# inputs        : data available for download from the  International Surface Geology project [http://certmapper.cr.usgs.gov/data/envision/index.html?widgets=geologymaps]; 
# outputs       : geotiff images projected in the "+proj=longlat +ellps=WGS84" system;
# remarks 1     : the shapefile was manually corrected for overlap and inconsistencies; 


library(RSAGA) 
library(rgdal)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"

# download files from server:
download.file("http://worldgrids.org/lib/exe/fetch.php?media=wd_geological_ages_shapefile.7z", destfile=paste(getwd(), "wd_geological_ages_shapefile.7z", sep="/"))
system("7za -e wd_geological_ages_shapefile.7z")
ogrInfo("wd_geological_ages.shp", "wd_geological_ages")
tbl <- read.dbf("wd_geological_ages.dbf")  ## Large DBF!
str(tbl)
tbl$dbf$GEA <- as.integer(tbl$dbf$reclass)
write.dbf(tbl, "wd_geological_ages.dbf")
levs <- levels(tbl$dbf$reclass)
levs.tbl <- data.frame(Number=1:length(levs), geological_age=levs)

## convert to raster (use the polygon ID):
myenv <- rsaga.env(path="C:/saga_vc")
rsaga.geoprocessor(lib="grid_gridding", module=0, param=list(USER_GRID="wd_geological_ages.sgrd", INPUT="wd_geological_ages.shp", FIELD=17, TARGET=0, LINE_TYPE=1, USER_SIZE=1/120, USER_XMIN=-180+1/120*0.5, USER_XMAX=180-1/120*0.5, USER_YMIN=-90+1/120*0.5, USER_YMAX=90-1/120*0.5), env=myenv) # takes cca 5 mins uses > 6GB!!
## export to geotiff:
#system(paste(gdal_translate, "wd_geological_ages.sdat GEAISG3a.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\""))
rsaga.geoprocessor(lib="io_gdal", module=2, param=list(GRIDS="wd_geological_ages.sgrd", FILE="wd_geological_ages.tif"), env=myenv)
## TH: GDAL warp with SAGA sdat file receives error "ERROR 3: Unable to seek to beginning of grid row."
system(paste(gdal_translate, "wd_geological_ages.tif GEAISG3a.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -a_nodata \"0\""))
GDALinfo("GEAISG3a.tif")

## Resample to other resolutions using the majority value!
# 2.5 km:
rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT="wd_geological_ages.sgrd", USER_GRID="GEAISG2a.sgrd", KEEP_TYPE=FALSE, SCALE_UP_METHOD=9, USER_SIZE=1/40, USER_XMIN=-180+1/40*0.5, USER_XMAX=180-1/40*0.5, USER_YMIN=-90+1/40*0.5, USER_YMAX=90-1/40*0.5), env=myenv) 
system(paste(gdal_translate, "GEAISG2a.sdat GEAISG2a.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -a_nodata \"0\""))
## 5 km:
rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT="GEAISG2a.sgrd", USER_GRID="GEAISG1a.sgrd", KEEP_TYPE=FALSE, SCALE_UP_METHOD=9, USER_SIZE=1/20, USER_XMIN=-180+1/20*0.5, USER_XMAX=180-1/20*0.5, USER_YMIN=-90+1/20*0.5, USER_YMAX=90-1/20*0.5), env=myenv)
system(paste(gdal_translate, "GEAISG1a.sdat GEAISG1a.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -a_nodata \"0\""))
## 20 km:
rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT="GEAISG1a.sgrd", USER_GRID="GEAISG0a.sgrd", KEEP_TYPE=FALSE, SCALE_UP_METHOD=9, USER_SIZE=1/5, USER_XMIN=-180+1/5*0.5, USER_XMAX=180-1/5*0.5, USER_YMIN=-90+1/5*0.5, USER_YMAX=90-1/5*0.5), env=myenv)
system(paste(gdal_translate, "GEAISG0a.sdat GEAISG0a.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -a_nodata \"0\""))

## compress:  
for(i in 0:3){
    outname = paste("GEAISG", i, 'a', sep="")
    system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), set.file.extension(outname, ".tif")))
    system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(outname, ".tif.gz"))
}

## make a PAL file:
library(RCurl)
# verify the certificate:
curl <- getCurlHandle()
options(RCurlOptions = list(capath = system.file("CurlSSL", "cacert.pem", package = "RCurl"), ssl.verifypeer = FALSE))
curlSetOpt(.opts = list(proxy = 'proxyserver:port'), curl = curl)

## legend entries:
cat(getURL("https://docs.google.com/spreadsheet/pub?key=0Ah5Ip0avaTkLdHhxTDVraGlZOU1maW8zeTQzRWFJd2c&single=true&gid=1&output=csv"), file = "legend.csv")
rgb.tbl <- read.csv("legend.csv")
str(rgb.tbl)
#rgb.tbl <- merge(levs.tbl, rgb.tbl, all.y=FALSE, all.x=TRUE)

## create a SAGA txt colour table
## convert to BGR codes:
BGR <- (rgb.tbl$B * 65536) + (rgb.tbl$G * 256) + rgb.tbl$R
## write a lookup table for SAGA GIS:
filename <- file("GEAISG.txt", "w", blocking=FALSE)
write("COLOR\tNAME\tDESCRIPTION\tMINIMUM\tMAXIMUM", filename)
for(i in 1:nrow(rgb.tbl)){
  write(paste(BGR[i], rgb.tbl[i,"geological_age"], paste("CL", i, sep=""), (rgb.tbl[i,"Number"]-1)+0.1, (rgb.tbl[i,"Number"])+0.1, sep="\t"), filename, append=TRUE)
}
close(filename)

## write a PAL file:
write.table(rgb.tbl[,c("geological_age","R","G","B")], "GEAISG.PAL", quote = FALSE,  col.names = FALSE)

## end of script;