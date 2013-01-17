# title         : getDMSP.R
# purpose       : download and resampling of Version 4 DMSP-OLS Nighttime Lights Time Series;
# reference     : [http://worldgrids.org/doku.php?id=wiki:ln1dms1]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, January 2013.
# inputs        : geotiff images (1km resolution) available for download from [http://ngdc.noaa.gov/dmsp/downloadV4composites.html];
# outputs       : maps projected in the "+proj=longlat +ellps=WGS84" system;
# remarks 1     : First download and install FWtools [http://fwtools.maptools.org];  

# ------------------------------------------------------------
# Initial settings and data download:
# ------------------------------------------------------------

library(R.utils)
library(utils)
library(rgdal)
library(RSAGA)
library(RCurl)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
gdalbuildvrt <- shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalbuildvrt.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"
workd <- normalizePath(getwd())

## get a list of files available on server:
URL <- "http://ngdc.noaa.gov/dmsp/data/web_data/v4composites/"
# list all maps available:
items <- strsplit(getURL(URL), "\n")[[1]]
# convert to a character vector:
zip.list <- items[grep(items, pattern=".tar")]
for(i in 1:length(zip.list)) {zip.list[i] <- strsplit(strsplit(zip.list[i], '.tar\">')[[1]][2], '</a></td>')[[1]][1]}
zip.list  # 31 image in total!
zip.list <- zip.list[!is.na(zip.list)]

## download files:
for(i in 1:length(zip.list)) {
 if(is.na(file.info(zip.list[i])$size)){
  download.file(paste(URL, zip.list[i], sep=""), destfile=paste(getwd(), "/", zip.list[i], sep=""), method="wget") 
 }
}

## untar:
for(i in 1:length(zip.list)) {
  tif.name = paste(strsplit(zip.list[i], ".tar")[[1]][1], "b_web.stable_lights.avg_vis.tif.gz", sep="")
  if(is.na(file.info(tif.name)$size)){
    try(untar(zip.list[i], files=tif.name))
  }
}
## un-gzip:
gz.name <- list.files(pattern=glob2rx("*.gz$"))
for(i in 1:length(gz.name)) {
  tif.name = paste(strsplit(gz.name[i], ".gz")[[1]][1])
  if(is.na(file.info(tif.name)$size)){
    gunzip(gz.name[i], overwrite=TRUE, remove=TRUE)
  }
}

## convert to SAGA GRIDS:
tif.lst <- list.files(pattern=glob2rx("*web.stable_lights.*.tif$"))
for(i in 1:length(tif.lst)) {
  rsaga.geoprocessor(lib="io_gdal", module=0, param=list(GRIDS=set.file.extension(tif.lst[i], ".sgrd"), FILES=tif.lst[i]))
  unlink(tif.lst[i])
}  

## run PCA in SAGA GIS:
sgrd.lst <- list.files(pattern=glob2rx("*web.stable_lights.*.sgrd$"))
rsaga.geoprocessor(lib="geostatistics_grid", module=8, param=list(GRIDS=paste(sgrd.lst, collapse=";", sep=""), PCA="LN1DMS3a.sgrd;LN2DMS3a.sgrd", METHOD=1, NFIRST=2)) ## TAKES >32GB memory!
## TH: results - PC1 93.85%, PC2 95.90%, PC3 96.68%

## mean value and sd:
rsaga.geoprocessor(lib="geostatistics_grid", module=4, param=list(GRIDS=paste(sgrd.lst, collapse=";", sep=""), MEAN="LNMDMS3a.sgrd", MIN="tmp.sgrd", MAX="tmp.sgrd", MAX="tmp.sgrd", VAR="tmp.sgrd", STDDEV="LNSDMS3a.sgrd", STDDEVLO="tmp.sgrd", STDDEVHI="tmp.sgrd"))

## resample to 1 km whole world:
out.lst <- c("LN1DMS3a", "LN2DMS3a", "LNMDMS3a")
## TH: this leads to many problems as the missing values need to be set manually in SAGA GIS, and the line on the top of the image has -99999 values which have to be basked out!
#mvFlag.lst <- c(1.838122, 0.088346, 0)
for(i in 1:length(out.lst)) {
  system(paste(gdalwarp, ' ', out.lst[i], '.sdat ', out.lst[i], '.tif -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""))
}
GDALinfo("LN1DMS3a.tif")

## resample to coarser resolutions:
# 20 km:
tif.lst <- list.files(pattern=glob2rx("LN*DMS3a.tif$"))
for(i in 1:length(tif.lst)){
  outname = gsub("3a", "2a", tif.lst[i])
  if(is.na(file.info(outname)$size)){ 
    system(paste(gdalwarp, ' ', tif.lst[i], ' ', outname, ' -r bilinear -te -180 -90 180 90 -tr ', 1/40,' ', 1/40, sep=""))
  }
  outname = gsub("3a", "1a", tif.lst[i])
  if(is.na(file.info(outname)$size)){ 
    system(paste(gdalwarp, ' ', tif.lst[i], ' ', outname, ' -r bilinear -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
  }  
  outname = gsub("3a", "0a", tif.lst[i])
  if(is.na(file.info(outname)$size)){ 
    system(paste(gdalwarp, ' ', tif.lst[i], ' ', outname, ' -r bilinear -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))
  }  
}

## compress and copy:
tif2.lst <- list.files(pattern=glob2rx("LN*DMS*a.tif$"))  
for(i in 1:length(tif2.lst)){
    system(paste("7za a", "-tgzip", set.file.extension(tif2.lst[i], ".tif.gz"), tif2.lst[i]))
    system(paste("xcopy", set.file.extension(tif2.lst[i], ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(tif2.lst[i], ".tif.gz"))
}

# end of script;
