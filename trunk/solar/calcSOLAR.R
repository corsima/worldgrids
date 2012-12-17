# title         : calcSOLAR.R
# purpose       : extraction of the potential solar insolation globally (mean value and standard deviation);
# reference     : [https://code.google.com/p/worldgrids/source/browse/]
# producer      : Prepared by M. Kilibarda and T. Hengl
# version       : 1
# inputs        : solar radiation derived for various days in SAGA GIS from globe DEM [http://worldgrids.org/doku.php?id=wiki:demsre3];
# outputs       : geotiff images projected in the "+proj=longlat +datum=WGS84" system;
# remarks 1     : First download and install FWtools [http://fwtools.maptools.org];

# -------------------------------------------
# Initial settings and data download:
# -------------------------------------------

library(RSAGA) 
library(rgdal)
library(GSIF)
library(maptools)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdalbuildvrt = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalbuildvrt.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"

## list of images:
dtif.lst <- list.files(pattern=glob2rx("*SRT3a.tif$"))
inf <- GDALinfo(dtif.lst[1])

## tiles:
b.l <- expand.grid(KEEP.OUT.ATTRS=FALSE, lonmin=seq(-180,0,by=180), latmin=seq(-90,0,by=90))
b.u <- expand.grid(KEEP.OUT.ATTRS=FALSE, lonmax=seq(0,180,by=180), latmax=seq(0,90,by=90))
btiles <- cbind(b.l, b.u)
str(btiles)

for(j in 1:nrow(btiles)){
  for(i in 1:length(dtif.lst)){
    outname <- strsplit(dtif.lst[i], ".tif")[[1]][1]
    tilename <- paste(outname, '_', j, '.sdat', sep="")
    # resample
    if(is.na(file.info(tilename)$size)){
      system(paste(gdalwarp, ' ', dtif.lst[i], ' ', tilename, ' -of \"SAGA\" -dstnodata \"255\" -r near -te ', btiles[j,"lonmin"] , ' ', btiles[j,"latmin"], ' ', btiles[j,"lonmax"] ,' ', btiles[j,"latmax"] ,' -tr ', 1/120, ' ', 1/120, sep=""))
    }
  }
  # derive derive mean, sd, min and max values:
  t.lst <- list.files(pattern=glob2rx(paste('*3a_', j, '.sgrd$', sep="")))
  if(is.na(file.info(paste("INMSRE3a_", j, ".sgrd", sep=""))$size)){
    rsaga.geoprocessor(lib="geostatistics_grid", module=4, param=list(GRIDS=paste(t.lst, collapse=";"), MEAN=paste("INMSRE3a_", j, ".sgrd", sep=""), MIN=paste("INLSRE3a_", j, ".sgrd", sep=""), MAX=paste("INHSRE3a_", j, ".sgrd", sep=""), STDDEV=paste("INSSRE3a_", j, ".sgrd", sep="")))
  }
}

# create a mosaic and compress files:
for(j in c("INMSRE", "INSSRE", "INLSRE", "INHSRE")){
  unlink("insol.vrt")
  system(paste(gdalbuildvrt, "insol.vrt", paste(j, "3a_", 1:nrow(btiles), ".sdat", sep="", collapse=" ")))
  # convert to geotif (1 km):
  if(is.na(file.info(paste(j, '3a.tif', sep=""))$size)){
  if(j == "INSSRE"){
    system(paste(gdalwarp, ' insol.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-99999\" ', j, '3a.tif', ' -r near -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""))
  } else {
    system(paste(gdalwarp, ' insol.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"255\" -ot Byte ', j, '3a.tif', ' -r near -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""))  
  }
  }
  # 2.5 km:
  if(is.na(file.info(paste(j, '2a.tif', sep=""))$size)){
  if(j == "INSSRE"){
    system(paste(gdalwarp, ' insol.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-99999\" ', j, '2a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/40,' ', 1/40, sep=""))
  } else {
  system(paste(gdalwarp, ' insol.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"255\" -ot Byte ', j, '2a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/40,' ', 1/40, sep=""))
  }
  }
  # 5 km:
  if(is.na(file.info(paste(j, '1a.tif', sep=""))$size)){
  if(j == "INSSRE"){
    system(paste(gdalwarp, ' insol.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-99999\" ', j, '1a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
  } else {
  system(paste(gdalwarp, ' insol.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"255\" -ot Byte ', j, '1a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
  }
  }
  # 20 km:
  if(is.na(file.info(paste(j, '0a.tif', sep=""))$size)){
  if(j == "INSSRE"){
    system(paste(gdalwarp, ' ', j, '1a.tif -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-99999\" ', j, '0a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))
  } else {
  system(paste(gdalwarp, ' ', j, '1a.tif -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"255\" -ot Byte ', j, '0a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))
  }
  }
}

## compress:  
for(j in c("INMSRE", "INSSRE", "INLSRE", "INHSRE")){
for(i in 0:3){
    outname = paste(j, i, 'a', sep="")
    system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), set.file.extension(outname, ".tif")))
    system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(outname, ".tif.gz"))
  }
}

# clean up temp files:
unlink(list.files(pattern=glob2rx("IN*SRE_*.*")))
unlink(list.files(pattern=glob2rx("I*SRT3a_*.*")))

# end of script;