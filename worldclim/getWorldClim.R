# title         : getWorldClim.R
# purpose       : download and export/resampling of BioClimatic parameters (www.worldclim.org);
# reference     : [http://worldgrids.org]
# producer      : Prepared by T. Hengl
# address       : In Wageningen, NL
# inputs        : images available at [http://biogeo.berkeley.edu/worldclim1_4/grid/cur/]; 
# outputs       : geotiff images projected in the "+proj=longlat +ellps=WGS84" system;
# remarks 1     : Description of the data available at [https://www.worldclim.org];
# remarks 2     : First download and install FWtools [http://fwtools.maptools.org]; 

# ------------------------------------------------------------
# Initial settings and data download:
# ------------------------------------------------------------

library(rgdal)
library(RSAGA)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
gdalbuildvrt = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalbuildvrt.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
si <- Sys.info()
crs <- "+proj=longlat +datum=WGS84"
outdir <- "G:/WORLDGRIDS/maps"

# download files from server:
download.file("http://biogeo.berkeley.edu/worldclim1_4/grid/cur/prec_30s_bil.zip", destfile=paste(getwd(), "prec_30s_bil.zip", sep="/"))
unzip(zipfile="prec_30s_bil.zip", exdir=getwd())
workd <- paste(gsub("/", "\\\\", getwd()), "\\", sep="")

## periods:
bil.name <- dir(path=getwd(), pattern=glob2rx("prec_*.bil"))
x <- as.numeric(sapply(bil.name, function(x){strsplit(strsplit(x, ".bil")[[1]][1], "_")[[1]][2]}))
monthtif.lst <- ifelse(x %in% c(11,12,1), "X1", ifelse(x %in% c(2,3,4), "X2", ifelse(x %in% c(5,6,7), "X3", "X4")))
monthtif.lst

## Convert to a SAGA GIS file:
for(i in 1:length(bil.name)){
  rsaga.geoprocessor(lib="io_gdal", module=0, param=list(GRIDS=set.file.extension(bil.name[i], ".sgrd"), FILES=bil.name[i]), show.output.on.console = FALSE)
}
unlink(bil.name)

for(k in levels(as.factor(monthtif.lst))){
    sel <- monthtif.lst %in% k
    if(is.na(file.info(paste("P", k, "WCL3a.tif", sep=""))$size)){
      try( rsaga.geoprocessor(lib="geostatistics_grid", module=4, param=list(GRIDS=paste(set.file.extension(bil.name[sel], ".sgrd"), collapse=";"), MEAN=paste("P", k, "WCL3a.sgrd", sep=""))) ) # MIN="tmp.sgrd", MAX="tmp.sgrd", STDDEV="tmp.sgrd"
      ## convert to a geotiff...
      unlink("tmp.tif")
      #system(paste(gdal_translate, " ", "P", k, "WCL3a.sdat tmp.tif -ot \"Int16\" -a_srs \"", crs, "\"", sep=""))
      rsaga.geoprocessor(lib="io_gdal", module=2, param=list(GRIDS=paste("P", k, "WCL3a.sgrd", sep=""), FILE="tmp.tif"), show.output.on.console = FALSE)
      system(paste(gdalwarp, " tmp.tif ", "P", k, "WCL3a.tif", " -dstnodata \"-9999\" -ot \"Int16\" -r near -te -180 -90 180 90 -tr ", 1/120, " ", 1/120, sep=""))
    }
}

for(k in levels(as.factor(monthtif.lst))){
  tifout <- paste("P", k, "WCL3a.tif", sep="")
  tifoutr <- paste("P", k, "WCL2a.tif", sep="")
  ## 2.5 km:
  if(!file.exists(tifoutr)){
    system(paste(gdalwarp, ' ', tifout, ' ', tifoutr, ' -r bilinear -dstnodata -32767 -te -180 -90 180 90 -tr ', 1/40,' ', 1/40, sep=""))
  }
  ## 5 km:
  tifoutr2 <- paste("P", k, "WCL1a.tif", sep="")
  if(!file.exists(tifoutr2)){
    system(paste(gdalwarp, ' ', tifout, ' ', tifoutr2, ' -r bilinear -dstnodata -32767 -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
  }
  ## 20 km:
  tifoutr3 <- paste("P", k, "WCL0a.tif", sep="")
  if(!file.exists(tifoutr3)){
    system(paste(gdalwarp, ' ', tifoutr2, ' ', tifoutr3, ' -r bilinear -dstnodata -32767 -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))
  }
}

## compress:
for(k in levels(as.factor(monthtif.lst))){  
  tifoutr <- paste("P", k, "WCL", sep="")
  for(i in 0:3){
    outname = paste(tifoutr, i, 'a', sep="")
    system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), set.file.extension(outname, ".tif")))
    system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(outname, ".tif.gz"))
  }
}

## clean up:
unlink(paste(set.file.extension(bil.name, ".***")))



# end of script;