# title         : getGLWD.R
# purpose       : download and resampling of Global Lakes and Wetlands Database;
# reference     : [http://worldgrids.org]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Aug 2012.
# inputs        : maps available for download from [http://www.worldwildlife.org/science/data/lakesandwetlands.html] and [http://wateriso.eas.purdue.edu/waterisotopes/pages/data_access/ArcGrids.html];
# outputs       : geotiff images projected in the "+proj=longlat +ellps=WGS84" system;
# remarks 1     : First download and install FWtools [http://fwtools.maptools.org];

# ------------------------------------------------------------
# Initial settings and data download:
# ------------------------------------------------------------

library(rgdal)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"

## download files from the server:
# download.file("http://www.worldwildlife.org/resources/media/glwd_31.zip", destfile=paste(getwd(), "glwd_31.zip", sep="/"))
# unzip(zipfile="glwd_31.zip", exdir=getwd())
# GDALinfo("glwd_31.img")

download.file("http://www.worldwildlife.org/science/data/WWFBinaryitem6599.zip", "WWFBinaryitem6599.zip")
unzip("WWFBinaryitem6599.zip")
GDALinfo("glwd_3")

# resample the 1 km, 2.5 km, 5 km and 20 km:
system(paste(gdalwarp, "glwd_3 -t_srs \"+proj=longlat +datum=WGS84\" GLWWWF3a.tif -ot \"Byte\" -r near -te -180 -90 180 90 -tr", 1/120, 1/120))
system(paste(gdalwarp, "glwd_3 -t_srs \"+proj=longlat +datum=WGS84\" GLWWWF2a.tif -ot \"Byte\" -r near -te -180 -90 180 90 -tr", 1/40, 1/40))
system(paste(gdalwarp, "glwd_3 -t_srs \"+proj=longlat +datum=WGS84\" GLWWWF1a.tif -ot \"Byte\" -r near -te -180 -90 180 90 -tr", 1/20, 1/20))
system(paste(gdalwarp, "glwd_3 -t_srs \"+proj=longlat +datum=WGS84\" GLWWWF0a.tif -ot \"Byte\" -r near -te -180 -90 180 90 -tr", 1/5, 1/5))

# Compress:
for(outname in c("GLWWWF3a.tif", "GLWWWF2a.tif", "GLWWWF1a.tif", "GLWWWF0a.tif")){
  if(is.na(file.info(paste(shortPathName(normalizePath(outdir)), paste(outname, "gz", sep="."), sep="\\"))$size)){
  system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), outname))
  system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
  unlink(set.file.extension(outname, ".tif.gz"))
}  # Compression takes > 15 mins
}


## MASK CLASSES OF INTEREST:

# read to SAGA sgrd format:
rsaga.geoprocessor(lib="io_gdal", module=0, param=list(GRIDS="glwd_31.sgrd", FILE="GLWWWF3a.tif")) 
# mask different classes:
rsaga.geoprocessor(lib="grid_calculus", module=1, param=list(INPUT="glwd_31.sgrd", RESULT="swamp.sgrd", FORMUL="ifelse(a=5,100,0)"))
rsaga.geoprocessor(lib="grid_calculus", module=1, param=list(INPUT="glwd_31.sgrd", RESULT="peatland.sgrd", FORMUL="ifelse(a=8,100,0)"))

for(j in c("swamp", "peatland")){
system(paste("C:\\PROGRA~2\\FWTOOL~1.7\\bin\\gdalwarp ", j, ".sdat -t_srs \"+proj=longlat +ellps=WGS84\" ", j, ".tif -r bilinear -te -180 -90 180 90 -tr 0.05 0.05 -ot \"Byte\" -dstnodata 255", sep=""))
}
 
# end of script;