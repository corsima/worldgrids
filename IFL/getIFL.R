# title         : getIFL.R
# purpose       : download and resampling of intact forest areas;
# reference     : [http://worldgrids.org]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Aug 2012.
# inputs        : shape file with borders of intact forests available for download from [http://www.intactforests.org/data.ifl.shp.html];
# outputs       : geotiff images projected in the "+proj=longlat +ellps=WGS84" system;
# remarks 1     : First download and install FWtools [http://fwtools.maptools.org];  


# ------------------------------------------------------------
# Initial settings and data download:
# ------------------------------------------------------------

library(rgdal)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_rasterize = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_rasterize.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"

# download files from the server:
download.file("http://www.intactforests.org/shp/world_ifl.zip", destfile=paste(getwd(), "world_ifl.zip", sep="/"))
unzip(zipfile="world_ifl.zip", exdir=getwd())
ogrInfo("world_ifl.shp", "world_ifl")  # this shape file is Large!

# download land mask map from the server:
download.file("http://worldgrids.org/lib/exe/fetch.php?media=lmbgsh3a.tif.gz", destfile=paste(getwd(), "lmbgsh3a.tif.gz", sep="/"))
system(paste("7za e lmbgsh3a.tif.gz"))
# mask out water bodies:
GDALinfo("LMBGSH3a.tif")
system(paste(gdal_translate, "LMBGSH3a.tif IFLGRE3a.tif -a_nodata 0"))
GDALinfo("IFLGRE3a.tif")

# system(paste(gdal_rasterize, "-l world_ifl -burn 1 world_ifl.shp IFLGRE3a.tif")) 
## takes too much time >40 mins!!!

# rasterize the shape file to 1 km resolution:
rsaga.geoprocessor(lib="grid_gridding", module=0, param=list(USER_GRID="mask.sgrd", INPUT="world_ifl.shp", FIELD=0, MULTIPLE=0, TARGET=0, LINE_TYPE=1, GRID_TYPE=0, USER_SIZE=1/120, USER_XMIN=-180+1/120*0.5, USER_XMAX=180-1/120*0.5, USER_YMIN=-90+1/120*0.5, USER_YMAX=90-1/120*0.5))
# convert to 0/1 mask
rsaga.geoprocessor("grid_calculus", 1, param=list(GRIDS="mask.sgrd", RESULT="mask_b.sgrd", FORMULA="eq(a,129)")) # takes time
## by default, SAGA puts everything into 32-bit Float, this is too big to handle:
rsaga.geoprocessor("grid_tools", 11, param=list(INPUT="mask_b.sgrd", OUTPUT="mask_b.sgrd", TYPE=1))
# invert numbers:
rsaga.geoprocessor("grid_calculus", 1, param=list(GRIDS="mask_b.sgrd", RESULT="mask_c.sgrd", FORMULA="eq(a,0)"))
rsaga.geoprocessor("grid_tools", 11, param=list(INPUT="mask_c.sgrd", OUTPUT="IFLGRE3a.sgrd", TYPE=1))
unlink("mask_b.***")
unlink("mask_c.***")
system(paste(gdal_translate, "IFLGRE3a.sdat IFLGRE3a.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\""))

# check validity:
GDALinfo("IFLGRE3a.tif")
# resample to 2.5, 5 and 20 km resolutions:
system(paste(gdalwarp, "IFLGRE3a.tif IFLGRE2a.tif -r bilinear -te -180 -90 180 90 -tr", 1/40, 1/40))
system(paste(gdalwarp, "IFLGRE3a.tif IFLGRE1a.tif -r bilinear -te -180 -90 180 90 -tr", 1/20, 1/20))
system(paste(gdalwarp, "IFLGRE3a.tif IFLGRE0a.tif -r bilinear -te -180 -90 180 90 -tr", 1/5, 1/5))

# Compress:
for(outname in c("IFLGRE0a.tif", "IFLGRE1a.tif", "IFLGRE2a.tif", "IFLGRE3a.tif")){
  if(is.na(file.info(paste(shortPathName(normalizePath(outdir)), paste(outname, "gz", sep="."), sep="\\"))$size)){
  system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), outname))
  system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
  unlink(set.file.extension(outname, ".tif.gz"))
}  # Compression takes > 15 mins
}


# end of script;
