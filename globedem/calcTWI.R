# title         : calcTWI.R
# purpose       : calculate SAGA Topographic Wetness Index derived using the DEMSRE3 (for land areas only);
# reference     : [http://worldgrids.org/doku.php?id=wiki:twisre3]
# producer      : Prepared by T. Hengl and Milan Kilibarda
# version       : 1
# inputs        : Global Relief Model based on SRTM 30+ and ETOPO DEM at 1/120 arcdeegre [http://worldgrids.org/maps/DEMSRE3a.tif.gz];
# outputs       : geotiff images projected in the "+proj=longlat +datum=WGS84" system;
# remarks 1     : First download and install FWtools [http://fwtools.maptools.org];
# remarks 2     : The calculation might miss some small islands;

library(RSAGA) 
library(rgdal)
library(raster)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
# download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
# unzip("7za920.zip")
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"

## download landmask and global DEM:
if(is.na(file.info("DEMSRE3a.tif")$size)|is.na(file.info("continents1km.zip")$size)){
  download.file("http://worldgrids.org/lib/exe/fetch.php?media=continents1km.zip", "continents1km.zip")
  system(paste("7za e continents1km.zip"))
  download.file("http://worldgrids.org/maps/DEMSRE3a.tif.gz", "DEMSRE3a.tif.gz")
  system(paste("7za e DEMSRE3a.tif.gz"))
}

load("continents.rda")
## tile and reproject land mask and DEM per continent:
for(j in 1:nrow(continents)){
  if(is.na(file.info(paste(continents[j,"conts"], '_DEMSRE3a.sdat', sep=""))$size)){
    system(paste(gdalwarp, ' DEMSRE3a.tif -t_srs \"', continents[j,"csy"], '\" ', continents[j,"conts"], '_DEMSRE3a.sdat -ot Int16 -of \"SAGA\" -r bilinear -te ', continents[j,"xmin"],' ', continents[j,"ymin"],' ', continents[j,"xmax"],' ', continents[j,"ymax"], ' -tr 1000 1000', sep=""))
  }
}

## for each tile, derive TWI:
for(j in 1:nrow(continents)){
  if(is.na(file.info(paste(continents[j,"conts"], "_TWISRE3a.sgrd", sep=""))$size)){
  rsaga.wetness.index(paste(continents[j,"conts"], "_DEMSRE3a.sgrd", sep=""), out.wetness.index=paste(continents[j,"conts"], "_TWISRE3a.sgrd", sep=""))
}
} ## This takes about 6 hours per area or almost 2 days of computing!!!

## mask out land areas and resample back to geographical coordinates:
for(j in 1:nrow(continents)){
  if(is.na(file.info(paste(continents[j,"conts"], '_TWISRE3a_ll.tif', sep=""))$size)){
  unlink("tmp.tif")
  x = readGDAL(paste(continents[j,"conts"], '_LMTGSH3a.tif', sep=""))
  x$TWI <- readGDAL(paste(continents[j,"conts"], '_TWISRE3a.sdat', sep=""))$band1
  x$TWI <- ifelse(x@data[,1]==1, x$TWI, NA)
  writeGDAL(x["TWI"], "tmp.tif", mvFlag=-99999)
  ## back-transform to geo coordinates:
  unlink(paste(continents[j,"conts"], '_TWISRE3a_ll.sdat', sep=""))
  if(j==5){
    system(paste(gdalwarp, ' tmp.tif ', continents[j,"conts"], '_TWISRE3a_ll.tif -r cubic -dstnodata \"-99999\" -s_srs \"', continents[j,"csy"],'\" -t_srs \"+proj=longlat +datum=WGS84\" -tr ', 1/120,' ', 1/120, sep=""))
  } else {
    system(paste(gdalwarp, ' tmp.tif ', continents[j,"conts"], '_TWISRE3a_ll.tif -r near -dstnodata \"-99999\" -s_srs \"', continents[j,"csy"],'\" -t_srs \"+proj=longlat +datum=WGS84\" -tr ', 1/120,' ', 1/120, sep=""))
  }
}
}

## convert to SAGA GIS format:
for(j in 1:nrow(continents)){
  system(paste(gdal_translate, ' ', continents[j,"conts"], '_TWISRE3a_ll.tif ', continents[j,"conts"], '_TWISRE3a_ll.sdat -a_nodata 0 -of \"SAGA\"', sep=""))
  unlink(paste(continents[j,"conts"], '_TWISRE3a_ll.tif', sep=""))
}

## merge maps: 
TWI.sgrd.list <- paste(continents$conts, "_TWISRE3a_ll.sgrd", sep="", collapse=";")
rsaga.geoprocessor(lib="grid_tools", module=3, param=list(GRIDS=TWI.sgrd.list, MERGED="TWISRE3a.sgrd", TYPE=7, INTERPOL=1, OVERLAP=0, MERGE_INFO_MESH_SIZE=1/120))

# convert to geotifs (1 km):
system(paste(gdalwarp, ' TWISRE3a.sgrd -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-99999\" TWISRE3a.tif -r near -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""))
# 2.5 km:
system(paste(gdalwarp, ' TWISRE3a.sgrd -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-99999\" TWISRE2a.tif -r bilinear -te -180 -90 180 90 -tr ', 1/40,' ', 1/40, sep=""))
# 5 km:
system(paste(gdalwarp, ' TWISRE3a.sgrd -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-99999\" TWISRE1a.tif -r bilinear -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
# 20 km:
system(paste(gdalwarp, ' TWISRE1a.tif -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-99999\" TWISRE0a.tif -r bilinear -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))


## compress:  
for(i in 0:3){
    outname = paste("TWISRE", i, 'a', sep="")
    system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), set.file.extension(outname, ".tif")))
    system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(outname, ".tif.gz"))
  }
}

# clean up temp files:
unlink(list.files(pattern=glob2rx("LA*MOD3a_*.*")))

# end of script;
