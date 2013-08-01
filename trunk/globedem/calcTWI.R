# title         : calcTWI.R
# purpose       : calculate SAGA Topographic Wetness Index derived using the DEMSRE3 (for land areas only);
# reference     : [http://worldgrids.org/doku.php?id=wiki:twisre3]
# producer      : Prepared by T. Hengl and Milan Kilibarda
# version       : 1
# inputs        : Global Relief Model based on SRTM 30+ and ETOPO DEM at 1/120 arcdeegre [http://worldgrids.org/maps/DEMSRE3a.tif.gz];
# outputs       : geotiff images projected in the "+proj=longlat +datum=WGS84" system;
# remarks 1     : First download and install FWtools [http://fwtools.maptools.org];
# remarks 2     : The calculation might miss some small islands;
# remarks 3     : SAGA version 2.0.8 was used;

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
  if(is.na(file.info(paste(continents[j,"conts"], '_TWISRE3a_ll.sdat', sep=""))$size)){
  unlink("tmp.tif")
  x = readGDAL(paste(continents[j,"conts"], '_LMTGSH3a.tif', sep=""))
  x$TWI <- readGDAL(paste(continents[j,"conts"], '_TWISRE3a.sdat', sep=""))$band1
  x$TWI <- ifelse(x@data[,1]==1, x$TWI, NA)
  writeGDAL(x["TWI"], "tmp.tif", mvFlag=-99999)
  ## back-transform to geo coordinates:
  unlink(paste(continents[j,"conts"], '_TWISRE3a_ll.sdat', sep=""))
  if(continents[j,"conts"]=="eu"){
    system(paste(gdalwarp, ' tmp.tif ', continents[j,"conts"], '_TWISRE3a_ll.tif -r cubic -dstnodata \"-99999\" -s_srs \"', continents[j,"csy"],'\" -t_srs \"+proj=longlat +datum=WGS84\" -tr ', 1/120,' ', 1/120, sep=""))
  } else {
    system(paste(gdalwarp, ' tmp.tif ', continents[j,"conts"], '_TWISRE3a_ll.tif -r near -dstnodata \"-99999\" -s_srs \"', continents[j,"csy"],'\" -t_srs \"+proj=longlat +datum=WGS84\" -tr ', 1/120,' ', 1/120, sep=""))
  }
}
}

## convert to SAGA GIS format:
for(j in 1:nrow(continents)){
  if(is.na(file.info(paste(continents[j,"conts"], '_TWISRE3a_ll.sdat', sep=""))$size)){
  system(paste(gdal_translate, ' ', continents[j,"conts"], '_TWISRE3a_ll.tif ', continents[j,"conts"], '_TWISRE3a_ll.sdat -a_nodata 0 -of \"SAGA\"', sep=""))
  #unlink(paste(continents[j,"conts"], '_TWISRE3a_ll.tif', sep=""))
  }
}

## merge maps: 
TWI.sgrd.list <- paste(continents$conts, "_TWISRE3a_ll.sgrd", sep="", collapse=";")
rsaga.geoprocessor(lib="grid_tools", module=3, param=list(GRIDS=TWI.sgrd.list, MERGED="TWISRE3a.sgrd", TYPE=7, INTERPOL=1, OVERLAP=0, MERGE_INFO_MESH_SIZE=1/120))

## filter missing pixels (due to reprojection problems)!
source('gridtiling.r')
system(paste(gdalwarp ,' -of SAGA -tr 1000 1000 -dstnodata -32767  -t_srs "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6378137.000000 +b=6356752.314245 +no_defs" /E:/TWI/TWISRE3a.sgrd',' ','TWISRE3a.sgrd',sep=""))
# make tiles for pixels filtering
gridtiling(file_path='twi.sgrd', # path to grid file in SAGA format, in this case DEM in sinusoidal projection, resolution 1 km                 
           tilename="twi", # prexif to be given to tile names                                                                                   
           tile_size=10000, # in cells                                                                                                          
           overlapping=100, # in cells                                                                                                          
           tiles_folder=paste(getwd(),'twi_tiles',sep='/'), # resulting folder                                                                  
           crs=CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6378137.000000 +b=6356752.314245 +no_defs"))                                           

# path to tiles
tiles1<- dir(path=paste(getwd(),"twi_tiles",sep='/'), pattern=glob2rx("*.sgrd"), full.names=T)   

# close gaps in each tile
for(j in 1:length(tiles1)) {      
  rsaga.geoprocessor(lib="grid_tools", 7,param=list(INPUT=tiles1[j], RESULT=tiles1[j], THRESHOLD=0.3 )  )
}

merged<-'twi_merg.sgrd' 

## merge maps: 
TWI.sgrd.list <- paste(tiles1, sep="", collapse=";")
rsaga.geoprocessor(lib="grid_tools", module=3, param=list(GRIDS=TWI.sgrd.list, MERGED=merged, TYPE=7, INTERPOL=1, OVERLAP=0, MERGE_INFO_MESH_SIZE=1/120))  

# assign projection to merged grid
rsaga.geoprocessor(lib="pj_proj4", module=0, 
                   param=list(GRIDS=merged,
                              CRS_PROJ4="+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6378137.000000 +b=6356752.314245 +no_defs"))
# Transfrom to WGS84             

rsaga.geoprocessor(lib="pj_proj4", module=7, 
                   param=list(SOURCE=merged,
                              SOURCE_PROJ="+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6378137.000000 +b=6356752.314245 +no_defs",
                              TARGET_PROJ="+proj=longlat +datum=WGS84",
                              INTERPOLATION='4',
                              TARGET_TYPE='0',
                              GET_USER_XMIN= '-180',
                              GET_USER_XMAX='180',
                              GET_USER_YMIN= '-90',
                              GET_USER_YMAX='90',
                              GET_USER_SIZE='0.008333',
                              GET_USER_GRID=merged  ))

# masking
rsaga.geoprocessor(lib="grid_tools", 24,param=list(GRID=merged, MASKED='TWISRE3a.sgrd', MASK='LMTGSH3a.sgrd' )  ) 
# reduce the storage size
rsaga.geoprocessor(lib="grid_calculus", 1, param=list(GRIDS='TWISRE3a.sgrd',          
                                                      RESULT='TWISRE3a.sgrd',
                                                      FORMULA='int((a-10)*10)'))   # FORMULA   ###############################

# Export Raster to GeoTIFF                                             
rsaga.geoprocessor(lib="io_gdal", 2, param=list(GRIDS='TWISRE3a.sgrd',          
                                                FILE='TWISRE3a_tmp.tif')) 


                                                                                         

# convert to geotifs (1 km):
system(paste(gdalwarp, ' TWISRE3a_tmp.tif -t_srs \"+proj=longlat +datum=WGS84\" -ot Byte -srcnodata \"-99999\" -dstnodata \"0\" TWISRE3a.tif -r near -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""))
# 2.5 km:
system(paste(gdalwarp, ' TWISRE3a_tmp.tif -t_srs \"+proj=longlat +datum=WGS84\" -ot Byte -srcnodata \"-99999\" -dstnodata \"0\" TWISRE2a.tif -r bilinear -te -180 -90 180 90 -tr ', 1/40,' ', 1/40, sep=""))
# 5 km:
system(paste(gdalwarp, ' TWISRE3a_tmp.tif -t_srs \"+proj=longlat +datum=WGS84\" -ot Byte -srcnodata \"-99999\" -dstnodata \"0\" TWISRE1a.tif -r bilinear -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
# 20 km:
system(paste(gdalwarp, ' TWISRE1a.tif -t_srs \"+proj=longlat +datum=WGS84\" -ot Byte -dstnodata \"-99999\" -dstnodata \"0\" TWISRE0a.tif -r bilinear -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))


## compress:  
for(i in 0:3){
    outname = paste("TWISRE", i, 'a', sep="")
    system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), set.file.extension(outname, ".tif")))
    system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(outname, ".tif.gz"))
  }

## clean up temp files:
#unlink(list.files(pattern=glob2rx("twi_*.*")))

# end of script;
