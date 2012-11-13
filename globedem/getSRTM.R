# title         : getSRTM.R
# purpose       : download and resampling SRTM DEM 30 plus / ETOPO DEM;
# reference     : [https://code.google.com/p/worldgrids/source/browse/]
# producer      : Prepared by T. Hengl
# version       : 1
# inputs        : maps publicaly available at [http://dds.cr.usgs.gov/srtm/version2_1/SRTM30/] [http://www.ngdc.noaa.gov/mgg/bathymetry/relief.html] and [ftp://ftp.ngdc.noaa.gov/mgg/global/relief/];
# outputs       : geotiff images projected in the "+proj=longlat +datum=WGS84" system;
# remarks 1     : First download and install FWtools [http://fwtools.maptools.org];
# remarks 2     : The resulting GlobeDEM is a combination (average) between SRMT 30+ and ETOPO DEM;
# remarks 3     : The srtm30_plus [ftp://topex.ucsd.edu/pub/srtm30_plus/] showed some inconsistancies and hence was finally removed from the scriopt 

# -------------------------------------------
# Initial settings and data download:
# -------------------------------------------

library(RSAGA) 
library(rgdal)
library(raster)
library(maptools)
library(RCurl)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
gdaldem = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdaldem.exe"))))
gdalbuildvrt = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalbuildvrt.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"

## tiling system:
p.l <- expand.grid(KEEP.OUT.ATTRS=FALSE, lonmin=seq(-180,120,by=60), latmin=seq(-90,45,by=45))
p.u <- expand.grid(KEEP.OUT.ATTRS=FALSE, lonmax=seq(-120,180,by=60), latmax=seq(-45,90,by=45))
tiles <- cbind(p.l, p.u)
str(tiles)

# convert to a polygon map:
tiles.poly <- tiles
tiles.poly$lon <- (tiles$lonmax - tiles$lonmin)/2 + tiles$lonmin 
tiles.poly$lat <- (tiles$latmax - tiles$latmin)/2 + tiles$latmin
tiles.pnt <- tiles.poly
tiles.pnt$label <- paste("P", 1:nrow(tiles.pnt), sep="")
coordinates(tiles.poly) <- ~ lon+lat
gridded(tiles.poly) <- TRUE
tiles.poly <- rasterToPolygons(raster(tiles.poly))
writePolyShape(tiles.poly, "tiles.shp")
coordinates(tiles.pnt) <- ~ lon+lat
writePointsShape(tiles.pnt["label"], "tiles_label.shp")

## ETOPO1 Global Relief Model
download.file("http://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/grid_registered/georeferenced_tiff/ETOPO1_Ice_g_geotiff.zip", destfile=paste(getwd(), "ETOPO1_Ice_g_geotiff.zip", sep="/"))
unzip("ETOPO1_Ice_g_geotiff.zip")
# split in blocks:
for(j in 1:nrow(tiles)){
  # resample
  if(is.na(file.info(paste('ETOPO1_Ice_g_', j, '.sdat', sep=""))$size)){
  system(paste(gdalwarp, ' ETOPO1_Ice_g.tif -ot Int16 -dstnodata -32767 -t_srs \"+proj=longlat +datum=WGS84\" ', paste('ETOPO1_Ice_g_', j, '.sdat', sep=""), ' -of \"SAGA\" -r bilinear -te ', tiles[j,"lonmin"] , ' ', tiles[j,"latmin"], ' ', tiles[j,"lonmax"] ,' ', tiles[j,"latmax"] ,' -tr ', 1/120, ' ', 1/120, sep=""))
  }
}

## GEBCO data set [https://www.bodc.ac.uk/data/online_delivery/gebco/gebco_08_grid/#sid]
# GDALinfo("gebco_08.nc")
rsaga.geoprocessor(lib="io_grid", 4, param=list(GRID="gebco_08.sgrd", FILE_DATA="gebco_08.nc", NX=43200, NY=21600, DXY=0.008333333333333333, XMIN=-182.546, YMIN=-90, DATA_TYPE=3, TOPDOWN=1, BYTEORDER_BIG=1, NODATA=65536, ZFACTOR=1, UNIT="meter"))
## TH: for some unknown reason the coordinates have been shifted in x dimension by -2.546 degrees (?!)
# split in blocks:
for(j in 1:nrow(tiles)){
  # resample
  if(is.na(file.info(paste('gebco_08_', j, '.sdat', sep=""))$size)){
  system(paste(gdalwarp, ' gebco_08.sdat -ot Int16 -dstnodata -32767 -t_srs \"+proj=longlat +datum=WGS84\" ', paste('gebco_08_', j, '.sdat', sep=""), ' -of \"SAGA\" -r bilinear -te ', tiles[j,"lonmin"] , ' ', tiles[j,"latmin"], ' ', tiles[j,"lonmax"] ,' ', tiles[j,"latmax"] ,' -tr ', 1/120, ' ', 1/120, sep=""))
  # mask areas that are not of interest:
  tmp <- readGDAL(paste('gebco_08_', j, '.sdat', sep=""))
  tmp$band1 <- ifelse(tmp$band1>0, NA, tmp$band1)
  writeGDAL(tmp[1], paste('gebco_08_', j, '.sdat', sep=""), 'SAGA', mvFlag=-32767, type='Int16')
  }
}

# download SRTM DEM (1 km):
## HIR: use the offical version from the USGS site? [http://dds.cr.usgs.gov/srtm/version2_1/SRTM30/]
download.file("ftp://ftp.ntsg.umt.edu/pub/data/SRTM30_1km/Merged/SRTM30_merge.hdr", destfile=paste(getwd(), "SRTM30_merge.hdr", sep="/"))
download.file("ftp://ftp.ntsg.umt.edu/pub/data/SRTM30_1km/Merged/SRTM30_merge.int16", destfile=paste(getwd(), "SRTM30_merge.int16", sep="/"))
GDALinfo("SRTM30_merge.int16")
system(paste(gdalwarp, ' SRTM30_merge.int16 -t_srs \"+proj=longlat +datum=WGS84\" SRTM30_merge.tif -r bilinear -te -180 -90 180 90 -tr ', 1/120, ' ', 1/120, sep=""))
# split in blocks:
for(j in 1:nrow(tiles)){
  # resample
  if(is.na(file.info(paste('SRTM30_merge_', j, '.sdat', sep=""))$size)){
  unlink(paste('SRTM30_merge_', j, '.tif', sep=""))
  system(paste(gdalwarp, ' SRTM30_merge.tif -ot Int16 -dstnodata -32767 -t_srs \"+proj=longlat +datum=WGS84\" ', paste('SRTM30_merge_', j, '.sdat', sep=""), ' -of \"SAGA\" -r bilinear -te ', tiles[j,"lonmin"] , ' ', tiles[j,"latmin"], ' ', tiles[j,"lonmax"] ,' ', tiles[j,"latmax"] ,' -tr ', 1/120, ' ', 1/120, sep=""))
  # mask areas that are not of interest:
  tmp <- readGDAL(paste('SRTM30_merge_', j, '.sdat', sep=""))
  tmp$band1 <- ifelse(tmp$band1==0, NA, tmp$band1)
  # oceans contain all empty pixels:
  if(!is.logical(tmp$band1)){
    tmp$band1 <- as.integer(tmp$band1)
  }
  writeGDAL(tmp[1], paste('SRTM30_merge_', j, '.sdat', sep=""), 'SAGA', mvFlag=-32767, type='Int16')
  unlink(paste('SRTM30_merge_', j, '.tif', sep=""))
  }
}
rm(tmp)

# merge the three DEMs:
for(j in 1:nrow(tiles)){
  if(is.na(file.info(paste("globedem_", j, ".sgrd", sep=""))$size)){
  rsaga.geoprocessor(lib="grid_tools", module=3, param=list(GRIDS=paste("SRTM30_merge_", j, ".sgrd;ETOPO1_Ice_g_", j, ".sgrd;gebco_08_", j, ".sgrd", sep=""), GRID_TARGET=paste("globedem_", j, ".sgrd", sep=""), MERGED=paste("globedem_", j, ".sgrd", sep=""), TYPE=4, INTERPOL=1, OVERLAP=1, MERGE_INFO_MESH_SIZE=1/120))
  }
}

# mosaic/average two DEMs and create a complete DEM
unlink("globedem.vrt")
system(paste(gdalbuildvrt, "globedem.vrt", paste("globedem_", 1:nrow(tiles), ".sdat", sep="", collapse=" ")))
# convert to geotif (1 km):
unlink("DEMSRE3a.tif")
system(paste(gdalwarp, "globedem.vrt -t_srs \"+proj=longlat +datum=WGS84\" DEMSRE3a.tif -r near -te -180 -90 180 90 -tr", 1/120, 1/120))
GDALinfo("DEMSRE3a.tif")
# 2.5 km:
system(paste(gdalwarp, "DEMSRE3a.tif DEMSRE2a.tif -r bilinear -te -180 -90 180 90 -tr", 1/40, 1/40))
# 5 km:
system(paste(gdalwarp, "DEMSRE3a.tif DEMSRE1a.tif -r bilinear -te -180 -90 180 90 -tr", 1/20, 1/20))
# 20 km:
system(paste(gdalwarp, "DEMSRE1a.tif DEMSRE0a.tif -r bilinear -te -180 -90 180 90 -tr", 1/5, 1/5))

# Tile the 2.5 km DEMs so they can be uploaded to WorldGrids.org:
b.l <- expand.grid(KEEP.OUT.ATTRS=FALSE, lonmin=seq(-180,0,by=180), latmin=seq(-90,0,by=90))
b.u <- expand.grid(KEEP.OUT.ATTRS=FALSE, lonmax=seq(0,180,by=180), latmax=seq(0,90,by=90))
btiles <- cbind(b.l, b.u)
for(j in 1:nrow(btiles)){
   outname = paste('DEMSRE2a_B', j, '.tif', sep="")
   system(paste(gdalwarp, ' globedem.vrt -t_srs \"+proj=longlat +datum=WGS84\" ', outname, ' -r bilinear -te ', btiles[j,"lonmin"] , ' ', btiles[j,"latmin"], ' ', btiles[j,"lonmax"] ,' ', btiles[j,"latmax"] ,' -tr ', 1/40, ' ', 1/40, sep=""))
   if(is.na(file.info(paste(shortPathName(normalizePath(outdir)), paste(outname, "gz", sep="."), sep="\\"))$size)){
    system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), outname))
    system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir))))
   }
   unlink(outname)
   unlink(set.file.extension(outname, ".tif.gz"))
}

# 1 km data:
for(j in 1:nrow(tiles)){
   outname = paste('DEMSRE3a_P', j, '.tif', sep="")
   system(paste(gdalwarp, ' globedem.vrt -t_srs \"+proj=longlat +datum=WGS84\" ', outname, ' -r near -te ', tiles[j,"lonmin"] , ' ', tiles[j,"latmin"], ' ', tiles[j,"lonmax"] ,' ', tiles[j,"latmax"] ,' -tr ', 1/120, ' ', 1/120, sep=""))
   if(is.na(file.info(paste(shortPathName(normalizePath(outdir)), paste(outname, "gz", sep="."), sep="\\"))$size)){
    # maximum compression possible -- THIS TAKES A BIT MORE TIME!:
    system(paste("7za a", "-tgzip -mx9", set.file.extension(outname, ".tif.gz"), outname))
    system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir))))
   }
   unlink(outname)
   unlink(set.file.extension(outname, ".tif.gz"))
}


# -------------------------------------------
# Derive basic DEM parameters:
# -------------------------------------------

for(i in 0:3){
  if(is.na(file.info(paste('SLPSRT', i, 'a.tif', sep=""))$size)){
    unlink("slope.tif")
    # derive slope map:
    system(paste(gdaldem, ' slope DEMSRE', i, 'a.tif slope.tif -s 111120 -p', sep=""))  # takes > 5 mins
    # round the numbers to one 100/255 percent:
    system(paste(gdal_translate, ' slope.tif SLPSRT', i, 'a.tif -ot Byte -scale 0 100 0 254', sep=""))
  }
}

# system(paste(gdaldem, "TRI DEMSRE3a.tif TRI1km.tif -s 111120"))
# system(paste(gdaldem, "TPI DEMSRE3a.tif TPI1km.tif -s 111120"))
# system(paste(gdaldem, "roughness DEMSRE3a.tif RGHSRE3a.tif -s 111120"))

# TO-DO: global SAGA TWI, solar insolation, and Valley depth



# -------------------------------------------
# Compress produced maps:
# -------------------------------------------

for(outname in c("DEMSRE0a.tif", "DEMSRE1a.tif", "DEMSRE2a.tif", "DEMSRE3a.tif", "SLPSRT0a.tif", "SLPSRT1a.tif", "SLPSRT2a.tif", "SLPSRT3a.tif")){
  if(is.na(file.info(paste(shortPathName(normalizePath(outdir)), paste(outname, "gz", sep="."), sep="\\"))$size)){
  system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), outname))
  system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
  unlink(set.file.extension(outname, ".tif.gz"))
}  # Compression takes > 15 mins
}

# Clean-up:
rm(tmp)
unlink("gebco_08_*.*")
unlink("ETOPO1_Ice_g_*.*")
unlink("SRTM30_merge_*.*")
unlink("slope5km.tif")
unlink("slope1km.tif")

# end script;