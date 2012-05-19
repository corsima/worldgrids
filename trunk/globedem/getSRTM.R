
# title         : getSRTM.R
# purpose       : download and resampling SRTM DEM 30 plus / ETOPO DEM;
# reference     : [https://code.google.com/p/worldgrids/source/browse/]
# producer      : Prepared by T. Hengl
# version       : 1
# inputs        : maps publicaly available at [ftp://topex.ucsd.edu/pub/srtm30_plus/srtm30/data/];
# outputs       : geotiff images projected in the "+proj=longlat +datum=WGS84" system;
# remarks 1     : First download and install FWtools [http://fwtools.maptools.org];
# remarks 2     : The resulting GlobeDEM is a combination (average) between SRMT 30+ and ETOPO DEM;  

# -------------------------------------------
# Initial settings and data download:
# -------------------------------------------

library(RSAGA) 
library(rgdal)
library(RCurl)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
gdaldem = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdaldem.exe"))))
gdalbuildvrt = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalbuildvrt.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
si <- Sys.info()
outdir <- "D:/WORLDGRIDS/maps"

# location of maps:
URL <- "ftp://topex.ucsd.edu/pub/srtm30_plus/srtm30/data/"
# list all tiles available:
items <- strsplit(getURL(URL), "\n")[[1]]
# convert to a character vector:
srtm.list <- items[grep(items, pattern=".srtm")]
srtm.list <- unlist(lapply(strsplit(srtm.list, ' '), function(x){x[length(x)]}))
srtm.list <- unlist(lapply(strsplit(srtm.list, '\r'), function(x){x[1]}))
str(srtm.list)

# download tiles from the server:
for(j in srtm.list){
  download.file(paste(URL, j, sep=""), destfile=paste(getwd(), j, sep="/")) 
}
# GDALinfo("e020n40.Bathymetry.srtm")
## Not a GDAL supported format!

# create a mosaic:
rsaga.geoprocessor(lib="io_grid", 9, param=list(GRID="srtmdem1km.sgrd", PATH=getwd(), XMIN=-180, XMAX=180, YMIN=-90, YMAX=90, TILE_PATH=paste(srtm.list, collate=";", sep="")))  # takes time and memory!
# GDALinfo("srtmdem1km.sdat")
## Split into two tiles otherwise SAGA GIS faces memory limit problems:
system(paste(gdalwarp, " srtmdem1km.sdat -t_srs \"+proj=longlat +datum=WGS84\" srtmdem1km_a.sdat -of \"SAGA\" -srcnodata \"-9999\" -dstnodata \"-32767\" -r bilinear -te -180 -90 0 90 -tr ", 1/120, " ", 1/120, sep=""))
system(paste(gdalwarp, " srtmdem1km.sdat -t_srs \"+proj=longlat +datum=WGS84\" srtmdem1km_b.sdat -of \"SAGA\" -srcnodata \"-9999\" -dstnodata \"-32767\" -r bilinear -te 0 -90 180 90 -tr ", 1/120, " ", 1/120, sep=""))

# download ETOPO DEM (contains also elevations for <-60 S):
download.file("ftp://ftp.ngdc.noaa.gov/mgg/global/relief/ETOPO1/bedrock/grid_registered/georeferenced_tiff/ETOPO1_Ice_g_geotiff.zip", destfile=paste(getwd(), "ETOPO1_Ice_g_geotiff.zip", sep="/"))
unzip(zipfile="ETOPO1_Ice_g_geotiff.zip", exdir=getwd())
GDALinfo("ETOPO1_Ice_g.tif")
system(paste(gdalwarp, " ETOPO1_Ice_g.tif -t_srs \"+proj=longlat +datum=WGS84\" ETOPO1_a.sdat -of \"SAGA\" -r bilinear -te 0 -90 180 90 -tr ", 1/120, " ", 1/120, sep=""))
system(paste(gdalwarp, " ETOPO1_Ice_g.tif -t_srs \"+proj=longlat +datum=WGS84\" ETOPO1_b.sdat -of \"SAGA\" -r bilinear -te -180 -90 0 90 -tr ", 1/120, " ", 1/120, sep=""))

# mosaic/average two DEMs and create a complete DEM:
rsaga.geoprocessor(lib="grid_tools", module=3, param=list(GRIDS="srtmdem1km_a.sgrd;ETOPO1_b.sgrd", GRID_TARGET="globedem_a.sgrd", MERGED="globedem_a.sgrd", TYPE=4, INTERPOL=1, OVERLAP=0, MERGE_INFO_MESH_SIZE=1/120))
rsaga.geoprocessor(lib="grid_tools", module=3, param=list(GRIDS="srtmdem1km_b.sgrd;ETOPO1_a.sgrd", GRID_TARGET="globedem_b.sgrd", MERGED="globedem_b.sgrd", TYPE=4, INTERPOL=1, OVERLAP=0, MERGE_INFO_MESH_SIZE=1/120))  ## memory consuming / takes ca 10 mins!

# convert to geotif (1 km):
system(paste(gdalbuildvrt, "globedem.vrt globedem_a.sdat globedem_b.sdat"))
system(paste(gdalwarp, "globedem.vrt -t_srs \"+proj=longlat +datum=WGS84\" DEMSRE3a.tif -r near -te -180 -90 180 90 -tr", 1/120, 1/120))
GDALinfo("DEMSRE3a.tif")

# convert to geotif (5 km):
system(paste(gdalwarp, "DEMSRE3a.tif DEMSRE1a.tif -r bilinear -te -180 -90 180 90 -tr", 6/120, 6/120))

# -------------------------------------------
# Derive basic DEM parameters:
# -------------------------------------------

system(paste(gdaldem, "slope DEMSRE3a.tif slope1km.tif -s 111120 -p"))  # takes > 5 mins
# round the numbers to one 100/255 percent:
system(paste(gdal_translate, "slope1km.tif SLPSRT3a.tif -ot Byte -scale 0 100 0 254"))
# system(paste(gdaldem, "TRI DEMSRE3a.tif TRI1km.tif -s 111120"))
# system(paste(gdaldem, "TPI DEMSRE3a.tif TPI1km.tif -s 111120"))
# system(paste(gdaldem, "roughness DEMSRE3a.tif RGHSRE3a.tif -s 111120"))
system(paste(gdaldem, "slope DEMSRE1a.tif slope5km.tif -s 111120 -p"))
# round the numbers to one 100/255 percent;
system(paste(gdal_translate, "slope5km.tif SLPSRT1a.tif -ot Byte -scale 0 100 0 254"))

# TO-DO: global SAGA TWI, solar insolation, and Valley depth


# -------------------------------------------
# Compress produced maps:
# -------------------------------------------

for(outname in c("DEMSRE1a.tif", "DEMSRE3a.tif", "SLPSRT1a.tif", "SLPSRT3a.tif")){
  system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), outname))
  system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
  unlink(set.file.extension(outname, ".tif.gz"))
}  # Compression takes > 15 mins

# Clean-up:
unlink("srtmdem1km_b.*")
unlink("srtmdem1km_a.*")
unlink("ETOPO1_a.*")
unlink("ETOPO1_b.*")
unlink("slope5km.tif")
unlink("slope1km.tif")
unlink("ETOPO1.*")

# end script;