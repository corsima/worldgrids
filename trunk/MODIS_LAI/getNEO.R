
# title         : getNEO.R
# purpose       : download and resampling of MODIS Leaf Area Index 0.1 arcdegree Monthly images;
# reference     : [https://code.google.com/p/worldgrids/source/browse/]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, June 2012.
# inputs        : images available at [http://neo.sci.gsfc.nasa.gov]; 
# outputs       : geotiff images projected in the "+proj=longlat +datum=WGS84" system;
# remarks 1     : Description of the data available at [https://lpdaac.usgs.gov/lpdaac/products/modis_products_table/leaf_area_index_fraction_of_photosynthetically_active_radiation/8_day_l4_global_1km/mcd15a2]; 
# remarks 2     : First download and install FWtools [http://fwtools.maptools.org];
# remarks 3     : NEO is not friendly for automated download (images need to be downloaded in several iterations);  

# ------------------------------------
# Initial settings and data download:
# ------------------------------------

library(rgdal)
library(RCurl)
library(maptools)
library(RSAGA)
library(XML)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
gdaldem = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdaldem.exe"))))
gdalbuildvrt = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalbuildvrt.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
si <- Sys.info()
outdir <- "D:/WORLDGRIDS/maps"

# get scene IDs and dates (begin):
ID.list <- as.list(rep(NA, 13))
Date.list <- as.list(rep(NA, 13))
for(i in 1:length(ID.list)){
 con <- url(paste("http://neo.sci.gsfc.nasa.gov/Search.html?minLat=&sse=&endDay=&maxLon=&ssn=&palette=rgb&ssw=&minLon=&sss=&startDay=&datasetId=MOD15A2_M_LAI&sourceDatasetId=&endYear=&format=TIFF&h=&size=full&w=&duration=day&maxLat=&pg=", i, "&startYear=&coverage=global&startMonth=&endMonth=", sep=""))
 tmp <- as.vector(scan(con, what="character"))
 ID.list[[i]] <- tmp[(which(tmp=="imageId")+2)[-1]]
 Date.list[[i]] <- as.Date(tmp[(which(tmp=="startDate")+2)[-1]])
 close(con)
# download images one by one: (open firefox and set the automatic download to the current directory)
for(j in 1:length(ID.list[[i]])){
  browseURL(paste("http://neo.sci.gsfc.nasa.gov/RenderData?si=", ID.list[[i]][j], "&cs=gs&format=TIFF", sep=""), browser="C:/Program Files (x86)/Mozilla Firefox/firefox.exe") 
} 
}

# download the missing images:
LAI.list <- dir(path=getwd(), pattern=glob2rx("*gs-167772161.0.TIFF"), full.names=FALSE)
LAI.name <- rep(NA, length(LAI.list)) 
for(i in 1:length(LAI.list)){
 LAI.name[[i]] <- strsplit(LAI.list[[i]], "gs-167772161.0.TIFF")[[1]]
}
for(i in 1:13){
for(j in 1:length(ID.list[[i]])){
  if(!(ID.list[[i]][j] %in% LAI.name)) {
  browseURL(paste("http://neo.sci.gsfc.nasa.gov/RenderData?si=", ID.list[[i]][j], "&cs=gs&format=TIFF", sep=""), browser="C:/Program Files (x86)/Mozilla Firefox/firefox.exe") 
  Sys.sleep(5)
}}}

# ------------------------------------
# Derivation of the mean and sd LAI for globe
# ------------------------------------

# import all images to SAGA format:
for(i in 1:13){
for(j in 1:length(ID.list[[i]])){
  fname <- paste("LAI_", gsub("-", "_", paste(Date.list[[i]][j])), sep="")
  system(paste(gdal_translate, paste(ID.list[[i]][j], "gs-167772161.0.TIFF", sep=""), set.file.extension(fname, ".sdat"), "-of \"SAGA\" -a_nodata 255"), show.output.on.console=FALSE)
}}

# derive mean and sd using all LAI images:
LAI.sgrd.list <- dir(path=getwd(), pattern=glob2rx("*.sgrd"), full.names=FALSE)   # 123 monthly maps
rsaga.geoprocessor(lib="geostatistics_grid", module=5, param=list(GRIDS=paste(LAI.sgrd.list, collapse=";"), MEAN="LAIm.sgrd", MIN="LAImin.sgrd", MAX="LAImax.sgrd", VAR="tmp.sgrd", STDDEV="LAIs.sgrd", STDDEVLO="tmp.sgrd", STDDEVHI="tmp.sgrd"))
# optional: derive mean LAI for each month
# Downscale to 5 km resolution:
unlink("LAIm_5km.sdat")
system(paste(gdalwarp, "LAIm.sdat -t_srs \"+proj=longlat +datum=WGS84\" LAIm_5km.sdat -srcnodata \"-99999\" -dstnodata \"-99999\" -of \"SAGA\" -r cubicspline -te -180 -90 180 90 -tr", 6/120, 6/120))
unlink("LAIs_5km.sdat")
system(paste(gdalwarp, "LAIs.sdat -t_srs \"+proj=longlat +datum=WGS84\" LAIs_5km.sdat -srcnodata \"-99999\" -dstnodata \"-99999\" -of \"SAGA\" -r cubicspline -te -180 -90 180 90 -tr", 6/120, 6/120))

# download mask map and fix some missing pixels:
download.file("http://spatial-analyst.net/worldmaps/landmask.zip", destfile=paste(getwd(), "landmask.zip", sep="/"))  # update with WorldGrids.org zip file!
unzip(zipfile="landmask.zip", exdir=getwd())
grids5km <- readGDAL("landmask.tif")
grids5km$LAIm <- readGDAL("LAIm_5km.sdat")$band1
grids5km$LAIs <- readGDAL("LAIs_5km.sdat")$band1
grids5km$LAImf <- round(ifelse(grids5km$band1==-32767, 0, ifelse(is.na(grids5km$LAIm), 0, grids5km$LAIm)), 1)
proj4string(grids5km) <- CRS("+proj=longlat +datum=WGS84")
# write to GeoTiff format:
writeGDAL(grids5km["LAImf"], "LAIMOD1a.tif", "GTiff", type="Byte")
writeGDAL(grids5km["LAIs"], "LASMOD1a.tif", "GTiff", type="Byte")
# soil mask:
grids5km$soilmask <- ifelse(grids5km$band1==-32767, 0, ifelse(grids5km$LAIm>0, 1, 0)) 
writeGDAL(grids5km["soilmask"], "SMKMOD1a.tif", "GTiff", type="Byte")
# check consistency:
GDALinfo("SMKMOD1a.tif")

for(outname in c("LAIMOD1a.tif", "LASMOD1a.tif", "SMKMOD1a.tif")){
  system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), outname))
  system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
  unlink(set.file.extension(outname, ".tif.gz"))
}

# Clean-up:
unlink("LAIs_5km.*")
rm(grids5km)

# end of script;