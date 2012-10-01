# title         : getGSHHS.R
# purpose       : land mask based on the GSHHS data;
# reference     : [http://worldgrids.org]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Aug 2012.
# inputs        : A Global Self-consistent, Hierarchical, High-resolution Shoreline Database (Shape file); 
# outputs       : geotiff images projected in the "+proj=longlat +ellps=WGS84" system and Google Maps system;
# remarks 1     : Data available at [http://www.ngdc.noaa.gov/mgg/shorelines/data/gshhs/version2.2.0/]; 

# ------------------------------------------------------------
# Initial settings and data download:
# ------------------------------------------------------------

library(RSAGA)
library(rgdal)
pixsize = 5000
g.csy <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"
g2.csy <- "+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=180.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs"
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"

# download files from server:
download.file("http://www.ngdc.noaa.gov/mgg/shorelines/data/gshhs/version2.2.0/GSHHS_shp_2.2.0.zip", destfile=paste(getwd(), "GSHHS_shp_2.2.0.zip", sep="/"))
download.file("ftp://ftp.soest.hawaii.edu/pwessel/gshhs/WDBII_shp_2.2.0.zip", "WDBII_shp_2.2.0.zip")
unzip(zipfile="GSHHS_shp_2.2.0.zip", exdir=getwd())
unzip(zipfile="WDBII_shp_2.2.0.zip", exdir=getwd())
ogrInfo("GSHHS_shp/h/GSHHS_h_L1.shp", "GSHHS_h_L1")

# even better choice is to use this map prepared by peterminton@yahool.com  [http://www.evs-islands.com/2007/11/data-global-land-mask-using-vectors.html]
# add a numeric value for classes:
lm.dbf <- read.dbf("Global GSHHS Land Mask.dbf")
lm.dbf$dbf$mask <- as.integer(lm.dbf$dbf$LAYER)
lm.levels <- levels(lm.dbf$dbf$LAYER)
lm.cols <- levels(lm.dbf$dbf$FILL_COLOR)
write.dbf(lm.dbf, "Global GSHHS Land Mask.dbf")
rm(lm.dbf)
ogrInfo("Global GSHHS Land Mask.shp", "Global GSHHS Land Mask")

# convert to raster (use the polygon ID):
rsaga.geoprocessor(lib="grid_gridding", module=0, param=list(USER_GRID="landmask.sgrd", INPUT="Global GSHHS Land Mask.shp", FIELD=21, MULTIPLE=0, TARGET=0, LINE_TYPE=1, GRID_TYPE=0, USER_SIZE=1/120, USER_XMIN=-180+1/120*0.5, USER_XMAX=180-1/120*0.5, USER_YMIN=-90+1/120*0.5, USER_YMAX=90-1/120*0.5)) # takes time!
# export to geotiff:
system(paste(gdal_translate, "landmask.sdat LMTGSH3a.tif -a_srs \"+proj=longlat +ellps=WGS84\" -ot \"Byte\""))
GDALinfo("LMTGSH3a.tif")
# resample to 1km, 2.5km, 5km and 20km:
system(paste(gdalwarp, "landmask.sdat LMTGSH1a.tif -s_srs \"+proj=longlat +ellps=WGS84\" -t_srs \"+proj=longlat +ellps=WGS84\" -ot \"Byte\" -r near -te -180 -90 180 90 -tr", 1/20, 1/20))
system(paste(gdalwarp, "landmask.sdat LMTGSH2a.tif -s_srs \"+proj=longlat +ellps=WGS84\" -t_srs \"+proj=longlat +ellps=WGS84\" -ot \"Byte\" -r near -te -180 -90 180 90 -tr", 1/40, 1/40))
system(paste(gdalwarp, "landmask.sdat LMTGSH0a.tif -s_srs \"+proj=longlat +ellps=WGS84\" -t_srs \"+proj=longlat +ellps=WGS84\" -ot \"Byte\" -r near -te -180 -90 180 90 -tr", 1/5, 1/5))

# ------------------------------------------------------------
# Prepare global land mask as a boolean map
# ------------------------------------------------------------

# convert to 0/1 map:
rsaga.geoprocessor("grid_calculus", 1, param=list(GRIDS="landmask.sgrd", RESULT="landmaskn.sgrd", FORMULA="lt(a,4)")) # takes time
## by default, SAGA puts everything into 32-bit Float, this is too big to handle:
rsaga.geoprocessor("grid_tools", 11, param=list(INPUT="landmaskn.sgrd", OUTPUT="landmaskn.sgrd", TYPE=1))

# resample to 1km, 2.5km, 5km and 20km:
system(paste(gdalwarp, "landmaskn.sdat LMBGSH3a.tif -s_srs \"+proj=longlat +ellps=WGS84\" -t_srs \"+proj=longlat +ellps=WGS84\" -r near -ot \"Byte\" -te -180 -90 180 90 -tr", 1/120, 1/120))
GDALinfo("LMBGSH3a.tif")
system(paste(gdalwarp, "landmaskn.sdat LMBGSH2a.tif -s_srs \"+proj=longlat +ellps=WGS84\" -t_srs \"+proj=longlat +ellps=WGS84\" -r bilinear -ot \"Float32\" -te -180 -90 180 90 -tr", 1/40, 1/40))  # TH: these maps now contain fine values indicating also proportion of land (>1x1 km in size) within 0.05 arcdegree blocks!
system(paste(gdalwarp, "landmaskn.sdat LMBGSH1a.tif -s_srs \"+proj=longlat +ellps=WGS84\" -t_srs \"+proj=longlat +ellps=WGS84\" -r bilinear -ot \"Float32\" -te -180 -90 180 90 -tr", 1/20, 1/20))
system(paste(gdalwarp, "landmaskn.sdat LMBGSH0a.tif -s_srs \"+proj=longlat +ellps=WGS84\" -t_srs \"+proj=longlat +ellps=WGS84\" -r bilinear -ot \"Float32\" -te -180 -90 180 90 -tr", 1/5, 1/5))

# Compress:
for(outname in c("LMTGSH3a.tif", "LMTGSH2a.tif", "LMTGSH1a.tif", "LMTGSH0a.tif", "LMBGSH3a.tif", "LMBGSH2a.tif", "LMBGSH1a.tif", "LMBGSH0a.tif")){
  if(is.na(file.info(paste(shortPathName(normalizePath(outdir)), paste(outname, "gz", sep="."), sep="\\"))$size)){
  system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), outname))
  system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
  unlink(set.file.extension(outname, ".tif.gz"))
}  # Compression takes > 15 mins
}


# ------------------------------------------------------------
# Global distance to continents map
# ------------------------------------------------------------

# reproject to Google Maps coordinate system
# (this needs to be done twice because globe is sphere!)
system(paste(gdalwarp, " landmask.sdat landmask5km_Lon0.sdat -of \"SAGA\" -s_srs \"+proj=longlat +ellps=WGS84\" -t_srs \"",g.csy , "\" -r near -ot \"Float32\" -te -20033925 -7705200 20036074 15544800 -tr 5000 5000", sep=""))
system(paste(gdalwarp, " landmask.sdat landmask5km_Lon180.sdat -of \"SAGA\" -s_srs \"+proj=longlat +ellps=WGS84\" -t_srs \"",g2.csy , "\" -r near -ot \"Float32\" -te -20033925 -7705200 20036074 15544800 -tr 5000 5000", sep=""))
# export to GeoTiff:
system(paste(gdal_translate, " landmask5km_Lon180.sdat land5kmB.tif -a_nodata 255 -a_srs \"",g2.csy ,"\" -ot \"Byte\"", sep=""))
system(paste(gdal_translate, " landmask5km_Lon0.sdat land5kmA.tif -a_nodata 255 -a_srs \"",g.csy ,"\" -ot \"Byte\"", sep=""))

# Convert raster map to polygons:
rsaga.geoprocessor(lib="shapes_grid", module=6, param=list(GRID="landmask5km_Lon180.sgrd", SHAPES="landmask5km_Lon180.shp", CLASS_ALL=1))
rsaga.geoprocessor(lib="shapes_grid", module=6, param=list(GRID="landmask5km_Lon0.sgrd", SHAPES="landmask5km_Lon0.shp", CLASS_ALL=1))

# mask out islands smaller than 100 x 100 km:
landmask <- readShapePoly("GSHHS_shp/h/GSHHS_h_L1.shp", repair=T, force_ring=T)
str(landmask, max.level=2)
# get the size of polygons:
landmask$area.pol <- sapply(slot(landmask, "polygons"), slot, "area")
# mask out the 'islands' - anything smaller than 1000x1000 km:
continents <- landmask[landmask$area.pol>10,]
proj4string(continents) <- CRS("+proj=longlat +ellps=WGS84")
writePolyShape(continents, "continents.shp")

# reproject:
rsaga.geoprocessor(lib="pj_proj4", 0, param=list(SOURCE_PROJ="\"+proj=longlat +datum=WGS84\"", TARGET_PROJ=paste('"', g.csy, '"', sep=""), SOURCE="continents.shp", TARGET="continents_Lon0.shp"))
# convert to raster:
rsaga.geoprocessor(lib="grid_gridding", module=0, param=list(GRID="continents_Lon0.sgrd", INPUT="continents_Lon0.shp", FIELD=2, TARGET_TYPE=0, LINE_TYPE=0, USER_CELL_SIZE=pixsize, USER_X_EXTENT_MIN=-20033925.36, USER_X_EXTENT_MAX=20036074.64, USER_Y_EXTENT_MIN=-7705200, USER_Y_EXTENT_MAX=15544800))
rsaga.geoprocessor(lib="grid_gridding", module=0, param=list(GRID="continents.sgrd", INPUT="continents.shp", FIELD=2, TARGET_TYPE=0, LINE_TYPE=0, USER_CELL_SIZE=0.05, USER_X_EXTENT_MIN=-179.975, USER_X_EXTENT_MAX=179.975, USER_Y_EXTENT_MIN=-89.975, USER_Y_EXTENT_MAX=89.975))
# fix the NA values: NODATA_VALUE	= 2
sgrd <- matrix((unlist(strsplit(readLines(file("continents.sgrd")), split="\t= "))), ncol=2, byrow=T)
sgrd
sgrd[12,2] <- "2"
write.table(sgrd, "continents.sgrd", quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t= ")
# derive buffer distance to continents:
## rsaga.geoprocessor(lib="grid_tools", module=10, param=list(SOURCE="continents5km.sgrd", DISTANCE="dcoast.sgrd", ALLOC="tmp.sgrd", BUFFER="tmp.sgrd", DIST=1e9, IVAL=pixsize))
# SAGA is too slow; ILWIS is more effcient:
system(paste(gdal_translate, " continents_Lon0.sdat continents_Lon0.mpr -of \"ILWIS\" -a_nodata 255 -ot \"Byte\"", sep=""))
system(paste(gdalwarp, " continents.sdat continents_Lon180.mpr -of \"ILWIS\" -s_srs \"+proj=longlat +ellps=WGS84\" -t_srs \"", g2.csy , "\" -r near -te -20033925 -7705200 20036074 15544800 -tr 5000 5000", sep=""))
# derived boolean maps:
shell(cmd=paste(ILWIS, " continents_Lon0B{dom=Bool.dom} = iff(continents_Lon0<1, ?, 1)", sep=""), wait=F)
shell(cmd=paste(ILWIS, " continents_Lon180B{dom=Bool.dom} = iff(continents_Lon180<1, ?, 1)", sep=""), wait=F)
# derive distance maps:
shell(cmd=paste(ILWIS, " dcoast_A.mpr{dom=value.dom;vr=0:100000000:1} = MapDistance(continents_Lon0B)", sep=""), wait=F)
shell(cmd=paste(ILWIS, " dcoast_B.mpr{dom=value.dom;vr=0:100000000:1} = MapDistance(continents_Lon180B)", sep=""), wait=F)
# combine the two maps:
system(paste(gdalwarp, " dcoast_B.mpr dcoast_Br.mpr -of \"ILWIS\" -s_srs \"", g2.csy , "\" -t_srs \"", g.csy , "\" -r bilinear -te -20033925 -7705200 20036074 15544800 -tr 5000 5000", sep=""))
shell(cmd=paste(ILWIS, " dist2con.mpr{dom=value.dom;vr=0:100000000:1} = min(dcoast_A, dcoast_Br)", sep=""), wait=F)
# resample to geographic coordinates:
system(paste(gdalwarp, " dist2con.sdat dist2con.tif -t_srs \"+proj=longlat +ellps=WGS84\" -s_srs \"", g.csy , "\" -r near -te -180 -90 180 90 -tr 0.05 0.05", sep=""))

# end of script;