# title         : getGSHHS.R
# purpose       : land mask based on the GSHHS data;
# reference     : [http://worldgrids.org]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Aug 2012.
# inputs        : A Global Self-consistent, Hierarchical, High-resolution Shoreline Database (Shape file); 
# outputs       : geotiff images projected in the "+proj=longlat +datum=WGS84" system and Google Maps system;
# remarks 1     : Data available at [http://www.ngdc.noaa.gov/mgg/shorelines/data/gshhs/version2.2.0/]; 

# ------------------------------------------------------------
# Initial settings and data download:
# ------------------------------------------------------------

library(RSAGA)
library(rgdal)
pixsize = 1000
## coordinate systems:
g.csy <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
g2.csy <- "+proj=robin +lon_0=120 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
g3.csy <- "+proj=robin +lon_0=-120 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
## another option - sinusoidal projection: "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
ogr2ogr = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/ogr2ogr.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
ILWIS <- "C:\\Ilwis3.4\\Ilwis30.exe -C"
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"

# download files from server:
download.file("http://www.ngdc.noaa.gov/mgg/shorelines/data/gshhs/version2.2.0/GSHHS_shp_2.2.0.zip", destfile=paste(getwd(), "GSHHS_shp_2.2.0.zip", sep="/"))
download.file("ftp://ftp.soest.hawaii.edu/pwessel/gshhs/WDBII_shp_2.2.0.zip", "WDBII_shp_2.2.0.zip")
unzip(zipfile="GSHHS_shp_2.2.0.zip", exdir=getwd())
unzip(zipfile="WDBII_shp_2.2.0.zip", exdir=getwd())
ogrInfo("GSHHS_shp/h/GSHHS_h_L1.shp", "GSHHS_h_L1")

# even better choice is to use this map prepared by peterminton@yahoo.com  [http://www.evs-islands.com/2007/11/data-global-land-mask-using-vectors.html]
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
system(paste(gdal_translate, "landmask.sdat LMTGSH3a.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\""))
GDALinfo("LMTGSH3a.tif")
# resample to 1km, 2.5km, 5km and 20km:
system(paste(gdalwarp, "landmask.sdat LMTGSH1a.tif -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -r near -te -180 -90 180 90 -tr", 1/20, 1/20))
system(paste(gdalwarp, "landmask.sdat LMTGSH2a.tif -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -r near -te -180 -90 180 90 -tr", 1/40, 1/40))
system(paste(gdalwarp, "landmask.sdat LMTGSH0a.tif -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -r near -te -180 -90 180 90 -tr", 1/5, 1/5))

# ------------------------------------------------------------
# Prepare global land mask as a boolean map
# ------------------------------------------------------------

# convert to 0/1 map:
rsaga.geoprocessor("grid_calculus", 1, param=list(GRIDS="landmask.sgrd", RESULT="landmaskn.sgrd", FORMULA="lt(a,4)")) # takes time
## by default, SAGA puts everything into 32-bit Float, this is too big to handle:
rsaga.geoprocessor("grid_tools", 11, param=list(INPUT="landmaskn.sgrd", OUTPUT="landmaskn.sgrd", TYPE=1))

# resample to 1km, 2.5km, 5km and 20km:
system(paste(gdalwarp, "landmaskn.sdat LMBGSH3a.tif -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"+proj=longlat +datum=WGS84\" -r near -ot \"Byte\" -te -180 -90 180 90 -tr", 1/120, 1/120))
GDALinfo("LMBGSH3a.tif")
system(paste(gdalwarp, "landmaskn.sdat LMBGSH2a.tif -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"+proj=longlat +datum=WGS84\" -r bilinear -ot \"Float32\" -te -180 -90 180 90 -tr", 1/40, 1/40))  # TH: these maps now contain fine values indicating also proportion of land (>1x1 km in size) within 0.05 arcdegree blocks!
system(paste(gdalwarp, "landmaskn.sdat LMBGSH1a.tif -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"+proj=longlat +datum=WGS84\" -r bilinear -ot \"Float32\" -te -180 -90 180 90 -tr", 1/20, 1/20))
system(paste(gdalwarp, "landmaskn.sdat LMBGSH0a.tif -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"+proj=longlat +datum=WGS84\" -r bilinear -ot \"Float32\" -te -180 -90 180 90 -tr", 1/5, 1/5))

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

## reproject to some equal areas system e.g. Sinusoidal coordinate system
## (this needs to be done three times to reduce effect of the projection system)
system(paste(gdalwarp, " LMBGSH3a.tif landmask_Lon0.mpr -of \"ILWIS\" -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"", g.csy , "\" -r near -tr 5000 5000", sep=""))
system(paste(gdalwarp, " LMBGSH3a.tif landmask_Lon120.mpr -of \"ILWIS\" -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"", g2.csy, "\" -r near -tr 5000 5000", sep=""))
system(paste(gdalwarp, " LMBGSH3a.tif landmask_LonM120.mpr -of \"ILWIS\" -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"", g3.csy, "\" -r near -tr 5000 5000", sep=""))

## derived boolean maps:
shell(cmd=paste(ILWIS, " landmask_Lon0f{dom=Bool.dom} = iff(landmask_Lon0=1, 1, ?)", sep=""), wait=F)
shell(cmd=paste(ILWIS, " landmask_Lon0w{dom=Bool.dom} = iff(landmask_Lon0=0, 1, ?)", sep=""), wait=F)
shell(cmd=paste(ILWIS, " landmask_Lon120f{dom=Bool.dom} = iff(landmask_Lon120=1, 1, ?)", sep=""), wait=F)
shell(cmd=paste(ILWIS, " landmask_Lon120w{dom=Bool.dom} = iff(landmask_Lon120=0, 1, ?)", sep=""), wait=F)
shell(cmd=paste(ILWIS, " landmask_LonM120f{dom=Bool.dom} = iff(landmask_LonM120=1, 1, ?)", sep=""), wait=F)
shell(cmd=paste(ILWIS, " landmask_LonM120w{dom=Bool.dom} = iff(landmask_LonM120=0, 1, ?)", sep=""), wait=F)

## derive distance maps:
shell(cmd=paste(ILWIS, " dcoast_Af.mpr{dom=value.dom;vr=0:100000000:1} = MapDistance(landmask_Lon0f)", sep=""), wait=F)
shell(cmd=paste(ILWIS, " dcoast_Aw.mpr{dom=value.dom;vr=0:100000000:1} = MapDistance(landmask_Lon0w)", sep=""), wait=F)
shell(cmd=paste(ILWIS, " dcoast_A.mpr{dom=value.dom;vr=-100000:100000:1} = (dcoast_Af - dcoast_Aw) / 1000", sep=""), wait=F)
shell(cmd=paste(ILWIS, " dcoast_Bf.mpr{dom=value.dom;vr=0:100000000:1} = MapDistance(landmask_Lon120f)", sep=""), wait=F)
shell(cmd=paste(ILWIS, " dcoast_Bw.mpr{dom=value.dom;vr=0:100000000:1} = MapDistance(landmask_Lon120w)", sep=""), wait=F)
shell(cmd=paste(ILWIS, " dcoast_B.mpr{dom=value.dom;vr=-100000:100000:1} = (dcoast_Bf - dcoast_Bw) / 1000", sep=""), wait=F)
shell(cmd=paste(ILWIS, " dcoast_Cf.mpr{dom=value.dom;vr=0:100000000:1} = MapDistance(landmask_LonM120f)", sep=""), wait=F)
shell(cmd=paste(ILWIS, " dcoast_Cw.mpr{dom=value.dom;vr=0:100000000:1} = MapDistance(landmask_LonM120w)", sep=""), wait=F)
shell(cmd=paste(ILWIS, " dcoast_C.mpr{dom=value.dom;vr=-100000:100000:1} = (dcoast_Cf - dcoast_Cw) / 1000", sep=""), wait=F)

## Resample back to geographic coords:
system(paste(gdalwarp, ' dcoast_A.mpr dcoast_A.sdat -of \"SAGA\" -s_srs \"', g.csy, '\" -t_srs \"+proj=longlat +datum=WGS84\" -r bilinear -tr 0.05 0.05 -te -180 -90 180 90', sep=""))
system(paste(gdalwarp, ' dcoast_B.mpr dcoast_B.sdat -of \"SAGA\" -s_srs \"', g2.csy, '\" -t_srs \"+proj=longlat +datum=WGS84\" -r bilinear -tr 0.05 0.05 -te 20 -90 180 90', sep=""))
system(paste(gdalwarp, ' dcoast_C.mpr dcoast_C.sdat -of \"SAGA\" -s_srs \"', g3.csy, '\" -t_srs \"+proj=longlat +datum=WGS84\" -r bilinear -tr 0.05 0.05 -te -180 -90 21 90', sep=""))


## derive minimum distance:
DIST.sgrd.list <- paste("dcoast_", LETTERS[1:3], ".sgrd", sep="", collapse=";")
rsaga.geoprocessor(lib="grid_tools", module=3, param=list(GRIDS=DIST.sgrd.list, MERGED="DICGSH1a.sgrd", TYPE=7, INTERPOL=1, OVERLAP=0, MERGE_INFO_MESH_SIZE=1))

## resample to geographic coordinates:
system(paste(gdal_translate, " DICGSH1a.sdat DICGSH1a.tif -a_srs \"+proj=longlat +datum=WGS84\"", sep=""))
## 20 km resolution:
system(paste(gdalwarp, ' DICGSH1a.sdat DICGSH0a.tif -r bilinear -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))

## compress and copy:
tif2.lst <- list.files(pattern=glob2rx("DICGSH*a.tif$"))  
for(i in 1:length(tif2.lst)){
    system(paste("7za a", "-tgzip", set.file.extension(tif2.lst[i], ".tif.gz"), tif2.lst[i]))
    system(paste("xcopy", set.file.extension(tif2.lst[i], ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(tif2.lst[i], ".tif.gz"))
}

# end of script;