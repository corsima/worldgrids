# title         : getGADM.R
# purpose       : download and resampling of Global Administrative DM data;
# reference     : [http://worldgrids.org]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Nov 2012.
# inputs        : data available for download from [http://www.gadm.org/world]; 
# outputs       : geotiff images projected in the "+proj=longlat +ellps=WGS84" system;
# remarks 1     : First download and install SAGA GIS [http://www.saga-gis.org] and FWtools [http://fwtools.maptools.org];  


# ------------------------------------------------------------
# Initial settings and data download:
# ------------------------------------------------------------

library(RSAGA) 
library(rgdal)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"

# download files from server:
download.file("http://www.gadm.org/data2/gadm_v2_shp.zip", destfile=paste(getwd(), "gadm_v2_shp.zip", sep="/"))
unzip(zipfile="gadm_v2_shp.zip", exdir=getwd())
ogrInfo("gadm2.shp", "gadm2")
## Get the country names (ISO standard)
gadm <- read.dbf("gadm2.dbf")  ## Large DBF!
str(gadm)
cnt = aggregate(gadm$dbf[,"ID_0"], by=list(gadm$dbf$ISO), FUN=mean)
str(cnt)

## make a PAL file:
rgb.tbl <- data.frame(R=NA, G=NA, B=NA)
library(XML)
tmp <- xmlTreeParse("cntgad.sprm", useInternalNodes = TRUE)
cols <- xmlRoot(tmp)[[11]]
col.lst <- sapply(sapply(xmlChildren(cols), function(x) x), xmlValue)
for(j in 1:length(col.lst)){
   rgb.tbl[j,c("R","G","B")] <- as.integer(substr(strsplit(col.lst[j], " ")[[1]], start=2, stop=4))
}
str(rgb.tbl)
## write a PAL file:
write.table(rgb.tbl, "cntgad.PAL", quote = FALSE,  col.names = FALSE)

## create a SAGA txt colour table
## convert to BGR codes:
BGR <- (rgb.tbl$B * 65536) + (rgb.tbl$G * 256) + rgb.tbl$R
## write a lookup table for SAGA GIS:
filename <- file("cntgad.txt", "w", blocking=FALSE)
write("COLOR\tNAME\tDESCRIPTION\tMINIMUM\tMAXIMUM", filename)
for(i in 1:nrow(cnt)){
  write(paste(BGR[i], cnt[i,"ISO"], paste("CL", i, sep=""), (i-1)+0.1, (i+1)+0.1, sep="\t"), filename, append=TRUE)
}
close(filename)

# convert to raster (use the polygon ID):
rsaga.geoprocessor(lib="grid_gridding", module=0, param=list(USER_GRID="countries.sgrd", INPUT="gadm2.shp", FIELD=1, TARGET=0, LINE_TYPE=1, USER_SIZE=1/120, USER_XMIN=-180+1/120*0.5, USER_XMAX=180-1/120*0.5, USER_YMIN=-90+1/120*0.5, USER_YMAX=90-1/120*0.5)) # takes cca 5 mins uses > 6GB!!
## export to geotiff:
# system(paste(gdal_translate, "countries.sdat CNTGAD3a.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\""))
rsaga.geoprocessor(lib="io_gdal", module=2, param=list(GRIDS="countries.sgrd", FILE="countries.tif"))
## TH: GDAL warp with SAGA sdat file receives error "ERROR 3: Unable to seek to beginning of grid row."
system(paste(gdal_translate, "countries.tif CNTGAD3a.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -a_nodata \"0\""))
# system(paste(gdalwarp, ' countries.tif CNTGAD3a.tif -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"255\" -ot Byte -r near -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""))
GDALinfo("CNTGAD3a.tif")

## Resample to other resolutions using the majority value!
# 2.5 km:
rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT="countries.sgrd", USER_GRID="CNTGAD2a.sgrd", KEEP_TYPE=FALSE, SCALE_UP_METHOD=9, USER_SIZE=1/40, USER_XMIN=-180+1/40*0.5, USER_XMAX=180-1/40*0.5, USER_YMIN=-90+1/40*0.5, USER_YMAX=90-1/40*0.5)) 
# system(paste(gdalwarp, ' countries.tif CNTGAD2a.tif -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"0\" -ot Byte -r near -te -180 -90 180 90 -tr ', 1/40,' ', 1/40, sep=""))
system(paste(gdal_translate, "CNTGAD2a.sdat CNTGAD2a.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -a_nodata \"0\""))
## 5 km:
rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT="CNTGAD2a.sgrd", USER_GRID="CNTGAD1a.sgrd", KEEP_TYPE=FALSE, SCALE_UP_METHOD=9, USER_SIZE=1/20, USER_XMIN=-180+1/20*0.5, USER_XMAX=180-1/20*0.5, USER_YMIN=-90+1/20*0.5, USER_YMAX=90-1/20*0.5))
# system(paste(gdalwarp, ' CNTGAD2a.tif CNTGAD1a.tif -t_srs \"+proj=longlat +datum=WGS84\" -ot Byte -r near -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
system(paste(gdal_translate, "CNTGAD1a.sdat CNTGAD1a.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -a_nodata \"0\""))
## 20 km:
rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT="CNTGAD1a.sgrd", USER_GRID="CNTGAD0a.sgrd", KEEP_TYPE=FALSE, SCALE_UP_METHOD=9, USER_SIZE=1/5, USER_XMIN=-180+1/5*0.5, USER_XMAX=180-1/5*0.5, USER_YMIN=-90+1/5*0.5, USER_YMAX=90-1/5*0.5))
# system(paste(gdalwarp, ' CNTGAD1a.tif CNTGAD0a.tif -t_srs \"+proj=longlat +datum=WGS84\" -ot Byte -r near -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))
system(paste(gdal_translate, "CNTGAD0a.sdat CNTGAD0a.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -a_nodata \"0\""))

## compress:  
for(i in 0:3){
    outname = paste("CNTGAD", i, 'a', sep="")
    system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), set.file.extension(outname, ".tif")))
    system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(outname, ".tif.gz"))
}

# clean up temp files:
unlink(list.files(pattern=glob2rx("countries.*")))

## prepare lines using the GADM:
# import shape(D:\Worldmaps\GADM\gadm1_lev0_lines.shp, D:\Worldmaps\GADM\gadm1_lev0_lines.ioc)
## Simplify lines (tunelling):
# gadm1_lev0_f.mps = SegmentMapTunneling(gadm1_lev0_lines,0.001000,yes)
# export Shapefile(gadm1_lev0_f.mps,worldborders_gadm)
# reproject map to Google Maps coordinate system:
# gadm_lev0_gc <- readShapeLines("worldborders_gadm_gc.shp")
# proj4string(gadm_lev0_gc) <- CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +wktext  +no_defs")
# save(gadm_lev0_gc, file="gadm_lev0_gc.RData")


# end of script;