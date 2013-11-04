# title         : getLandForms.R
# purpose       : download and resampling of Global landform classification maps;
# reference     : [http://worldgrids.org]
# producer      : Prepared by T. Hengl
# address       : In Wageningen, NL, Nov 2012.
# inputs        : SCALA physiographic map of the world by Clemens Eisank / Iwahashi and Richard J. Pike (2007) map available for download from [http://gisstar.gsi.go.jp/terrain/front_page.htm]; ; 
# outputs       : geotiff images projected in the "+proj=longlat +ellps=WGS84" system;
# remarks 1     : First download and install SAGA GIS [http://www.saga-gis.org] and FWtools [http://fwtools.maptools.org]; 

library(RSAGA) 
library(rgdal)
library(foreign)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"
## Robinson projection system:
# [http://spatialreference.org/ref/esri/54030/]
prj.1 <- "+proj=robin +lon_0=120 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
prj.2 <- "+proj=robin +lon_0=40 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
prj.3 <- "+proj=robin +lon_0=-80 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
prj.l <- list(prj.1, prj.2, prj.3)

## Automated object-based classification of topography from SRTM data [http://dx.doi.org/10.1016/j.geomorph.2011.12.001]
## Data can be downloaded from [http://zgis202.plus.sbg.ac.at/LandformClassification/default.aspx]
#system(paste(gdalwarp, "L3_physiographic.tif", "-t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"255\" L3POBI3a.tif -r near -te -180 -90 180 90 -tr", 1/120, 1/120))

## updated map, now as a shapefile (3 Nov 2013):
tbl <- read.dbf("L3_physioReg.dbf")  ## Large DBF!
str(tbl)
tbl$dbf$L3POBI3a <- as.integer(tbl$dbf$Class_name)
write.dbf(tbl, "L3_physioReg.dbf")
levs <- levels(tbl$dbf$Class_name)

## convert to raster (use the output class):
myenv <- rsaga.env(path="C:/saga_vc")
rsaga.geoprocessor(lib="grid_gridding", module=0, param=list(USER_GRID="L3_physioReg.sgrd", INPUT="L3_physioReg.shp", FIELD=18, TARGET=0, LINE_TYPE=1, USER_SIZE=1/120, USER_XMIN=-180+1/120*0.5, USER_XMAX=180-1/120*0.5, USER_YMIN=-90+1/120*0.5, USER_YMAX=90-1/120*0.5), env=myenv) # takes cca 5 mins uses > 4GB!!
## export to geotiff:
rsaga.geoprocessor(lib="io_gdal", module=2, param=list(GRIDS="L3_physioReg.sgrd", FILE="L3_physioReg.tif"), env=myenv)
system(paste(gdal_translate, "L3_physioReg.tif L3POBI3b.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -a_nodata \"0\""))
## check validity:
GDALinfo("L3POBI3b.tif")
unlink("L3_physioReg.tif")

# resample to 2.5, 5 and 20 km resolutions:
system(paste(gdalwarp, "L3POBI3b.tif L3POBI2b.tif -dstnodata \"0\" -r near -te -180 -90 180 90 -tr", 1/40, 1/40))
system(paste(gdalwarp, "L3POBI3b.tif L3POBI1b.tif -dstnodata \"0\" -r near -te -180 -90 180 90 -tr", 1/20, 1/20))
system(paste(gdalwarp, "L3POBI3b.tif L3POBI0b.tif -dstnodata \"0\" -r near -te -180 -90 180 90 -tr", 1/5, 1/5))

# Compress:
for(outname in c("L3POBI0b.tif", "L3POBI1b.tif", "L3POBI2b.tif", "L3POBI3b.tif")){
  if(is.na(file.info(paste(shortPathName(normalizePath(outdir)), paste(outname, "gz", sep="."), sep="\\"))$size)){
  system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), outname))
  system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
  unlink(set.file.extension(outname, ".tif.gz"))
}  # Compression takes > 15 mins
}

## make a PAL file:
rgb.tbl <- data.frame(x=c(1,3,2,5,8,6,7,4), Group.1=c("flat plains", "high mountains", "high hills", "low mountains", "tablelands", "rough low hills", "smooth low hills", "irregular plains"), R=c(115,206,201,226,232,245,244,188), G=c(180,206,144,176,181,219,237,206), B=c(116,206,124,185,128,144,159,130))
rgb.tbl <- rgb.tbl[order(rgb.tbl$x),] 

## create a SAGA txt colour table
## convert to BGR codes:
BGR <- (rgb.tbl$B * 65536) + (rgb.tbl$G * 256) + rgb.tbl$R
## write a lookup table for SAGA GIS:
filename <- file("L3POBI3b.txt", "w", blocking=FALSE)
write("COLOR\tNAME\tDESCRIPTION\tMINIMUM\tMAXIMUM", filename)
for(i in 1:nrow(rgb.tbl)){
  write(paste(BGR[i], rgb.tbl[i,"Group.1"], paste("CL", i, sep=""), (rgb.tbl[i,"x"]-1)+0.1, rgb.tbl[i,"x"]+0.1, sep="\t"), filename, append=TRUE)
}
close(filename)

## write a PAL file:
write.table(rgb.tbl[,c("Group.1","R","G","B")], "L3POBI3b.PAL", quote = FALSE,  col.names = FALSE)


## Iwahashi and Richard J. Pike (2007) map:
## download files from server:
z.lst <- c("Asia_Oceania", "America", "Europe_Africa") 
for(j in 1:length(z.lst)){
  if(is.na(file.info(paste(z.lst[j], "zip", sep="."))$size)){
    download.file(paste("http://gisstar.gsi.go.jp/terrain/original/", z.lst[j], ".zip", sep=""), destfile=paste(getwd(), paste(z.lst[j], "zip", sep="."), sep="/"))
  }
  unzip(zipfile=paste(z.lst[j], "zip", sep="."), exdir=getwd())
}
GDALinfo("America/Robinson80W_N_S_America/america_class")
GDALinfo("Europe_Africa/Robinson40E_Europe_Africa/europe_class")

## three tiles:
t.lst <- c("Asia_Oceania/Robinson120E_Asia_Oceania/asia_class", "Europe_Africa/Robinson40E_Europe_Africa/europe_class", "America/Robinson80W_N_S_America/america_class")

## We only need two blocks:
system(paste(gdal_translate, t.lst[2], "-a_srs", paste('\"', prj.l[[2]], '\"', sep=""), 'tmp.tif'))
unlink(paste('IWLSRE3a_', 2, '.tif', sep=""))
system(paste(gdalwarp, "tmp.tif", "-t_srs \"+proj=longlat +datum=WGS84\"", paste('IWLSRE3a_', 2, '.tif', sep=""), " -r near -te -30 -90 180 90 -tr", 1/120, 1/120))


system(paste(gdalwarp, t.lst[2], "-s_srs", paste('\"', prj.l[[2]], '\"', sep=""), "-t_srs \"+proj=longlat +datum=WGS84\"", paste('IWLSRE3a_', 2, '.tif', sep=""), " -r near -te -30 -90 180 90 -tr", 1/120, 1/120))
unlink(paste('IWLSRE3a_', 3, '.tif', sep=""))
system(paste(gdalwarp, t.lst[3], "-s_srs", paste('\"', prj.l[[3]], '\"', sep=""), "-t_srs \"+proj=longlat +datum=WGS84\"", paste('IWLSRE3a_', 3, '.tif', sep=""), " -r near -te -180 -90 -30 90 -tr", 1/120, 1/120))


# check validity:
GDALinfo("IWLSRE3a_1.tif")
# resample to 2.5, 5 and 20 km resolutions:
system(paste(gdalwarp, "IWLSRE3a.tif IWLSRE2a.tif -r near -te -180 -90 180 90 -tr", 1/40, 1/40))
system(paste(gdalwarp, "IWLSRE3a.tif IWLSRE1a.tif -r near -te -180 -90 180 90 -tr", 1/20, 1/20))
system(paste(gdalwarp, "IWLSRE3a.tif IWLSRE0a.tif -r near -te -180 -90 180 90 -tr", 1/5, 1/5))

# Compress:
for(outname in c("IWLSRE0a.tif", "IWLSRE1a.tif", "IWLSRE2a.tif", "IWLSRE3a.tif")){
  if(is.na(file.info(paste(shortPathName(normalizePath(outdir)), paste(outname, "gz", sep="."), sep="\\"))$size)){
  system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), outname))
  system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
  unlink(set.file.extension(outname, ".tif.gz"))
}  # Compression takes > 15 mins
}


# end of script;