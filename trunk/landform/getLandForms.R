# title         : getLandForms.R
# purpose       : download and resampling of Global landform classification maps;
# reference     : [http://worldgrids.org]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Nov 2012.
# inputs        : SCALA physiographic map of the world / Iwahashi and Richard J. Pike (2007) map available for download from [http://gisstar.gsi.go.jp/terrain/front_page.htm]; ; 
# outputs       : geotiff images projected in the "+proj=longlat +ellps=WGS84" system;
# remarks 1     : First download and install SAGA GIS [http://www.saga-gis.org] and FWtools [http://fwtools.maptools.org]; 

library(RSAGA) 
library(rgdal)
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
system(paste(gdalwarp, "L3_physiographic.tif", "-t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"255\" L3POBI3a.tif -r near -te -180 -90 180 90 -tr", 1/120, 1/120))
# check validity:
GDALinfo("L3POBI3a.tif")
# resample to 2.5, 5 and 20 km resolutions:
system(paste(gdalwarp, "L3POBI3a.tif L3POBI2a.tif -dstnodata \"255\" -r near -te -180 -90 180 90 -tr", 1/40, 1/40))
system(paste(gdalwarp, "L3POBI3a.tif L3POBI1a.tif -dstnodata \"255\" -r bilinear -te -180 -90 180 90 -tr", 1/20, 1/20))
system(paste(gdalwarp, "L3POBI3a.tif L3POBI0a.tif -dstnodata \"255\" -r bilinear -te -180 -90 180 90 -tr", 1/5, 1/5))

# Compress:
for(outname in c("L3POBI0a.tif", "L3POBI1a.tif", "L3POBI2a.tif", "L3POBI3a.tif")){
  if(is.na(file.info(paste(shortPathName(normalizePath(outdir)), paste(outname, "gz", sep="."), sep="\\"))$size)){
  system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), outname))
  system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
  unlink(set.file.extension(outname, ".tif.gz"))
}  # Compression takes > 15 mins
}

## make a PAL file:
cnt = data.frame(Group.1=c("flat plains", "high mountains", "high hills", "low mountains", "tablelands", "rough low hills", "smooth low hills", "irregular plains"), x=c(1:8))
rgb.tbl <- data.frame(R=c(115,206,201,226,232,245,244,188), G=c(180,206,144,176,181,219,237,206), B=c(116,206,124,185,128,144,159,130))
## write a PAL file:
write.table(rgb.tbl, "l3pobi.PAL", quote = FALSE,  col.names = FALSE)

## create a SAGA txt colour table
## convert to BGR codes:
BGR <- (rgb.tbl$B * 65536) + (rgb.tbl$G * 256) + rgb.tbl$R
## write a lookup table for SAGA GIS:
filename <- file("l3pobi.txt", "w", blocking=FALSE)
write("COLOR\tNAME\tDESCRIPTION\tMINIMUM\tMAXIMUM", filename)
for(i in 1:nrow(cnt)){
  write(paste(BGR[i], cnt[i,"Group.1"], paste("CL", i, sep=""), (cnt[i,"x"]-1)+0.1, cnt[i,"x"]+0.1, sep="\t"), filename, append=TRUE)
}
close(filename)



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