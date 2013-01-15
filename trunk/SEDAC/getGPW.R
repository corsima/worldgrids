# title         : getGPW.R
# purpose       : derivation of PCs using the gridded population density maps from Socioeconomic data and application centre (SEDAC), Columbia University;
# reference     : [http://worldgrids.org/doku.php?id=wiki:pd1gpw1]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Jan 2013.
# inputs        : maps publicaly available at [http://sedac.ciesin.columbia.edu/data/collection/gpw-v3];
# outputs       : geotiff images projected in the "+proj=longlat +ellps=WGS84" system;
# remarks 1     : First download and install FWtools [http://fwtools.maptools.org];

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


## download all maps in BIL format:
for(year in c("90", "95", "00", "05", "10", "15")){
system(paste('"c:/Program Files/Mozilla Firefox/firefox.exe"', 'http://sedac.ciesin.columbia.edu/gpw/global.jsp?file=gpwv3&data=pdens&type=bil&resolut=25&year=', year, sep=""), wait = FALSE)
}

## unzip and reproject files:
zip.lst <- dir(path=getwd(), pattern=glob2rx("*_bil_25.zip"), full.names=FALSE, recursive=FALSE)
for(i in 1:length(zip.lst)){  unzip(zip.lst[i], exdir=getwd()) }
## reproject to 5 km resolution and add poles:
bil.lst <- dir(path=getwd(), pattern=glob2rx("*ag.bil"), full.names=FALSE, recursive=FALSE)
for(i in 1:length(bil.lst)){
  tif.name <- set.file.extension(bil.lst[i], ".tif")
  if(is.na(file.info(tif.name)$size)){
    system(paste(gdalwarp, bil.lst[i], tif.name, '-ot \"Float32\" -r bilinear -te -180 -90 180 90 -tr', 1/20, 1/20))
  }
}
## clean up:
unlink(dir(path=getwd(), pattern=glob2rx("*.bil"), full.names=FALSE, recursive=FALSE))
unlink(dir(path=getwd(), pattern=glob2rx("*.hdr"), full.names=FALSE, recursive=FALSE))
unlink(dir(path=getwd(), pattern=glob2rx("*.blw"), full.names=FALSE, recursive=FALSE))
unlink(dir(path=getwd(), pattern=glob2rx("*.stx"), full.names=FALSE, recursive=FALSE))

# ------------------------------------------------------------
# Derivation of PC1-2:
# ------------------------------------------------------------

tif.lst <- dir(path=getwd(), pattern=glob2rx("*ag.tif"), full.names=FALSE, recursive=FALSE)
glds <- readGDAL(tif.lst[1])
for(i in 2:length(tif.lst)){ glds@data[,tif.lst[i]] <- readGDAL(tif.lst[i])$band1 }
names(glds) <- tif.lst
proj4string(glds) = "+proj=longlat +datum=WGS84"

## filter the missing values:
x <- glds@data
#x[is.na(x)] <- 0
for(i in 1:ncol(x)){ x[,i] <- log1p(x[,i]) } 
summary(x)

formulaString <- ~ glds00ag.tif + glds05ag.tif + glds10ag.tif + glds90ag.tif + glds95ag.tif + glds15ag.tif
pcs <- prcomp(formula=formulaString, x)
str(pcs)
## copy values: 
glds@data[,"PC1"] <- pcs$x[,1]
glds@data[,"PC2"] <- pcs$x[,2]
glds@data[,"M"] <- rowSums(glds@data[,1:6])/6

# ------------------------------------------------------------
# Resampling and export to geotiff:
# ------------------------------------------------------------

writeGDAL(glds["PC1"], "PD1GPW1a.tif", driver="GTiff", mvFlag=-99999)
writeGDAL(glds["PC2"], "PD2GPW1a.tif", driver="GTiff", mvFlag=-99999)
writeGDAL(glds["M"], "PDMGPW1a.tif", driver="GTiff", mvFlag=-99999)
GDALinfo("PDMGPW1a.tif")

## 20 km resolution:
out.lst <- c("PD1GPW", "PD2GPW", "PDMGPW")
for(j in out.lst){ system(paste(gdalwarp, ' ', j, '1a.tif', ' ', j, '0a.tif -r bilinear -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep="")) }

## compress and copy:
tif2.lst <- list.files(pattern=glob2rx("PD*a.tif$"))  
for(i in 1:length(tif2.lst)){
    system(paste("7za a", "-tgzip", set.file.extension(tif2.lst[i], ".tif.gz"), tif2.lst[i]))
    system(paste("xcopy", set.file.extension(tif2.lst[i], ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(tif2.lst[i], ".tif.gz"))
}

# end of script;