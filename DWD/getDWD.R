# title         : getDWD.R
# purpose       : downscale precipitation images (GPCC) to 5 km resolution using cubic splines;
# reference     : [http://worldgrids.org/doku.php?id=wiki:p01dwd1]
# producer      : Prepared by T. Hengl
# version       : 1
# inputs        : Global monthly precipitation images at 0.25 degrees [ftp://ftp-anon.dwd.de/pub/data/gpcc/html/gpcc_normals_download.html];
# outputs       : geotiff images projected in the "+proj=longlat +datum=WGS84" system;
# remarks 1     : First download and install FWtools [http://fwtools.maptools.org];
# remarks 2     : The calculation might miss some small islands;

library(RSAGA) 
library(rgdal)
library(raster)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"

## download landmask at 5 km resolution:
if(is.na(file.info("lmtgsh1a.tif.gz")$size)){
  download.file("http://worldgrids.org/lib/exe/fetch.php?media=lmtgsh1a.tif.gz", "lmtgsh1a.tif.gz")
  system(paste("7za e lmtgsh1a.tif.gz"))
}

## download precipitation images:
if(is.na(file.info("gpcc_precipitation_normals_rr_version_2011_0_25_degree.zip")$size)){
  download.file("ftp://ftp-anon.dwd.de/pub/data/gpcc/gpcc_normals_v2011/gpcc_precipitation_normals_rr_version_2011_0_25_degree.zip", "gpcc_precipitation_normals_rr_version_2011_0_25_degree.zip")
  system(paste("7za e gpcc_precipitation_normals_rr_version_2011_0_25_degree.zip"))
}

## read rasters to R:
pr <- read.table("gpcc_precipitation_normals_rr_version_2011_0_25_degree", header=FALSE, skip=11)
str(pr)
unlink("gpcc_precipitation_normals_rr_version_2011_0_25_degree")
## mask out missing values:
for(i in 1:ncol(pr)){ pr[,i] <- ifelse(pr[,i]<0, NA, pr[,i]) }
names(pr) <- c(paste("P0", 1:9, "DWD", sep=""), paste("P", 10:12, "DWD", sep=""), "PTADWD")
grid025 <- expand.grid(lon=seq(-180+.25/2, 180-.25/2, by=.25), lat=seq(-90+.25/2, 90-.25/2, by=.25))
gridded(grid025) <- ~lon+lat
proj4string(grid025) = "+proj=longlat +datum=WGS84"
grid025 <- SpatialPixelsDataFrame(as.matrix(data.frame(lon=grid025@coords[,1], lat=grid025@coords[nrow(grid025@coords):1,2])), data=pr, proj4string=grid025@proj4string, grid=grid025@grid)
#spplot(grid025[1]) 

# convert to geotifs (5 km):
grid005 <- readGDAL("LMTGSH1a.tif")
grid005$band1 <- ifelse(grid005$band1==129, NA, 1)
grid005 <- as(grid005, "SpatialPixelsDataFrame") ## takes ca 5 mins!

for(i in 1:(ncol(grid025)-1)){
  if(is.na(file.info(paste(names(grid025)[i], "1a.tif", sep=""))$size)){  
    writeGDAL(grid025[i], paste(names(grid025)[i], "X.tif", sep=""), "GTiff", mvFlag=-9999, type="Int16")
    unlink("tmp.tif")
    system(paste(gdalwarp, ' ', paste(names(grid025)[i], "X.tif", sep=""), ' tmp.tif -dstnodata \"-9999\" -ot \"Int16\" -r cubicspline -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
    ## land mask:
    grid005@data[,names(grid025)[i]] <- readGDAL("tmp.tif")$band1[grid005@grid.index]
    writeGDAL(grid005[names(grid025)[i]], paste(names(grid025)[i], "1a.tif", sep=""), "GTiff", mvFlag=-9999, type="Int16")
  }
}

## resample to coarser resolutions:
# 20 km:
tif.lst <- list.files(pattern=glob2rx("P*DWD1a.tif$"))
for(i in 1:length(tif.lst)){
  outname = gsub("1a", "0a", tif.lst[i])
  if(is.na(file.info(outname)$size)){ 
    system(paste(gdalwarp, ' ', tif.lst[i], ' -dstnodata \"-9999\" ', outname, ' -r bilinear -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))
  }
}

## compress and copy:
tif2.lst <- list.files(pattern=glob2rx("P*DWD*a.tif$"))  
for(i in 1:length(tif2.lst)){
    system(paste("7za a", "-tgzip", set.file.extension(tif2.lst[i], ".tif.gz"), tif2.lst[i]))
    system(paste("xcopy", set.file.extension(tif2.lst[i], ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(tif2.lst[i], ".tif.gz"))
}

# clean up temp files:
unlink(list.files(pattern=glob2rx("P*DWDX.tif$")))

# end of script;
