# title         : makeSoilMask.R
# purpose       : Global soil mask at multiple resolutions
# reference     : Worldgrids.org;
# producer      : Prepared by T. Hengl
# address       : In Wageningen, NL.
# inputs        : Global Raster Water Mask, MODIS land cover map and Leaf Area Index map (see WorldGrids.org -> layers);
# outputs       : geotiff images projected in the WGS84 coordinate system;
# remarks 1     : Methodology explained at http://gsif.isric.org;

# -------------------------------------
# Initial settings and data download:
# -------------------------------------

library(rgdal)
library(raster)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
gdalbuildvrt = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalbuildvrt.exe"))))
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"
## Get 7z:
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")

## download maps:
g.lst <- c("G11IGB3", "LAMMOD3", "WMKMOD3")
for(j in 1:length(g.lst)){
  outname = paste(g.lst[j], "a.tif.gz", sep="")
  outname.tif = paste(g.lst[j], "a.tif", sep="")
  if(is.na(file.info(outname.tif)$size)){
    if(is.na(file.info(outname)$size)){
      download.file(paste("http://worldgrids.org/lib/exe/fetch.php?media=", outname, sep=""), outname)
    }
    system(paste("7za e", outname))
    unlink(outname)
  } 
}

## These rasters are HUGE, hence we run raster calcs per 24 tiles...

# -------------------------------------
# Tiling system:
# -------------------------------------

tif.lst <- list.files(pattern=glob2rx("*3a.tif$"))
info <- GDALinfo(tif.lst[1])
t.l <- expand.grid(KEEP.OUT.ATTRS=FALSE, lonmin=seq(-180,120,by=60), latmin=seq(-90,45,by=45))
t.u <- expand.grid(KEEP.OUT.ATTRS=FALSE, lonmax=seq(-120,180,by=60), latmax=seq(-45,90,by=45))
tiles <- cbind(t.l, t.u)
tiles$offset.y <- round(info[["rows"]]*(90-tiles$latmax)/180)
tiles$offset.x <- info[["columns"]] + round(info[["columns"]]*(tiles$lonmin-180)/360)
tiles$region.dim.y <- round(info[["rows"]]*(tiles$latmax-tiles$latmin)/180)
tiles$region.dim.x <- round(info[["columns"]]*(tiles$lonmax-tiles$lonmin)/360)

# -------------------------------------
# Process per tile:
# -------------------------------------

for(i in 1:nrow(tiles)){
  if(is.na(file.info(paste("SMKISR3a_", i, ".tif", sep=""))$size)){
  ## import data:
  x <- readGDAL(tif.lst[1], offset=c(tiles$offset.y[i], tiles$offset.x[i]), region.dim=c(tiles$region.dim.y[i], tiles$region.dim.x[i]))
  for(j in 2:length(tif.lst)){
    x@data[,tif.lst[j]] <- readGDAL(tif.lst[j], offset=c(tiles$offset.y[i], tiles$offset.x[i]), region.dim=c(tiles$region.dim.y[i], tiles$region.dim.x[i]))$band1
  }
  names(x) <- tif.lst
  ## derive soil mask values:
  x$smask <- ifelse(x$G11IGB3a.tif==16&is.na(x$LAMMOD3a.tif), 3, ifelse(x$G11IGB3a.tif==13, 2, ifelse(x$LAMMOD3a.tif>0&!x$LAMMOD3a.tif>100|!(x$G11IGB3a.tif==16|x$G11IGB3a.tif==13|x$G11IGB3a.tif==0|x$G11IGB3a.tif==15)&x$WMKMOD3a.tif<60, 1, NA))) ## takes 1-2 mins per tile!
  ## write a geotif:
  writeGDAL(x["smask"], paste("SMKISR3a_", i, ".tif", sep=""), type="Byte", mvFlag=0)
}
}

## Create a mosaic:
lst.mos <- list.files(pattern=glob2rx("SMKISR3a_*.tif"))
lst.mos
unlink("my_liste.txt")
cat(lst.mos, sep="\n", file="my_liste.txt", append=TRUE)
unlink("glc.vrt")
system(paste(gdalbuildvrt, "-input_file_list my_liste.txt glc.vrt"))
## resample:
unlink("SMKISR3a.tif")
system(paste(gdalwarp, ' glc.vrt SMKISR3a.tif -t_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -r near -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""))

# -------------------------------------
# Export and save:
# -------------------------------------

rgb.tbl <- data.frame(Category=c("soils with vegetation cover","urban areas","bare soil areas"), Number=1:3, R=c(50,242,194), G=c(252,7,194), B=c(23,7,194))
## create a SAGA txt colour table
## convert to BGR codes:
BGR <- (rgb.tbl$B * 65536) + (rgb.tbl$G * 256) + rgb.tbl$R
## write a lookup table for SAGA GIS:
filename <- file("SMKISR.txt", "w", blocking=FALSE)
write("COLOR\tNAME\tDESCRIPTION\tMINIMUM\tMAXIMUM", filename)
for(i in 1:nrow(rgb.tbl)){
  write(paste(BGR[i], rgb.tbl[i,"Category"], paste("CL", i, sep=""), (rgb.tbl[i,"Number"]-1)+0.1, (rgb.tbl[i,"Number"])+0.1, sep="\t"), filename, append=TRUE)
}
close(filename)

## write a PAL file:
write.table(rgb.tbl[,c("Category","R","G","B")], "SMKISR.PAL", quote = FALSE,  col.names = FALSE)

## resample:
system(paste(gdalwarp, ' SMKISR3a.tif SMKISR2a.tif -r near -ot \"Byte\" -te -180 -90 180 90 -tr ', 1/40,' ', 1/40, sep=""))
unlink("WMKMOD1a.tif")
system(paste(gdalwarp, ' SMKISR3a.tif SMKISR1a.tif -r near -ot \"Byte\" -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
unlink("WMKMOD0a.tif")
system(paste(gdalwarp, ' SMKISR3a.tif SMKISR0a.tif -r near -ot \"Byte\" -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))
unlink("wmask.tif")
system(paste(gdalwarp, ' SMKISR3a.tif smask.tif -r near -ot \"Byte\" -te -180 -90 180 90 -tr ', 1,' ', 1, sep=""))


## Zip the output files:
for(i in 0:3){
  outname = paste("SMKISR", i, 'a', sep="")
  system(paste("7za a", "-tgzip -mx=9", set.file.extension(outname, ".tif.gz"), set.file.extension(outname, ".tif"))) 
  system(paste("xcopy", normalizePath(set.file.extension(outname, ".tif.gz")), shortPathName(normalizePath(outdir))))
  unlink(set.file.extension(outname, ".tif.gz"))
}

## end of script;