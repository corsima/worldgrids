# title         : processMOD13A3.R
# purpose       : processing of the MODIS 1 km resolution monthly EVI images;
# reference     : [http://worldgrids.org]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Aug 2012.
# inputs        : 1 km geotifs in the format "EVI_Day_1km_*_*_*.tif" (global composites produced using FWTools);
# outputs       : geotiff images projected in the "+proj=longlat +ellps=WGS84" system;
# remarks 1     : First download and install FWtools [http://fwtools.maptools.org];

library(rgdal)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
gdalbuildvrt = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalbuildvrt.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"

## tiling system (24 tiles):
p.l <- expand.grid(KEEP.OUT.ATTRS=FALSE, lonmin=seq(-180,120,by=60), latmin=seq(-90,45,by=45))
p.u <- expand.grid(KEEP.OUT.ATTRS=FALSE, lonmax=seq(-120,180,by=60), latmax=seq(-45,90,by=45))
tiles <- cbind(p.l, p.u)
str(tiles)

## list files:
dtif.lst <- list.files(pattern=glob2rx("EVI_1km_*_*_*.tif$"))
GDALinfo(dtif.lst[1])

## EVI images:
for(j in 1:nrow(tiles)){
  for(i in 1:length(dtif.lst)){
    outname <- strsplit(dtif.lst[i], ".tif")[[1]][1]
    tilename <- paste(outname, '_', j, '.sdat', sep="")
    # resample
    if(is.na(file.info(tilename)$size)){
      system(paste(gdalwarp, ' ', dtif.lst[i], ' ', tilename, ' -of \"SAGA\" -dstnodata \"-3000\" -r near -te ', tiles[j,"lonmin"] , ' ', tiles[j,"latmin"], ' ', tiles[j,"lonmax"] ,' ', tiles[j,"latmax"] ,' -tr ', 1/120, ' ', 1/120, sep=""))
    }
  }
  # derive derive mean, sd, min and max values:
  t.lst <- list.files(pattern=glob2rx(paste('EVI_1km_*_*_*_', j, '.sgrd$', sep="")))
  if(is.na(file.info(paste("EVMMOD3a_", j, ".sgrd", sep=""))$size)){
  rsaga.geoprocessor(lib="geostatistics_grid", module=4, param=list(GRIDS=paste(t.lst, collapse=";"), MEAN=paste("EVMMOD3a_", j, ".sgrd", sep=""), MIN=paste("EVLMOD3a_", j, ".sgrd", sep=""), MAX=paste("EVHMOD3a_", j, ".sgrd", sep=""), STDDEV=paste("EVSMOD3a_", j, ".sgrd", sep="")))
  # rsaga.geoprocessor(lib="geostatistics_grid", module=8, param=list(GRIDS=paste(t.lst, collapse=";"), PCA=paste("EV_MOD3a_", j, "_", sep=""), METHOD=1, NFIRST=3))
  }
}

# create a mosaic and compress files:
for(j in c("EVMMOD", "EVSMOD", "EVLMOD", "EVHMOD")){
  unlink("mod13a3.vrt")
  system(paste(gdalbuildvrt, "mod13a3.vrt", paste(j, "3a_", 1:nrow(tiles), ".sdat", sep="", collapse=" ")))
  # convert to geotif (1 km):
  if(is.na(file.info(paste(j, '3a.tif', sep=""))$size)){
  system(paste(gdalwarp, ' mod13a3.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-3000\" -ot Int16 ', j, '3a.tif', ' -r near -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""))
  }
  # 2.5 km:
  if(is.na(file.info(paste(j, '2a.tif', sep=""))$size)){
  system(paste(gdalwarp, ' mod13a3.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-3000\" -ot Int16 ', j, '2a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/40,' ', 1/40, sep=""))
  }
  # 5 km:
  if(is.na(file.info(paste(j, '1a.tif', sep=""))$size)){
  system(paste(gdalwarp, ' mod13a3.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-3000\" -ot Int16 ', j, '1a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
  }
  # 20 km:
  system(paste(gdalwarp, ' ', j, '1a.tif', ' -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-3000\" -ot Int16 ', j, '0a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))
}

## compress:  
for(j in c("EVMMOD", "EVSMOD", "EVLMOD", "EVHMOD")){
for(i in 0:3){
    outname = paste(j, i, 'a', sep="")
    system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), set.file.extension(outname, ".tif")))
    system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(outname, ".tif.gz"))
  }
}

# clean up temp files:
unlink(list.files(pattern=glob2rx("EVI_1km_*_*_*_*.*")))
unlink(list.files(pattern=glob2rx("EV*MOD3a_*.*")))

# end of script;