# title         : processMOD11A2.R
# purpose       : processing of the MODIS 1 km resolution composites (intermadiate products);
# reference     : [http://worldgrids.org]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Aug 2012.
# inputs        : 1 km geotifs in the format "LST_Day_1km_*_*_*.tif" (global composites produced using FWTools);
# outputs       : geotiff images projected in the "+proj=longlat +ellps=WGS84" system;
# remarks 1     : First download and install FWtools [http://fwtools.maptools.org];

library(rgdal)
library(RSAGA)
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
dtif.lst <- list.files(pattern=glob2rx("LST_Day_1km_*_*_*.tif$"))
ntif.lst <- list.files(pattern=glob2rx("LST_Night_1km_*_*_*.tif$"))
x <- as.numeric(sapply(dtif.lst, function(x){strsplit(x, "_")[[1]][5]}))
monthtif.lst <- ifelse(x %in% c(12,1), "X1", ifelse(x %in% c(2,3), "X2", ifelse(x %in% c(4,5), "X3", ifelse(x %in% c(6,7), "X4", ifelse(x %in% c(8,9), "X5", "X6")))))

download.file("http://worldgrids.org/lib/exe/fetch.php?media=smkisr3a.tif.gz", destfile=paste(getwd(), "/", "smkisr3a.tif.gz", sep=""))
system("7za e smkisr3a.tif.gz")

## land mask:
for(j in 1:nrow(tiles)){
  landtile <- paste("SMKISR3a_T", j, ".sdat", sep="")
  if(is.na(file.info(landtile)$size)){
    system(paste(gdalwarp, ' SMKISR3a.tif ', landtile, ' -of \"SAGA\" -srcnodata 0 -dstnodata \"-32767\" -ot Int16 -r near -te ', tiles[j,"lonmin"] , ' ', tiles[j,"latmin"], ' ', tiles[j,"lonmax"] ,' ', tiles[j,"latmax"] ,' -tr ', 1/120, ' ', 1/120, sep=""))
  }
}

## DAY TIME - derive principle components for each tile:
for(j in 1:nrow(tiles)){
  for(i in 1:length(dtif.lst)){
    outname <- strsplit(dtif.lst[i], ".tif")[[1]][1]
    tilename <- paste(outname, '_', j, '.sdat', sep="")
    # resample
    if(is.na(file.info(tilename)$size)){
      system(paste(gdalwarp, ' ', dtif.lst[i], ' ', tilename, ' -of \"SAGA\" -dstnodata \"-32767\" -r near -te ', tiles[j,"lonmin"] , ' ', tiles[j,"latmin"], ' ', tiles[j,"lonmax"] ,' ', tiles[j,"latmax"] ,' -tr ', 1/120, ' ', 1/120, sep=""))
    }
  }
  # derive mean, sd, min and max values:
  t.lst <- list.files(pattern=glob2rx(paste('LST_Day_1km_*_*_*_', j, '.sgrd$', sep="")))
  if(is.na(file.info(paste("TDMMOD3a_", j, ".sgrd", sep=""))$size)){
  rsaga.geoprocessor(lib="geostatistics_grid", module=4, param=list(GRIDS=paste(t.lst, collapse=";"), MEAN=paste("TDMMOD3a_", j, ".sgrd", sep=""), MIN=paste("TDLMOD3a_", j, ".sgrd", sep=""), MAX=paste("TDHMOD3a_", j, ".sgrd", sep=""), STDDEV=paste("TDSMOD3a_", j, ".sgrd", sep="")))
  }
}

## Mean temperature for 6 seasons:
for(j in 1:nrow(tiles)){
  for(i in 1:length(dtif.lst)){
    outname <- strsplit(dtif.lst[i], ".tif")[[1]][1]
    tilename <- paste(outname, '_', j, '.sdat', sep="")
    # resample
    if(is.na(file.info(tilename)$size)){
      system(paste(gdalwarp, ' ', dtif.lst[i], ' ', tilename, ' -of \"SAGA\" -dstnodata \"-32767\" -r near -te ', tiles[j,"lonmin"] , ' ', tiles[j,"latmin"], ' ', tiles[j,"lonmax"] ,' ', tiles[j,"latmax"] ,' -tr ', 1/120, ' ', 1/120, sep=""))
    }
  }
  ## derive mean value per month:
  t.lst <- list.files(pattern=glob2rx(paste('LST_Day_1km_*_*_*_', j, '.sgrd$', sep="")))
  for(k in levels(as.factor(monthtif.lst))){
    sel <- monthtif.lst %in% k
    if(is.na(file.info(paste("T", k, "MOD3a_", j, ".sgrd", sep=""))$size)){
      rsaga.geoprocessor(lib="geostatistics_grid", module=4, param=list(GRIDS=paste(t.lst[sel], collapse=";"), MEAN=paste("T", k, "MOD3a_", j, "_f.sgrd", sep=""))) ## , MIN="tmp.sgrd", MAX="tmp.sgrd", STDDEV="tmp.sgrd"
      ## filter out all missing pixels...
      rsaga.geoprocessor(lib="grid_tools", module=7, param=list(INPUT=paste("T", k, "MOD3a_", j, "_f.sgrd", sep=""), MASK=paste("SMKISR3a_T", j, ".sgrd", sep=""), RESULT=paste("T", k, "MOD3a_", j, ".sgrd", sep=""), THRESHOLD=.1))
    }
  }
}

## create a mosaic and compress files:
for(j in c("TDMMOD", "TDSMOD", "TDLMOD", "TDHMOD", paste("TX", 1:6, "MOD", sep=""))){
  unlink("mod11a2.vrt")
  system(paste(gdalbuildvrt, "mod11a2.vrt", paste(j, "3a_", 1:nrow(tiles), ".sdat", sep="", collapse=" ")))
  # convert to geotif (1 km):
  if(is.na(file.info(paste(j, '3a.tif', sep=""))$size)){
  system(paste(gdalwarp, ' mod11a2.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-32767\" -ot Int16 ', j, '3a.tif', ' -r near -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""))
  }
  # 2.5 km:
  if(is.na(file.info(paste(j, '2a.tif', sep=""))$size)){
  system(paste(gdalwarp, ' mod11a2.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-32767\" -ot Int16 ', j, '2a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/40,' ', 1/40, sep=""))
  }
  # 5 km:
  if(is.na(file.info(paste(j, '1a.tif', sep=""))$size)){
  system(paste(gdalwarp, ' mod11a2.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-32767\" -ot Int16 ', j, '1a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
  }
  # 20 km:
  system(paste(gdalwarp, ' ', j, '1a.tif', ' -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-32767\" -ot Int16 ', j, '0a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep="")) 
}

## compress:  
for(j in c("TDMMOD", "TDSMOD", "TDLMOD", "TDHMOD", paste("TX", 1:6, "MOD", sep=""))){
for(i in 0:3){
    outname = paste(j, i, 'a', sep="")
    system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), set.file.extension(outname, ".tif")))
    system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(outname, ".tif.gz"))
  }
}

# clean up temp files:
unlink(list.files(pattern=glob2rx("LST_Day_1km_*_*_*_*.*")))
unlink(list.files(pattern=glob2rx("TD*MOD3a_*.*")))
unlink(list.files(pattern=glob2rx("TX*MOD3a_*_f.*")))
unlink(list.files(pattern=glob2rx("TX*MOD3a_*.*")))

## list files:
ntif.lst <- list.files(pattern=glob2rx("LST_Night_1km_*_*_*.tif$"))

## NIGHT TIME - derive stats for each tile:
for(j in 1:nrow(tiles)){
  for(i in 1:length(ntif.lst)){
    outname <- strsplit(ntif.lst[i], ".tif")[[1]][1]
    tilename <- paste(outname, '_', j, '.sdat', sep="")
    # resample
    if(is.na(file.info(tilename)$size)){
      system(paste(gdalwarp, ' ', ntif.lst[i], ' ', tilename, ' -of \"SAGA\" -dstnodata \"-32767\" -r near -te ', tiles[j,"lonmin"] , ' ', tiles[j,"latmin"], ' ', tiles[j,"lonmax"] ,' ', tiles[j,"latmax"] ,' -tr ', 1/120, ' ', 1/120, sep=""))
    }
  }
  # derive mean, sd, min and max values:
  t.lst <- list.files(pattern=glob2rx(paste('LST_Night_1km_*_*_*_', j, '.sgrd$', sep="")))
  if(is.na(file.info(paste("TNMMOD3a_", j, ".sgrd", sep=""))$size)){
  rsaga.geoprocessor(lib="geostatistics_grid", module=4, param=list(GRIDS=paste(t.lst, collapse=";"), MEAN=paste("TNMMOD3a_", j, ".sgrd", sep=""), MIN=paste("TNLMOD3a_", j, ".sgrd", sep=""), MAX=paste("TNHMOD3a_", j, ".sgrd", sep=""), STDDEV=paste("TNSMOD3a_", j, ".sgrd", sep="")))
  }
}

# create a mosaic and compress files:
for(j in c("TNMMOD", "TNSMOD", "TNLMOD", "TNHMOD")){
  unlink("mod11a2.vrt")
  system(paste(gdalbuildvrt, "mod11a2.vrt", paste(j, "3a_", 1:nrow(tiles), ".sdat", sep="", collapse=" ")))
  # convert to geotif (1 km):
  if(is.na(file.info(paste(j, '3a.tif', sep=""))$size)){
  system(paste(gdalwarp, ' mod11a2.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-32767\" -ot Int16 ', j, '3a.tif', ' -r near -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""))
  }
  # 2.5 km:
  if(is.na(file.info(paste(j, '2a.tif', sep=""))$size)){
  system(paste(gdalwarp, ' mod11a2.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-32767\" -ot Int16 ', j, '2a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/40,' ', 1/40, sep=""))
  }
  # 5 km:
  if(is.na(file.info(paste(j, '1a.tif', sep=""))$size)){
  system(paste(gdalwarp, ' mod11a2.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-32767\" -ot Int16 ', j, '1a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
  }
  # 20 km:
  system(paste(gdalwarp, ' ', j, '1a.tif', ' -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-32767\" -ot Int16 ', j, '0a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))
}

## compress:  
for(j in c("TNMMOD", "TNSMOD", "TNLMOD", "TNHMOD")){
for(i in 0:3){
    outname = paste(j, i, 'a', sep="")
    system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), set.file.extension(outname, ".tif")))
    system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(outname, ".tif.gz"))
  }
}

# clean up temp files:
unlink(list.files(pattern=glob2rx("LST_Night_1km_*_*_*_*.*")))
unlink(list.files(pattern=glob2rx("TN*MOD3a_*.*")))

# end of script;