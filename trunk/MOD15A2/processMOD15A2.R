# title         : processMOD15A2.R
# purpose       : processing of the MODIS 1 km resolution 8-day LAI images (years 2001 and 2011);
# reference     : [http://worldgrids.org]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Nov 2012.
# inputs        : 1 km geotifs in the format "LAI_Day_1km_*_*_*.tif" (global composites produced using FWTools);
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
dtif.lst <- list.files(pattern=glob2rx("LAI_1km_*_*_*.tif$"))
GDALinfo(dtif.lst[1])

## LAI images:
for(j in 1:nrow(tiles)){
  for(i in 1:length(dtif.lst)){
    outname <- strsplit(dtif.lst[i], ".tif")[[1]][1]
    tilename <- paste(outname, '_', j, '.sdat', sep="")
    # resample
    if(is.na(file.info(tilename)$size)){
      system(paste(gdalwarp, ' ', dtif.lst[i], ' ', tilename, ' -of \"SAGA\" -dstnodata \"255\" -r near -te ', tiles[j,"lonmin"] , ' ', tiles[j,"latmin"], ' ', tiles[j,"lonmax"] ,' ', tiles[j,"latmax"] ,' -tr ', 1/120, ' ', 1/120, sep=""))
      ## mask out values outside 0-100 range:
      x = raster(tilename)
      x = calc(x, function(a){ifelse(a > 100, NA, a)})
      ## TH: if min(x) = max(x) then the writeRaster function fails;
      if(!is.na(minValue(x))&!is.na(maxValue(x))){
        writeRaster(x, tilename, format="SAGA", overwrite=TRUE, datatype="INT1U") 
      } else {
        unlink(tilename)
      }
    }
  }
  # derive derive mean, sd, min and max values:
  t.lst <- list.files(pattern=glob2rx(paste('LAI_1km_*_*_*_', j, '.sgrd$', sep="")))
  if(is.na(file.info(paste("LAMMOD3a_", j, ".sgrd", sep=""))$size)){
    rsaga.geoprocessor(lib="geostatistics_grid", module=4, param=list(GRIDS=paste(t.lst, collapse=";"), MEAN=paste("LAMMOD3a_", j, ".sgrd", sep=""), MIN=paste("LALMOD3a_", j, ".sgrd", sep=""), MAX=paste("LAHMOD3a_", j, ".sgrd", sep=""), STDDEV=paste("LASMOD3a_", j, ".sgrd", sep="")))
  }
}

# create a mosaic and compress files:
for(j in c("LAMMOD", "LASMOD", "LALMOD", "LAHMOD")){
  unlink("mod15a2.vrt")
  system(paste(gdalbuildvrt, "mod15a2.vrt", paste(j, "3a_", 1:nrow(tiles), ".sdat", sep="", collapse=" ")))
  # convert to geotif (1 km):
  if(is.na(file.info(paste(j, '3a.tif', sep=""))$size)){
  system(paste(gdalwarp, ' mod15a2.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-3000\" -ot Int16 ', j, '3a.tif', ' -r near -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""))
  }
  # 2.5 km:
  if(is.na(file.info(paste(j, '2a.tif', sep=""))$size)){
  system(paste(gdalwarp, ' mod15a2.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-3000\" -ot Int16 ', j, '2a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/40,' ', 1/40, sep=""))
  }
  # 5 km:
  if(is.na(file.info(paste(j, '1a.tif', sep=""))$size)){
  system(paste(gdalwarp, ' mod15a2.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-3000\" -ot Int16 ', j, '1a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
  }
  # 20 km:
  system(paste(gdalwarp, ' ', j, '1a.tif', ' -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-3000\" -ot Int16 ', j, '0a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))
}

## compress:  
for(j in c("LAMMOD", "LASMOD", "LALMOD", "LAHMOD")){
for(i in 0:3){
    outname = paste(j, i, 'a', sep="")
    system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), set.file.extension(outname, ".tif")))
    system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(outname, ".tif.gz"))
  }
}

# clean up temp files:
unlink(list.files(pattern=glob2rx("LAI_1km_*_*_*_*.*")))
unlink(list.files(pattern=glob2rx("LA*MOD3a_*.*")))

# end of script;