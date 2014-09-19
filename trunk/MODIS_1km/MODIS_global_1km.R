# title         : MODIS_global_1km.R
# purpose       : Processing MODIS images 1 km
# reference     : Worldgrids.org;
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Feb 2012.
# inputs        : MOD11A2 and MOD13A2 for years 2001 and 2011 (downloaded from [ftp://e4ftl01.cr.usgs.gov]);
# outputs       : geotiff images projected in the WGS84 coordinate system;
# remarks 1     : This script is available from www.worldgrids.org;
# remarks 2     : the whole world is 648 MODIS tiles in the Sinusoidal projection;
# remarks 3     : this code is Windows OS specific;

## 1. GENERAL SETTINGS:

library(rgdal)
library(RSAGA)
myenv <- rsaga.env(path="C:/saga_vc")
# library(modis)
fw.path <- utils::readRegistry("SOFTWARE\\WOW6432NODE\\FWTools")$Install_Dir
gdalwarp <- shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate <- shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
gdalbuildvrt <- shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalbuildvrt.exe"))))
gdalinfo <- shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalinfo.exe"))))
workd <- normalizePath(getwd()) 
outdir <- "Z:/1km"
LAI.dir = "Z:/WORLDGRIDS/MOD15A2"
GLC.dir = "Z:/WORLDGRIDS/MOD12Q1"
EVI.dir = "Z:/WORLDGRIDS/MOD13A3"
# MODIS projection system:
modis.prj <- "+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6371007.181 +b=6371007.181 +units=m +no_defs"
crs <- "+proj=longlat +datum=WGS84"

## Get 7z:
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")

## ----------------------------------
## 2. LIST LOCAL HDF FILES:
## ----------------------------------

# list products:
pl <- list.dirs(path=getwd())
pl2 <- list.dirs()
# mask out all directories not of interest:
pl <- pl[grep(pl2, pattern="MOD*")]
pl <- pl[sapply((sapply(pl, strsplit, "/")), length) > (length(strsplit(getwd(), "/")[[1]])+1)]
## list all files within each date:
fl <- NULL
for(j in 1:length(pl)){
   fl[[j]] <- list.files(pl[j], pattern="*.hdf$")
}
# save temp data:
save.image(".RData")

## ----------------------------------
## 3. MOSAICK AND RESAMPLE:
## ----------------------------------

## LST images: 
lst.sel <- grep(pl, pattern="MOD11A2")
lst.sel <- lst.sel[sapply(pl[lst.sel], function(x){length(strsplit(x, "/")[[1]])>4})]

for(j in lst.sel){    # for each date;

 ## output file name:
 outname <- strsplit(pl[j], "/")[[1]]
 ## check if numeric or date?
 if(!is.na(as.numeric(outname[length(outname)]))){
   f_outname <- format(as.Date(paste(outname[length(outname)-1], outname[length(outname)], sep="-"), format="%Y-%j"), "%Y_%m_%d")
   doutname <- paste('LST_Day_1km_', f_outname, '.tif', sep="")
   noutname <- paste('LST_Night_1km_', f_outname, '.tif', sep="") 
 } else {
  doutname <- paste('LST_Day_1km_', gsub("\\.", "_", outname[length(outname)]),'.tif', sep="")
   noutname <- paste('LST_Night_1km_', gsub("\\.", "_", outname[length(outname)]),'.tif', sep="") 
 } 

 ## check if the file exists already: 
 if(is.na(file.info(doutname)$size)|is.na(file.info(noutname)$size)){

  hddir <- normalizePath(pl[j])
  for(i in 1:length(fl[[j]])){  # for each HDF file:
  # remove all corrupted HDF files:
  inf <- system(paste(gdalinfo, ' ', hddir, '\\', fl[[j]][i], sep=""), intern = TRUE)
  if(any(unlist(strsplit(inf, " ")) %in% "ERROR")){
     unlink(paste(hddir, '\\', fl[[j]][i], sep=""))
  } 
  ## generate a name:
  LSTd.name <- paste(gsub("\\.", "_", strsplit(fl[[j]][i], "005")[[1]][1]), "LST_Day_1km", sep="") 
  LSTn.name <- paste(gsub("\\.", "_", strsplit(fl[[j]][i], "005")[[1]][1]), "LST_Night_1km", sep="")

  ## extract Day time LST:
  if(is.na(file.info(paste(LSTd.name, ".tif", sep=""))$size)){
     try(system(paste(gdal_translate, ' HDF4_EOS:EOS_GRID:\"', hddir, '\\', fl[[j]][i], '\":MODIS_Grid_8Day_1km_LST:LST_Day_1km -a_nodata \"0\" ', LSTd.name, '.tif', sep=""), show.output.on.console = FALSE))
  }
  # extract night time LST:
  if(is.na(file.info(paste(LSTn.name, ".tif", sep=""))$size)){
     try(system(paste(gdal_translate, ' HDF4_EOS:EOS_GRID:\"', hddir, '\\', fl[[j]][i], '\":MODIS_Grid_8Day_1km_LST:LST_Night_1km -a_nodata \"0\" ', LSTn.name, '.tif', sep=""), show.output.on.console = FALSE))
  }
 } ## end extraction of bands (648 tiles) - it takes cca 20 MINS per date!!;

 ## list all individual images:
 tifd.lst <- list.files(pattern="*LST_Day_1km.tif$")

 ## DAY TIME IMAGES - generate a VRT (virtual mosaic):
 system(paste(gdalbuildvrt, "Day.vrt *LST_Day_1km.tif"))
 #system(paste(gdalwarp, ' Day.vrt ', doutname, ' -s_srs \"', modis.prj, '\" -t_srs \"', crs, '\" -r bilinear -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""))   ## this would be the easiest, but there is no need for such a precision!
 ## Create 4 mosaicks
 system(paste(gdalwarp, ' Day.vrt Ablock.tif -s_srs \"', modis.prj, '\" -t_srs \"', crs, '\" -r bilinear -te -180 -90 0 0 -tr ', 1/120,' ', 1/120, sep=""), show.output.on.console = FALSE)
 system(paste(gdalwarp, ' Day.vrt Bblock.tif -s_srs \"', modis.prj, '\" -t_srs \"', crs, '\" -r bilinear -te 0 -90 180 0 -tr ', 1/120,' ', 1/120, sep=""), show.output.on.console = FALSE)
 system(paste(gdalwarp, ' Day.vrt Cblock.tif -s_srs \"', modis.prj, '\" -t_srs \"', crs, '\" -r bilinear -te -180 0 0 90 -tr ', 1/120,' ', 1/120, sep=""), show.output.on.console = FALSE)
 system(paste(gdalwarp, ' Day.vrt Dblock.tif -s_srs \"', modis.prj, '\" -t_srs \"', crs, '\" -r bilinear -te 0 0 180 90 -tr ', 1/120,' ', 1/120, sep=""), show.output.on.console = FALSE)  ##  3-5 MINS
 ## for each block:
 for(k in c("Ablock","Bblock","Cblock","Dblock")){
   # Convert to a SAGA GIS file:
   rsaga.geoprocessor(lib="io_gdal", module=0, param=list(GRIDS=set.file.extension(k, ".sgrd"), FILES=set.file.extension(k, ".tif")), show.output.on.console = FALSE)
   ## fix the NA values: NODATA_VALUE = 0
   sgrd <- matrix((unlist(strsplit(readLines(file(set.file.extension(k, ".sgrd"))), split="\t= "))), ncol=2, byrow=T)
   sgrd[12,2] <- "0"
   filec <- file(set.file.extension(k, ".sgrd"))
   write.table(sgrd, file=filec, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t= ")
   ## Convert to degree celsius and round the numbers:
   rsaga.geoprocessor(lib="grid_calculus", module=1, param=list(GRIDS=set.file.extension(k, ".sgrd"), RESULT=paste(k, "f.sgrd", sep="_"), FORMULA="a*0.02-273.15"), show.output.on.console = FALSE)
   ## Convert to a geotiff:
   system(paste(gdal_translate, ' ', paste(k, "f.sdat", sep="_"), ' ', paste(k, "f.tif", sep="_"), ' -ot \"Int16\" -a_srs \"', crs, '\" ', sep=""), show.output.on.console = FALSE)
 }
 ## GENERATE A MOSAICK WITH PROPER CRS:
 system(paste(gdalbuildvrt, "Day_f.vrt *block_f.tif"))
 system(paste(gdalwarp, ' Day_f.vrt ', doutname, ' -s_srs \"', crs, '\" -t_srs \"', crs, '\" -ot \"Int16\"  -srcnodata \"-32768\" -dstnodata \"-32768\" -r near -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""), show.output.on.console = FALSE)   ## 5 MINS
  ## recycle temporary file:
 unlink("Day.vrt");  unlink("Day_f.vrt")
 unlink(list.files(pattern="*block*")); unlink(tifd.lst) 
 ## Zip the output files:
 system(paste("7za a", "-tgzip", set.file.extension(doutname, ".tif.gz"), doutname)) 
 system(paste("xcopy", set.file.extension(doutname, ".tif.gz"), shortPathName(normalizePath(outdir))))
 unlink(set.file.extension(doutname, ".tif.gz"))

 ## list all individual images (night):
 tifn.lst <- list.files(pattern="*LST_Night_1km.tif$")
 
 ## NIGHT TIME IMAGES - generate a VRT (virtual mosaic):
 system(paste(gdalbuildvrt, "Night.vrt *LST_Night_1km.tif"))
 ## Create 4 mosaicks
 system(paste(gdalwarp, ' Night.vrt Ablock.tif -s_srs \"', modis.prj, '\" -t_srs \"', crs, '\" -r bilinear -te -180 -90 0 0 -tr ', 1/120,' ', 1/120, sep=""), show.output.on.console = FALSE)
 system(paste(gdalwarp, ' Night.vrt Bblock.tif -s_srs \"', modis.prj, '\" -t_srs \"', crs, '\" -r bilinear -te 0 -90 180 0 -tr ', 1/120,' ', 1/120, sep=""), show.output.on.console = FALSE)
 system(paste(gdalwarp, ' Night.vrt Cblock.tif -s_srs \"', modis.prj, '\" -t_srs \"', crs, '\" -r bilinear -te -180 0 0 90 -tr ', 1/120,' ', 1/120, sep=""), show.output.on.console = FALSE)
 system(paste(gdalwarp, ' Night.vrt Dblock.tif -s_srs \"', modis.prj, '\" -t_srs \"', crs, '\" -r bilinear -te 0 0 180 90 -tr ', 1/120,' ', 1/120, sep=""), show.output.on.console = FALSE)  ##  3-5 MINS
 ## for each block:
 for(k in c("Ablock","Bblock","Cblock","Dblock")){
   ## Convert to a SAGA GIS file:
   rsaga.geoprocessor(lib="io_gdal", module=0, param=list(GRIDS=set.file.extension(k, ".sgrd"), FILES=set.file.extension(k, ".tif")), show.output.on.console = FALSE)
   ## fix the NA values: NODATA_VALUE = 0
   sgrd <- matrix((unlist(strsplit(readLines(file(set.file.extension(k, ".sgrd"))), split="\t= "))), ncol=2, byrow=T)
   sgrd[12,2] <- "0"
   filec <- file(set.file.extension(k, ".sgrd"))
   write.table(sgrd, file=filec, quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t= ")
   ## Convert to degree celsius and round the numbers:
   rsaga.geoprocessor(lib="grid_calculus", module=1, param=list(GRIDS=set.file.extension(k, ".sgrd"), RESULT=paste(k, "f.sgrd", sep="_"), FORMULA="a*0.02-273.15"), show.output.on.console = FALSE)
   ## Convert to a geotiff:
   system(paste(gdal_translate, ' ', paste(k, "f.sdat", sep="_"), ' ', paste(k, "f.tif", sep="_"), ' -ot \"Int16\" -a_srs \"', crs, '\" ', sep=""))
 }
 ## GENERATE A MOSAICK WITH PROPER CRS:
 system(paste(gdalbuildvrt, "Night_f.vrt *block_f.tif"))
 system(paste(gdalwarp, ' Night_f.vrt ', noutname, ' -s_srs \"', crs, '\" -t_srs \"', crs, '\" -ot \"Int16\"  -srcnodata \"-32768\" -dstnodata \"-32768\" -r near -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""), show.output.on.console = FALSE)   ## 5 MINS
 ## recycle temporary file:
 unlink("Night.vrt");  unlink("Night_f.vrt")
 unlink(list.files(pattern="*block*")); unlink(tifn.lst) 
 ## Zip the output files:
 system(paste("7za a", "-tgzip", set.file.extension(noutname, ".tif.gz"), noutname)) 
 system(paste("xcopy", set.file.extension(noutname, ".tif.gz"), shortPathName(normalizePath(outdir))))
 unlink(set.file.extension(noutname, ".tif.gz"))

 }
}


## MODIS/Terra Vegetation Indices Monthly L3 Global 1km SIN Grid: 
evi.sel <- grep(pl, pattern="MOD13A3")

for(j in evi.sel){    # for each date;

 # check if the file exists already:
 outname <- strsplit(pl[j], "/")[[1]]
 doutname <- paste('EVI_1km_', gsub("\\.", "_", outname[length(outname)]),'.tif', sep="") 
 
 if(is.na(file.info(doutname)$size)){
 hddir <- normalizePath(pl[j])

  for(i in 1:length(fl[[j]])){  # for each HDF file:
  # remove all corrupted HDF files:
  inf <- system(paste(gdalinfo, ' ', hddir, '\\', fl[[j]][i], sep=""), intern = TRUE)
  if(any(unlist(strsplit(inf, " ")) %in% "ERROR")){
     unlink(paste(hddir, '\\', fl[[j]][i], sep=""))
  } 
  # generate a name:
  EVI.name <- paste(gsub("\\.", "_", strsplit(fl[[j]][i], "005")[[1]][1]), "EVI_1km", sep="") 

  # extract EVI band:
  if(is.na(file.info(paste(EVI.name, ".tif", sep=""))$size)){
     try(system(paste(gdal_translate, ' HDF4_EOS:EOS_GRID:\"', hddir, '\\', fl[[j]][i], '\":MOD_Grid_monthly_1km_VI:\"1 km monthly EVI\" -a_nodata \"0\" ', EVI.name, '.tif', sep=""), show.output.on.console = FALSE))
  }
 } # end extraction of bands (648 tiles) - it takes cca 20 MINS per date!!;
                              
 # list all individual images:
 tif.lst <- list.files(pattern="*EVI_1km.tif$")

 # generate a VRT (virtual mosaic):
 system(paste(gdalbuildvrt, "EVI.vrt *EVI_1km.tif"))
 # Create a mosaick:
 system(paste(gdalwarp, ' EVI.vrt ', doutname, ' -s_srs \"', modis.prj, '\" -t_srs \"', crs, '\" -ot \"Int16\"  -srcnodata \"-3000\" -dstnodata \"-3000\" -r bilinear -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""), show.output.on.console = FALSE)  ## 5 MINS
 # recycle temporary file:
 unlink("EVI.vrt")
 unlink(tif.lst) 
 # Zip the output files:
 system(paste("7za a", "-tgzip", set.file.extension(doutname, ".tif.gz"), doutname)) 
 system(paste("xcopy", set.file.extension(doutname, ".tif.gz"), shortPathName(normalizePath(outdir))))
 unlink(set.file.extension(doutname, ".tif.gz"))

 }
}

## MODIS Leaf Area Index: 
lai.sel <- grep(pl, pattern="MOD15A2")

for(j in lai.sel){    # for each date;

 # check if the file exists already:
 outname <- strsplit(pl[j], "/")[[1]]
 doutname <- paste(LAI.dir, '/', 'LAI_1km_', gsub("\\.", "_", outname[length(outname)]),'.tif', sep="") 
 
 if(is.na(file.info(doutname)$size)){
 hddir <- normalizePath(pl[j])

  for(i in 1:length(fl[[j]])){  # for each HDF file:
  # remove all corrupted HDF files:
  inf <- system(paste(gdalinfo, ' ', hddir, '\\', fl[[j]][i], sep=""), intern = TRUE)
  if(any(unlist(strsplit(inf, " ")) %in% "ERROR")){
     unlink(paste(hddir, '\\', fl[[j]][i], sep=""))
  } 
  # generate a name:
  LAI.name <- paste(gsub("\\.", "_", strsplit(fl[[j]][i], "005")[[1]][1]), "LAI_1km", sep="") 

  # extract LAI band:
  if(is.na(file.info(paste(LAI.dir, "/", LAI.name, ".tif", sep=""))$size)){
     try(system(paste(gdal_translate, ' HDF4_EOS:EOS_GRID:\"', hddir, '\\', fl[[j]][i], '\":MOD_Grid_MOD15A2:Lai_1km -a_nodata 255 ', LAI.name, '.tif', sep=""), show.output.on.console = FALSE))
  }
 } # end extraction of bands (648 tiles) - it takes cca 20 MINS per date!!;
                              
 # list all individual images:
 tif.lst <- list.files(pattern="*LAI_1km.tif$")

 # generate a VRT (virtual mosaic):
 system(paste(gdalbuildvrt, "LAI.vrt *LAI_1km.tif"))
 # Create a mosaick:
 system(paste(gdalwarp, ' LAI.vrt ', doutname, ' -s_srs \"', modis.prj, '\" -t_srs \"', crs, '\" -srcnodata 255 -dstnodata 255 -r near -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""), show.output.on.console = FALSE)  ## 5 MINS
 # recycle temporary file:
 unlink("LAI.vrt")
 unlink(tif.lst) 
 # Zip the output files:
 system(paste("7za a", "-tgzip -mx=9", set.file.extension(doutname, ".tif.gz"), doutname)) 
 system(paste("xcopy", normalizePath(set.file.extension(doutname, ".tif.gz")), shortPathName(normalizePath(outdir))))
 unlink(set.file.extension(doutname, ".tif.gz"))

 }
}


## MODIS Land cover images: 
glc.sel <- grep(pl, pattern="MOD12Q1")

for(j in glc.sel){    # for each date;

 # check if the file exists already:
 outname <- strsplit(pl[j], "/")[[1]]
 doutname <- paste(GLC.dir, '/', 'GLC_1km_', gsub("\\.", "_", outname[length(outname)]),'.tif', sep="") 
 
 if(is.na(file.info(doutname)$size)){
 hddir <- normalizePath(pl[j])

  for(i in 1:length(fl[[j]])){  # for each HDF file:
  # remove all corrupted HDF files:
  inf <- system(paste(gdalinfo, ' ', hddir, '\\', fl[[j]][i], sep=""), intern = TRUE)
  if(any(unlist(strsplit(inf, " ")) %in% "ERROR")){
     unlink(paste(hddir, '\\', fl[[j]][i], sep=""))
  } 
  # generate a name:
  GLC.name <- paste(gsub("\\.", "_", strsplit(fl[[j]][i], "\\.004\\.")[[1]][1]), "_GLC_1km", sep="") 

  # extract IGBP band:
  if(is.na(file.info(paste(GLC.dir, "/", GLC.name, ".tif", sep=""))$size)){
     try(system(paste(gdal_translate, ' HDF4_EOS:EOS_GRID:\"', hddir, '\\', fl[[j]][i], '\":MOD12Q1:Land_Cover_Type_1 -a_nodata 255 ', GLC.name, '.tif', sep=""), show.output.on.console = FALSE))
  }
 } # end extraction of bands (317 tiles);
                              
 # list all individual images:
 tif.lst <- list.files(pattern=glob2rx("*_GLC_1km.tif$"))
 unlink("my_liste.txt")
 cat(tif.lst, file="my_liste.txt", fill=TRUE, append=TRUE)
 ## generate a VRT (virtual mosaic):
 system(paste(gdalbuildvrt, "-input_file_list my_liste.txt GLC.vrt"))
 ## Create a mosaick:
 system(paste(gdalwarp, ' GLC.vrt ', doutname, ' -s_srs \"', modis.prj, '\" -t_srs \"', crs, '\" -srcnodata 255 -dstnodata 255 -r near -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""), show.output.on.console = FALSE)  ## 5 MINS
  
 ## recycle temporary file:
 unlink("GLC.vrt")
 unlink(tif.lst) 

 }
}

## ----------------------------------
## 4. OTHER PROCESSING:
## ----------------------------------

# OPTIONAL: reproject to some local coordinate system:


# end of the script;
