# title         : processMOD12Q1.R
# purpose       : processing of the MODIS Land Cover product (MOD12Q1:Land_Cover_Type_1);
# reference     : [http://worldgrids.org/doku.php?id=wiki:layers#modis_products]
# producer      : Prepared by T. Hengl
# addresss      : In Wageningen, NL.
# inputs        : data available for download via FTP from [http://ladsweb.nascom.nasa.gov/];
# outputs       : geotiff images projected in the "+proj=longlat +datum=WGS84" system;
# remarks 1     : First download and install FWtools [http://fwtools.maptools.org]; 

# -------------------------------------
# Initial settings and data download:
# -------------------------------------

library(rgdal)
library(RSAGA)
myenv <- rsaga.env(path="C:/saga_vc")
library(XML)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
gdalbuildvrt = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalbuildvrt.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"

## list all geotifs:
tif.lst <- list.files(pattern=glob2rx("GLC_1km*.tif$"))
tif.lst

## resample:
for(k in 1:length(tif.lst)){
  ki <- as.numeric(strsplit(tif.lst[k], "_")[[1]][3])-2000
  cln1 <- ifelse(ki<10, as.character(paste("0", ki, sep="")), as.character(ki))
  file.copy(tif.lst[k], paste('G',cln1,'IGB3a.tif', sep=""))
  # 2.5 km
  outn = paste('G',cln1,'IGB2a.tif', sep="")
  if(is.na(file.info(outn)$size)){
    system(paste(gdalwarp, ' ', tif.lst[k], ' ', outn, ' -srcnodata 255 -dstnodata 255 -r near -te -180 -90 180 90 -tr ', 1/40,' ', 1/40, sep=""))
  }
  # 5 km:
  outn = paste('G',cln1,'IGB1a.tif', sep="")
  if(is.na(file.info(outn)$size)){
    system(paste(gdalwarp, ' ', tif.lst[k], ' ', outn, ' -srcnodata 255 -dstnodata 255 -r near -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
  }
  # 20 km:
  outn = paste('G',cln1,'IGB0a.tif', sep="")
  if(is.na(file.info(outn)$size)){
    system(paste(gdalwarp, ' ', tif.lst[k], ' ', outn, ' -srcnodata 255 -dstnodata 255 -r near -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))
  }
}
## compress:
for(k in 1:length(tif.lst)){
  ki <- as.numeric(strsplit(tif.lst[k], "_")[[1]][3])-2000
  cln1 <- ifelse(ki<10, as.character(paste("0", ki, sep="")), as.character(ki))  
  for(i in 0:3){
    outname = paste('G',cln1,'IGB', i, 'a', sep="")
    system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), set.file.extension(outname, ".tif")))
    system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(outname, ".tif.gz"))
}}


## tiling system:
t.l <- expand.grid(KEEP.OUT.ATTRS=FALSE, lonmin=seq(-180,120,by=60), latmin=seq(-90,45,by=45))
t.u <- expand.grid(KEEP.OUT.ATTRS=FALSE, lonmax=seq(-120,180,by=60), latmax=seq(-45,90,by=45))
tiles <- cbind(t.l, t.u)

## make a PAL file:
library(RCurl)
# verify the certificate:
curl <- getCurlHandle()
options(RCurlOptions = list(capath = system.file("CurlSSL", "cacert.pem", package = "RCurl"), ssl.verifypeer = FALSE))
curlSetOpt(.opts = list(proxy = 'proxyserver:port'), curl = curl)

## legend entries:
cat(getURL("https://docs.google.com/spreadsheet/pub?key=0Ah5Ip0avaTkLdHhxTDVraGlZOU1maW8zeTQzRWFJd2c&single=true&gid=0&output=csv"), file = "legend.csv")
rgb.tbl <- read.csv("legend.csv")
str(rgb.tbl)
#rgb.tbl <- merge(levs.tbl, rgb.tbl, all.y=FALSE, all.x=TRUE)

## create a SAGA txt colour table
## convert to BGR codes:
BGR <- (rgb.tbl$B * 65536) + (rgb.tbl$G * 256) + rgb.tbl$R
## write a lookup table for SAGA GIS:
filename <- file("GXXIGB.txt", "w", blocking=FALSE)
write("COLOR\tNAME\tDESCRIPTION\tMINIMUM\tMAXIMUM", filename)
for(i in 1:nrow(rgb.tbl)){
  write(paste(BGR[i], rgb.tbl[i,"Category"], paste("CL", i, sep=""), (rgb.tbl[i,"Number"]-1)+0.1, (rgb.tbl[i,"Number"])+0.1, sep="\t"), filename, append=TRUE)
}
close(filename)

## write a PAL file:
write.table(rgb.tbl[,c("Category","R","G","B")], "GXXIGB.PAL", quote = FALSE,  col.names = FALSE)

## levels of special interest:
## Snow and ice -- 15
## Barren or sparsely vegetated -- 16

## split world in 24 blocks:
# THIS TAKES >5 mins
for(k in 1:length(tif.lst)){
  dout <- strsplit(tif.lst[k], "_")[[1]][3]
  for(j in 1:nrow(tiles)){
    if(is.na(file.info(paste('glc_', dout, '_', j, '.sdat', sep=""))$size)){  
      ## tile:
      system(paste(gdalwarp, ' ', tif.lst[k], ' -t_srs \"+proj=longlat +datum=WGS84\" glc_', dout, "_", j, '.sdat -of \"SAGA\" -ot \"Byte\" -r near -te ', tiles[j,"lonmin"] , ' ', tiles[j,"latmin"], ' ', tiles[j,"lonmax"] ,' ', tiles[j,"latmax"] , sep=""), show.output.on.console = FALSE)
    }
  }
}

## mask values per each class:
for(i in 1:16){
  cln <- ifelse(i<10, as.character(paste("0", i, sep="")), as.character(i))
  if(is.na(file.info(paste('L', cln, 'IGB3a.tif', sep=""))$size)){
  for(j in 1:nrow(tiles)){
    for(k in 1:length(tif.lst)){
    dout <- strsplit(tif.lst[k], "_")[[1]][3]
      if(is.na(file.info(paste('glc_', dout, "_", j, '_', i, '.sdat', sep=""))$size)){
        rsaga.geoprocessor(lib="grid_calculus", module=1, param=list(GRIDS=paste('glc_', dout, "_", j, '.sgrd', sep=""), RESULT=paste('glc_', dout, "_", j, '_', i, '.sgrd', sep=""), FORMULA=paste("ifelse(a=", i,",100,0)", sep="")), show.output.on.console = FALSE, env=myenv)
      }
    }
    ## list files:
    lst <- list.files(pattern=glob2rx(paste('glc_*', '_', j, '_', i, '.sgrd', sep="")))
    ## derive mean value for all years:
    if(is.na(file.info(paste("L",cln,"IGB3a_", j, ".sgrd", sep=""))$size)){
      rsaga.geoprocessor(lib="geostatistics_grid", module=4, param=list(GRIDS=paste(lst, collapse=";"), MEAN=paste("L",cln,"IGB3a_", j, ".sgrd", sep=""), MIN="tmp.sgrd", MAX="tmp.sgrd", STDDEV="tmp.sgrd"), show.output.on.console = FALSE, env=myenv)
    } 
  }
  ## create a mosaic:
  lst.mos <- list.files(pattern=glob2rx(paste("L",cln,"IGB3a_*.sgrd", sep="")))
  rsaga.geoprocessor(lib="grid_tools", module=3, param=list(GRIDS=paste(lst.mos, collapse=";"), TYPE=1, INTERPOL=0, OVERLAP=0, TARGET=0, USER_XMIN=-180+(1/120)/2, USER_XMAX=180-(1/120)/2, USER_YMIN=-90+(1/120)/2, USER_YMAX=90-(1/120)/2, USER_SIZE=1/120, USER_GRID="mos.sgrd"), env=myenv, show.output.on.console = FALSE)
  ## resample:
  system(paste(gdal_translate, ' mos.sdat ', 'L', cln, 'IGB3a.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -a_nodata \"255\"', sep=""))
}
}

l.lst <- list.files(pattern=glob2rx("L*IGB3a.tif$"))
l.lst

## resample:
for(k in 1:length(l.lst)){
  cln1 <- ifelse(k<10, as.character(paste("0", k, sep="")), as.character(k))
  # 2.5 km
  outn = paste('L',cln1,'IGB2a.tif', sep="")
  if(is.na(file.info(outn)$size)){
    system(paste(gdalwarp, ' ', l.lst[k], ' ', outn, ' -r bilinear -te -180 -90 180 90 -tr ', 1/40,' ', 1/40, sep=""))
  }
  # 5 km:
  outn = paste('L',cln1,'IGB1a.tif', sep="")
  if(is.na(file.info(outn)$size)){
    system(paste(gdalwarp, ' ', l.lst[k], ' ', outn, ' -r bilinear -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
  }
  # 20 km:
  outn = paste('L',cln1,'IGB0a.tif', sep="")
  if(is.na(file.info(outn)$size)){
    system(paste(gdalwarp, ' ', l.lst[k], ' ', outn, ' -r bilinear -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))
  }
}
## compress:
for(k in 1:length(l.lst)){
  cln1 <- ifelse(k<10, as.character(paste("0", k, sep="")), as.character(k))  
  for(i in 0:3){
    outname = paste('L',cln1,'IGB', i, 'a', sep="")
    system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), set.file.extension(outname, ".tif")))
    system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(outname, ".tif.gz"))
}}


## clean up:
unlink("glc_20**_**_**.***")
unlink("L**IGB3a_**.***")
unlink("glc_20**_**.***")
unlink("L**IGB3a_*.***")  
unlink("mos.***")
unlink("*.xml")

## end of script;