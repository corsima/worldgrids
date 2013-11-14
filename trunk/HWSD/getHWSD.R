# title         : getHWSD.R
# purpose       : download and resampling of the Harmonized World Soil DB 7 March, 2012 - Version 1.21  (FAO Soil map of the World);
# reference     : [https://code.google.com/p/worldgrids/source/browse/]
# producer      : Prepared by T. Hengl
# version       : 1
# address       : In Wageningen, NL.
# inputs        : data available for download from [http://www.iiasa.ac.at/Research/LUC/External-World-soil-database/HTML/]; 
# outputs       : geotiff images projected in the "+proj=longlat +ellps=WGS84" system;
# remarks 1     : First download and install SAGA GIS [http://www.saga-gis.org] and FWtools [http://fwtools.maptools.org];  
# remarks 2     : Merging <5 km grids and table values requires >4GB and high performance computing;
# remarks 3     : Detailed description of the database available at [http://www.iiasa.ac.at/Research/LUC/External-World-soil-database/HWSD_Documentation.pdf]; the output maps are available from [http://worldgrids.org/doku.php?id=wiki:layers#harmonized_world_soil_database_images]  


# ------------------------------------------------------------
# Initial settings and data download:
# ------------------------------------------------------------

library(RODBC)
library(rgdal)
library(GSIF)
library(RSAGA)
library(snowfall)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdalbuildvrt <- shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalbuildvrt.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
outdir <- "G:/WORLDGRIDS/maps"

# donwload the raster image and the DB:
download.file("http://webarchive.iiasa.ac.at/Research/LUC/External-World-soil-database/HWSD_Data/HWSD_RASTER.zip", destfile=paste(getwd(), "HWSD_RASTER.zip", sep="/"))
unzip(zipfile="HWSD_RASTER.zip", exdir=getwd()) # 1.5 GB image!!
download.file("http://webarchive.iiasa.ac.at/Research/LUC/External-World-soil-database/HWSD_Data/HWSD.mdb", destfile=paste(getwd(), "HWSD.mdb", sep="/"))  
download.file("http://webarchive.iiasa.ac.at/Research/LUC/External-World-soil-database/HWSD_Data/HWSD_META.mdb", destfile=paste(getwd(), "HWSD_META.mdb", sep="/"))


## download soil mask:
download.file("http://worldgrids.org/lib/exe/fetch.php?media=smkisr3a.tif.gz", destfile=paste(getwd(), "/", "smkisr3a.tif.gz", sep=""))
system("7za e smkisr3a.tif.gz")
GDALinfo("hwsd.bil")

## tiling system:
grd <- GDALinfo("SMKISR3a.tif")
tiles <- getSpatialTiles(grd, block.x=30, block.y=30)
str(tiles)

unlink("hwsd.tif")
for(j in 1:nrow(tiles)){
  ## tile HWSD map:
  HWSDtile <- paste("hwsd_T", j, ".sdat", sep="")
  if(is.na(file.info(HWSDtile)$size)){
    system(paste(gdalwarp, ' hwsd.bil ', HWSDtile, ' -s_srs \"+proj=longlat +ellps=WGS84\" -t_srs \"+proj=longlat +ellps=WGS84\" -of \"SAGA\" -srcnodata 0 -dstnodata 65535 -r near -te ', tiles[j,"xl"] , ' ', tiles[j,"yl"], ' ', tiles[j,"xu"] ,' ', tiles[j,"yu"] ,' -tr ', 1/120, ' ', 1/120, sep=""))
  }
}

## resample also to 5 km resolution:
system(paste(gdalwarp, " hwsd.bil -t_srs \"+proj=longlat +ellps=WGS84\" hwsd.tif -ot \"Byte\" -dstnodata 0 -r near -te -180 -90 180 90 -tr 0.05 0.05"))
unlink("hwsd.bil")

## import the HWSD database to R:
cHWSD <- odbcConnect(dsn="HWSD")
cHWSD_md <- odbcConnect(dsn="HWSD_META")
tbl.list <- sqlTables(cHWSD)$TABLE_NAME
# get soil classes:
HWSD.SMU <- sqlQuery(cHWSD, query="SELECT ID, MU_GLOBAL, SU_SYMBOL FROM HWSD_SMU") # soil mapping units and soil types (37);
str(HWSD.SMU)

# export the reclassification table:
HWSD.SMU$minimum <- round(HWSD.SMU$MU_GLOBAL-0.5, 1)
HWSD.SMU$maximum <- round(HWSD.SMU$MU_GLOBAL+0.5, 1)
HWSD.SMU$MU_GLOBALc <- as.factor(HWSD.SMU$MU_GLOBAL)
HWSD.SMU$newclass <- as.integer(HWSD.SMU$SU_SYMBOL)
write.table(HWSD.SMU[,c("MU_GLOBAL", "newclass")], "class_SMU.txt", sep = "\t", quote=FALSE, row.names=FALSE)

## get all attributes:
HWSD_DATA <- sqlFetch(cHWSD, "HWSD_DATA", na.strings=c("NA", "<NA>", ""))
str(HWSD_DATA)
## test some summaries:
#summary(HWSD_DATA$T_OC); summary(HWSD_DATA$T_CACO3); summary(HWSD_DATA$T_CASO4); summary(HWSD_DATA$T_CASO4)

## Fetch metadata
HWSD_DATA.mt <- sqlColumns(cHWSD, "HWSD_DATA")  
# description tables:
mt.list <- tbl.list[grep(tbl.list, pattern="^D_")]
HWSD.mt_list <- NULL
for(j in 1:length(mt.list)){ HWSD.mt_list[[j]] <- sqlFetch(cHWSD, mt.list[j]) }

# get metadata:
HWSD_META <- sqlFetch(cHWSD_md, "HWSD_METADATA")
str(HWSD_META)

# Soil properties of interest:
sv.list <- c("MU_GLOBAL", names(HWSD_DATA)[-c(1:13)])
# reformat values:
HWSD_DATA$MU_GLOBAL <- as.factor(HWSD_DATA$MU_GLOBAL)
HWSD_DATA$T_TEXTURE <- as.factor(HWSD_DATA$T_TEXTURE)
HWSD_DATA$T_USDA_TEX_CLASS <- as.factor(HWSD_DATA$T_USDA_TEX_CLASS)
HWSD_DATA$S_USDA_TEX_CLASS <- as.factor(HWSD_DATA$S_USDA_TEX_CLASS)
HWSD_DATA$ADD_PROP <- as.factor(HWSD_DATA$ADD_PROP)
HWSD_DATA$DRAINAGE <- as.factor(HWSD_DATA$DRAINAGE)
HWSD_DATA$PHASE1 <- as.factor(HWSD_DATA$PHASE1)
HWSD_DATA$PHASE2 <- as.factor(HWSD_DATA$PHASE2)
HWSD_DATA$SWR <- ifelse(HWSD_DATA$SWR==0, NA, HWSD_DATA$SWR)
HWSD_DATA$ROOTS <- ifelse(HWSD_DATA$ROOTS==0, NA, HWSD_DATA$ROOTS)
HWSD_DATA$IL <- ifelse(HWSD_DATA$IL==0, NA, HWSD_DATA$IL)
HWSD_DATA$x <- 1/HWSD_DATA$SHARE
# find the dominant class (highest percentage):
rnk <- aggregate(. ~ MU_GLOBAL, HWSD_DATA[order(HWSD_DATA$MU_GLOBAL), c("MU_GLOBAL", "x")], FUN=rank, ties.method="first", na.action = na.pass)
HWSD_DATA[order(HWSD_DATA$MU_GLOBAL),"x"] <- unlist(rnk$x)
str(HWSD_DATA)
## 48148 values for 16327 MU_GLOBAL IDs, so the values need to be aggregated based on "SHARE";

# ------------------------------------------------------------
# generate percentage cover per each WRB soil type:
# ------------------------------------------------------------

## clean up classes:
#write.csv(data.frame(number=1:length(levels(HWSD_DATA$SU_SYM74)), names=levels(HWSD_DATA$SU_SYM74)), "SU_SYM74.csv")
#write.csv(data.frame(number=1:length(levels(HWSD_DATA$SU_SYM90)), names=levels(HWSD_DATA$SU_SYM90)), "SU_SYM90.csv")
library(RCurl)
# verify the certificate:
curl <- getCurlHandle()
options(RCurlOptions = list(capath = system.file("CurlSSL", "cacert.pem", package = "RCurl"), ssl.verifypeer = FALSE))
curlSetOpt(.opts = list(proxy = 'proxyserver:port'), curl = curl)

## HWSDB classification system:
cat(getURL("https://docs.google.com/spreadsheet/pub?key=0Ah5Ip0avaTkLdG54OU1STllDZnM2TEk5WXhvS21JdlE&single=true&gid=0&output=csv"), file = "WRBtbl.csv")
WRB.tbl <- read.csv("WRBtbl.csv", na.strings = c("NA",""))

cat(getURL("https://docs.google.com/spreadsheet/pub?key=0Ah5Ip0avaTkLdG54OU1STllDZnM2TEk5WXhvS21JdlE&single=true&gid=5&output=csv"), file = "cleanup_SU_SYM74.csv")
cleanup_SU_SYM74 <- read.csv("cleanup_SU_SYM74.csv", na.strings = c("NA",""))
str(cleanup_SU_SYM74)

cat(getURL("https://docs.google.com/spreadsheet/pub?key=0Ah5Ip0avaTkLdG54OU1STllDZnM2TEk5WXhvS21JdlE&single=true&gid=6&output=csv"), file = "cleanup_SU_SYM90.csv")
cleanup_SU_SYM90 <- read.csv("cleanup_SU_SYM90.csv", na.strings = c("NA",""))
str(cleanup_SU_SYM90)

WRB1 <- merge(HWSD_DATA[,c("MU_GLOBAL","SU_SYM74","SHARE")], cleanup_SU_SYM74[,c("names","Group")], by.x="SU_SYM74", by.y="names", all=FALSE)
WRB2 <- merge(HWSD_DATA[,c("MU_GLOBAL","SU_SYM90","SHARE")], cleanup_SU_SYM90[,c("SYM90","Group")], by.x="SU_SYM90", by.y="SYM90", all=FALSE)
WRB <- rbind(WRB1[,c("MU_GLOBAL","SHARE","Group")], WRB2[,c("MU_GLOBAL","SHARE","Group")])
## some 100 MU_GLOBAL do not have a match... this can be ingored

## merge per WRB group per tile:
wrapper.HWSDtile <- function(i){
  HWSDtile <- paste("hwsd_T", i, ".sdat", sep="")
  HWSDRaster <- readGDAL(HWSDtile)
  names(HWSDRaster) <- "MU_GLOBAL"
  ## skip if empty set...
 if(!all(is.na(HWSDRaster$MU_GLOBAL))){
  HWSDRaster <- as(HWSDRaster, "SpatialPixelsDataFrame")
  HWSDRaster$MU_GLOBAL <- as.factor(HWSDRaster$MU_GLOBAL)
  HWSDRaster$index <- as.factor(HWSDRaster@grid.index)
  ## takes time and RAM!!
  gc()
  ## check coverage for each class:
  for(k in levels(WRB$Group)){
    HWSDouttile <- paste("hwsd_T", i, "_", k, ".tif", sep="")   
    if(!file.exists(HWSDouttile)){
      HWSDRaster$SHARE <- NULL
      sel <- WRB$Group == k & !is.na(WRB$Group)
      x <- merge(x=HWSDRaster@data, y=WRB[sel,], all.x=TRUE, all.y=FALSE, by="MU_GLOBAL", sort=FALSE)
      gc()
      if(all(is.na(x$SHARE))){
        HWSDRaster$SHARE <- 0
      } else {
        ## Some classes appear in the same MU_GLOBAL hence we need to aggregate!!
        xa <- aggregate(x$SHARE, by=list(x$index), sum)
        gc()
        HWSDRaster@data[as.integer(xa$Group.1),"SHARE"] <- xa$x
        HWSDRaster$SHARE <- ifelse(is.na(HWSDRaster$SHARE), 0, HWSDRaster$SHARE)
      } 
      writeGDAL(HWSDRaster["SHARE"], HWSDouttile, type="Byte", mvFlag=255, "GTiff")
    }
  }
 }
}

sfInit(parallel=TRUE, cpus=6)
sfLibrary(rgdal)
sfLibrary(sp)

## export all objects that are used in this function:
sfExport("tiles", "WRB")
system.time(x <- sfLapply(1:nrow(tiles), wrapper.HWSDtile))
sfStop()

## Create mosaics:
for(k in levels(WRB$Group)){
  tifout <- paste("G", WRB.tbl$Code[which(WRB.tbl$Group == k)], "HWS3a.tif", sep="")
  if(!file.exists(tifout)){
    ## list all tiles of specific class:
    #tif.lst <- list.files(pattern=paste("hwsd_*_", k,".tif$"))
    ## create a mosaic:
    unlink("hwsd.vrt")
    system(paste(gdalbuildvrt, " hwsd.vrt ", " hwsd_*_", k, ".tif", sep=""))
    system(paste(gdalwarp, " hwsd.vrt ", tifout, " -ot \"Byte\" -dstnodata 255 -r near -te -180 -90 180 90 -tr ", 1/120," ", 1/120, sep=""), show.output.on.console = FALSE)
  }
}
## check if everything is OK:
GDALinfo("GACHWS3a.tif")

## resample to various resolutions
## create a mosaic and compress files:
for(k in levels(WRB$Group)){
  tifout <- paste("G", WRB.tbl$Code[which(WRB.tbl$Group == k)], "HWS3a.tif", sep="")
  ## 2.5 km:
  tifoutr <- paste("G", WRB.tbl$Code[which(WRB.tbl$Group == k)], "HWS2a.tif", sep="")
  if(!file.exists(tifoutr)){
    system(paste(gdalwarp, ' ', tifout, ' ', tifoutr, ' -dstnodata 255 -r bilinear -te -180 -90 180 90 -tr ', 1/40,' ', 1/40, sep=""))
  }
  ## 5 km:
  tifoutr <- paste("G", WRB.tbl$Code[which(WRB.tbl$Group == k)], "HWS1a.tif", sep="")
  if(!file.exists(tifoutr)){
    system(paste(gdalwarp, ' ', tifout, ' ', tifoutr, ' -dstnodata 255 -r bilinear -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
  }
  ## 20 km:
  tifoutr <- paste("G", WRB.tbl$Code[which(WRB.tbl$Group == k)], "HWS0a.tif", sep="")  
  if(!file.exists(tifoutr)){
    system(paste(gdalwarp, ' ', tifout, ' ', tifoutr, ' -dstnodata 255 -r bilinear -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))
  }
}

## compress:
for(k in levels(WRB$Group)){  
  tifoutr <- paste("G", WRB.tbl$Code[which(WRB.tbl$Group == k)], "HWS", sep="")
  for(i in 0:3){
    outname = paste(tifoutr, i, 'a', sep="")
    system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), set.file.extension(outname, ".tif")))
    system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(outname, ".tif.gz"))
  }
}

HWSDRaster <- readGDAL("hwsd.tif")
# mask the water bodies (the NA values are actually "65534"):
HWSDRaster$band1 <- ifelse(HWSDRaster$band1==0, NA, HWSDRaster$band1)
# summary(HWSDRaster)
names(HWSDRaster) <- "MU_GLOBAL"
HWSDRaster <- as(HWSDRaster, "SpatialPixelsDataFrame")
gc()

# convert to table data:
HWSDRaster.pnt <- data.frame(HWSDRaster)
gc()
HWSDRaster.pnt$MU_GLOBAL <- as.factor(HWSDRaster.pnt$MU_GLOBAL)
HWSDRaster.SMU <- merge(x=HWSD.SMU, y=HWSDRaster.pnt, all.y=TRUE, all.x=FALSE, by="MU_GLOBAL")  # takes time and RAM!
# HWSDRaster.pntf <- subset(HWSDRaster.pnt, !is.na(HWSDRaster.pnt$x)&!is.na(HWSDRaster.pnt$y))
coordinates(HWSDRaster.SMU) <- ~x+y
gridded(HWSDRaster.SMU) <- TRUE
# fullgrid(HWSDRaster.SMU) <- TRUE
HWSDRaster.SMU$SU <- as.integer(HWSDRaster.SMU$SU_SYMBOL)
HWSDRaster.SMU$MU <- as.integer(paste(HWSDRaster.SMU$MU_GLOBAL))
writeGDAL(HWSDRaster.SMU["SU"], "HWSD_SMU.tif", "GTiff")

## aggregate per MU (a weighted average):
## HIR: It would be also usefull to derive the dominant class?
HWSD_M <- data.frame(MU_GLOBAL=levels(HWSD_DATA$MU_GLOBAL))
for(j in 2:length(sv.list)){
   if(is.numeric(HWSD_DATA[,sv.list[j]])){
   sm <- aggregate(. ~ MU_GLOBAL, data=data.frame(x = HWSD_DATA$SHARE/100 * HWSD_DATA[,sv.list[j]], MU_GLOBAL=HWSD_DATA$MU_GLOBAL), sum, na.action=na.omit)
   sb <- aggregate(. ~ MU_GLOBAL, data=data.frame(x = HWSD_DATA$SHARE/100 * as.numeric(as.logical(HWSD_DATA[,sv.list[j]])), MU_GLOBAL=HWSD_DATA$MU_GLOBAL), sum, na.action=na.omit)  ## TH: this one should always be 100% but I did not check;
   smb <- join(x=data.frame(MU_GLOBAL=HWSD_M$MU_GLOBAL, i=row.names(HWSD_M)), y=sm, by="MU_GLOBAL")[,"x"] / join(x=data.frame(MU_GLOBAL=HWSD_M$MU_GLOBAL, i=row.names(HWSD_M)), y=sb, by="MU_GLOBAL")[,"x"]
  # somewhere the product is null and needs to replaced (devision by 0 leads to NaN)
  HWSD_M[,sv.list[j]] <- ifelse(is.nan(smb), 0, smb)
  } else {
    if(is.factor(HWSD_DATA[,sv.list[j]])){
      sf <- HWSD_DATA[which(HWSD_DATA$x==1),c("MU_GLOBAL", sv.list[j])]
      HWSD_M[,sv.list[j]] <- join(x=HWSD_M, y=sf, by="MU_GLOBAL")[,sv.list[j]]
  }
  }
}  ## takes few minutes!
# round up the ordinary vars again:
HWSD_M$SWR <- as.integer(round(HWSD_M$SWR, 0))
HWSD_M$ROOTS <- as.integer(round(HWSD_M$ROOTS, 0))
HWSD_M$IL <- as.integer(round(HWSD_M$IL, 0))
HWSD_M$AWC_CLASS <- as.integer(round(HWSD_M$AWC_CLASS, 0))

## merge MU_GLOBAL and soil property of interest
## TH: This takes increadible 8GB of RAM!
gc()
HWSDRaster.DATA <- merge(x=HWSD_M[,sv.list], y=HWSDRaster.pnt, all.y=TRUE, all.x=FALSE, by="MU_GLOBAL")
gc()
# convert to a spatial layer:
coordinates(HWSDRaster.DATA) <- ~x+y
gridded(HWSDRaster.DATA) <- TRUE
proj4string(HWSDRaster.DATA) <- CRS("+proj=longlat +ellps=WGS84")
# str(HWSDRaster.DATA@data)
gc()

# write all layers to geotif:
for(j in 2:length(sv.list)){
  if(is.na(file.info(paste("HWSD_", sv.list[j],".tif", sep=""))$size)){
  if(is.numeric(HWSDRaster.DATA@data[,j])){
  # check the numeric precision:
  rn <- range(HWSDRaster.DATA@data[,j], na.rm=TRUE)
  cn <- length(levels(as.factor(HWSDRaster.DATA@data[,j])))
  # write to geotif:
  if(rn[1]>=0 & rn[2]<=255 & cn<255){
  writeGDAL(HWSDRaster.DATA[j], paste("HWSD_", sv.list[j],".tif", sep=""), type="Byte", mvFlag=255)
  }
  else { if(class(HWSDRaster.DATA@data[,j])=="integer") {
  writeGDAL(HWSDRaster.DATA[j], paste("HWSD_", sv.list[j],".tif", sep=""), type="Int16", mvFlag="-99999")
  }
  else {
  writeGDAL(HWSDRaster.DATA[j], paste("HWSD_", sv.list[j],".tif", sep=""), mvFlag="-99999")
  } }
  } else {
    if(is.factor(HWSDRaster.DATA@data[,j])){
     HWSDRaster.DATA$img <- as.integer(HWSDRaster.DATA@data[,j])
     gc()
     writeGDAL(HWSDRaster.DATA["img"], paste("HWSD_", sv.list[j],".tif", sep=""), type="Byte", mvFlag=255)
    }
  }
  }
}
# This will export all soil parameters as separate GeoTiffs;

## Resample to 5 km grid:
outname.lst <- paste(c("MUG", "TTX", "DRA", "RDP", "AWC", "PS1", "PS2", "ROT", "ILA", "SWR", "APR", "TGR", "TSN", "TSL", "TCL", "TUX", "TBR", "TBD", "TOC", "TPH", "TCC", "TCS", "TBS", "TTE", "TCA", "TCS", "TES", "TEC", "SGR", "SSN", "SSL", "SCL", "SUX", "SBR", "SBD", "SOC", "SPH", "SCC", "SCS", "SBS", "STE", "SCA", "SCS", "SES", "SEC"), "HWS", sep="")
# View(data.frame(x=sv.list, y=outname.lst))
for(j in 2:length(sv.list)){
   # needs to get the missing value manually!!
   inf <- GDALinfo(paste('HWSD_', sv.list[j], '.tif', sep=""))
   if(is.na(file.info(paste(outname.lst[j], '1a.tif', sep=""))$size)){ 
   # 5 km:
   system(paste(gdalwarp, ' HWSD_', sv.list[j], '.tif -t_srs \"+proj=longlat +ellps=WGS84\" ', outname.lst[j], '1a.tif -srcnodata ', attr(inf, "df")$NoDataValue ,' -dstnodata ', attr(inf, "df")$NoDataValue ,' -r near -te -180 -90 180 90 -tr 0.05 0.05', sep=""))
   }
   if(is.na(file.info(paste(outname.lst[j], '0a.tif', sep=""))$size)){
   # 20 km:
   if(is.factor(HWSDRaster.DATA@data[,j])|is.integer(HWSDRaster.DATA@data[,j])){
    system(paste(gdalwarp, ' HWSD_', sv.list[j], '.tif -t_srs \"+proj=longlat +ellps=WGS84\" ', outname.lst[j], '0a.tif -srcnodata ', attr(inf, "df")$NoDataValue ,' -dstnodata ', attr(inf, "df")$NoDataValue ,' -r near -te -180 -90 180 90 -tr 0.2 0.2', sep=""))
   } else {
    system(paste(gdalwarp, ' HWSD_', sv.list[j], '.tif -t_srs \"+proj=longlat +ellps=WGS84\" ', outname.lst[j], '0a.tif -srcnodata ', attr(inf, "df")$NoDataValue ,' -dstnodata ', attr(inf, "df")$NoDataValue ,' -r bilinear -te -180 -90 180 90 -tr 0.2 0.2', sep=""))
   }
   }
   # unlink(paste('HWSD_', sv.list[j], '.tif', sep=""))
}
GDALinfo("SWRHWS1a.tif") # must be Byte with 255 mask
GDALinfo("TPHHWS1a.tif") # must be Float32 with -99999 mask


## Compress and copy:
for(outname in c(paste(outname.lst[-1], "1a.tif", sep=""), paste(outname.lst[-1], "0a.tif", sep=""))){
  if(is.na(file.info(paste(shortPathName(normalizePath(outdir)), paste(outname, "gz", sep="."), sep="\\"))$size)){
  system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), outname))
  system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
  unlink(set.file.extension(outname, ".tif.gz"))
}  # Compression takes > 15 mins
}




## write attribute table and metadata:
write.table(HWSD_DATA[,-1], "hwsdmu.csv", row.names=FALSE, sep=";", quote=FALSE)
write.table(HWSD_DATA.mt[,-c(1:3)], "hwsd_.csv", row.names=FALSE, sep=";", quote=FALSE) # [,c("COLUMN_NAME","REMARKS")]
# description tables for each soil property:
for(j in 1:length(mt.list)){
  write.table(HWSD.mt_list[[j]][,c("CODE","VALUE")], paste("HWSD_", strsplit(mt.list[j], "D_")[[1]][2], ".csv", sep=""), row.names=FALSE, sep=";", quote=FALSE)
}

### Some numeric variables need to be manually formated:
selc <- grep(mt.list, pattern="D_ROOTS")
x <- HWSD.mt_list[[selc]]$VALUE
HWSD.mt_list[[selc]]$VALUE <- ifelse(x==">80", 90, ifelse(x=="0-20", 10, ifelse(x=="0-80", 40, ifelse(x=="20-40", 30, ifelse(x=="40-60", 50, ifelse(x=="60-80", 70, ifelse(x=="None", 0, NA)))))))
selc <- grep(mt.list, pattern="D_IL")
x <- HWSD.mt_list[[selc]]$VALUE
HWSD.mt_list[[selc]]$VALUE <- ifelse(x=="< 40", 20, ifelse(x=="> 150", 180, ifelse(x=="40-80", 60, ifelse(x=="80-150", 115, ifelse(x=="None", 0, NA)))))

# convert to actual values:
for(k in 1:length(mt.list)){
    selc <- grep(names(HWSD@data), pattern=strsplit(mt.list[k], "^D_")[[1]][2])
    if(length(nzchar(selc))>0){
      for(p in 1:length(selc)){
      spr <- data.frame(CODE=HWSD@data[,selc[p]])
      sval <- merge(x=spr, y=HWSD.mt_list[[k]], by="CODE", all.x=TRUE, all.y=FALSE)
      ## This takes time!
      HWSD@data[,selc[p]] <- sval$VALUE
}}}

# attach complete names:     
for(j in 1:length(names(HWSD@data))){
    # add complete var name:
    selr <- which(HWSD_DATA.mt$COLUMN_NAME %in% names(HWSD@data)[j])
    attr(HWSD@data[,j], "remarks") <- HWSD_DATA.mt[selr,"REMARKS"]
    # this takes time!     
}
str(HWSD@data)

HWSD.url <- "http://www.iiasa.ac.at/Research/LUC/External-World-soil-database/HWSD_Documentation.pdf"
### Statistics for soil properties:
nv <- length(names(HWSD))-1
HWSD_attr <- data.frame(varname=names(HWSD)[-1], fullname=rep(NA, nv), description=rep(NA, nv), href=rep(HWSD.url, nv), units=rep(NA, nv), minval=rep(0, nv), mminval=rep(NA, nv), avgval=rep(NA, nv), sdval=rep(NA, nv), mmaxval=rep(NA, nv), priority=rep(4, nv), maxval=rep(100, nv))
for(j in 1:nv){
     selr <- which(HWSD_META$FIELD %in% names(HWSD@data)[j])
     HWSD_attr$units[j] <- paste(HWSD_META[selr,"UNIT"])
     HWSD_attr$fullname[j] <- attr(HWSD@data[,j+1], "remarks")
     if(class(HWSD@data[,j+1])=="factor"){
         HWSD_attr$avgval[j] <- names(summary(HWSD@data[,j+1])[1]) 
         HWSD_attr$units[j] <- paste(1:length(levels(HWSD@data[,j+1])), ":", levels(HWSD@data[,j+1]), collapse="; ")
         }
     else { if(class(HWSD@data[,j+1])=="numeric"|class(HWSD@data[,j+1])=="integer"|class(HWSD@data[,j+1])=="atomic") {
         xr <- signif(range(HWSD@data[,j+1], na.rm=TRUE), 3)
         xm <- signif(mean(HWSD@data[,j+1], na.rm=TRUE), 3)
         xs <- signif(sd(HWSD@data[,j+1], na.rm=TRUE), 3)
         HWSD_attr$mminval[j] <- xr[1]
         HWSD_attr$mmaxval[j] <- xr[2]         
         HWSD_attr$avgval[j] <- xm
         HWSD_attr$sdval[j] <- xs                  
     }}
}

## this can now be written to a table:
str(HWSD_attr)
write.csv(HWSD_attr, "HWSD_attr.csv")

# Save the R image:
save(HWSD, file="HWSD.RData")
# 4,241,944 pixels
# http://globalsoilmap.net/data/HWSD.RData





# end of script;