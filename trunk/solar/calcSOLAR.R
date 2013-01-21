# title         : calcSOLAR.R
# purpose       : extraction of the potential solar insolation globally (mean value and standard deviation);
# reference     : [https://code.google.com/p/worldgrids/source/browse/] 
# reference     : [http://worldgrids.org/doku.php?id=wiki:inmsre3"]
# producer      : Prepared by M. Kilibarda and T. Hengl
# version       : 1
# inputs        : solar radiation derived for various days in SAGA GIS from globe DEM [http://worldgrids.org/doku.php?id=wiki:demsre3];
# outputs       : geotiff images projected in the "+proj=longlat +datum=WGS84" system;
# remarks 1     : First download and install FWtools [http://fwtools.maptools.org]; This script takes very long time, depends on pc, but more then 15 days.

# -------------------------------------------
# Initial settings and data download:
# -------------------------------------------

library(sp)
library(RSAGA) 
library(rgdal)
library(GSIF)
library(maptools)
# split strings in names and values parts
require(stringr)
# Multicore processing
library(doMC)

fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdalbuildvrt = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalbuildvrt.exe"))))
gdal_merge.py = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_merge.py"))))

download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"


registerDoMC(7) # 7 cores processing in parallel
getDoParWorkers()

# list of dates for calculation 

list.dates<-as.Date(c("2011-01-01","2011-01-09","2011-01-17","2011-01-25","2011-02-02","2011-02-10","2011-02-18","2011-02-26","2011-03-06",
                      "2011-03-14","2011-03-22","2011-03-30","2011-04-07","2011-04-15","2011-04-23","2011-05-01","2011-05-09","2011-05-17",
                      "2011-05-25","2011-06-02","2011-06-10","2011-06-18","2011-06-26","2011-07-04","2011-07-12","2011-07-20","2011-07-28",
                      "2011-08-05","2011-08-13","2011-08-21","2011-08-29","2011-09-06","2011-09-14","2011-09-22","2011-09-30","2011-10-08",
                      "2011-10-16","2011-10-24","2011-11-01","2011-11-09","2011-11-17","2011-11-25","2011-12-03","2011-12-11","2011-12-19","2011-12-27"))

# the first result contains days from 2010-12-28 to 2011-01-04 ("2010-12-28" "2010-12-29" "2010-12-30" "2010-12-31" "2011-01-01" "2011-01-02" "2011-01-03" "2011-01-04" )
list.dates<-data.frame(modis_dates=list.dates,date_a=list.dates-4,date_b=list.dates+3)
# preparation for SAGA 
list.dates$day_a= as.character( as.numeric(format(list.dates$date_a, "%d") ) -1 )
list.dates$day_b= as.character(as.numeric(format(list.dates$date_b, "%d") ) -1)
list.dates$mon_a= as.character(as.numeric(format(list.dates$date_a, "%m") ) -1)
list.dates$mon_b= as.character(as.numeric(format(list.dates$date_b, "%m") ) -1 )

#########################  Function for grid tiling with overlap ###################################################

    gridtiling<-function(file_path="demsre3.sgrd", # path to grid file in SAGA format 
                         tilename="tile", # prexif to be given to tile names
                         tile_size=1000, # in cells 
                         overlapping=50, # in cells
                         tiles_folder=paste(getwd(),'tiles',sep='/'), # resulting folder
                         crs=CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6378137.000000 +b=6356752.314245 +no_defs")  ) {
      
      require(sp)  
      require(rgdal)
      require(RSAGA)
      require(stringr)
      
      hd <- readLines(file_path, n=-1) # read grid header 
      dir.create(tiles_folder)
      
      # split strings in names and values parts
      require(stringr)
      hd <- unlist(str_split(hd, "\t= "))
      
      # make a dataframe from the header data
      header <- data.frame(Parameter=hd[seq(1,length(hd),2)],
                           Value= hd[seq(2,length(hd),2)]  ,  stringsAsFactors=FALSE)
      
      cel=as.numeric(header[header$Parameter== 'CELLSIZE',2] ) # velicina celije
      
      # rsaga.get.usage("shapes_grid",10) 
      # Grid System Extent 
      rsaga.geoprocessor(lib="shapes_grid", module=10, param=list(PARAMETERS_GRID_SYSTEM_NX=as.numeric(header[header$Parameter== 'CELLCOUNT_X',2] ),
                                                                  PARAMETERS_GRID_SYSTEM_NY=as.numeric(header[header$Parameter== 'CELLCOUNT_Y',2] ), 
                                                                  PARAMETERS_GRID_SYSTEM_X=as.numeric(header[header$Parameter== 'POSITION_XMIN',2] ), 
                                                                  PARAMETERS_GRID_SYSTEM_Y=as.numeric(header[header$Parameter== 'POSITION_YMIN',2] ),
                                                                  PARAMETERS_GRID_SYSTEM_D=as.numeric(header[header$Parameter== 'CELLSIZE',2] ),
                                                                  SHAPES='extent.shp' ,CELLS=0))
      
      extent= readOGR(".", "extent")
      
      bb=t(extent@bbox)
      
      bb[1,]=bb[1,]-overlapping*cel
      bb[2,]=bb[2,]+overlapping*cel
      
      step=tile_size*cel+overlapping*cel
      nlon=ceiling(diff(bb[,1]) / step )
      nlat=ceiling(diff(bb[,2]) / step )
      
      
      p.l <- expand.grid(KEEP.OUT.ATTRS=FALSE, lonmin=seq(bb[1,1],bb[1,1]+(nlon-1)*step,   by=step), latmin=seq(bb[1,2],bb[1,2]+(nlat-1)*step  ,by=step)) - overlapping*cel    
      p.u <- expand.grid(KEEP.OUT.ATTRS=FALSE, lonmax=seq(bb[1,1]+step, bb[1,1]+nlon*step   ,by=step), latmax=seq(bb[1,2]+step   ,bb[1,2]+nlat*step,  by=step))  +overlapping*cel  
      ptiles <- cbind(p.l, p.u)
      
      poligoni=as.list(rep(NA,length(ptiles$lonmin)))
      
      for(i in 1:length(ptiles$lonmin)) {
        
        x <- rbind(ptiles$lonmin[i], ptiles$lonmax[i],ptiles$lonmax[i],ptiles$lonmin[i],ptiles$lonmin[i])
        y <- rbind( ptiles$latmin[i], ptiles$latmin[i], ptiles$latmax[i], ptiles$latmax[i],ptiles$latmin[i])
        
        poligoni[[i]]<-  Polygons( list(Polygon(cbind(x,y))) ,i)  
      }
      
      poll<-SpatialPolygons(poligoni)
      
      
      tile<-SpatialPolygonsDataFrame(poll, data.frame(tile=paste('tile',1:length(ptiles$lonmin),sep="")), match.ID = F)
      
      proj4string(tile)<- crs
      
      wd=getwd()
      proj4string(tile)<- crs
      writeOGR(tile, dsn=getwd(), layer='tile', driver='ESRI Shapefile', overwrite_layer=TRUE)
      
      # rsaga.get.usage("shapes_grid",2)
      # Grid Statistics for Polygons 
      rsaga.geoprocessor(lib="shapes_grid", module=2, param=list(GRIDS=paste(getwd(),'/',file_path,sep=''), 
                                                                 POLYGONS="tile.shp", RESULT="tile_stat.shp",MAX=TRUE,VAR=FALSE, STDDEV=FALSE,QUANTILE=FALSE))
      tile_stat <- readOGR(".", "tile_stat")
      
      # summary(tile_stat@data)
      # Remove polygons (tiles) with nodata 
      tile_stat<-tile_stat[!is.na(tile_stat@data[,2]),]
      tile_stat$TILE<-1:length(tile_stat$TILE)
      writeOGR(tile_stat, dsn=getwd(), layer='tile_o', driver='ESRI Shapefile', overwrite_layer=TRUE)
      
      
      for(i in 1:length(tile_stat@data[,2]) ) {
        
        tile1<-tile_stat[i,]  
        # tile1 <- spTransform(tile1, crs)  #bbox is ok now
        writeOGR(tile1, dsn=getwd(), layer='tile1', driver='ESRI Shapefile', overwrite_layer=TRUE)
        
        # rsaga.get.usage("shapes_grid",2)
        # Grid Statistics for Polygons
        rsaga.geoprocessor(lib="shapes_grid", module=7, param=list(INPUT=paste(getwd(),'/',file_path,sep=''), POLYGONS="tile1.shp",
                                                                   OUTPUT=paste(tiles_folder,"/",tilename,i,".sgrd",sep="")  ))
        # fix nodata values coused by croping grid of integer type
        aa=readGDAL(paste(tiles_folder,"/",tilename,i,'.sdat',sep=""))
        
        aa$band1=ifelse(aa$band1== as.numeric(header[header$Parameter== 'NODATA_VALUE',2] ),NA,aa$band1)
        aa$band1=ifelse(aa$band1==-32767 ,NA,aa$band1)
        
        writeGDAL(aa,paste(tiles_folder,"/",tilename,i,'.sdat',sep=""), "SAGA",mvFlag=-99999)
        
        # crop to data
        rsaga.geoprocessor(lib="grid_tools", 17,param=list(INPUT=paste(tiles_folder,"/",tilename,i,".sgrd",sep=""),  OUTPUT=paste(tiles_folder,"/",tilename,i,".sgrd",sep="") )  )
      }
    }# end of function

#########################  End of tiling function ###########################################################################

# setwd("/home/kili/Solar")



gridtiling(file_path="demsre3.sgrd", # path to grid file in SAGA format, in this case DEM in sinusoidal projection, resolution 1 km 
           tilename="dem", # prexif to be given to tile names
           tile_size=1000, # in cells 
           overlapping=50, # in cells
           tiles_folder=paste(getwd(),'demsre3_tiles',sep='/'), # resulting folder
           crs=CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6378137.000000 +b=6356752.314245 +no_defs"))

gridtiling(file_path="lat.sgrd", # path to grid file in SAGA format, in this case grid of latitudes in sinusoidal projection, resolution 1 km
           tilename="lat", # prexif to be given to tile names
           tile_size=1000, # in cells 
           overlapping=50, # in cells
           tiles_folder=paste(getwd(),'lat_tiles',sep='/'), # resulting folder
           crs=CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6378137.000000 +b=6356752.314245 +no_defs"))
gridtiling(file_path="lon.sgrd", # path to grid file in SAGA format, in this case grid of longitudes in sinusoidal projection, resolution 1 km
           tilename="lon", # prexif to be given to tile names
           tile_size=1000, # in cells 
           overlapping=50, # in cells
           tiles_folder=paste(getwd(),'lon_tiles',sep='/'), # resulting folder
           crs=CRS("+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6378137.000000 +b=6356752.314245 +no_defs"))

demsre3o<-dir(path=paste(getwd(),"demsre3_tiles",sep='/'), pattern=glob2rx("*.sgrd"), full.names=T)
lato<-dir(path=paste(getwd(),"lat_tiles",sep='/'), pattern=glob2rx("*.sgrd"), full.names=T)
lono<-dir(path=paste(getwd(),"lon_tiles",sep='/'), pattern=glob2rx("*.sgrd"), full.names=T)

tot_nam<-gsub(".sgrd","tot.sgrd", dir(paste(path=getwd(),"/demsre3_tiles",sep=""), pattern=glob2rx("*.sgrd"), full.names=F) )

# parallel processing
foreach(j=1:46) %dopar%{
  
  br_mozaika=j
  datum=gsub("-","_", as.character(list.dates$modis_dates[br_mozaika]) )
  
  
  tot_fol=paste(getwd(),"/total",datum,sep="")
  dif_fol=paste(getwd(),"/diffuse",datum,sep="")
  dir_fol=paste(getwd(),"/direct",datum,sep="")
  
  dir.create(tot_fol)
  dir.create(dif_fol)
  dir.create(dir_fol)
  
  for(i in 1:length(demsre3o))  {
    
    
    
    rsaga.geoprocessor(lib="ta_lighting", 2, param=list(GRD_DEM=paste(getwd(),"/demsre3_tiles/dem",i,'.sgrd',sep=""),
                                                        GRD_LAT=paste(getwd(),"/lat_tiles/lat",i,'.sgrd',sep=""),
                                                        GRD_LON=paste(getwd(),"/lon_tiles/lon",i,'.sgrd',sep=""),
                                                        GRD_DIRECT=paste(dir_fol,'dir.sgrd',sep="/"),
                                                        GRD_DIFFUS=paste(dif_fol,'dif.sgrd',sep="/"),
                                                        GRD_TOTAL=paste(tot_fol,'/tot',i,'.sgrd',sep=""), 
                                                        PERIOD="2", DHOUR=6, 
                                                        DDAYS=1, 
                                                        DAY_A=list.dates$day_a[br_mozaika],
                                                        MON_A=list.dates$mon_a[br_mozaika],
                                                        DAY_B=list.dates$day_b[br_mozaika],
                                                        MON_B=list.dates$mon_b[br_mozaika]) ,intern=F )
    
  } # end i
  
} # end j (foreach)


ime=c('I01' , 'I02'  ,'I03',  'I04',  'I05','I06'  ,'I07',	'I08',	'I09',	'I10',paste("I",11:46,sep=""))
ime=paste(ime,'SRT3a.tif',sep='')
wd=getwd()           

for(j in 1:46){
  
  br_mozaika=j
  datum=gsub("-","_", as.character(list.dates$modis_dates[br_mozaika]) )
  
  setwd(wd)
  
  tot_fol=paste(getwd(),"/total",datum,sep="")
  dif_fol=paste(getwd(),"/diffuse",datum,sep="")
  dir_fol=paste(getwd(),"/direct",datum,sep="")
  
  setwd(tot_fol)
  fajlovi<-dir(path=tot_fol, pattern=glob2rx("*.sdat"), full.names=F) #  calc. radiation all tiles 
  merged=paste(" ",wd,"/inss",datum,".sdat",sep="") #  1. mosaic
  
  ff =paste("\"",fajlovi,"\"" ,sep="", collapse=" ") 
  system(paste(gdal_merge.py, "  -of saga -ot Byte -o ",merged," -n -99999 -a_nodata -99999 ",ff, sep=""))  # merge all files
  
  # rsaga.get.usage("pj_proj4",0)
  # Set Coordinate Reference System              
  rsaga.geoprocessor(lib="pj_proj4", module=0, 
                     param=list(GRIDS=paste(wd,"/inss",datum,".sgrd",sep=""),
                                CRS_PROJ4="+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6378137.000000 +b=6356752.314245 +no_defs"))
  
  # rsaga.get.usage("pj_proj4",7)
  # Transfrom to WGS84             
  
  rsaga.geoprocessor(lib="pj_proj4", module=7, 
                     param=list(SOURCE=paste(wd,"/inss",datum,".sgrd",sep=""),
                                SOURCE_PROJ="+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6378137.000000 +b=6356752.314245 +no_defs",
                                TARGET_PROJ="+proj=longlat +datum=WGS84",
                                INTERPOLATION='0',
                                TARGET_TYPE='0',
                                GET_USER_XMIN= '-180',
                                GET_USER_XMAX='180',
                                GET_USER_YMIN= '-90',
                                GET_USER_YMAX='90',
                                GET_USER_SIZE='0.008333',
                                GET_USER_GRID=paste(wd,"/inss",datum,".sgrd",sep="")
                     ))
  
  # rsaga.get.usage("io_gdal",2)   
  # Export Raster to GeoTIFF 
  rsaga.geoprocessor(lib="io_gdal", 2, param=list(GRIDS=paste(wd,"/inss",datum,".sgrd",sep=""), 
                                                  FILE=paste(wd,ime[j],sep="/") ))
  
  
} # end of j # end of solar radiation creation, 46 8days tiff images are calculated

## list of images:
dtif.lst <- list.files(pattern=glob2rx("*SRT3a.tif$"))
inf <- GDALinfo(dtif.lst[1])

## tiles:
b.l <- expand.grid(KEEP.OUT.ATTRS=FALSE, lonmin=seq(-180,0,by=180), latmin=seq(-90,0,by=90))
b.u <- expand.grid(KEEP.OUT.ATTRS=FALSE, lonmax=seq(0,180,by=180), latmax=seq(0,90,by=90))
btiles <- cbind(b.l, b.u)
str(btiles)

for(j in 1:nrow(btiles)){
  for(i in 1:length(dtif.lst)){
    outname <- strsplit(dtif.lst[i], ".tif")[[1]][1]
    tilename <- paste(outname, '_', j, '.sdat', sep="")
    # resample
    if(is.na(file.info(tilename)$size)){
      system(paste(gdalwarp, ' ', dtif.lst[i], ' ', tilename, ' -of \"SAGA\" -dstnodata \"255\" -r near -te ', btiles[j,"lonmin"] , ' ', btiles[j,"latmin"], ' ', btiles[j,"lonmax"] ,' ', btiles[j,"latmax"] ,' -tr ', 1/120, ' ', 1/120, sep=""))
    }
  }
  # derive derive mean, sd, min and max values:
  t.lst <- list.files(pattern=glob2rx(paste('*3a_', j, '.sgrd$', sep="")))
  if(is.na(file.info(paste("INMSRE3a_", j, ".sgrd", sep=""))$size)){
    rsaga.geoprocessor(lib="geostatistics_grid", module=4, param=list(GRIDS=paste(t.lst, collapse=";"), MEAN=paste("INMSRE3a_", j, ".sgrd", sep=""), MIN=paste("INLSRE3a_", j, ".sgrd", sep=""), MAX=paste("INHSRE3a_", j, ".sgrd", sep=""), STDDEV=paste("INSSRE3a_", j, ".sgrd", sep="")))
  }
}

# create a mosaic and compress files:
for(j in c("INMSRE", "INSSRE", "INLSRE", "INHSRE")){
  unlink("insol.vrt")
  system(paste(gdalbuildvrt, "insol.vrt", paste(j, "3a_", 1:nrow(btiles), ".sdat", sep="", collapse=" ")))
  # convert to geotif (1 km):
  if(is.na(file.info(paste(j, '3a.tif', sep=""))$size)){
  if(j == "INSSRE"){
    system(paste(gdalwarp, ' insol.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-99999\" ', j, '3a.tif', ' -r near -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""))
  } else {
    system(paste(gdalwarp, ' insol.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"255\" -ot Byte ', j, '3a.tif', ' -r near -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""))  
  }
  }
  # 2.5 km:
  if(is.na(file.info(paste(j, '2a.tif', sep=""))$size)){
  if(j == "INSSRE"){
    system(paste(gdalwarp, ' insol.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-99999\" ', j, '2a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/40,' ', 1/40, sep=""))
  } else {
  system(paste(gdalwarp, ' insol.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"255\" -ot Byte ', j, '2a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/40,' ', 1/40, sep=""))
  }
  }
  # 5 km:
  if(is.na(file.info(paste(j, '1a.tif', sep=""))$size)){
  if(j == "INSSRE"){
    system(paste(gdalwarp, ' insol.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-99999\" ', j, '1a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
  } else {
  system(paste(gdalwarp, ' insol.vrt -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"255\" -ot Byte ', j, '1a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
  }
  }
  # 20 km:
  if(is.na(file.info(paste(j, '0a.tif', sep=""))$size)){
  if(j == "INSSRE"){
    system(paste(gdalwarp, ' ', j, '1a.tif -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"-99999\" ', j, '0a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))
  } else {
  system(paste(gdalwarp, ' ', j, '1a.tif -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"255\" -ot Byte ', j, '0a.tif', ' -r bilinear -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))
  }
  }
}

## compress:  
for(j in c("INMSRE", "INSSRE", "INLSRE", "INHSRE")){
for(i in 0:3){
    outname = paste(j, i, 'a', sep="")
    system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), set.file.extension(outname, ".tif")))
    system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(outname, ".tif.gz"))
  }
}

# clean up temp files:
unlink(list.files(pattern=glob2rx("inss*.*")))
unlink(list.files(pattern=glob2rx("IN*SRE_*.*")))
unlink(list.files(pattern=glob2rx("I*SRT3a_*.*")))

# clean up temp folders:
unlink(list.dirs(pattern=glob2rx("dif*")), recursive=T)
unlink(list.dirs(pattern=glob2rx("dir*")), recursive=T)
unlink(list.dirs(pattern=glob2rx("tot*")), recursive=T)
# end of script;