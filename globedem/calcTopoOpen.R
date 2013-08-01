# title         : calcTopoOpen.R
# purpose       : calculate SAGA Topographic openness derived using the DEMSRE3 (for land areas only);
# reference     : [http://worldgrids.org/doku.php?id=wiki:twisre3]
# producer      : Prepared by Milan Kilibarda
# version       : 1
# inputs        : Global Relief Model based on SRTM 30+ and ETOPO DEM at 1/120 arcdeegre [http://worldgrids.org/maps/DEMSRE3a.tif.gz];
# outputs       : geotiff images projected in the "+proj=longlat +datum=WGS84" system;
# remarks 1     : First download and install FWtools [http://fwtools.maptools.org];
# remarks 2     : The calculation might miss some small islands;
# remarks 3     : SAGA Version: 2.1.0 was used;

library(RSAGA) 
library(rgdal)
library(raster)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
# download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
# unzip("7za920.zip")
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"

# clip grid with polygon rsaga.get.usage("shapes_grid",7) 2.0.8
pols<- dir(path="E:/TopoOpenes/tilesCont", pattern=glob2rx("*shp"), full.names=T)  
grds=gsub("*.shp","_DEMSRE3a.sgrd", pols)
saga210=rsaga.env(path="C:/saga_vc")
for(i in 1:length(pols)  ) {
rsaga.geoprocessor(lib="shapes_grid", module=7, 
                   param=list(OUTPUT=grds[i] ,
                              INPUT="DEMSRE3a.sgrd",
                              POLYGONS=pols[i]  ) )  }



saga208path="C:/Program Files/R/R-3.0.0/library/RSAGA/saga_vc"
# saga 210  rsaga.env(path="C:/saga_vc")
rsaga.get.usage("ta_lighting",5, env=saga210)
OPN_p=gsub("*_DEMSRE3a.sgrd","_OPN_p.sgrd", grds)
OPN_n=gsub("*_DEMSRE3a.sgrd","_OPN_n.sgrd", grds)


# for(i in 1:length(grds)  ) {
#   rsaga.geoprocessor(lib="ta_lighting", module=5, env=saga210,
#                      param=list(DEM=grds[i] ,
#                                 POS=OPN_p[i],
#                                 NEG=pols[i]  ) )  }

OPN_p=pols<- dir(path="E:/TopoOpenes/tilesCont/TOPN_pos", pattern=glob2rx("*sgrd"), full.names=T)
## merge maps: 
TOP.sgrd.list <- paste(OPN_p, sep="", collapse=";")
rsaga.geoprocessor(lib="grid_tools", module=3, param=list(GRIDS=TOP.sgrd.list, MERGED="OPISRE3a.sgrd", TYPE=7, INTERPOL=1, OVERLAP=0, MERGE_INFO_MESH_SIZE=1/120))  

# assign projection to merged grid
rsaga.geoprocessor(lib="pj_proj4", module=0, 
                   param=list(GRIDS="OPISRE3a.sgrd",
                              CRS_PROJ4="+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6378137.000000 +b=6356752.314245 +no_defs"))
# Transfrom to WGS84             

rsaga.geoprocessor(lib="pj_proj4", module=7, 
                   param=list(SOURCE="OPISRE3a.sgrd",
                              SOURCE_PROJ="+proj=sinu +lon_0=0 +x_0=0 +y_0=0 +a=6378137.000000 +b=6356752.314245 +no_defs",
                              TARGET_PROJ="+proj=longlat +datum=WGS84",
                              INTERPOLATION='4',
                              TARGET_TYPE='0',
                              GET_USER_XMIN= '-180',
                              GET_USER_XMAX='180',
                              GET_USER_YMIN= '-90',
                              GET_USER_YMAX='90',
                              GET_USER_SIZE='0.008333',
                              GET_USER_GRID="OPISRE3a.sgrd"  ))

# reduce the storage size
rsaga.geoprocessor(lib="grid_calculus", 1, param=list(GRIDS='OPISRE3a.sgrd',          
                                                      RESULT='OPISRE3a.sgrd',
                                                      FORMULA='int(a*1000)'))   # FORMULA

# Export Raster to GeoTIFF                                             
rsaga.geoprocessor(lib="io_gdal", 2, param=list(GRIDS='OPISRE3a.sgrd',          
                                                FILE='OPISRE3a_tmp.tif'))


# convert to geotifs (1 km):
system(paste(gdalwarp, ' OPISRE3a_tmp.tif -t_srs \"+proj=longlat +datum=WGS84\" -ot UInt16 -srcnodata \"-99999\" -dstnodata \"0\" OPISRE3a.tif -r near -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""))
# 2.5 km:
system(paste(gdalwarp, ' OPISRE3a_tmp.tif -t_srs \"+proj=longlat +datum=WGS84\" -ot UInt16 -srcnodata \"-99999\" -dstnodata \"0\" OPISRE2a.tif -r bilinear -te -180 -90 180 90 -tr ', 1/40,' ', 1/40, sep=""))
# 5 km:
system(paste(gdalwarp, ' OPISRE3a_tmp.tif -t_srs \"+proj=longlat +datum=WGS84\" -ot UInt16 -srcnodata \"-99999\" -dstnodata \"0\" OPISRE1a.tif -r bilinear -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
# 20 km:
system(paste(gdalwarp, ' OPISRE1a.tif -t_srs \"+proj=longlat +datum=WGS84\" -ot UInt16 -dstnodata \"-99999\" -dstnodata \"0\" OPISRE0a.tif -r bilinear -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))

## compress:  
for(i in 0:3){
  outname = paste("OPISRE", i, 'a', sep="")
  system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), set.file.extension(outname, ".tif")))
  system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
  unlink(set.file.extension(outname, ".tif.gz"))
}

