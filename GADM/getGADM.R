# title         : getGADM.R
# purpose       : download and resampling of Global Administrative DM data;
# reference     : [http://worldgrids.org]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, Nov 2012.
# inputs        : data available for download from [http://www.gadm.org/world]; 
# outputs       : geotiff images projected in the "+proj=longlat +ellps=WGS84" system;
# remarks 1     : First download and install SAGA GIS [http://www.saga-gis.org] and FWtools [http://fwtools.maptools.org]; 


library(RSAGA) 
library(rgdal)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"

# download files from server:
download.file("http://www.gadm.org/data2/gadm_v2_shp.zip", destfile=paste(getwd(), "gadm_v2_shp.zip", sep="/"))
unzip(zipfile="gadm_v2_shp.zip", exdir=getwd())
ogrInfo("gadm2.shp", "gadm2")
## Get the country names (ISO standard)
gadm <- read.dbf("gadm2.dbf")  ## Large DBF!
str(gadm)
gadm$dbf[which(gadm$dbf$NAME_0=="Angola")[1],"ISO"]
cnt = aggregate(gadm$dbf[,"ID_0"], by=list(gadm$dbf$ISO), FUN=mean)
str(cnt)
ISO.country = cnt
names(ISO.country) <- c("ISO", "id")
save(ISO.country, file="ISO.country.rda")

## make a PAL file:
rgb.tbl <- data.frame(R=NA, G=NA, B=NA)
library(XML)
tmp <- xmlTreeParse("cntgad.sprm", useInternalNodes = TRUE)
cols <- xmlRoot(tmp)[[11]]
col.lst <- sapply(sapply(xmlChildren(cols), function(x) x), xmlValue)
for(j in 1:length(col.lst)){
   rgb.tbl[j,c("R","G","B")] <- as.integer(substr(strsplit(col.lst[j], " ")[[1]], start=2, stop=4))
}
str(rgb.tbl)
## write a PAL file:
write.table(rgb.tbl, "cntgad.PAL", quote = FALSE,  col.names = FALSE)

## create a SAGA txt colour table
## convert to BGR codes:
BGR <- (rgb.tbl$B * 65536) + (rgb.tbl$G * 256) + rgb.tbl$R
## write a lookup table for SAGA GIS:
filename <- file("cntgad.txt", "w", blocking=FALSE)
write("COLOR\tNAME\tDESCRIPTION\tMINIMUM\tMAXIMUM", filename)
for(i in 1:nrow(cnt)){
  write(paste(BGR[i], cnt[i,"Group.1"], paste("CL", i, sep=""), (cnt[i,"x"]-1)+0.1, (cnt[i,"x"]+1)+0.1, sep="\t"), filename, append=TRUE)
}
close(filename)

## convert to raster (use the polygon ID):
rsaga.geoprocessor(lib="grid_gridding", module=0, param=list(USER_GRID="countries.sgrd", INPUT="gadm2.shp", FIELD=1, TARGET=0, LINE_TYPE=1, USER_SIZE=1/120, USER_XMIN=-180+1/120*0.5, USER_XMAX=180-1/120*0.5, USER_YMIN=-90+1/120*0.5, USER_YMAX=90-1/120*0.5)) # takes cca 5 mins uses > 6GB!!
## export to geotiff:
# system(paste(gdal_translate, "countries.sdat CNTGAD3a.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\""))
rsaga.geoprocessor(lib="io_gdal", module=2, param=list(GRIDS="countries.sgrd", FILE="countries.tif"))
## TH: GDAL warp with SAGA sdat file receives error "ERROR 3: Unable to seek to beginning of grid row."
system(paste(gdal_translate, "countries.tif CNTGAD3a.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -a_nodata \"0\""))
# system(paste(gdalwarp, ' countries.tif CNTGAD3a.tif -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"255\" -ot Byte -r near -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""))
GDALinfo("CNTGAD3a.tif")

## Resample to other resolutions using the majority value!
# 2.5 km:
rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT="countries.sgrd", USER_GRID="CNTGAD2a.sgrd", KEEP_TYPE=FALSE, SCALE_UP_METHOD=9, USER_SIZE=1/40, USER_XMIN=-180+1/40*0.5, USER_XMAX=180-1/40*0.5, USER_YMIN=-90+1/40*0.5, USER_YMAX=90-1/40*0.5)) 
# system(paste(gdalwarp, ' countries.tif CNTGAD2a.tif -t_srs \"+proj=longlat +datum=WGS84\" -dstnodata \"0\" -ot Byte -r near -te -180 -90 180 90 -tr ', 1/40,' ', 1/40, sep=""))
system(paste(gdal_translate, "CNTGAD2a.sdat CNTGAD2a.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -a_nodata \"0\""))
## 5 km:
rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT="CNTGAD2a.sgrd", USER_GRID="CNTGAD1a.sgrd", KEEP_TYPE=FALSE, SCALE_UP_METHOD=9, USER_SIZE=1/20, USER_XMIN=-180+1/20*0.5, USER_XMAX=180-1/20*0.5, USER_YMIN=-90+1/20*0.5, USER_YMAX=90-1/20*0.5))
# system(paste(gdalwarp, ' CNTGAD2a.tif CNTGAD1a.tif -t_srs \"+proj=longlat +datum=WGS84\" -ot Byte -r near -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
system(paste(gdal_translate, "CNTGAD1a.sdat CNTGAD1a.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -a_nodata \"0\""))
## 20 km:
rsaga.geoprocessor(lib="grid_tools", module=0, param=list(INPUT="CNTGAD1a.sgrd", USER_GRID="CNTGAD0a.sgrd", KEEP_TYPE=FALSE, SCALE_UP_METHOD=9, USER_SIZE=1/5, USER_XMIN=-180+1/5*0.5, USER_XMAX=180-1/5*0.5, USER_YMIN=-90+1/5*0.5, USER_YMAX=90-1/5*0.5))
# system(paste(gdalwarp, ' CNTGAD1a.tif CNTGAD0a.tif -t_srs \"+proj=longlat +datum=WGS84\" -ot Byte -r near -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))
system(paste(gdal_translate, "CNTGAD0a.sdat CNTGAD0a.tif -a_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -a_nodata \"0\""))

## compress:  
for(i in 0:3){
    outname = paste("CNTGAD", i, 'a', sep="")
    system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), set.file.extension(outname, ".tif")))
    system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(outname, ".tif.gz"))
}

## clean up temp files:
unlink(list.files(pattern=glob2rx("countries.*")))

## Continents:
download.file("http://www.naturalearthdata.com/download/10m/cultural/ne_10m_admin_0_countries.zip", "ne_10m_admin_0_countries.zip") 
unzip("ne_10m_admin_0_countries.zip")
## import to R:
admin.wrld = readShapePoly("ne_10m_admin_0_countries.shp")
proj4string(admin.wrld) <- "+proj=latlong +datum=WGS84"
admin.wrld.l <- as(admin.wrld, "SpatialLinesDataFrame")

## continents and land mask:
download.file("http://worldgrids.org/lib/exe/fetch.php?media=lmtgsh3a.tif.gz", "lmtgsh3a.tif.gz")
system(paste("7za e lmtgsh3a.tif.gz"))

## Australia and New Zealand: 
au.csy = "+proj=aea +lat_1=-18 +lat_2=-36 +lat_0=0 +lon_0=132 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
## (Africa) proj4: 
af.csy = "+proj=laea +lat_0=5 +lon_0=20 +x_0=0 +y_0=0 +units=m +ellps=WGS84 +datum=WGS84"
## South Asia:
sas.csy = "+proj=lcc +lat_1=7 +lat_2=-32 +lat_0=-15 +lon_0=125 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
## North Asia:
nas.csy = "+proj=lcc +lat_1=15 +lat_2=65 +lat_0=30 +lon_0=95 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
## Europe:
eu.csy = "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"
## North America: 
na.csy = "+proj=aea +lat_1=29.5 +lat_2=45.5 +lat_0=23 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
## South/Central America:
sa.csy = "+proj=aea +lat_1=-5 +lat_2=-42 +lat_0=-32 +lon_0=-60 +x_0=0 +y_0=0 +ellps=aust_SA +units=m +no_defs"
## Antartica:
ant.csy = "+proj=stere +lat_0=-90 +lat_ts=-71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
## Artic:
art.csy = "+proj=stere +lat_0=90 +lat_ts=71 +lon_0=0 +k=1 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
## North Pacific ocean (Hawaii):
npo.csy = "+proj=aea +lat_1=8 +lat_2=18 +lat_0=13 +lon_0=-157 +x_0=0 +y_0=0 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
## North Atlantic ocean (Azores):
nao.csy = "+proj=utm +zone=26 +ellps=intl +towgs84=-104,167,-38,0,0,0,0 +units=m +no_defs"
## South Pacific Ocean:
poly.csy = "+proj=aea +lat_1=-20 +lat_2=-60 +lat_0=-42 +lon_0=-155 +x_0=0 +y_0=0 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
## Antartic land:
anta.csy = "+proj=utm +zone=42 +south +ellps=intl +towgs84=145,-187,103,0,0,0,0 +units=m +no_defs"

## bounding boxes:
conts = c("au", "af", "sas", "nas", "eu", "na", "sa", "ant", "art", "npo", "nao", "poly", "anta")
xmin = c(-5013000, -5090000, -10733000, -4191000, 2426000, -4792000, -4013000, -2963000, -2963000, -2404000, 89000, -2357000, -1958000)
xmax = c(5946000, 4470000, 1454000, 4096000, 7294000, 3353000, 3289000, 2983000, 2983000, 451000, 768000, 4439000, 1029000)
ymin = c(-6420000, -4616000, 961000, 316000, 1428000, -1856000, -3380000, -2708000, -2708000, -1392000, 4041000, 1044000, 3963000)
ymax = c(600000, 4004000, 9400000, 5966000, 5447000, 5753000, 5780000, 2524000, 2524000, 1950000, 4426000, 4246000, 4951000)
csy = c(au.csy, af.csy, sas.csy, nas.csy, eu.csy, na.csy, sa.csy, ant.csy, art.csy, npo.csy, nao.csy, poly.csy, anta.csy)
continents <- data.frame(conts, xmin, xmax, ymin, ymax, csy)
save(continents, file="continents.rda")

## tile and reproject land mask and DEM per continent:
for(j in 1:nrow(continents)){
  if(is.na(file.info(paste(continents[j,"conts"], '_LMTGSH3a.tif', sep=""))$size)){
    system(paste(gdalwarp, ' LMTGSH3a.tif -t_srs \"', continents[j,"csy"], '\" ', continents[j,"conts"], '_LMTGSH3a.tif -ot Byte -r near -te ', continents[j,"xmin"],' ', continents[j,"ymin"],' ', continents[j,"xmax"],' ', continents[j,"ymax"], ' -tr 1000 1000', sep=""))
  }
}

## Some pixels close to longitudes 180 / -180, and latitudes 90 / -90 get messed up and need to be fixed manually:
#land.art <- readGDAL("art_lmtgsh3a.mpr")
#names(land.art) = "mask"
#land.art$mask <- ifelse(land.art$mask==129, NA, land.art$mask)
#save(land.art, file="land.art.rda", compress="xz")
#writeGDAL(land.art, "art_LMTGSH3a.tif", type="Byte", mvFlag=129)

## The output maps are at:
#http://worldgrids.org/lib/exe/fetch.php?media=continents1km.zip

# end of script;