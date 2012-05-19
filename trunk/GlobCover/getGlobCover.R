
# title         : getGlobCover.R
# purpose       : download and resampling of GlobCover Land Cover version V2.2 image;
# reference     : [https://code.google.com/p/worldgrids/source/browse/]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, June 2012.
# inputs        : data available for download from [http://ionia1.esrin.esa.int];
# outputs       : geotiff images projected in the "+proj=longlat +datum=WGS84" system;
# remarks 1     : First download and install FWtools [http://fwtools.maptools.org]; 
# remarks 2     : The GlobCover Land Cover product is the highest resolution (300 meters) Global Land Cover product produced until 2011 - but it is not advisable to work with global resolution of 300 m (e.g. to run GIS analysis on a standard PC);

# -------------------------------------
# Initial settings and data download:
# -------------------------------------

library(rgdal)
library(RSAGA)
library(gdata)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
gdaldem = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdaldem.exe"))))
gdalbuildvrt = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalbuildvrt.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
si <- Sys.info()
outdir <- "D:/WORLDGRIDS/maps"

# download files from the server:
download.file("ftp://guestglobcover:Pir8Eefo@us-ext-nas.eo.esa.int/global/Globcover_200412_200606_V2.2_Global.zip", destfile=paste(getwd(), "Globcover.zip", sep="/"))
unzip(zipfile="Globcover.zip", exdir=getwd())
GDALinfo("GLOBCOVER_200412_200606_V2.2_Global_CLA.tif")

# --------------------------------------
# Extract maps / classes separately:
# --------------------------------------

class.lst <- read.dbf("GLOBCOVER_200412_200606_V2.2_Global_CLA.tif.vat.dbf")
perl <- gdata:::findPerl("perl")
class.names <- read.xls("Globcover_Legend.xls", perl=perl)
class.names$Label

# resample to 1 km resolution / create indicator maps:
system(paste(gdalwarp, "GLOBCOVER_200412_200606_V2.2_Global_CLA.tif -t_srs \"+proj=longlat +datum=WGS84\" GLCESA3a.tif -r near -te -180 -90 180 90 -tr", 1/120, 1/120))
GDALinfo("GLCESA3a.tif")

# split world in 4 blocks:
system(paste(gdalwarp, 'GLOBCOVER_200412_200606_V2.2_Global_CLA.tif -t_srs \"+proj=longlat +datum=WGS84\" GLOBCOVER_a.tif -ot \"Byte\" -r near -te -180 -90 0 0 -tr', 1/240, 1/240))
rsaga.geoprocessor(lib="io_gdal", 0, param=list(GRIDS="GLOBCOVER_a.sgrd", FILES="GLOBCOVER_a.tif"))
unlink("GLOBCOVER_a.tif")
system(paste(gdalwarp, 'GLOBCOVER_200412_200606_V2.2_Global_CLA.tif -t_srs \"+proj=longlat +datum=WGS84\" GLOBCOVER_b.tif -ot \"Byte\" -r near -te -180 0 0 90 -tr', 1/240, 1/240))
rsaga.geoprocessor(lib="io_gdal", 0, param=list(GRIDS="GLOBCOVER_b.sgrd", FILES="GLOBCOVER_b.tif"))
unlink("GLOBCOVER_b.tif")
system(paste(gdalwarp, 'GLOBCOVER_200412_200606_V2.2_Global_CLA.tif -t_srs \"+proj=longlat +datum=WGS84\" GLOBCOVER_c.tif -ot \"Byte\" -r near -te 0 0 180 90 -tr', 1/240, 1/240))
rsaga.geoprocessor(lib="io_gdal", 0, param=list(GRIDS="GLOBCOVER_c.sgrd", FILES="GLOBCOVER_c.tif"))
unlink("GLOBCOVER_c.tif")
system(paste(gdalwarp, 'GLOBCOVER_200412_200606_V2.2_Global_CLA.tif -t_srs \"+proj=longlat +datum=WGS84\" GLOBCOVER_d.tif -ot \"Byte\" -r near -te 0 -90 180 0 -tr', 1/240, 1/240))
rsaga.geoprocessor(lib="io_gdal", 0, param=list(GRIDS="GLOBCOVER_d.sgrd", FILES="GLOBCOVER_d.tif"))
unlink("GLOBCOVER_d.tif")


# For classes "14" and "200" derive coverage at 1 km grid:
for(i in c(2,20)){
    for(j in c("a","b","c","d")){
      rsaga.geoprocessor(lib="grid_calculus", module=1, param=list(GRIDS=paste("GLOBCOVER_", j, ".sgrd", sep=""), RESULT=paste("class_", j, sep=""), FORMULA=paste("ifelse(a=", class.lst$dbf$Value[i],",100,0)", sep="")))
    }
    # crate a mosaick:
    system(paste(gdalbuildvrt, "globedem.vrt class_a.tif class_b.tif class_c.tif class_d.tif"))    
    # resample each block to 1 km resolution:
    if(i==2){
      system(paste(gdalwarp, "globedem.vrt -t_srs \"+proj=longlat +datum=WGS84\" CRPESA3a.tif -ot \"Byte\" -r bilinear -te -180 -90 180 90 -tr", 1/120, 1/120))    
    }
    if(i==20){
      system(paste(gdalwarp, "globedem.vrt -t_srs \"+proj=longlat +datum=WGS84\" BARESA3a.tif -ot \"Byte\" -r bilinear -te -180 -90 180 90 -tr", 1/120, 1/120))    
    }
}
GDALinfo("GLCESA3a.tif")

# ------------------------------------------------------------
# obtain GLC2000 (slightly outdated):
# ------------------------------------------------------------

# donwload a copy of the Global Land Cover 2000 map:
download.file("http://bioval.jrc.ec.europa.eu/glc2000/products/glc2000_v1_1_Tiff.zip", destfile=paste(getwd(), "glc2000_v1_1_Tif.zip", sep="/"))
unzip(zipfile="glc2000_v1_1_Tiff.zip", exdir=getwd())
GDALinfo("Tiff/glc2000_v1_1.tif") # projection system missing!
# resample from 1 km to 5 km grid:
system(paste(gdalwarp, 'Tiff/glc2000_v1_1.tif -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"+proj=longlat +datum=WGS84\" GLCJRC1a.tif -dstnodata 23 -r near -te -180 -90 180 90 -tr', 6/120, 6/120))
GDALinfo("GLCJRC1a.tif")

# -------------------------------------------
# Compress produced maps:
# -------------------------------------------

for(outname in c("GLCJRC1a.tif", "GLCESA3a.tif")){
  system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), outname))
  system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
  unlink(set.file.extension(outname, ".tif.gz"))
}  # Compression takes > 15 mins

# Clean-up:
unlink("class_*.tif")
unlink("GLOBCOVER_*.tif")

# end of script;

