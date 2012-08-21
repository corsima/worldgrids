# title         : getGlobCover.R
# purpose       : download and resampling of GlobCover Land Cover version V2.2 image;
# reference     : [http://worldgrids.org/doku.php?id=wiki:glcesa3]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, August 2012.
# inputs        : data available for download from [http://ionia1.esrin.esa.int];
# outputs       : geotiff images projected in the "+proj=longlat +datum=WGS84" system;
# remarks 1     : First download and install FWtools [http://fwtools.maptools.org]; 
# remarks 2     : The GlobCover Land Cover product is the highest resolution (300 meters) Global Land Cover product produced until 2011 but the image has 64800 rows and 129600 columns;

# -------------------------------------
# Initial settings and data download:
# -------------------------------------

library(rgdal)
library(RSAGA)
library(XML)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
gdaldem = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdaldem.exe"))))
gdalbuildvrt = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalbuildvrt.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"

# download files from the server:
download.file("ftp://guestglobcover:Pir8Eefo@us-ext-nas.eo.esa.int/global/Globcover_200412_200606_V2.2_Global.zip", destfile=paste(getwd(), "Globcover.zip", sep="/"))
unzip(zipfile="Globcover.zip", exdir=getwd())
GDALinfo("GLOBCOVER_200412_200606_V2.2_Global_CLA.tif")

# --------------------------------------
# Extract maps:
# --------------------------------------

class.lst <- read.dbf("GLOBCOVER_200412_200606_V2.2_Global_CLA.tif.vat.dbf")
library(gdata)
# perl <- gdata:::findPerl("perl")
perl = "C:/Perl64/bin/perl.exe"
class.names <- read.xls("Globcover_Legend.xls", perl=perl)
class.names$Label

## resample to 1, 2.5 and 5 km resolution:
system(paste(gdalwarp, "GLOBCOVER_200412_200606_V2.2_Global_CLA.tif -t_srs \"+proj=longlat +datum=WGS84\" GLCESA3a.tif -r near -te -180 -90 180 90 -tr", 1/120, 1/120))
GDALinfo("GLCESA3a.tif")
system(paste(gdalwarp, "GLOBCOVER_200412_200606_V2.2_Global_CLA.tif -t_srs \"+proj=longlat +datum=WGS84\" GLCESA3a.tif -r near -te -180 -90 180 90 -tr", 1/120, 1/120))
system(paste(gdalwarp, "GLOBCOVER_200412_200606_V2.2_Global_CLA.tif -t_srs \"+proj=longlat +datum=WGS84\" GLCESA2a.tif -r near -te -180 -90 180 90 -tr", 1/40, 1/40))
system(paste(gdalwarp, "GLOBCOVER_200412_200606_V2.2_Global_CLA.tif -t_srs \"+proj=longlat +datum=WGS84\" GLCESA1a.tif -r near -te -180 -90 180 90 -tr", 1/20, 1/20))
system(paste(gdalwarp, "GLOBCOVER_200412_200606_V2.2_Global_CLA.tif -t_srs \"+proj=longlat +datum=WGS84\" GLCESA0a.tif -r near -te -180 -90 180 90 -tr", 1/5, 1/5))

## Compress and copy:
for(outname in c("GLCESA3a.tif", "GLCESA2a.tif", "GLCESA1a.tif", "GLCESA0a.tif")){
  if(is.na(file.info(paste(shortPathName(normalizePath(outdir)), paste(outname, "gz", sep="."), sep="\\"))$size)){
  system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), outname))
  system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
  unlink(set.file.extension(outname, ".tif.gz"))
}  # Compression takes > 15 mins
}


# --------------------------------------
# Create indicator maps:
# --------------------------------------

## tiling system:
t.l <- expand.grid(KEEP.OUT.ATTRS=FALSE, lonmin=seq(-180,120,by=60), latmin=seq(-90,45,by=45))
t.u <- expand.grid(KEEP.OUT.ATTRS=FALSE, lonmax=seq(-120,180,by=60), latmax=seq(-45,90,by=45))
tiles <- cbind(t.l, t.u)

## split world in 24 blocks:
# THIS TAKES >2 hrs
for(j in 1:nrow(tiles)){
  if(is.na(file.info(paste('glc_', j, '.sdat', sep=""))$size)){  
  # resample to 500 m resolution:
  system(paste(gdalwarp, ' GLOBCOVER_200412_200606_V2.2_Global_CLA.tif -t_srs \"+proj=longlat +datum=WGS84\" glc_', j, '.sdat -of \"SAGA\" -ot \"Byte\" -r near -te ', tiles[j,"lonmin"] , ' ', tiles[j,"latmin"], ' ', tiles[j,"lonmax"] ,' ', tiles[j,"latmax"] ,' -tr ', 1/240, ' ', 1/240, sep=""))
  }
  # mask values per each class:
  for(i in 1:length(levels(as.factor(class.lst$dbf$Value)))){
    if(is.na(file.info(paste('glc_', j, '_', i, 'r.sdat', sep=""))$size)){
    rsaga.geoprocessor(lib="grid_calculus", module=1, param=list(GRIDS=paste('glc_', j, '.sgrd', sep=""), RESULT=paste('glc_', j, '_', i, '.sgrd', sep=""), FORMULA=paste("ifelse(a=", class.lst$dbf$Value[i],",100,0)", sep="")), show.output.on.console = FALSE)
    # resample to 1 km resolution:
    system(paste(gdalwarp, ' glc_', j, '_', i, '.sdat -t_srs \"+proj=longlat +datum=WGS84\" glc_', j, '_', i, 'r.sdat -of \"SAGA\" -ot \"Byte\" -r bilinear -te ', tiles[j,"lonmin"] , ' ', tiles[j,"latmin"], ' ', tiles[j,"lonmax"] ,' ', tiles[j,"latmax"] ,' -tr ', 1/120, ' ', 1/120, sep=""), show.output.on.console = FALSE)
    # clean up:
    unlink(paste('glc_', j, '_', i, '.sdat', sep=""))  
    unlink(paste('glc_', j, '_', i, '.sgrd', sep=""))  
    unlink(paste('glc_', j, '_', i, '.mgrd', sep=""))  
  }
  }
}
# clean up class tiles:
for(j in 1:nrow(tiles)){ unlink(paste('glc_', j, '.****', sep="")) }

## glue all tiles together (per class):  
for(i in 1:length(levels(as.factor(class.lst$dbf$Value)))){
  if(i < 10){ out <- paste("G0", i, "ESA3a", sep="") } else { out <- paste("G", i, "ESA3a", sep="") }
  if(is.na(file.info(paste(out, '.tif', sep=""))$size)){
  # list files:
  lst <- list.files(pattern=glob2rx(paste('glc_*_', i, 'r.sdat', sep="")))
  # create a mosaic:
  system(paste(gdalbuildvrt, "glc.vrt", paste(lst, collapse=" ")))
  # resample:
  system(paste(gdalwarp, ' glc.vrt ', out, '.tif -t_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -r near -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""))    
  }
}

for(i in 1:length(levels(as.factor(class.lst$dbf$Value)))){
  ## 2.5 and 5.6 km:
  if(i < 10){ out1 <- paste("G0", i, "ESA0a", sep="") } else { out1 <- paste("G", i, "ESA0a", sep="") } 
  if(i < 10){ out2 <- paste("G0", i, "ESA1a", sep="") } else { out2 <- paste("G", i, "ESA1a", sep="") }
  if(i < 10){ out3 <- paste("G0", i, "ESA2a", sep="") } else { out3 <- paste("G", i, "ESA2a", sep="") }    
  lst <- list.files(pattern=glob2rx(paste('glc_*_', i, 'r.sdat', sep="")))
  # create a mosaic:
  system(paste(gdalbuildvrt, "glc.vrt", paste(lst, collapse=" ")))
  # resample:
  if(is.na(file.info(paste(out2, '.tif', sep=""))$size)){
  system(paste(gdalwarp, ' glc.vrt ', out2, '.tif -t_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -r bilinear -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
  }
  if(is.na(file.info(paste(out3, '.tif', sep=""))$size)){
  system(paste(gdalwarp, ' glc.vrt ', out3, '.tif -t_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -r bilinear -te -180 -90 180 90 -tr ', 1/40,' ', 1/40, sep=""))
  }
  if(is.na(file.info(paste(out1, '.tif', sep=""))$size)){
  system(paste(gdalwarp, ' glc.vrt ', out1, '.tif -t_srs \"+proj=longlat +datum=WGS84\" -ot \"Byte\" -r bilinear -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))
  }
}
unlink("glc.vrt")

## Compress and copy:
lst <- list.files(pattern=glob2rx('G**ESA*a.tif$'))
for(outname in lst){
  if(is.na(file.info(paste(shortPathName(normalizePath(outdir)), paste(outname, "gz", sep="."), sep="\\"))$size)){
  system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), outname))
  system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
  unlink(set.file.extension(outname, ".tif.gz"))
}  # Compression takes > 15 mins
}

## clean up temp files
unlink(list.files(pattern=glob2rx(paste('glc_*_', i, 'r.sdat', sep=""))))
unlink(list.files(pattern=glob2rx(paste('glc_*_', i, 'r.sgrd', sep=""))))


# ------------------------------------------------------------
# GLC2000 image (slightly outdated):
# ------------------------------------------------------------

# donwload a copy of the Global Land Cover 2000 map:
download.file("http://bioval.jrc.ec.europa.eu/glc2000/products/glc2000_v1_1_Tiff.zip", destfile=paste(getwd(), "glc2000_v1_1_Tif.zip", sep="/"))
unzip(zipfile="glc2000_v1_1_Tiff.zip", exdir=getwd())
GDALinfo("Tiff/glc2000_v1_1.tif") # projection system missing!
# resample from 1, 2.5 km and 5 km grid:
system(paste(gdalwarp, 'Tiff/glc2000_v1_1.tif -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"+proj=longlat +datum=WGS84\" GLCJRC0a.tif -dstnodata 23 -r near -te -180 -90 180 90 -tr', 1/5, 1/5))
system(paste(gdalwarp, 'Tiff/glc2000_v1_1.tif -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"+proj=longlat +datum=WGS84\" GLCJRC1a.tif -dstnodata 23 -r near -te -180 -90 180 90 -tr', 1/20, 1/20))
system(paste(gdalwarp, 'Tiff/glc2000_v1_1.tif -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"+proj=longlat +datum=WGS84\" GLCJRC2a.tif -dstnodata 23 -r near -te -180 -90 180 90 -tr', 1/40, 1/40))
system(paste(gdalwarp, 'Tiff/glc2000_v1_1.tif -s_srs \"+proj=longlat +datum=WGS84\" -t_srs \"+proj=longlat +datum=WGS84\" GLCJRC3a.tif -dstnodata 23 -r near -te -180 -90 180 90 -tr', 1/120, 1/120))
GDALinfo("GLCJRC1a.tif")
## Compress and copy:
lst <- list.files(pattern=glob2rx('GLCJRC*a.tif$'))
for(outname in lst){
  if(is.na(file.info(paste(shortPathName(normalizePath(outdir)), paste(outname, "gz", sep="."), sep="\\"))$size)){
  system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), outname))
  system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
  unlink(set.file.extension(outname, ".tif.gz"))
}  # Compression takes > 15 mins
}

# ------------------------------------
# create the SLD files:
# ------------------------------------

# read the legend:
glc_leg <- read.xls("Tiff/Global_Legend.xls", perl=perl)
str(glc_leg)
Format_Information_Content = "GTiff"
obj.name = "GLCJRC1a.tif"
sld.file = "GLCJRC1a.sld"
Citation_title = "Global Land Cover 2000"
ColorMap_type = "intervals"
pal <- rgb(red=glc_leg$Red, green=glc_leg$Green, blue=glc_leg$Blue) 
label <- glc_leg$CLASSNAMES 
bounds <- glc_leg$VALUE
opacity = c(rep(1, length(label)-1), 0)  
# last class full transparency!

l1 = newXMLNode("StyledLayerDescriptor", attrs=c(version="1.0.0"), namespaceDefinitions=c("xsi:schemaLocation"="http://www.opengis.net/sld StyledLayerDescriptor.xsd", "sld"="http://www.opengis.net/sld", "ogc"="http://www.opengis.net/ogc", "gml"="http://www.opengis.net/gml"))
l2 <- newXMLNode("NamedLayer", parent = l1)
l3 <- newXMLNode("Name", paste(Citation_title, "(", Format_Information_Content, ")", sep=""), parent = l2)
l3b <- newXMLNode("UserStyle", parent = l2)
l4 <- newXMLNode("Title", paste(obj.name, "style", sep="_"), parent = l3b)
l4b <- newXMLNode("FeatureTypeStyle", parent = l3b)
l5 <- newXMLNode("Rule", parent = l4b)
l6 <- newXMLNode("RasterSymbolizer", parent = l5)
l7 <- newXMLNode("ColorMap", attrs=c(type=ColorMap_type), parent = l6)
txt <- sprintf('<ColorMapEntry color="#%s" quantity="%.2f" label="%s" opacity="%.1f"/>', pal, bounds, label, opacity)
parseXMLAndAdd(txt, l7)
saveXML(l1, sld.file)

system(paste("xcopy", sld.file, shortPathName(normalizePath("D:/WORLDGRIDS/sld"))))

# end of script;

