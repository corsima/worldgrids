
# title         : getGSMaP.R
# purpose       : download and resampling of the Global Satellite Mapping of Precipitation (GSMaP) products;
# reference     : [https://code.google.com/p/worldgrids/]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, June 2012.
# inputs        : images available at [http://sharaku.eorc.jaxa.jp/GSMaP_crest/]; 
# outputs       : geotiff images projected in the "+proj=longlat +datum=WGS84" system;
# remarks 1     : Description of the data set available at [ftp://hokusai.eorc.jaxa.jp/pub/gsmap_crest/gsmap_dataformat.pdf]; 
# remarks 2     : First download and install FWtools [http://fwtools.maptools.org];

# ------------------------------------
# Initial settings and data download:
# ------------------------------------

library(rgdal)
library(R.utils)
library(plotKML)
require(RColorBrewer)

fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
gdal_translate = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdal_translate.exe"))))
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"

download.file("ftp://hokusai.eorc.jaxa.jp/pub/gsmap_crest/MVK+/monthly/MVK+_2003-2006.tar.gz", destfile=paste(getwd(), "MVK+_2003-2006.tar.gz", sep="/"), mode='wb', method='wget')
# unzip:
gunzip("MVK+_2003-2006.tar.gz")
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
system("7za e -ttar MVK+_2003-2006.tar")   # 49 files


# list all binary files in dir:
pr.list <- dir(path=getwd(), pattern=glob2rx("*.v484"))
for(j in 1:length(pr.list)){
  pr.name <-  paste("PRE", substr(pr.list[j], nchar(pr.list[j])-25, nchar(pr.list[j])-20), sep="")
  # convert to SAGA format:
  rsaga.geoprocessor(lib="io_grid", module=4, param=list(GRID=set.file.extension(pr.name, ".sgrd"), FILE_DATA=pr.list[j], NX=3600, NY=1200, DXY=0.1, XMIN=0.05, YMIN=-59.95, NODATA=-999, DATA_OFFSET=0, LINE_OFFSET=0, DATA_TYPE=6, BYTEORDER_BIG=0, TOPDOWN=1))
}

# ------------------------------------
# Derive mean and sd:
# ------------------------------------

prg.list <- dir(path=getwd(), pattern=glob2rx("*.sgrd"))
rsaga.geoprocessor(lib="geostatistics_grid", module=4, param=list(GRIDS=paste(prg.list, collapse=";"), MEAN="PREm.sgrd", MIN="tmp.sgrd", MAX="tmp.sgrd", VAR="tmp.sgrd", STDDEV="PREs.sgrd", STDDEVLO="tmp.sgrd", STDDEVHI="tmp.sgrd"))

# read to R:
grids <- readGDAL("PREm.sdat")
grids$m <- grids$band1*24*30    # convert to mm/month:
grids$s <- readGDAL("PREs.sdat")$band1*100
grids$band1 <- NULL
grids.pnt <- data.frame(grids)
# fix coordinates:
grids.pnt$x <- ifelse(grids.pnt$x > 180, grids.pnt$x-360, grids.pnt$x)
coordinates(grids.pnt) <- ~x+y   # takes cca 10 mins;
gridded(grids.pnt) <- TRUE
fullgrid(grids.pnt) <- TRUE
proj4string(grids.pnt) <- CRS("+proj=longlat +datum=WGS84")
writeGDAL(grids.pnt[1], "PRECm.sdat", "SAGA", mvFlag=0)
writeGDAL(grids.pnt[2], "PRECs.sdat", "SAGA", mvFlag=0)
GDALinfo("PRECm.sdat")
# resample to 5 km grid:
system(paste(gdalwarp, "PRECm.sdat -t_srs \"+proj=longlat +datum=WGS84\" PRECm5km.sdat -of \"SAGA\" -srcnodata 0 -dstnodata 0 -r bilinear -te -180 -90 180 90 -tr 0.05 0.05"))

# Obtain the PREC map from worldclim and average the two rainfall maps:
download.file("http://biogeo.berkeley.edu/worldclim1_4/grid/cur/bio_2-5m_bil.zip", destfile=paste(getwd(), "bio_2-5m_bil.zip", sep="/"))
unzip(zipfile="bio_2-5m_bil.zip", exdir=getwd())
system(paste(gdalwarp, "bio12.bil -t_srs \"+proj=longlat +datum=WGS84\" biocl12.sdat -of \"SAGA\" -dstnodata 55537 -r bilinear -te -180 -90 180 90 -tr 0.05 0.05"))
# mask the missing values:
rsaga.geoprocessor(lib="grid_calculus", module=1, param=list(INPUT="biocl12.sgrd", RESULT="biocl12f.sgrd", FORMUL="ifelse(a=55537,65535,a)")) 
rsaga.geoprocessor(lib="grid_calculus", module=1, param=list(INPUT="PRECm.sgrd", RESULT="PREC.sgrd", FORMUL="a*12")) # annual rainfall;
# derive mean value:
rsaga.geoprocessor(lib="geostatistics_grid", module=4, param=list(GRIDS="PREC.sgrd;biocl12f.sgrd", MEAN="PREGSM1a.sgrd", MIN="tmp.sgrd", MAX="tmp.sgrd", VAR="tmp.sgrd", STDDEV="tmp.sgrd", STDDEVLO="tmp.sgrd", STDDEVHI="tmp.sgrd"))

# convert to Geotif:
system(paste(gdal_translate, "PREGSM1a.sdat PREGSM1a.tif -a_srs \"+proj=longlat +datum=WGS84\""))
GDALinfo("PREGSM1a.tif")
system(paste(gdalwarp, "PREGSM1a.sdat -t_srs \"+proj=longlat +datum=WGS84\" PREGSM0a.tif -srcnodata 0 -dstnodata 0 -r bilinear -te -180 -90 180 90 -tr 0.2 0.2"))

# ------------------------------------
# Compress produced maps:
# ------------------------------------

for(outname in c("PREGSM1a.tif", "PREGSM0a.tif")){
  if(is.na(file.info(paste(shortPathName(normalizePath(outdir)), paste(outname, "gz", sep="."), sep="\\"))$size)){
    system(paste("7za a", "-tgzip", set.file.extension(outname, ".tif.gz"), outname))
    system(paste("xcopy", set.file.extension(outname, ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(outname, ".tif.gz"))
  }
}  # Compression takes > 15 mins

# Clean-up:
rm(grids.pnt)
rm(grids)

# ------------------------------------
# create the SLD file:
# ------------------------------------

Format_Information_Content = "GTiff"
obj.name = "PREGSM1a.tif"
sld.file = "PREGSM1a.sld"
Citation_title = "Total annual precipitation"
ColorMap_type = "intervals"
opacity = 1
data(SAGA_pal)
pal <- colorRamp(SAGA_pal[["SG_COLORS_WHITE_BLUE"]], space = "rgb", interpolate = "linear") 
bounds <- round(seq(0, 1734.5, length.out=100), 0)
obj = data.frame(ranks=1:100, color=rgb(pal(scales::rescale(bounds))/ 255))

l1 = newXMLNode("StyledLayerDescriptor", attrs=c(version="1.0.0"), namespaceDefinitions=c("xsi:schemaLocation"="http://www.opengis.net/sld StyledLayerDescriptor.xsd", "sld"="http://www.opengis.net/sld", "ogc"="http://www.opengis.net/ogc", "gml"="http://www.opengis.net/gml"))
l2 <- newXMLNode("NamedLayer", parent = l1)
l3 <- newXMLNode("Name", paste(Citation_title, "(", Format_Information_Content, ")", sep=""), parent = l2)
l3b <- newXMLNode("UserStyle", parent = l2)
l4 <- newXMLNode("Title", paste(obj.name, "style", sep="_"), parent = l3b)
l4b <- newXMLNode("FeatureTypeStyle", parent = l3b)
l5 <- newXMLNode("Rule", parent = l4b)
l6 <- newXMLNode("RasterSymbolizer", parent = l5)
l7 <- newXMLNode("ColorMap", attrs=c(type=ColorMap_type), parent = l6)
txt <- sprintf('<ColorMapEntry color="#%s" quantity="%.2f" label="%s" opacity="%.1f"/>', obj$color, bounds, bounds, rep(opacity, length(obj$color)))
parseXMLAndAdd(txt, l7)
saveXML(l1, sld.file)

system(paste("xcopy", "PREGSM1a.sld", shortPathName(normalizePath("G:/WORLDGRIDS/sld"))))

# end of script; 
