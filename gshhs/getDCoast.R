
# title         : getDCoast.R
# purpose       : derivation of the distance to coast line map;
# reference     : [http://spatial-analyst.net/wiki/index.php?title=Global_datasets]
# producer      : Prepared by T. Hengl
# last update   : In Amsterdam, NL, March 2010.
# inputs        : land mask map previously derived;
# outputs       : geotiff images projected in the "+proj=longlat +ellps=WGS84" system;
# remarks 1     : First download and install FWtools [http://fwtools.maptools.org] and ILWIS GIS;
# remarks 2     : Distances are expressed in arcdegrees, which is not really optimal (distortions on the North pole are significant);   


# ------------------------------------------------------------
# Initial settings and data download:
# ------------------------------------------------------------

library(rgdal)
ILWIS <- "C:\\Progra~1\\N52\\Ilwis35\\IlwisClient.exe -C"

# download image from the server:
download.file("http://spatial-analyst.net/worldmaps/landmask.zip", destfile=paste(getwd(), "landmask.zip", sep="/"))
unzip(zipfile="landmask.zip", exdir=getwd())
# import to ILWIS:
shell(cmd=paste(ILWIS, " import tiff(landmask.tif, landmask)", sep=""), wait=F)

# derived boolean maps:
shell(cmd=paste(ILWIS, " landmasko{dom=Bool.dom} = iff(landmask>-31349, ?, 1)", sep=""), wait=F)
shell(cmd=paste(ILWIS, " landmaskw{dom=Bool.dom} = iff(landmask>-31349, 1, ?)", sep=""), wait=F)

# derive distance maps:
shell(cmd=paste(ILWIS, " dcoast_o.mpr{dom=value.dom;vr=0.00:100.00:0.01} = MapDistance(landmasko)", sep=""), wait=F)
shell(cmd=paste(ILWIS, " dcoast_w.mpr{dom=value.dom;vr=0.00:100.00:0.01} = MapDistance(landmaskw)", sep=""), wait=F)

# combine the two maps:
shell(cmd=paste(ILWIS, " dcoast{dom=value.dom;vr=-100.00:100.00:0.01} = -dcoast_w+dcoast_o", sep=""), wait=F)

# convert to a tiff image;
system(paste("C:\\PROGRA~1\\FWTOOL~1.7\\bin\\gdal_translate dcoast.mpr -a_srs \"+proj=longlat +ellps=WGS84\" dcoast.tif -a_nodata -9999", sep=""))

# end of script;
