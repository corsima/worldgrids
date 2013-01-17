# title         : getBioVal.R
# purpose       : download and resampling of the global map of Accessibility;
# reference     : [http://worldgrids.org/doku.php?id=wiki:gacgem1]
# producer      : Prepared by T. Hengl
# last update   : In Wageningen, NL, January 2013.
# inputs        : map available for download from [http://bioval.jrc.ec.europa.eu/products/gam/download.htm];
# outputs       : geotiff images projected in the "+proj=longlat +ellps=WGS84" system;
# remarks 1     : First download and install FWtools [http://fwtools.maptools.org] and ILWIS GIS;  

library(rgdal)
fw.path = utils::readRegistry("SOFTWARE\\WOW6432Node\\FWTools")$Install_Dir
gdalwarp = shQuote(shortPathName(normalizePath(file.path(fw.path, "bin/gdalwarp.exe"))))
download.file("http://downloads.sourceforge.net/project/sevenzip/7-Zip/9.20/7za920.zip", destfile=paste(getwd(), "/", "7za920.zip", sep="")) 
unzip("7za920.zip")
si <- Sys.info()
outdir <- "G:/WORLDGRIDS/maps"

download.file("http://bioval.jrc.ec.europa.eu/products/gam/download/access_50k.zip", "access_50k.zip")
unzip("access_50k.zip")
GDALinfo("access_50k/acc_50k")

## resample to coarser resolutions:
system(paste(gdalwarp, ' access_50k/acc_50k GACGEM3a.tif -dstnodata \"-9999\" -r near -te -180 -90 180 90 -tr ', 1/120,' ', 1/120, sep=""))
system(paste(gdalwarp, ' access_50k/acc_50k GACGEM2a.tif -dstnodata \"-9999\" -r bilinear -te -180 -90 180 90 -tr ', 1/40,' ', 1/40, sep=""))
system(paste(gdalwarp, ' access_50k/acc_50k GACGEM1a.tif -dstnodata \"-9999\" -r bilinear -te -180 -90 180 90 -tr ', 1/20,' ', 1/20, sep=""))
system(paste(gdalwarp, ' access_50k/acc_50k GACGEM0a.tif -dstnodata \"-9999\" -r bilinear -te -180 -90 180 90 -tr ', 1/5,' ', 1/5, sep=""))

## compress and copy:
tif.lst <- list.files(pattern=glob2rx("GACGEM*a.tif$"))  
for(i in 1:length(tif.lst)){
    system(paste("7za a", "-tgzip", set.file.extension(tif.lst[i], ".tif.gz"), tif.lst[i]))
    system(paste("xcopy", set.file.extension(tif.lst[i], ".tif.gz"), shortPathName(normalizePath(outdir)))) 
    unlink(set.file.extension(tif.lst[i], ".tif.gz"))
}

# end of script;