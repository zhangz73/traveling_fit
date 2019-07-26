library(raster)
library(ggplot2)
library(rgdal)
library(maptools)
library(gpclib)
library(sp)

map_cor_ad <- function(map_data, shapefile, vec_cor){
  ad_sp <- spTransform(raster_uganda_ad2,
            CRS("+proj=longlat +datum=WGS84"))
  cor_sp <- SpatialPointsDataFrame(map_data[,vec_cor], mapdata)
  crs(cor_sp) <- "+proj=longlat +datum=WGS84"
  o <- over(cor_sp, ad_sp)
  ret <- data.frame(cbind(data.frame(map_data[,vec_cor]), data.frame(o)))
  return(ret)
}

ret <- map_cor_ad(mapdata, raster_uganda_ad2, c(4,5))