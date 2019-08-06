library(raster)
library(ggplot2)
library(rgdal)
library(maptools)
library(gpclib)
library(sp)
library(parallel)

no_cores <- detectCores() - 1

map_cor_ad <- function(map_data, shapefile, vec_cor){
  ad_sp <- spTransform(shapefile,
            CRS("+proj=longlat +datum=WGS84"))
  cor_sp <- SpatialPointsDataFrame(map_data[,vec_cor], map_data)
  crs(cor_sp) <- "+proj=longlat +datum=WGS84"
  
  n <- 100
  parts <- split(x = 1:nrow(cor_sp), f = cut(1:nrow(cor_sp), n))
  cl <- makeCluster(no_cores, type = "FORK")
  
  o <- parLapply(cl = cl, X = 1:n, fun = function(x) {
    ov <- over(cor_sp[parts[[x]],], ad_sp)
    gc()
    return(ov)
  })
  
  ret <- data.frame(cbind(data.frame(map_data), data.frame(o)))
  return(ret)
}
