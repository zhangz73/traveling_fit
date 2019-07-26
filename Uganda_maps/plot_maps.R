library(raster)
library(ggplot2)
library(rgdal)
library(maptools)
library(gpclib)
library(sp)

uganda_admin1 <- readOGR("J:\\WORK\\11_geospatial\\06_original shapefiles\\GAUL_admin\\admin1\\g2015_2014_1","g2015_2014_1")
mapdata <- read.csv("UGA_MACRO_data.csv", header = T)

raster_ad1 <- shapefile("J:\\WORK\\11_geospatial\\06_original shapefiles\\GAUL_admin\\admin1\\g2015_2014_1\\g2015_2014_1.shp")
raster_uganda_ad1 <- raster_ad1[raster_ad1$ADM0_NAME=="Uganda",]
raster_uganda_ad1_shp <- broom::tidy(raster_uganda_ad1, 
                                     region="ADM1_NAME")

raster_ad2 <- shapefile("J:\\WORK\\11_geospatial\\06_original shapefiles\\GAUL_admin\\admin2\\g2015_2014_2\\g2015_2014_2.shp")
raster_uganda_ad2 <- raster_ad2[raster_ad2$ADM0_NAME=="Uganda",]

raster_uganda_ad2_shp <- broom::tidy(raster_uganda_ad2, 
                                     region="ADM2_NAME")

p1 <- ggplot()
p2 <- p1 + geom_polygon(data = spTransform(raster_uganda_ad1,
      CRS("+proj=longlat +datum=WGS84")), aes(x = long, y = lat),
                        fill = NA, color = "black")
p3 <- p2 + geom_point(data = mapdata, aes(x = x, y = y, fill = pop),
                      shape = 22, size = 2) +
        scale_fill_gradient(low = "yellow", high = "red", trans="log10")
p3
p4 <- p3 + geom_point(data = raster_uganda, aes(x = long, y = lat, 
                          color = factor(raster_uganda$ADM1_NAME)),
                      size = 3, shape = 0)
p4



p1 <- ggplot()
p2 <- p1 + geom_polygon(data = spTransform(raster_uganda_ad1,
      CRS("+proj=longlat +datum=WGS84")), aes(x = long, y = lat),
                        fill = NA, color = 'Black')
p3 <- p2 + geom_polygon(data = raster_uganda_ad1_shp, 
                        aes(x = long, y = lat, fill = id, group=id))
p3



p1 <- ggplot()
p2 <- p1 + geom_polygon(data = spTransform(raster_uganda_ad2,
      CRS("+proj=longlat +datum=WGS84")), aes(x = long, y = lat),
      fill = NA, color = 'Black')
p3 <- p2 + geom_polygon(data = raster_uganda_ad2_shp, 
                      aes(x = long, y = lat, fill = id, group=id))
p3