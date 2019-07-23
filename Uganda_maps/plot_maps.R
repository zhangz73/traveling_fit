library(raster)
library(ggplot2)
library(rgdal)

uganda_admin1 <- readOGR("J:\\WORK\\11_geospatial\\06_original shapefiles\\GAUL_admin\\admin1\\g2015_2014_1","g2015_2014_1")
mapdata <- read.csv("UGA_MACRO_data.csv", header = T)

raster_ad1 <- shapefile("J:\\WORK\\11_geospatial\\06_original shapefiles\\GAUL_admin\\admin1\\g2015_2014_1\\g2015_2014_1.shp")
raster_uganda <- raster_ad1[raster_ad1$ADM0_NAME=="Uganda",]

p1 <- ggplot()
p2 <- p1 + geom_polygon(data = spTransform(raster_uganda,
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