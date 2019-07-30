library(igraph)
library(maps)
library(geosphere)
library(ggplot2)
library(rgdal)
library(maptools)
library(gpclib)
library(raster)
library(data.table)

source("map_coor_ad.R")

plot_w <- function(url_to_shp, w_df){
  area <- shapefile(url_to_shp)
  area_sp <- spTransform(area, CRS("+proj=longlat +datum=WGS84"))
  dt <- data.table(w_df)
  small_df <- dt[,.(lon = mean(s_x), lat = mean(s_y)), by = sad2]

  p1 <- ggplot()
  p2 <- p1 + geom_polygon(data = area_sp, aes(x = long, y = lat),
                          fill = NA, color = "black")
  p3 <- p2 + geom_point(data = small_df, aes(x = lon, y = lat,
                          fill = sad2, group = sad2), shape = 22,
                        size = 5)

  col.1 <- adjustcolor("orange red", alpha=0.4)
  col.2 <- adjustcolor("orange", alpha=0.4)
  edge.pal <- colorRampPalette(c(col.1, col.2), alpha = TRUE) 
  edge.col <- edge.pal(100)
  
  p4 <- p3
  for(i in 1:nrow(w_df)){
    row <- w_df[i,]
    
    arc <- gcIntermediate(c(row$s_x, row$s_y), c(row$d_x, row$d_y),
                          n = 1000, addStartEnd = T)
    edge.ind <- round(100 * row$w / max(w_df$w))
    
    if(edge.ind > 0){
      arc <- data.frame(arc)
      p4 <- p4 + geom_line(data = arc, aes(x = lon, y = lat),
                          col = edge.col[edge.ind],
                         lwd = edge.ind / 30,
                         arrow = arrow(length = 
                            unit(0.25, "inches")))
    }
  }
  p4 <- p4 + xlab('long') + ylab('lat')
  
  print(p4)
}

## Plot Commuting Flow For Bioko Island
map_data <- read.csv("Bioko/BI_survey_data.csv", header = T)
map_data <- data.table(map_data)
loc_coor <- map_data[, by = ad2, .(x = mean(lon), y = mean(lat))]
w_df_id <- c()
locs <- c("Baney", "Malabo", "Peri", "Luba", "Riaba", "Moka", "Ureka")
for(i in 1:nrow(map_data)){
  row <- map_data[i,]
  w_ban <- row$ti_ban
  w_mal <- row$ti_mal / 2
  w_per <- row$ti_mal / 2
  w_lub <- row$ti_lub
  w_ria <- row$ti_ria
  w_mok <- row$ti_mok
  w_ure <- row$ti_ure
  ws <- c(w_ban, w_mal, w_per, w_lub, w_ria, w_mok, w_ure)
  for(j in 1:length(locs)){
    curr <- data.frame(sad2 = row$ad2, dad2 = locs[j], w = ws)
    w_df_id <- data.frame(rbind(data.frame(w_df_id), curr))
  }
}
cors <- c()
for(i in 1:nrow(w_df_id)){
  row <- w_df_id[i,]
  c1 <- loc_coor[ad2==row$sad2, .(x, y)]
  c2 <- loc_coor[ad2==row$dad2, .(x, y)]
  curr <- data.frame(s_x = c1$x, s_y = c1$y, d_x = c2$x, d_y = c2$y)
  cors <- data.frame(rbind(data.frame(cors), curr))
}
w_df0 <- data.frame(cbind(data.frame(w_df_id), data.frame(cors)))
w_df0 <- data.table(w_df0)
w_df <- w_df0[sad2 != dad2, .(w = sum(w), s_x = mean(s_x), 
                  s_y = mean(s_y), d_x = mean(d_x), 
                  d_y = mean(d_y)), by = .(sad2, dad2)]


## Plot Commuting Flow For Uganda
uga_dt <- fread("uga_pred_with_ad2.csv", sep = ",", header = T)
areaId_ad2 <- unique(uga_dt[, .(areaId = source_areaId, ad2 = ad2)])

uga_rel_dt <- fread("uga_rel_data.csv", sep = ",", header = T)



shp_dir <- "J:\\WORK\\11_geospatial\\06_original shapefiles\\GAUL_admin\\admin2\\g2015_2014_2\\g2015_2014_2.shp"
plot_w(shp_dir, w_df)
