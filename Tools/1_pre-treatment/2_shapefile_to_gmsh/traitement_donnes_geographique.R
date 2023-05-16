require(rgdal)
library("leaflet")
library("raster")




res = c(19143,2082, 7311,19107,
        19080,
        19141,
        18896,
        18685,
        18761,
        18782,
        18832,
        18538,
        18593,
        18348,
        18305,
        18201,
        18219,
        18426,
        2661,
        2222,
        2196,
        2266,
        1666,
        1444,
        1512,
        1262,
        1196,
        1077,
        1141,
        948	,
        1316,
        19789,
        19386,
        19321,
        2051,
        508	,
        19549,
        19860,
        19903,
        20053,
        20161,
        20184,
        1587,
        1407,
        1082,
        525,
        657	,
        469	,
        293	,
        19163,
        2245,
        1167,
        1782,
        928	,
        727	,
        214	,
        359	,
        20085, 1)

# ===============================
# source data 
# ==============================
shape <- readOGR( "/home/livillenave/Téléchargements/contour.shp")
contour = rgeos::gUnaryUnion(shape)
# writeOGR(contour,  "/home/livillenave/Documents/data/AUDE/true/full_lit_mineur.shp", 
#          driver = "ESRI Shapefile")
laisse_crue1 =readOGR( "/home/livillenave/Documents/data/DONNEES/prepared_data/laisses_crue/TOPODOC/2018_OCT_15_AUDE_EZI_S.shp")
laisse_crue2 =readOGR( "/home/livillenave/Documents/data/DONNEES/prepared_data/laisses_crue/TOPODOC/2018_OCT_15_ORBI_EZI_S.shp")
laisse_crue3 =readOGR( "/home/livillenave/Documents/data/DONNEES/prepared_data/laisses_crue/CEREMA/20180CT15AUDE_EZI_S.shp")
laisse_crue4 =readOGR( "/home/livillenave/Documents/data/DONNEES/prepared_data/laisses_crue/AGERIN/2018OCT15AUDE_EZI_S.shp")
laisse_crue5 =readOGR( "/home/livillenave/Documents/data/DONNEES/prepared_data/laisses_crue/IMSRN/2018OCT15AUDE_EZI_S.shp")
whole_catchment = readRDS("/home/livillenave/Documents/data/CONFIG_500_lag/Y1612010/contour.rds")



library("sp")
library("rgdal")
coords = contour@polygons[[1]]@Polygons[[1]]@coords
sp_line <- SpatialLines(list(Lines(list(Line(coords)), ID=1)))
# set coordinate reference system with SpatialLines(..., proj4string=CRS(...))
# e.g. CRS("+proj=longlat +datum=WGS84")
sp_line_df <- SpatialLinesDataFrame(sp_line, data=data.frame(ID=1))



# ===============================
# prepare
# ==============================

laisse_crue <- union(laisse_crue1, laisse_crue2)
laisse_crue <- union(laisse_crue, laisse_crue3)
laisse_crue <- union(laisse_crue, laisse_crue4)
laisse_crue <- union(laisse_crue, laisse_crue5)
bc <- shape[shape@data$ID %in% res,]

writeOGR(sp_line_df, "/home/livillenave/Documents/data/AUDE/true/", layer="full_lit_mineur", driver="ESRI Shapefile")

writeOGR(laisse_crue, "/home/livillenave/Documents/data/AUDE/true/", layer="laisse_crue", driver="ESRI Shapefile")

# ===============================
# project
# ==============================
shape2 = spTransform(shape, sp::CRS("+init=epsg:4326"))
contour2 = spTransform(contour, sp::CRS("+init=epsg:4326"))
laisse_crue  <-  spTransform(laisse_crue, sp::CRS("+init=epsg:4326"))
bc2 <-  spTransform(bc, sp::CRS("+init=epsg:4326"))
whole_catchment2 <- spTransform(whole_catchment, sp::CRS("+init=epsg:4326"))

z = read.csv("/home/livillenave/Documents/data/CONFIG_500/Y1612010/lilian_info.csv")


z$calibration_gauge = z$code_bv %in% c("Y1135010",
                                       "Y1232010",
                                       "Y1314010",
                                       "Y1364010",
                                       "Y1605050")


z$calibration_gauge <- ifelse(z$calibration_gauge, "red", "blue")

coordinates(z) <- c("xl93_model", "yl93_model")
crs(z)  <- CRS("+init=epsg:2154") 
z2 <- sp::spTransform(x = z, CRSobj = CRS("+proj=longlat +datum=WGS84") )
#z2 <- as.data.frame(z2)






t = readRDS("/home/livillenave/Documents/data/CONFIG_500_lag/Y1612010/coupled_files/500/compare/treshold_500.rds")
# domain <- projectRaster(from = t$domain, crs = CRS("+proj=longlat +datum=WGS84"))

domain <- sp::spTransform(x = t$shape, CRSobj = CRS("+proj=longlat +datum=WGS84") )
plot(domain)                   


r = t$domain
values(r) <- NA
values(r)[t$id_connections$id_inflowed_cells] <- 0
values(r)[t$id_connections$id_main_inflows] <- 1
values(r)[t$id_connections$id_lat_inflows] <- 2
r <- projectRaster(from = r, crs = CRS("+proj=longlat +datum=WGS84"), method = "ngb")



main_inflow = as.data.frame(xyFromCell(t$domain, cell = t$id_connections$id_main_inflows))
coordinates(main_inflow) <- c("x", "y")
crs(main_inflow) <- CRS("+init=epsg:2154") 
main_inflow <- sp::spTransform(x = main_inflow, CRSobj = CRS("+proj=longlat +datum=WGS84") )



m <- leaflet() %>%
  addTiles() %>%  # Add default OpenStreetMap map tiles  
  addPolygons(data = whole_catchment2, color = "black", weight = 5,
              opacity = 1, stroke = TRUE, fill = FALSE,
  group = "whole catchment") %>%
  addPolygons(data = shape2, color = "black", weight = 0.5, stroke = TRUE, fill = FALSE,
              group = "mesh") %>%
  addPolygons(data = contour2, color = "red", stroke=TRUE, weight=1,
              group = "contour_hydrau" ) %>%
  addPolygons(data = laisse_crue, color = "blue", stroke=TRUE, weight=1, opacity = 0.5,
              group = "laisse_crue" ) %>%
  addCircleMarkers(data = z2, color = ~calibration_gauge, fillColor = ~calibration_gauge,
                   popup = paste0("Code BH : " , z2@data$code_bv),
                   radius = 2, 
                   stroke = TRUE, weight = 4, opacity = 1,
                   group = "stations"
  )%>%
  addPolygons(data = domain, color = "yellow", stroke=FALSE, weight=1,
              group = "contour_hydro" )%>%
  addCircleMarkers(data = main_inflow, color = "purple", fillColor = "purple",
                   radius = 2, 
                   stroke = TRUE, weight = 4, opacity = 1,
                   group = "main_inflow") %>% 
  addPolygons(data = bc2, color = "green", weight = 5,opacity = 1,
              group = "bc") 

# on embellit la carte 
m <- m %>%
  addScaleBar(
    position = c("bottomright"),
    options = scaleBarOptions()
  )

# m on ajoute l'interactivité
m <- m %>%
  addLayersControl(
    baseGroups = c("whole catchment"),
    overlayGroups = c("contour_hydrau",  "mesh", "laisse_crue", "stations", "contour_hydro", "main_inflow", "bc"),
    options = layersControlOptions(collapsed = FALSE)
  )

m <- m %>% hideGroup(group = c("contour_hydrau",  "mesh", "laisse_crue", "stations", "contour_hydro", "main_inflow", "bc"))
m




###############################################
###############################################

require(osmdata)
require(sf)
library(lwgeom)
require(ggplot2)

coord = c(2.252356,48.784685,2.443243,48.946392)
city = "Paris"
country = "France"
EPSG = 4326 

# obtain coordinates for ggplot
bbx <- getbb(paste0(city, ", ", country))

bbx[1,1] = coord[1]
bbx[2,1] = coord[2]

bbx[1,2] = coord[3]
bbx[2,2] = coord[4]

q0 <- opq(bbox = coord) 


# Defining a buffer
buffer <- 0
p_big <- rbind(c(bbx[1] - buffer, bbx[2] - buffer),
               c(bbx[1] - buffer, bbx[4] + buffer),
               c(bbx[3] + buffer, bbx[4] + buffer),
               c(bbx[3] + buffer, bbx[2] - buffer),
               c(bbx[1] - buffer, bbx[2] - buffer))

# Putting the coordinates into a squared polygon object
pol <- st_polygon(list(p_big)) %>% st_geometry

# Providing the SRID (here, unprojected lon/lat)
st_crs(pol) <- EPSG



# Get water data
water <- opq(bbox = st_bbox(pol)) %>%
  add_osm_feature(key = "water") %>%
  osmdata_sf()

# I thought this would be enough but when I plot it, it was incomplete
# So I need to get more specific query

# Get river data
river <- opq(bbox = st_bbox(pol)) %>%
  add_osm_feature(key = 'waterway', value = "river") %>%
  osmdata_sf()


ggplot() + 
  geom_sf(data=water$osm_polygons, inherit.aes = FALSE, lwd=0, fill = "green", color = "green") + 
  geom_sf(data=water$osm_multipolygons, inherit.aes = FALSE, lwd=0, fill = "blue", color = "blue") +
  geom_sf(data=river$osm_multilines, inherit.aes = FALSE, lwd=0, fill = "red", color = "red") + 
  coord_sf(xlim = c(min(p_big[,1]), max(p_big[,1])), ylim = c(min(p_big[,2]), max(p_big[,2])), expand = FALSE)






q <- opq(bbox = c(2.183743, 42.84516, 2.957069, 43.26482))
q <- opq(bbox = c(2.5, 4.9, 2.6, 43))


river <- q %>%
  add_osm_feature(key = 'waterway' )%>%
  osmdata_sp()


river$osm_lines[!is.na(river$osm_lines@data$width),]
plot(contour2, col = "red")
plot(river$osm_lines, add = TRUE)
plot(river$osm_lines[!is.na(river$osm_lines@data$width),], col = "green", add =TRUE)

which(!is.na(river$osm_lines@data$width))
