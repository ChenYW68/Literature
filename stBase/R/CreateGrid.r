CreateGrid <- function(data, max.edge = c(0.6, 0.8),
                        offset = c(0.7, 0.3),
                        cutoff = 2, col = "black", size = 1)
{
  sampel.coords = as.matrix(data[, c("LON", "LAT")])

  global.coords =as.matrix(Map_BTH[, 1:2])

  mesh <- inla.mesh.2d(
    #boundary = boundary
     loc.domain = global.coords
    , loc = sampel.coords[, 1:2]
    , max.edge = max.edge # max. edge length
    , offset = offset# the inner and outer part
    , cutoff = cutoff
  )
  
  data <- as.data.frame(data)
  coordinates(data) = ~ LON + LAT

  p <- ggplot() + gg(mesh) + #geom_sf(col = "black") +
    ggtitle(paste("Vertices: ", mesh$n)) +
    geom_polygon(data = Map_BTH,
                 aes(x = long,y = lat, group = group),
                 colour = 'black',
                 fill = NA) + 
    # coord_sf(datum = st_crs(5880)) +
     gg(data, col = col, size = size) +
    theme_bw() + #coord_fixed() +
    labs(x =  "longitude", y = "latitude") +
    theme(axis.text = element_text(size = 8, colour = "black")
          , axis.title = element_text(size = 14, colour = "black")
          # , legend.title = element_text(size = 12, colour = "black")
          # , legend.text = element_text(size = 8, colour = "black")
    )
  grid.coords <- mesh[["loc"]][, 1:2]
  N.BAUs <- nrow(grid.coords)
  N.BAUs
  ###############################################################
  ###############################################################
  BAUs.Dist <- as.matrix(dist(coordinates(grid.coords)))
  ###############################################################
  ###############################################################
  grid.coords <- as.data.frame(grid.coords)
  colnames(grid.coords) <- c("LON", "LAT")

  grid.coords$ID = 1:nrow(grid.coords)

  grid.coords <- spCoords.transform(grid.coords, method = 2)

#   LCC <- "+proj=lcc +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0
# +y_0=0 +ellps=GRS80 +datum=NAD83 +units=km +no_defs"
#   Coord <- grid.coords[, c("LON", "LAT")] %>% setDF() %>% as.matrix()
#   # transforming latitude and longitude using the Lambert conformal projection
#   sp.coord <- SpatialPoints(coords=Coord, proj4string=CRS("+proj=longlat +ellps=WGS84"))
#   
#   sp.coord.trans <- spTransform(sp.coord,  CRS(LCC))
#   
#   Co <- as.data.frame(sp.coord.trans) %>% setnames(c("LON", "LAT"), c("LON_X", "LAT_Y"))
#   
#   grid.coords$LON_X <- Co$LON_X
#   grid.coords$LAT_Y <- Co$LAT_Y
  
  setcolorder(grid.coords, c(3, 1:2, 4:5))
  adjacent.matrix <- as.matrix(mesh[["graph"]][["vv"]])
  grid <- list(grid.coords = grid.coords,
               adjacent.matrix = adjacent.matrix,
               mesh = mesh,
               plot.grid = p)
  return(grid)
}
