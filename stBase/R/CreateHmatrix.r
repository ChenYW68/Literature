CreateHmatrix <- function(grid = NULL, site,
                           method = c("indictor", "INLA", "Gaussian"),
                           threshold = 1e-1, factor = 1,
                           distance = T, cs = 0.5,
                           normal.constant = 1e3,
                           distance.method = "geodetic:km")
{
  if(is.null(grid))
  {
    cat("The parameter grid should be null!")
  }
  grid_coords <- grid$grid.coords
  ###########################################################################
  grid_coords <- grid_coords %>% setorderv("ID", 1)
  N.BAUs = nrow(grid_coords)
  if(distance.method == "geodetic:km"){
    BAUs.Dist <- fields::rdist((grid_coords[, c( "LON_X","LAT_Y")]))/normal.constant
    Hdist <- fields::rdist(site[, c("LON_X", "LAT_Y")],
                   grid_coords[, c("LON_X", "LAT_Y")])/normal.constant
  }else{
    BAUs.Dist <- fields::rdist((grid_coords[, c( "LON","LAT")])) 
    Hdist <- fields::rdist(site[, c("LON", "LAT")], grid_coords[, c("LON", "LAT")])
  }
  
  ###########################################################################
  #                                 create Q
  ###########################################################################
  if(distance)
  {
    Adj.Mat = -grid$adjacent.matrix*BAUs.Dist
  }else{
    Adj.Mat = -grid$adjacent.matrix
  }
  for(i in 1:N.BAUs)
  {
    Adj.Mat[i, i] = -sum(Adj.Mat[i, ])
  }
  colnames(Adj.Mat) <- row.names(Adj.Mat) <- grid$grid.coords$ID
  ###########################################################################
  site <- site %>% setorder(SITEID)
  SiteId <- site$SITEID
  ###########################################################################
  #                                 Create H
  ###########################################################################
  if(method == "indictor")
  {
    # Temp1 <- Temp2 <- grid_coords
    # dist <- iDist(grid_coords[, c( "LON_X","LAT_Y")],
    #               grid_coords[, c( "LON_X","LAT_Y")])/1e3
    # h0 = factor
    # H <- matrix(h0, N.BAUs, N.BAUs)
    # for(s in 1:N.BAUs)
    # {
    #   # if(s %% 1000 ==0)
    #   # {cat("s = ", s, "\n")}
    #   j.index1 <- which(dist[s, ] > threshold)
    #   # j.index2 <- which(dist[s, ] <= threshold)
    #   H[s, j.index1] <- 0
    #   # a <- dist[s, j.index2]
    #   # dist[s, s] <- min(a[a>0])/2
    #   # H[s, j.index2] = (1/dist[s, j.index2])/sum(1/dist[s, j.index2])
    # }
    
    
    ###########################################################################
    #                                create  Hs
    ###########################################################################
    Hs <- matrix(factor, nrow = length(SiteId), ncol = N.BAUs)
    # threshold <- 1e5
    # if(distance.method == "geodetic:km"){
    # dist <- iDist(site[, c("LON_X", "LAT_Y")], 
    #               grid_coords[, c("LON_X", "LAT_Y")])/normal.constant
    # }else{
    #   dist <- iDist(site[, c("LON", "LAT")], 
    #                 grid_coords[, c("LON", "LAT")])/normal.constant
    # }
    if(!is.null(cs)){
      threshold <- max(Hdist)*cs
    }
    for(s in 1:length(SiteId))
    {
      j.index1 <- which(Hdist[s, ] > threshold)
      j.index2 <- which(Hdist[s, ] <= threshold)
      #set 0, if distance > threshold
      Hs[s, j.index1] <- 0
      # a <- dist[s, j.index2]
      # dist[s, which(dist[s, ] == 0)] <- min(a[a>0])/2
      # Hs[s, j.index2] = (1/dist[s, j.index2])/sum(1/dist[s, j.index2])
    }
  }else if(method == "INLA"){
    grid.coords <- as.matrix(grid_coords[, c( "LON", "LAT")])
    # H <- as.matrix(inla.spde.make.A(grid$mesh, loc = grid.coords[, 1:2]))
    sampel.coords = as.matrix(site[, c("LON", "LAT")])
    Hs <- as.matrix(inla.spde.make.A(grid$mesh, 
                                     loc = sampel.coords[, 1:2]))
  }else{
    # H = exp(-BAUs.Dist^2/(cs*max(BAUs.Dist)))
    # # Hdist <- iDist(site[, c("LON", "LAT")], grid_coords[, c("LON", "LAT")])
    # Hs = exp(-Hdist^2/(cs*max(Hdist)))
    # 
    # H = Wendland(BAUs.Dist
    #              , theta = max(BAUs.Dist)*cs
    #              , dimension = 1, k =1)
    Hs = Wendland(Hdist
                  , theta = max(Hdist)*cs
                  , dimension = 1, k =1)
  }
  # Hs = grid$Hs
  rownames(Hs) <- SiteId
  
  ###########################################################################
  Data_Str <- list(#H = H
    # ,
    Hs = Hs
    , G = Adj.Mat
    , BAUs.Dist = BAUs.Dist
    , Hs.Dist = Hdist
    , N.BAUs = N.BAUs
    , threshold = threshold)
  ###########################################################################
  return(Data_Str)
}
