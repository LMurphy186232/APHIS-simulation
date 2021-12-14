###############################################################################
# This file contains the code and parameter values for calculating tree risk
# according to the NY - Long Island model, designated "model 13". The model is
# a function of distance to sources, DBH, number of Acer neighbors within 30 m,
# the mean distance to the forest edges in each of the 8 cardinal and 
# intercardinal directions (in m; if the tree isn't in a forest landcover, this
# is 0).
#
# Source pressure is rasterized, calculated for 100m hexagonal grid cells across
# the landscape. A tree will take the source pressure value of the grid cell
# that it's in.
#
# model13 <- function(dsn, delta, delta2, betas, alpha, gamma, 
#                     b1, b2, c1, c2, e1, e2, mu) {
#  # 2 dispersal kernels 
#  d1 <- delta1 * exp((dsn/1000)*(distances^delta))
#  
#  disp <- betas[beta_ind] * d1
#  sp <- rowSums(disp)
# 
#  source_term <-  1/(1 + (sp/alpha)^gamma)
#  dbh_term <- 1/(1 + (targets$dbh/b1)^b2)                                         
#  den_term <- 1/(1 + (targets$acer_in_30m/c1)^c2)                                  
#  dis_term <- 1/(1 + (targets$mean_noforestdist/e1)^e2)                           
#  prob <- mu * source_term * dbh_term * den_term * dis_term            
#  prob
#}
###############################################################################

ny.par <- list(
  dsn   = -2.494357/1000.0,
  delta = 1.056781,
  betas = c(0.1582599, 0.1005197, 0.005741494, 0.03213881, 0.006482716),
  gamma = -0.5183086,
  alpha = 6.911606,
  b1    = 98.10091,
  b2    = 1.08E-07,
  c1    = 0.1420186,
  c2    = 0.5222963,
  e1    = 9248.323,
  e2    = 3.97579,
  mu    = 0.7401554)

# Distance to look for neighbor Acers - 30 m, but in feet here
acer_dist <- 30 * 3.28084

# Max distance for source trees, in feet
maxdist = 5280





#-----------------------------------------------------------------------------#
# This will do setup for the model. It sets up the hexagonal grid system.
#-----------------------------------------------------------------------------#
do_setup <- function() {
  # Cell size of grid to calculate source pressure. In model 13 this is 100 m
  cellsize <- 100 * 3.28084 # convert to feet
  
  tree_dat$uid <- 1:nrow(tree_dat)
  
  #---------------------------------------------------------------------------#
  # We need to create the hexagonal grid that sp is based on.
  # Create a rectangle spatial polygon bounding all points
  #---------------------------------------------------------------------------#
  xcol = which(names(tree_dat) == "x.ft")
  ycol = which(names(tree_dat) == "y.ft")
  trees <<- SpatialPointsDataFrame(tree_dat[,c(xcol, ycol)], tree_dat)
  
  minX <- min(trees$x.ft) - (cellsize*2)
  maxX <- max(trees$x.ft) + (cellsize*2)
  minY <- min(trees$y.ft) - (cellsize*2)
  maxY <- max(trees$y.ft) + (cellsize*2)
  coords <- rbind(c(minX, minY), 
                  c(minX, maxY), 
                  c(maxX, maxY),
                  c(maxX, minY),
                  c(minX, minY))
  
  # Make a Polygon
  p <- Polygon(coords, hole = F)
  
  # Make a Polygons - separate class
  poly_list <- list(Polygons(list(p), ID = 1))
  
  # Make our SpatialPolygons object
  sp_poly <- SpatialPolygons(poly_list)                                                             
  
  #-----------------------------------------------------------------------------#
  # Hex a grid over the top. 
  #-----------------------------------------------------------------------------#
  hex_pts <-spsample(sp_poly,type="hexagonal",cellsize=cellsize)
  hex_polys <- HexPoints2SpatialPolygons(hex_pts)
  
  
  # Get the centroids of the hex polys
  cell_centers <- gCentroid(hex_polys, byid = T)
  
  #-----------------------------------------------------------------------------#
  # Extract the polygon for each tree
  #-----------------------------------------------------------------------------#
  # x will be a vector of indexes of hex_polys that each tree falls into
  x <- over(trees, hex_polys)
  cell_id_for_tree <- names(hex_polys)[x]
  trees$cell <<- cell_id_for_tree
  
  #-----------------------------------------------------------------------------#
  # Filter and package cells data
  #-----------------------------------------------------------------------------#
  # Remove hex polys with no trees
  cell_centers <- cell_centers[unique(x),]
  hex_polys <- hex_polys[unique(x)]
  
  #save(cell_centers, hex_polys, file = "GIS data.Rdata")
  
  cells <<- data.frame(ID = names(hex_polys), 
                      x = cell_centers@coords[,1], 
                      y = cell_centers@coords[,2])
  cells <<- cells[cells$ID %in% trees$cell,]
  rm(hex_polys, cell_centers, hex_pts, maxX, maxY, minX, minY, coords,
     x, xcol, ycol, poly_list, sp_poly, cell_id_for_tree)
  
  # Clear the infestations and pick trees randomly to be the new infested
  if (sum(trees$infested) == 0) {
    trees$infested <<- 0
    trees$infested[sample(1:nrow(trees), 1)] <<- 1
  }
  
  # Some new fields we'll want
  if (length(which(names(trees) == "year_infested")) == 0) {
    trees$year_infested <<- NA
    trees$year_infested[trees$infested == 1] <<- start_year - 1
  }
  if (length(which(names(trees) == "year_detected")) == 0) {
    trees$year_detected <<- NA
  }
  if (length(which(names(trees) == "being_surveyed")) == 0) {
    trees$being_surveyed <<- FALSE
  }
}
#-----------------------------------------------------------------------------#




#-----------------------------------------------------------------------------#
# Calculate tree risk according to model 13
#-----------------------------------------------------------------------------#
get_tree_risk <- function() {
  #---------------------------------------------------------------------------#
  # Update the number of Acer neighbors in 30 meters of each tree
  #---------------------------------------------------------------------------#
  # Divide up the area into chunks, to process faster. 
  # Add a buffer to be absolutely sure we get all trees.
  xchunks <- c(seq(from = min(trees$x.ft), to = max(trees$x.ft), by = 500), max(trees$x.ft))
  ychunks <- c(seq(from = min(trees$y.ft), to = max(trees$y.ft), by = 500), max(trees$y.ft))
  
  acer_in_30m  <- rep(NA, nrow(trees))
  
  count <- 0
  for (i in 2:length(xchunks)) {
    
    # Two sets of trees: core trees in the chunk, and potential neighbors 
    # in the buffer
    xmin <- xchunks[i-1]
    xmax <- xchunks[i]  
    
    xminbuff <- xmin - 120
    xmaxbuff <- xmax + 120
    
    for (j in 2:length(ychunks)) {
      ymin <- ychunks[j-1]
      ymax <- ychunks[j]  
      
      yminbuff <- ymin - 120
      ymaxbuff <- ymax + 120
      
      # Get the trees within the core area area
      chunk_trees <- subset(trees, x.ft >= xmin & 
                              x.ft <= xmax & 
                              y.ft >= ymin &
                              y.ft <= ymax)
      
      # Get the potential neighbors within the buffer
      buffer_trees <- subset(trees, x.ft >= xminbuff & 
                               x.ft <= xmaxbuff & 
                               y.ft >= yminbuff &
                               y.ft <= ymaxbuff)
      
      
      if (nrow(chunk_trees) > 0) {
        # Process each tree
        for (k in 1:nrow(chunk_trees)) {
          
          # Get distance to all neighbors
          dists <- sqrt((chunk_trees$x.ft[k] - buffer_trees$x.ft)^2 + 
                          (chunk_trees$y.ft[k] - buffer_trees$y.ft)^2)
          
          # Add up who's in the radius - subtracting 1 for self
          chunk_trees$acer_in_30m [k] = sum(dists <= acer_dist) - 1
        }
        
        # Add back to master dataset
        ind <- match(chunk_trees$uid, trees$uid)
        acer_in_30m [ind] <- chunk_trees$acer_in_30m
      }
    }
  }
  #---------------------------------------------------------------------------#
  
  
  
  #---------------------------------------------------------------------------#
  # Get the tree risk
  #---------------------------------------------------------------------------#
  # Get sources: omit those that are still pre-emergence
  sources <- subset(trees, infested == 1 & 
                    (year - year_infested) > lag)
  bind <- pmin(5, year - sources$year_infested - lag)
  prob <- rep(0, nrow(trees))
  if (nrow(sources) > 0) {
    cells$sp <- 0
    for (i in 1:nrow(cells)) {
      dv <- sqrt((cells$x[i] - sources$x.ft)^2 + 
                 (cells$y[i] - sources$y.ft)^2)
      x <- which(dv <= maxdist)
      if (length(x) > 0) {
        #stop("NN")
        dd1 <- exp((ny.par$dsn)*dv[x]^ny.par$delta)
        disp <- ny.par$betas[bind[x]] * dd1
        if (any(is.nan(disp))) stop("NOPE")
        cells$sp[i] <- sum(disp, na.rm = T)
      }
    }
    ind <- match(trees$cell, cells$ID)
    source_term <- ifelse(cells$sp[ind] == 0, 1, 
                          1/(1 + (cells$sp[ind]/ny.par$alpha)^ny.par$gamma))
    dbh_term    <- 1/(1 + (trees$dbh/ny.par$b1)^ny.par$b2)
    den_term    <- 1/(1 + (acer_in_30m/ny.par$c1)^ny.par$c2)
    dis_term    <- 1/(1 + (trees$mean_noforestdist/ny.par$e1)^ny.par$e2)
    prob <- ny.par$mu * source_term * dbh_term * den_term * dis_term  
    
    # Make sure trees not near a source get no risk
    no_risk <- which(cells$sp == 0)
    if (length(no_risk) > 0) {
      prob[trees$cell %in% cells$ID[no_risk]] <- 0
    }
  }
  prob
}