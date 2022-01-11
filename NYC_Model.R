###############################################################################
# This file contains the code and parameter values for calculating tree risk
# according to the NYC model, designated "model 11". The model is
# a function of distance and direction to sources, DBH, number of Acer 
# neighbors within 30 m, the mean distance to the forest edges in each of the 
# 8 cardinal and intercardinal directions (in m; if the tree isn't in a forest
# landcover, this is 0).
#
# Source pressure is rasterized, calculated for 100m hexagonal grid cells across
# the landscape. A tree will take the source pressure value of the grid cell
# that it's in. 
#
# The original model fit a series of beta values, for relative year of 
# detection. This informed detection probability more than infectiousness.
# When simulating spread from active sources, we use the single largest beta
# value.
#
# model11 <- function(dsn, betas, gamma, mu, alpha, dirs, delta1, 
#                     delta2, dsn, dof, dbf, a, b, c) {
#  d1 <- delta[1] * exp((dsn/1000)*(distances))
#  d2 <- delta[2] * exp(-0.5*(log( distances/dof)/dbf)^2)
#  disp <- betas[beta_ind] * dirs[dir_ind] * (d1 + d2)
#  sp <- rowSums(disp) # This is on a cell basis
#
#  source_term <-  1/(1 + (sp/alpha)^gamma)
#  dbh_term <- 1/(1 + (targets$dbh/b1)^b2)                                         
#  den_term <- 1/(1 + (targets$acer_in_30m/c1)^c2)                                  
#  dis_term <- 1/(1 + (targets$mean_noforestdist/e1)^e2)                           
#  prob <- mu * source_term * dbh_term * den_term * dis_term            
#  prob
#}
###############################################################################

nyc.par <- list(
  dsn	= -5.245118,
  delta1 =	0.0003467984,
  delta2 = 0.04536906,
  dof = 15837.44,
  dbf = 0.05532833,
  beta = 0.6132349,
  dirs = c(1,
           0.2359833,
           0.006759266,
           1.232661e-33,
           0.0355077,
           0.002069478,
           0.001262272,
           0.004827997,
           0.002,
           0.02251909,
           0.7987717,
           0.001849527),
  gamma = -0.3628943,
  alpha = 0.001007975,
  b1 = 30.87094,
  b2 = -0.7562917,
  c1 = 823.9778,
  c2 = 9.975651,
  e1 = 9946.723, 
  e2 = 3.889605e-05,
  mu = 0.4206378)

if (linear_unit == "feet") {
  
  # Distance to look for neighbor Acers - 30 m, but in feet here
  acer_dist <- 30 * 3.28084
  
  # Max distance for source trees, in feet
  maxdist = 5280
  
  # Cell size of grid to calculate source pressure. In model 13 this is 100 m
  cellsize <- 100 * 3.28084 # convert to feet
  
} else if (linear_unit == "meters") {
  
  # Distance to look for neighbor Acers - 30 m
  acer_dist <- 30
  
  # Max distance for source trees, feet to meters
  maxdist = 5280 * 0.3048
  
  # Cell size of grid to calculate source pressure. In model 13 this is 100 m
  cellsize <- 100
  
} else {
  stop(paste0("\"", linear_unit, "\" is an invalid choice for linear units.\""))
}



if (!any(names(tree_dat) == "x")) stop("Missing field \"x\" in tree dataset.")
if (!any(names(tree_dat) == "y")) stop("Missing field \"y\" in tree dataset.")
if (!any(names(tree_dat) == "mean_noforestdist")) stop("Missing field \"mean_noforestdist\" in tree dataset.")





#-----------------------------------------------------------------------------#
# This will do setup for the model. It sets up the hexagonal grid system.
#-----------------------------------------------------------------------------#
do_setup <- function() {
  
  
  tree_dat$uid <- 1:nrow(tree_dat)
  
  #---------------------------------------------------------------------------#
  # We need to create the hexagonal grid that sp is based on.
  # Create a rectangle spatial polygon bounding all points
  #---------------------------------------------------------------------------#
  xcol = which(names(tree_dat) == "x")
  ycol = which(names(tree_dat) == "y")
  trees <<- SpatialPointsDataFrame(tree_dat[,c(xcol, ycol)], tree_dat)
  
  minX <- min(trees$x) - (cellsize*2)
  maxX <- max(trees$x) + (cellsize*2)
  minY <- min(trees$y) - (cellsize*2)
  maxY <- max(trees$y) + (cellsize*2)
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
  if (length(which(names(trees) == "infested")) == 0) {
    trees$infested <<- 0 
  }
  if (sum(trees$infested) == 0) {
    trees$infested <<- 0
    trees$infested[sample(1:nrow(trees), num_initial_outbreak)] <<- 1
  }
  
  # Some new fields we'll want
  if (length(which(names(trees) == "year_infested")) == 0) {
    trees$year_infested <<- NA
    trees$year_infested[trees$infested == 1] <<- start_year - 1
  }
  if (length(which(names(trees) == "year_removed")) == 0) {
    trees$year_removed <<- NA
  }
  if (length(which(names(trees) == "being_surveyed")) == 0) {
    trees$being_surveyed <<- FALSE
  }
}
#-----------------------------------------------------------------------------#




#-----------------------------------------------------------------------------#
# Calculate tree risk according to model 11
#-----------------------------------------------------------------------------#
get_tree_risk <- function() {
  #---------------------------------------------------------------------------#
  # Update the number of Acer neighbors in 30 meters of each tree
  #---------------------------------------------------------------------------#
  # Divide up the area into chunks, to process faster. 
  # Add a buffer to be absolutely sure we get all trees.
  xchunks <- c(seq(from = min(trees$x), to = max(trees$x), by = 500), max(trees$x))
  ychunks <- c(seq(from = min(trees$y), to = max(trees$y), by = 500), max(trees$y))
  
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
      chunk_trees <- subset(trees, x >= xmin & 
                              x <= xmax & 
                              y >= ymin &
                              y <= ymax)
      
      # Get the potential neighbors within the buffer
      buffer_trees <- subset(trees, x >= xminbuff & 
                               x <= xmaxbuff & 
                               y >= yminbuff &
                               y <= ymaxbuff)
      
      
      if (nrow(chunk_trees) > 0) {
        # Process each tree
        for (k in 1:nrow(chunk_trees)) {
          
          # Get distance to all neighbors
          dists <- sqrt((chunk_trees$x[k] - buffer_trees$x)^2 + 
                        (chunk_trees$y[k] - buffer_trees$y)^2)
          
          # Add up who's in the radius - subtracting 1 for self
          chunk_trees$acer_in_30m [k] = sum(dists <= acer_dist &
                                            is.na(buffer_trees$year_removed)) - 1
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
  chunk <- pi*2 / 12
  # Get sources: omit those that are still pre-emergence
  sources <- subset(trees, infested == 1 & 
                    (year - year_infested) > lag)
  prob <- rep(0, nrow(trees))
  if (nrow(sources) > 0) {
    cells$sp <- 0
    for (i in 1:nrow(cells)) {
      
      # Calculate distance to sources
      dv <- sqrt((cells$x[i] - sources$x)^2 + 
                 (cells$y[i] - sources$y)^2)
      
      # Select which sources are within the effective range
      x <- which(dv <= maxdist)
      if (length(x) > 0) {
        
        # Determine direction to the sources
        X <- sources$x[x] - cells$x[i]
        Y <- sources$y[x] - cells$y[i]
        azi <- rep(NA, length(x))
        
        #Calculate the azimuth - correct for quadrant of "to" relative to "from"
        #(counting clockwise from upper right)
        ## first quadrant
        azi <- ifelse (Y > 0 & X >= 0, atan(X/Y), azi)
        
        # second quadrant
        azi <- ifelse (X > 0 & Y <= 0, (pi / 2.0) + atan((-1.0 * Y)/X), azi)
        
        # third quadrant
        azi <- ifelse (X <= 0 & Y < 0, pi + atan((-1.0 * X) / (-1.0 * Y)), azi)
        
        # fourth quadrant
        azi <- ifelse (X < 0 & Y >= 0, ((1.5 * pi) + atan(Y / (-1.0 * X))), azi)
        
        dir1 <- (trunc((azi/chunk), 0)) + 1
        
        dd1 <- nyc.par$delta1 * exp((nyc.par$dsn/1000)*dv[x])
        dd2 <- nyc.par$delta2 * exp(-0.5*(log( dv[x]/nyc.par$dof)/nyc.par$dbf)^2)
        disp <- nyc.par$beta * nyc.par$dirs[dir1] * (dd1 + dd2)
        if (any(is.nan(disp))) stop("NOPE")
        cells$sp[i] <- sum(disp, na.rm = T)
      }
    }
    ind <- match(trees$cell, cells$ID)
    source_term <- ifelse(cells$sp[ind] == 0, 1, 
                          1/(1 + (cells$sp[ind]/nyc.par$alpha)^nyc.par$gamma))
    dbh_term    <- 1/(1 + (trees$dbh/nyc.par$b1)^nyc.par$b2)
    den_term    <- 1/(1 + (acer_in_30m/nyc.par$c1)^nyc.par$c2)
    dis_term    <- 1/(1 + (trees$mean_noforestdist/nyc.par$e1)^nyc.par$e2)
    prob <- nyc.par$mu * source_term * dbh_term * den_term * dis_term  
    
    # Make sure trees not near a source get no risk
    no_risk <- which(cells$sp == 0)
    if (length(no_risk) > 0) {
      prob[trees$cell %in% cells$ID[no_risk]] <- 0
    }
  }
  prob
}