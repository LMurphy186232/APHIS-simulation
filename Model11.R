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

bakdir <- getwd()
setwd(source_directory)
sourceCpp("countacers.cc")
setwd(bakdir)

# Check to make sure our "par" has everything in it that we need
if (any(!c("dsn", "delta1", "delta2", "dof", "dbf", "beta", "dirs",
          "gamma", "alpha", "b1", "b2", "c1", "c2", "e1",
          "e2", "mu") %in% names(par))) stop("Par missing required parameters")

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
  xcol = which(names(tree_dat) == "x")
  ycol = which(names(tree_dat) == "y")
  trees <<- SpatialPointsDataFrame(tree_dat[,c(xcol, ycol)], tree_dat)
  
  createHexGridOnTrees()
  
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

  # If this is the first year - have to count them all
  if (year == start_year) {
    trees$acer_in_30m  <<- rep(NA, nrow(trees))
    live <- which(is.na(trees$year_removed))
    trees$acer_in_30m[live] <<- count_acer_neighbors(trees$x[live], trees$y[live], acer_dist)
  } else {
    # If nobody died last year, there's no updating to do
    if (any(trees$year_removed == (year-1), na.rm=T)) {
      dead <- which(trees$year_removed == (year-1))
      live <- which(is.na(trees$year_removed))
      vv <- update_acer_neighbors(trees$x[live], trees$y[live],
                                  trees$x[dead], trees$y[dead],
                                  acer_dist)
      trees$acer_in_30m[live] <<- trees$acer_in_30m[live] - vv
          
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
      dv <- sqrt((cells$X[i] - sources$x)^2 + 
                 (cells$Y[i] - sources$y)^2)
      x <- which(dv <= maxdist)
      if (length(x) > 0) {
        
        # Determine direction to the sources
        X <- sources$x[x] - cells$X[i]
        Y <- sources$y[x] - cells$Y[i]
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
        
        dd1 <- par$delta1 * exp((par$dsn)*dv[x])
        dd2 <- par$delta2 * exp(-0.5*(log( dv[x]/par$dof)/par$dbf)^2)
        disp <- par$beta * par$dirs[dir1] * (dd1 + dd2)
        if (any(is.nan(disp))) stop("NOPE")
        cells$sp[i] <- sum(disp, na.rm = T)
      }
    }
    ind <- match(trees$cell, cells$ID)
    source_term <- ifelse(cells$sp[ind] == 0, 1, 
                          1/(1 + (cells$sp[ind]/par$alpha)^par$gamma))
    dbh_term    <- 1/(1 + (trees$dbh/par$b1)^par$b2)
    den_term    <- 1/(1 + (trees$acer_in_30m/par$c1)^par$c2)
    dis_term    <- 1/(1 + (trees$mean_noforestdist/par$e1)^par$e2)
    prob <- par$mu * source_term * dbh_term * den_term * dis_term  
    
    # Make sure trees not near a source get no risk
    no_risk <- which(cells$sp == 0)
    if (length(no_risk) > 0) {
      prob[trees$cell %in% cells$ID[no_risk]] <- 0
    }
  }
  prob
}