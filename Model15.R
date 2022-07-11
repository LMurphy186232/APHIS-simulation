###############################################################################
# This file contains the code and parameter values for calculating tree risk
# according to the SC model, designated "model 15". The model is
# a function of distance and direction to sources, DBH, and number of Acer 
# neighbors within 30 m.
#
# Source pressure is NOT rasterized. 
#
# The original model fit a series of beta values, for relative year of 
# detection. This informed detection probability more than infectiousness.
# When simulating spread from active sources, we use the single largest beta
# value.
#
# model15 <- function(dsn, delta, betas, alpha, gamma, 
#                     b1, b2, c1, c2, mu) {
#  # 1 dispersal kernel 
#  d1 <- exp((dsn/1000)*(distances^delta))
#  
#  disp <- betas[beta_ind] * d1
#  sp <- rowSums(disp)
# 
#  source_term <-  1/(1 + (sp/alpha)^gamma)
#  dbh_term <- 1/(1 + (targets$dbh/b1)^b2)                                         
#  den_term <- 1/(1 + (targets$acer_in_30m/c1)^c2)                                  
#  prob <- mu * source_term * dbh_term * den_term
#  prob
#}
###############################################################################

bakdir <- getwd()
setwd(source_directory)
sourceCpp("countacers.cc")
setwd(bakdir)

# Check to make sure our "par" has everything in it that we need
if (any(!c("dsn", "delta", "beta", "gamma", "alpha", "b1", "b2",
           "c1", "c2", "mu") %in% names(par))) stop("Par missing required parameters")

if (linear_unit == "feet") {
  
  # Distance to look for neighbor Acers - 30 m, but in feet here
  acer_dist <- 30 * 3.28084
  
  # Max distance for source trees, in feet
  maxdist = 5280
  
} else if (linear_unit == "meters") {
  
  # Distance to look for neighbor Acers - 30 m
  acer_dist <- 30
  
  # Max distance for source trees, feet to meters
  maxdist = 5280 * 0.3048
  
} else {
  stop(paste0("\"", linear_unit, "\" is an invalid choice for linear units.\""))
}



if (!any(names(tree_dat) == "x")) stop("Missing field \"x\" in tree dataset.")
if (!any(names(tree_dat) == "y")) stop("Missing field \"y\" in tree dataset.")



#-----------------------------------------------------------------------------#
# This will do setup for the model. It sets up the hexagonal grid system.
#-----------------------------------------------------------------------------#
do_setup <- function() {
  
  
  tree_dat$uid <- 1:nrow(tree_dat)
  xcol = which(names(tree_dat) == "x")
  ycol = which(names(tree_dat) == "y")
  trees <<- SpatialPointsDataFrame(tree_dat[,c(xcol, ycol)], tree_dat)
  
  #createHexGridOnTrees()
  
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
# Calculate tree risk according to model 15
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
  # Get sources: omit those that are still pre-emergence
  sources <- subset(trees, infested == 1 & 
                   (year - year_infested) > lag)
  prob <- rep(0, nrow(trees))
  sp   <- rep(0, nrow(trees))
  if (nrow(sources) > 0) {
    for (i in 1:nrow(trees)) {
      dv <- sqrt((trees$x[i] - sources$x)^2 + 
                 (trees$y[i] - sources$y)^2)
      x <- which(dv <= maxdist & dv > 0)
      if (length(x) > 0) {
        dd1 <- exp((par$dsn)*dv[x]^par$delta)
        disp <- par$beta * dd1
        if (any(is.nan(disp))) stop("NOPE")
        sp[i] <- sum(disp, na.rm = T)
      }
    }
    source_term <- ifelse(sp == 0, 1, 
                          1/(1 + (sp/par$alpha)^par$gamma))
    dbh_term    <- 1/(1 + (trees$dbh/par$b1)^par$b2)
    den_term    <- 1/(1 + (trees$acer_in_30m/par$c1)^par$c2)
    prob <- par$mu * source_term * dbh_term * den_term
    
    # Make sure trees not near a source get no risk
    no_risk <- which(sp == 0)
    if (length(no_risk) > 0) {
      prob[no_risk] <- 0
    }
  }
  prob
}