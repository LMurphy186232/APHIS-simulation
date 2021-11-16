###############################################################################
# Long Island spread simulation
# 
# This script simulates ALB spread for Long Island based on model fitting done
# by Lora Murphy and Charles Canham on behalf of APHIS. These model results are
# from the model designated "model 13". 
# 
# Inputs:
# A data frame of trees (all should be Acers), with the following fields: 
# (Any others included will be ignored)
# -- infested: integer, 0 being uninfested, 1 being infested. The dataset 
#    should have 1s for the trees that constitute the initial outbreak.
# -- x.ft, y.ft: location coordinates. Point of origin doesn't matter, this
#    will only be used to calculate distances between points. The linear 
#    unit for these should be in feet.
# -- mean_noforestdist: the mean distance to the forest edges in each of the
#    8 cardinal and intercardinal directions, in m. If the tree isn't in 
#    a forest landcover, this is 0.
# -- dbh: tree's DBH in cm
###############################################################################

library(sp)
library(rgeos)

#set.seed(111)
working_directory <- "C:/users/lora/documents/Projects/APHIS2/NY/Simulation"

#-----------------------------------------------------------------------------#
# Run settings:
#-----------------------------------------------------------------------------#
# Year to start simulation at - this will count from here. It doesn't matter
# what it is as long as it's an integer.
start_year <- 1

# How long to run simulation
sim_length <- 10

#-----------------------------------------------------------------------------#
# Biology settings:
#-----------------------------------------------------------------------------#
# Lag period between oviposition and emergence, in years
lag <- 1

dsn   <- -2.494357/1000.0
delta <- 1.056781
betas <- c(0.1582599, 0.1005197, 0.005741494, 0.03213881, 0.006482716)
gamma <- -0.5183086
alpha <- 6.911606
b1    <- 98.10091
b2    <- 1.08E-07
c1    <- 0.1420186
c2    <- 0.5222963
e1    <- 9248.323
e2    <- 3.97579
mu    <- 0.7401554

# Distance to look for neighbor Acers - 30 m, but in feet here
acer_dist <- 30 * 3.28084

#-----------------------------------------------------------------------------#
# Management settings:
#-----------------------------------------------------------------------------#
# Probability that a post-emergence tree in a unit being surveyed is detected
# and removed
prob_surv_post_detect <- 0.9

# Probability that a pre-emergence tree in a unit being surveyed is detected
# and removed
prob_surv_prem_detect <- 0.7

# Probability that a post-emergence tree in a unit NOT being surveyed is
# detected. This is also what will be used to find the initial tree in an
# infestation
prob_nosurv_post_detect <- 0.05

# Probability that a pre-emergence tree NOT in a survey unit will be detected.
# This is also what will be used to find the initial tree in an infestation.
prob_nosurv_prem_detect <- 0

# Budget. The max number of trees that the budget will allow to be removed 
# per year. This is the simplistic beginning, assuming all trees cost the same.
# This could be elaborated into a cost function and total budget if warranted.
max_trees_removed_per_year <- 10000000
#-----------------------------------------------------------------------------#













###############################################################################
#-----------------------------------------------------------------------------#
# START WORK
#-----------------------------------------------------------------------------#
###############################################################################

# Cell size of grid to calculate source pressure. In model 13 this is 100 m
cellsize <- 100 * 3.28084 # convert to feet

#-----------------------------------------------------------------------------#
# Load and prepare data
#-----------------------------------------------------------------------------#
setwd(working_directory)
tree_dat <- readRDS("LI_trees_simcell100.rds")
tree_dat$uid <- 1:nrow(tree_dat)

#-----------------------------------------------------------------------------#
# We need to create the hexagonal grid that sp is based on.
# Create a rectangle spatial polygon bounding all points
#-----------------------------------------------------------------------------#
xcol = which(names(tree_dat) == "x.ft")
ycol = which(names(tree_dat) == "y.ft")
trees <- SpatialPointsDataFrame(tree_dat[,c(xcol, ycol)], tree_dat)

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
trees$cell <- cell_id_for_tree

#-----------------------------------------------------------------------------#
# Filter and package cells data
#-----------------------------------------------------------------------------#
# Remove hex polys with no trees
cell_centers <- cell_centers[unique(x),]
hex_polys <- hex_polys[unique(x)]

#save(cell_centers, hex_polys, file = "GIS data.Rdata")

cells <- data.frame(ID = names(hex_polys), 
                    x = cell_centers@coords[,1], 
                    y = cell_centers@coords[,2])
cells <- cells[cells$ID %in% trees$cell,]
rm(hex_polys, cell_centers, hex_pts, maxX, maxY, minX, minY, coords,
   x, xcol, ycol, poly_list, sp_poly, cell_id_for_tree)

# Clear the infestations and pick trees randomly to be the new infested
trees$infested <- 0
trees$infested[sample(1:nrow(trees), 20)] <- 1

# Some new fields we'll want
trees$year_infested <- NA
trees$year_infested[trees$infested == 1] <- start_year - 1
trees$year_detected <- NA
trees$being_surveyed <- FALSE


maxdist = 5280

# This list will hold the survey areas each time step
surveys_done = list()
surveys_done[[10]] <- NA
#-----------------------------------------------------------------------------#

# This is useful for debugging but a terrible idea otherwise. If this is 
# false, the number of nearby Acer neighbors will not be counted each 
# timestep. (It will always be counted in timestep 1.)
update_acers <- T


###############################################################################
#-----------------------------------------------------------------------------#
# Simulation loop
#-----------------------------------------------------------------------------#
###############################################################################
removed_trees <- NULL
end_year <- (start_year-1) + sim_length
pb=winProgressBar(min=start_year, max=end_year)
yearcount <- 1
for (year in start_year:end_year) {

  if (update_acers || year == start_year) {
  #---------------------------------------------------------------------------#
  # Update the number of Acer neighbors in 30 meters of each tree
  #---------------------------------------------------------------------------#
  # Divide up the area into chunks, to process faster. 
  # Add a buffer to be absolutely sure we get all trees.
  xchunks <- c(seq(from = min(trees$x.ft), to = max(trees$x.ft), by = 500), max(trees$x.ft))
  ychunks <- c(seq(from = min(trees$y.ft), to = max(trees$y.ft), by = 500), max(trees$y.ft))
  
  trees$acer_in_30m  <- NA
  
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
        trees$acer_in_30m [ind] <- chunk_trees$acer_in_30m
      }
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
  if (nrow(sources) > 0) {
    cells$sp <- 0
    for (i in 1:nrow(cells)) {
      dv <- sqrt((cells$x[i] - sources$x.ft)^2 + 
                 (cells$y[i] - sources$y.ft)^2)
      x <- which(dv <= maxdist)
      if (length(x) > 0) {
        dd1 <- exp((dsn)*dv[x]^delta)
        disp <- betas[1] * dd1
        if (any(is.nan(disp))) stop("NOPE")
        cells$sp[i] <- sum(disp, na.rm = T)
      }
    }
    ind <- match(trees$cell, cells$ID)
    source_term <- ifelse(cells$sp[ind] == 0, 1, 
                          1/(1 + (cells$sp[ind]/alpha)^gamma))
    dbh_term    <- 1/(1 + (trees$dbh/b1)^b2)
    den_term    <- 1/(1 + (trees$acer_in_30m/c1)^c2)
    dis_term    <- 1/(1 + (trees$mean_noforestdist/e1)^e2)
    prob <- mu * source_term * dbh_term * den_term * dis_term  
    
    # Make sure trees not near a source get no risk
    no_risk <- which(cells$sp == 0)
    if (length(no_risk) > 0) {
      prob[trees$cell %in% cells$ID[no_risk]] <- 0
    }
    
    #-----------------------------------------------------------------------#
    # Choose trees to become newly infested based on probability
    #-----------------------------------------------------------------------#
    new_infested <- which(runif(n = nrow(trees)) < prob)
    trees$year_infested[new_infested] <- year
    trees$infested[new_infested] <- 1
    
  }
  #-------------------------------------------------------------------------#
  
  
  
  
  
  #-------------------------------------------------------------------------#
  # Management: what areas are being surveyed? Look at last year's 
  # detections, and survey a mile and a half outside of each
  #-------------------------------------------------------------------------#
  surveys_done[[yearcount]] <- NA
  trees$being_surveyed <- FALSE
  if (yearcount > 1) {
    detected_last_year <- which(removed_trees$year_removed == year-1)
    
    if (length(detected_last_year) > 0) {
      
      # Create the shapefile of survey regions based on the locations of
      # the detected trees
      xcol = which(names(removed_trees) == "x.ft")
      ycol = which(names(removed_trees) == "y.ft")
      points <- SpatialPoints(removed_trees[detected_last_year,c(xcol, ycol)])
      if (is.na(points)) stop("BAD")
      if (length(points) == 0) stop("BAD")
      surveys_done[[yearcount]] <- gBuffer(points, width=7290)
    }
  }
  # Identify which trees are being surveyed this year
  if (!is.null(surveys_done[[yearcount]]) && 
      !is.na(surveys_done[[yearcount]])) {
    x <- over(trees, surveys_done[[yearcount]])
    trees$being_surveyed <- !is.na(x)
  }
  saveRDS(surveys_done, "surveys_done.RDS")
  #-------------------------------------------------------------------------#
  
  
  
  #-------------------------------------------------------------------------#
  # Management: prepare to track the number of trees removed, in case we
  # have a budget
  #-------------------------------------------------------------------------#
  trees_left_in_budget <- max_trees_removed_per_year
  
  
  
  #-------------------------------------------------------------------------#
  # Management: figure out which post-emergence trees in surveyed areas
  # get detected. These will be the first to be removed if budget is limited
  #-------------------------------------------------------------------------#
  maybe_to_remove <- !is.na(trees$year_infested)      &
    year - trees$year_infested > lag &
    trees$being_surveyed
  
  if (sum(maybe_to_remove) > 0) {
    to_remove <- maybe_to_remove & 
      runif(n=nrow(trees)) < prob_surv_post_detect
    
    # Check to make sure we're within budget. If we can't remove everything,
    # pick which ones to remove, random selection
    if (sum(to_remove) > trees_left_in_budget) {
      x <- which(to_remove)
      x <- sample(x, trees_left_in_budget)
      to_remove <- rep(FALSE, length(to_remove))
      to_remove[x] <- TRUE
    }
    
    if (sum(to_remove) > 0) {
      this_year_removed_trees <- trees[which(to_remove),]
      this_year_removed_trees$year_removed <- year
      if (is.null(removed_trees)) {
        removed_trees <- this_year_removed_trees
      } else {
        removed_trees <- rbind(removed_trees, this_year_removed_trees)
      }
      
      # Track what we've removed so far, in case there's a budget
      trees_left_in_budget <- trees_left_in_budget -
                              nrow(this_year_removed_trees)
      rm(this_year_removed_trees)
      
      trees <- trees[!to_remove,]
    }
  }
  #-------------------------------------------------------------------------#
  
  
  #-------------------------------------------------------------------------#
  # Management: figure out which pre-emergence trees in surveyed areas
  # get detected. I'm assuming these will be second choice for removal
  #-------------------------------------------------------------------------#
  maybe_to_remove <- !is.na(trees$year_infested)       &
                     year - trees$year_infested <= lag &
                     trees$being_surveyed
  
  if (sum(maybe_to_remove) > 0 && trees_left_in_budget > 0) {
    
    to_remove <- maybe_to_remove & 
                 runif(n=nrow(trees)) < prob_surv_prem_detect
    
    # Stay within budget. If we can't remove everything, pick which ones to 
    # remove, random selection
    if (sum(to_remove) > trees_left_in_budget) {
      x <- which(to_remove)
      x <- sample(x, trees_left_in_budget)
      to_remove <- rep(FALSE, length(to_remove))
      to_remove[x] <- TRUE
    }
    if (sum(to_remove) > 0) {
      this_year_removed_trees <- trees[which(to_remove),]
      this_year_removed_trees$year_removed <- year
      if (is.null(removed_trees)) {
        removed_trees <- this_year_removed_trees
      } else {
        removed_trees <- rbind(removed_trees, this_year_removed_trees)
      }
      
      # Track what we've removed so far, in case there's a budget
      trees_left_in_budget <- trees_left_in_budget -
        nrow(this_year_removed_trees)
      rm(this_year_removed_trees)
      
      trees <- trees[!to_remove,]
    }
  }
  #-------------------------------------------------------------------------#
  
  
  
  #-------------------------------------------------------------------------#
  # Management: figure out which post-emergence trees NOT in surveyed areas
  # get detected
  #-------------------------------------------------------------------------#
  maybe_to_remove <- !is.na(trees$year_infested)      &
    year - trees$year_infested > lag &
    !trees$being_surveyed
  
  if (sum(maybe_to_remove) > 0 && trees_left_in_budget > 0) {
    
    to_remove <- maybe_to_remove & 
      runif(n=nrow(trees)) < prob_nosurv_post_detect
    
    # Stay within budget. If we can't remove everything, pick which ones to 
    # remove, random selection
    if (sum(to_remove) > trees_left_in_budget) {
      x <- which(to_remove)
      x <- sample(x, trees_left_in_budget)
      to_remove <- rep(FALSE, length(to_remove))
      to_remove[x] <- TRUE
    }
    if (sum(to_remove) > 0) {
      this_year_removed_trees <- trees[which(to_remove),]
      this_year_removed_trees$year_removed <- year
      if (is.null(removed_trees)) {
        removed_trees <- this_year_removed_trees
      } else {
        removed_trees <- rbind(removed_trees, this_year_removed_trees)
      }
      
      # Track what we've removed so far, in case there's a budget
      trees_left_in_budget <- trees_left_in_budget -
        nrow(this_year_removed_trees)
      rm(this_year_removed_trees)
      
      trees <- trees[!to_remove,]
    }
  }
  #-------------------------------------------------------------------------#
  
  
  #-------------------------------------------------------------------------#
  # Management: figure out which pre-emergence trees NOT in surveyed areas
  # get detected
  #-------------------------------------------------------------------------#
  maybe_to_remove <- !is.na(trees$year_infested)       &
                     year - trees$year_infested <= lag &
                     !trees$being_surveyed
  
  if (sum(maybe_to_remove) > 0 && trees_left_in_budget > 0) {
    
    to_remove <- maybe_to_remove & 
                 runif(n=nrow(trees)) < prob_nosurv_prem_detect
    
    # Stay within budget. If we can't remove everything, pick which ones to 
    # remove, random selection
    if (sum(to_remove) > trees_left_in_budget) {
      x <- which(to_remove)
      x <- sample(x, trees_left_in_budget)
      to_remove <- rep(FALSE, length(to_remove))
      to_remove[x] <- TRUE
    }
    if (sum(to_remove) > 0) {
      this_year_removed_trees <- trees[which(to_remove),]
      this_year_removed_trees$year_removed <- year
      if (is.null(removed_trees)) {
        removed_trees <- this_year_removed_trees
      } else {
        removed_trees <- rbind(removed_trees, this_year_removed_trees)
      }
      
      # Track what we've removed so far, in case there's a budget
      trees_left_in_budget <- trees_left_in_budget -
        nrow(this_year_removed_trees)
      rm(this_year_removed_trees)
      
      trees <- trees[!to_remove,]
    }
  }
  #-------------------------------------------------------------------------#
  
  
  
  
  
  
  
  yearcount <- yearcount + 1
  setWinProgressBar(pb, title = paste0("year ", year, " complete"), value = year)
}
close(pb)




save(trees, removed_trees, file="Results_surveys.Rdata")

#load("Results_surveys.Rdata")
#surveys_done <- readRDS("surveys_done.RDS")
windows()
for (year in start_year:end_year) {
  plot(trees$x.ft, trees$y.ft, pch=20, col="black",
       main = paste("Year", year))
  x = which(trees$infested == 1 & trees$year_infested <= year)
  if (length(x) > 0) {
    points(trees$x.ft[x], trees$y.ft[x], pch=20, col="red")
  }
  x = which(removed_trees$year_removed  >= year & 
            removed_trees$year_infested <= year)
  if (length(x) > 0) {
    points(trees$x.ft[x], trees$y.ft[x], pch=20, col="red")
  }
  x = which(removed_trees$year_removed <= year)
  if (length(x) > 0) {
    points(trees$x.ft[x], trees$y.ft[x], pch=20, col="gray")
  }
  if (!is.null(surveys_done[[year]]) && !is.na(surveys_done[[year]])) 
    plot(surveys_done[[year]],add=T, col="blue")
  savePlot(paste0("run3year", year), type="png")
}

###############################################################################
testing <- F
if (testing) {
  sources <- subset(trees, infested == 1)
  # Do a probability loop by hand to check the calculations
  cells$sp <- 0
  for (i in 1:nrow(cells)) {
    dv <- sqrt((cells$x.ft[i] - sources$x.ft)^2 + 
               (cells$y.ft[i] - sources$y.ft)^2)
    x <- which(dv <= maxdist)
    if (length(x) > 0) {
      dd1 <- exp((dsn)*dv[x]^delta)
      disp <- betas[1] * dd1
      if (any(is.nan(disp))) stop("NOPE")
      cells$sp[i] <- sum(disp, na.rm = T)
    }
  }
  ind <- match(trees$cell, cells$ID)
  source_term <- 1/(1 + (cells$sp[ind]/alpha)^gamma)
  dbh_term    <- 1/(1 + (trees$dbh/b1)^b2)
  den_term    <- 1/(1 + (trees$acer_in_30m/c1)^c2)
  dis_term    <- 1/(1 + (trees$mean_noforestdist/e1)^e2)
  prob <- mu * source_term * dbh_term * den_term * dis_term     
  
  
  
  
  trees <- trees[1:40,]
  
  ind <- (match(trees$cell, cells$ID)-1)
  prob2 <- sim_step(cells$x.ft, cells$y.ft, 
                    sources$x.ft, sources$y.ft, 
                    rep(1, nrow(sources)), 
                    trees$dbh,
                    trees$acer_in_30m,
                    trees$mean_noforestdist,
                    betas, 
                    ind,
                    dsn, delta,
                    gamma, alpha,
                    b1, b2, c1, c2,
                    e1, e2, mu, maxdist)
}
