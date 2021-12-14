###############################################################################
# Runs the simulation.
###############################################################################
setwd(source_directory)
source("LI_Model.R")

#tree_dat <- readRDS("SimTrees.rds")
trees <- NULL
do_setup()
rm(tree_dat)

# This list will hold the survey areas each time step
surveys_done = list()
surveys_done[[sim_length]] <- NA
#-----------------------------------------------------------------------------#

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
  
  
  #-----------------------------------------------------------------------#
  # Choose trees to become newly infested based on probability
  #-----------------------------------------------------------------------#
  prob <- get_tree_risk()
  new_infested <- which(runif(n = nrow(trees)) < prob & 
                          trees$infested == 0)
  trees$year_infested[new_infested] <- year
  trees$infested[new_infested] <- 1
  
  
  #-------------------------------------------------------------------------#
  
  
  
  
  
  #-------------------------------------------------------------------------#
  # Management: what areas are being surveyed? Look at last year's 
  # detections, and survey a mile and a half outside of each
  #-------------------------------------------------------------------------#
  surveys_done[[yearcount]] <- NA
  trees$being_surveyed <- FALSE
  if (yearcount > 1) {
    detected_last_year <- which(removed_trees$year_removed == year-1)
    
    # If we have just a ton, limit it to 500 to avoid killing gBuffer
    if (length(detected_last_year) > 500) {
      detected_last_year <- detected_last_year[sample(1:length(detected_last_year), 500)]
    }
    
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


#-----------------------------------------------------------------------------#
# Graphing
#-----------------------------------------------------------------------------#
#load("Results_surveys.Rdata")
#surveys_done <- readRDS("surveys_done.RDS")
windows()
for (year in start_year:end_year) {
  
  # Collect some statistics
  total_infested <- length(which(trees$year_infested <= year)) +
    length(which(removed_trees$year_infested <= year))
  newly_infested <- length(which(trees$year_infested == year)) +
    length(which(removed_trees$year_infested == year))
  total_removed <- 0
  newly_removed <- 0
  if (!is.null(removed_trees)) {
    total_removed <- length(which(removed_trees$year_removed <= year))
    newly_removed <- length(which(removed_trees$year_removed == year))  
  }
  
  plot(trees$x.ft, trees$y.ft, pch=20, col="black",
       main = paste("Year", year,
                    "\nTotal infested:", total_infested, 
                    "Newly infested:", newly_infested,
                    "\nTotal removed:" , total_removed,
                    "Newly removed:" , newly_removed),
       xlab="", ylab="", xaxt="n", yaxt="n")
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
    plot(surveys_done[[year]],add=T, border="blue")
  savePlot(paste0("run1year", year), type="png")
}
