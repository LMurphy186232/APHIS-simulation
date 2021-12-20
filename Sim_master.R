###############################################################################
# Runs the simulation.
###############################################################################
# Do a bit of double-checking...
if (model != "LI") stop(paste0("\"", model, "\" is an invalid model choice."))

setwd(source_directory)
source("LI_Model.R")


trees <- NULL
do_setup()
rm(tree_dat)

# This list will hold the survey areas each time step
surveys_done = list()
surveys_done[[sim_length]] <- NA
#-----------------------------------------------------------------------------#
setwd(output_directory)
###############################################################################
#-----------------------------------------------------------------------------#
# Simulation loop
#-----------------------------------------------------------------------------#
###############################################################################
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
    detected_last_year <- which(trees$year_removed == year-1)
    
    # If we have just a ton, limit it to 500 to avoid killing gBuffer
    if (length(detected_last_year) > 500) {
      detected_last_year <- detected_last_year[sample(1:length(detected_last_year), 500)]
    }
    
    if (length(detected_last_year) > 0) {
      
      # Create the shapefile of survey regions based on the locations of
      # the detected trees
      xcol = which(names(trees) == "x")
      ycol = which(names(trees) == "y")
      points <- SpatialPoints(trees[detected_last_year,c(xcol, ycol)])
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
  #saveRDS(surveys_done, "surveys_done.RDS")
  #-------------------------------------------------------------------------#
  
  
  
  #-------------------------------------------------------------------------#
  # Management: prepare to track the number of trees removed, in case we
  # have a budget
  #-------------------------------------------------------------------------#
  if (year_to_begin_management >= year) {
    trees_left_in_budget <- max_trees_removed_per_year
    
    
    
    #-------------------------------------------------------------------------#
    # Management: figure out which post-emergence trees in surveyed areas
    # get detected. These will be the first to be removed if budget is limited
    #-------------------------------------------------------------------------#
    maybe_to_remove <- !is.na(trees$year_infested)      &
                       is.na(trees$year_removed)        &
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
        trees$year_removed[which(to_remove)] <- year
        
        # Track what we've removed so far, in case there's a budget
        trees_left_in_budget <- trees_left_in_budget - sum(to_remove)
      }
    }
    #-------------------------------------------------------------------------#
    
    
    #-------------------------------------------------------------------------#
    # Management: figure out which pre-emergence trees in surveyed areas
    # get detected. I'm assuming these will be second choice for removal
    #-------------------------------------------------------------------------#
    maybe_to_remove <- !is.na(trees$year_infested)       &
                       is.na(trees$year_removed)         &
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
        trees$year_removed[which(to_remove)] <- year
        
        # Track what we've removed so far, in case there's a budget
        trees_left_in_budget <- trees_left_in_budget - sum(to_remove)
      }
    }
    #-------------------------------------------------------------------------#
    
    
    
    #-------------------------------------------------------------------------#
    # Management: figure out which post-emergence trees NOT in surveyed areas
    # get detected
    #-------------------------------------------------------------------------#
    maybe_to_remove <- !is.na(trees$year_infested)      &
                       is.na(trees$year_removed)        &
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
        trees$year_removed[which(to_remove)] <- year
        
        # Track what we've removed so far, in case there's a budget
        trees_left_in_budget <- trees_left_in_budget - sum(to_remove)
      }
    }
    #-------------------------------------------------------------------------#
    
    
    #-------------------------------------------------------------------------#
    # Management: figure out which pre-emergence trees NOT in surveyed areas
    # get detected
    #-------------------------------------------------------------------------#
    maybe_to_remove <- !is.na(trees$year_infested)       &
                       is.na(trees$year_removed)         &
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
        trees$year_removed[which(to_remove)] <- year
        
        # Track what we've removed so far, in case there's a budget
        trees_left_in_budget <- trees_left_in_budget - sum(to_remove)
      }
    }
    #-------------------------------------------------------------------------#

  }

  yearcount <- yearcount + 1
  setWinProgressBar(pb, title = paste0("year ", year, " complete"), value = year)
}
close(pb)


trees$being_surveyed <- NULL
saveRDS(trees, file=paste0(output_root, "_trees.Rdata"))
saveRDS(surveys_done, file=paste0(output_root, "_surveys.Rdata"))