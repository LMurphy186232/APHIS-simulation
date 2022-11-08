###############################################################################
# Runs the simulation. IF YOU WANT TO RUN THE SIMULATION, THIS FILE ISN'T THE
# ONE YOU USE: Open "Simulation.R" instead.
###############################################################################
# Do a bit of double-checking...
if (!model %in% c(11, 13, 15)) stop(paste0("\"", model, 
                                     "\" is an invalid model choice."))

# Set random seed if desired
if (!is.null(random_seed)) set.seed(random_seed)

setwd(source_directory)
if (model == 13) {
  source("Model13.R")
} 
if (model == 11) {
  source("Model11.R")
}
if (model == 15) {
  source("Model15.R")
}
source("Helper_functions.R")

trees <- NULL
do_setup()
rm(tree_dat)

# Determine whether the surveys are prescribed or not
surveys_are_prescribed <- F
if (exists("prescribed_survey")) {
  surveys_are_prescribed <- !is.null(prescribed_survey)
}

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
pb=winProgressBar(min=(start_year-1), max=end_year)
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
  if (!surveys_are_prescribed) {
    
    #-------------------------------------------------------------------------#
    # Surveys are not prescribed by the user. Make decisions about whether 
    # there will be surveys and where, based on last year's detections.
    #-------------------------------------------------------------------------#
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
        proj4string(points) <- proj4string(trees)
        if (length(points) == 0) stop("BAD")
        ss <- st_buffer(st_as_sf(points), dist=survey_radius)
        surveys_done[[yearcount]] <- as(st_geometry(ss), "Spatial")
        
      }
    }
    # Identify which trees are being surveyed this year
    if (!is.null(surveys_done[[yearcount]]) && 
        mode(surveys_done[[yearcount]]) == "S4") {
      x <- over(trees, surveys_done[[yearcount]])
      
      #x <- over(surveys_done[[yearcount]], trees)
      trees$being_surveyed <- !is.na(x)
    }
  } else {
    #-------------------------------------------------------------------------#
    # Surveys are prescribed. Choose trees based on which units, if any, are
    # scheduled to be surveyed.
    #-------------------------------------------------------------------------#
    x <- which(prescribed_survey$year == year)
    if (length(x) > 0) {
      trees$being_surveyed <- trees$unit %in% prescribed_survey$unit[x]
      bstart = 1
      # Make an approximate geometry of the survey area(s)
      for (xi in 1:length(x)) {
        # Create the shapefile of survey regions based on the locations of
        # the detected trees
        xcol = which(names(trees) == "x")
        ycol = which(names(trees) == "y")
        points <- SpatialPoints(trees[trees$being_surveyed,c(xcol, ycol)])
        
        # Make a buffer of tree locations
        buf1 <-st_buffer(st_as_sf(points), dist=100)
        # Dissolve them together
        buf2 <- st_union(buf1)
        surveys_done[[yearcount]] <- as(st_geometry(buf2), "Spatial")
      }
    }
  }
  #saveRDS(surveys_done, "surveys_done.RDS")
  #---------------------------------------------------------------------------#
  
  
  
  #---------------------------------------------------------------------------#
  # Management: prepare to track the number of trees removed, in case we
  # have a budget
  #---------------------------------------------------------------------------#
  if (year_to_begin_management <= year) {
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
  
  #---------------------------------------------------------------------------#
  # Grow trees
  #---------------------------------------------------------------------------#
  trees$dbh <- ifelse(is.na(trees$year_removed), 
                      trees$dbh + trees$dbh * growth_slope + growth_intercept, 
                      trees$dbh)
  
  yearcount <- yearcount + 1
  setWinProgressBar(pb, title = paste0("year ", year, " complete"), value = year)
}
close(pb)


trees$being_surveyed <- NULL
saveRDS(trees, file=paste0(output_root, "_trees.rds"))
saveRDS(surveys_done, file=paste0(output_root, "_surveys.rds"))
