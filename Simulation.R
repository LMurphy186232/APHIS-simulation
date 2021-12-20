###############################################################################
# ALB spread simulation
# 
# This script simulates ALB spread.
###############################################################################

library(sp)
library(rgeos)

#-----------------------------------------------------------------------------#
# Input data: a dataframe of trees that are potential ALB victims. All are 
# assumed to be Acers. The dataframe should have the following fields, with
# these exact case-sensitive names (there can be other fields, they'll be 
# ignored):
# -- "x", "y": location coordinates. Point of origin doesn't matter, this
#    will only be used to calculate distances between points. SPECIFY THE 
#    LINEAR UNITS BELOW (feet, meters).
# -- "mean_noforestdist": the mean distance to the forest edges in each of the
#    8 cardinal and intercardinal directions, in m. If the tree isn't in 
#    a forest landcover, this is 0.
# -- "dbh": tree's DBH in cm
# 
# OPTIONAL:
# -- "infested": integer, 0 being uninfested, 1 being infested. If this field
#    is included and there is at least 1 tree with a 1 status, the trees 
#    marked as 1 will be the seed of the outbreak. If either this field is not
#    included or all the values are 0, a select number of trees will be
#    randomly assigned to be the seed of the outbreak (number settable below).
#    
# Output: a dataframe, like trees, with additional fields for dates of 
# infestation and removal.
#-----------------------------------------------------------------------------#


# Location of the source R script files
source_directory <- "C:/users/lora/documents/Projects/APHIS2/Simulation/APHIS-Simulation"

# Where to put output
output_directory <- "C:/users/lora/documents/Projects/APHIS2/NY/Simulation/Test"
# Output name. This will be used as the root for various files that may be
# produced
output_root <- "testrun"

#-----------------------------------------------------------------------------#
# Load input tree file according to desired method; let it be called 
# "tree_dat"
#-----------------------------------------------------------------------------#
data_directory <- "C:/users/lora/documents/Projects/APHIS2/Simulation"
setwd(data_directory)
tree_dat <- readRDS("LI_trees.rds")
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
# Run settings:
#-----------------------------------------------------------------------------#
# Year to start simulation at - this will count from here. It doesn't matter
# what it is as long as it's an integer.
start_year <- 1

# How long, in years, to run simulation
sim_length <- 10


# Linear units for coordinates: one of "feet" or "meters"
linear_unit <- "feet"

#-----------------------------------------------------------------------------#
# Biology settings:
#-----------------------------------------------------------------------------#
# Lag period between oviposition and emergence, in years
lag <- 1

# Model desired: Currently only "LI" supported
model <- "LI"


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
max_trees_removed_per_year <- 100000000

# Year to begin any management. Allows for a spin-up period. Management is not
# guaranteed to start this year - depending on probabilities, the 
# infestation may yet go unnoticed for a while longer.
year_to_begin_management <- 1

# Number of trees to randomly select for infestation to start the outbreak
# in the first year. THIS IS IGNORED IF THERE ARE ANY TREES ALREADY INFESTED
# IN THE INPUT DATASET.
num_initial_outbreak <- 0
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
# Run the model
#-----------------------------------------------------------------------------#
setwd(source_directory)
source("Sim_master.R")


#-----------------------------------------------------------------------------#
# Graphics
#-----------------------------------------------------------------------------#
make_gif <- T
if (make_gif) {
  library(dplyr)
  library(raster)
  library(magick)
  library(tidyverse)
}

#-----------------------------------------------------------------------------#
# Graphing
#-----------------------------------------------------------------------------#
windows()
for (year in start_year:end_year) {
  
  # Collect some statistics
  total_infested <- length(which(trees$year_infested <= year))
  newly_infested <- length(which(trees$year_infested == year))
  total_removed <- 0
  newly_removed <- 0
  total_removed <- length(which(trees$year_removed <= year))
  newly_removed <- length(which(trees$year_removed == year))  
  
  
  plot(trees$x, trees$y, pch=20, col="black",
       main = paste("Year", year,
                    "\nTotal infested:", total_infested, 
                    "Newly infested:", newly_infested,
                    "\nTotal removed:" , total_removed,
                    "Newly removed:" , newly_removed),
       xlab="", ylab="", xaxt="n", yaxt="n")
  x = which(trees$infested == 1 & trees$year_infested <= year)
  if (length(x) > 0) {
    points(trees$x[x], trees$y[x], pch=20, col="red")
  }
  x = which(trees$year_removed  >= year & 
            trees$year_infested <= year)
  if (length(x) > 0) {
    points(trees$x[x], trees$y[x], pch=20, col="pink")
  }
  x = which(trees$year_removed <= year)
  if (length(x) > 0) {
    points(trees$x[x], trees$y[x], pch=20, col="gray")
  }
  if (!is.null(surveys_done[[year]]) && !is.na(surveys_done[[year]])) 
    plot(surveys_done[[year]],add=T, border="blue")
  savePlot(paste0("runyear", formatC(year, width=2, flag=0)), type="png")
}


list.files(pattern = "png$", full.names = TRUE) %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps= 2,loop = 1) %>% # animates, can opt for number of loops
  image_write(paste0(output_root,".gif")) # write to current dir

