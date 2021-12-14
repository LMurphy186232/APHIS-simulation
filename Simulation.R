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

working_directory <- "C:/users/lora/documents/Projects/APHIS2/NY/Simulation"
source_directory <- "C:/users/lora/documents/Projects/APHIS2/Simulation/APHIS-Simulation"

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
#-----------------------------------------------------------------------------#







#-----------------------------------------------------------------------------#
# Load and prepare data
#-----------------------------------------------------------------------------#
setwd(working_directory)
tree_dat <- readRDS("LI_trees_simcell100.rds")

#-----------------------------------------------------------------------------#
# Run the model
#-----------------------------------------------------------------------------#
setwd(source_directory)
source("Sim_master.R")