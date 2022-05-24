###############################################################################
# Use this file as a basis for simulations. You can save a copy for each
# simulation you would like to run, to preserve settings.
###############################################################################
library(sp)
library(sf)
library(Rcpp)


# Location of the source R script files
source_directory <- "~/APHIS-Simulation"

# Where to put output
output_directory <- ""
# Output name. This will be used as the root for various files that may be
# produced
output_root <- "testrun"

#-----------------------------------------------------------------------------#
# Load input tree file according to desired method; call it "tree_dat". Make
# it has no NAs in DBH, x, and y. These should all be trees of the genus Acer.
#-----------------------------------------------------------------------------#
tree_dat <- readRDS("LI_trees.rds")
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
# Run settings:
#-----------------------------------------------------------------------------#
# Year to start simulation at. It doesn't matter what it is as long as it's an 
# integer. This can be a counter or an actual year if you prefer.
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

# Model desired: Currently 11 and 13 supported
model <- 13

# Settings for linear growth model for trees. 
# DBH(t+1) = DBH(t) + (slope * DBH(t) + intercept)
# Tree growth slope:
growth_slope <- 0.01

# Tree growth intercept:
growth_intercept <- 0


#-----------------------------------------------------------------------------#
# Management settings:
#-----------------------------------------------------------------------------#
# Probability that a post-emergence tree in a unit being surveyed is detected
# and removed, 0 to 1
prob_surv_post_detect <- 0.9

# Probability that a pre-emergence tree in a unit being surveyed is detected
# and removed, 0 to 1
prob_surv_prem_detect <- 0.7

# Probability that a post-emergence tree in a unit NOT being surveyed is
# detected, 0 to 1. 
prob_nosurv_post_detect <- 0.05

# Probability that a pre-emergence tree in a unit NOT being surveyed will be 
# detected.
prob_nosurv_prem_detect <- 0

# Budget. The max number of trees that the budget will allow to be removed 
# per year. 
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
# Model parameters
# These should be in a list called "par". Several sets of values can be found
# in the file "Parameters.R". Example:
# 
# Long Island, NY (best model is model 13):
par <- list(
 dsn   = -2.494357/1000.0, #source pressure (sp)
 delta = 1.056781,         #source pressure (sp)
 beta = 0.1582599,         #source pressure (sp)
 gamma = -0.5183086,       #source term
 alpha = 6.911606,         #source term
 b1    = 98.10091,         #dbh term
 b2    = 1.08E-07,         #dbh term
 c1    = 0.1420186,        #Acer neighbor density term
 c2    = 0.5222963,        #Acer neighbor density term
 e1    = 9248.323,         #landcover term
 e2    = 3.97579,          #landcover term
 mu    = 0.7401554)        #upper limit of infestation prob 
#-----------------------------------------------------------------------------#
 


#-----------------------------------------------------------------------------#
# Run the model
#-----------------------------------------------------------------------------#
setwd(source_directory)
source("Sim_master.R")

#-----------------------------------------------------------------------------#
# Output graphics. Some basic functions are included with the script below.
#-----------------------------------------------------------------------------#
setwd(source_directory)
source("Sim_graphics.R")
makeMaps(start_year, end_year, paste0(output_directory, "/", output_root, "map"))
makeLineGraphs(start_year, end_year, paste0(output_directory, "/", output_root, "lines"))

