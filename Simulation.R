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
# Load input tree file according to desired method; let it be called 
# "tree_dat"
#-----------------------------------------------------------------------------#
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

# Model desired: Currently 11 and 13 supported
model <- 13

# Tree growth slope:
growth_slope <- 0.01
# Tree growth intercept:
growth_intercept <- 0


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
# Model parameters
# These should be in a list called "par". Default values for different 
# regions are included here. 
# 
# Long Island, NY (model 13):
# par <- list(
# dsn   = -2.494357/1000.0, #source pressure (sp)
# delta = 1.056781,         #source pressure (sp)
# beta = 0.1582599,         #source pressure (sp)
# gamma = -0.5183086,       #source term
# alpha = 6.911606,         #source term
# b1    = 98.10091,         #dbh term
# b2    = 1.08E-07,         #dbh term
# c1    = 0.1420186,        #Acer neighbor density term
# c2    = 0.5222963,        #Acer neighbor density term
# e1    = 9248.323,         #landcover term
# e2    = 3.97579,          #landcover term
# mu    = 0.7401554)        #upper limit of infestation prob
# 
# NYC (model 11):
# par <- list(
# dsn	= -5.245118/1000.0, #d1 in source pressure (sp)
# delta1 =	0.0003467984, #d1 in source pressure (sp)
# delta2 = 0.04536906,    #d2 in source pressure (sp)
# dof = 15837.44,         #d2 in source pressure (sp)
# dbf = 0.05532833,       #d2 in source pressure (sp)
# beta = 0.6132349,       #source pressure (sp)
# dirs = c(1,             #sources at 0-30 degrees
#          0.2359833,     #sources at 30-60 degrees 
#          0.006759266,   #sources at 60-90 degrees
#          1.232661e-33,  #sources at 90-120 degrees
#          0.0355077,     #sources at 120-150 degrees
#          0.002069478,   #sources at 150-180 degrees
#          0.001262272,   #sources at 180-210 degrees
#          0.004827997,   #sources at 210-240 degrees
#          0.002,         #sources at 240-270 degrees
#          0.02251909,    #sources at 270-300 degrees
#          0.7987717,     #sources at 300-330 degrees
#          0.001849527),  #sources at 330-360 degrees
# gamma = -0.3628943,     #source term
# alpha = 0.001007975,    #source term
# b1 = 30.87094,          #dbh term
# b2 = -0.7562917,        #dbh term
# c1 = 823.9778,          #Acer neighbor density term
# c2 = 9.975651,          #Acer neighbor density term
# e1 = 9946.723,          #landcover term
# e2 = 3.889605e-05,      #landcover term
# mu = 0.4206378)         #upper limit of infestation prob
# 
# Worcester, MA (model 11)
# par <- list(
# dsn = -11.59742,
# delta1 = 0.08356635,
# delta2 = 1.35E-16,
# dof = 1760.66,
# dbf = 0.1692018,
# beta = 1.04E-09,
# dirs = c(5.51E-69,
#          1.19E-57,
#          5.10E-55,
#          5.04E-49,
#          1.36E-24,
#          4.07E-26,
#          5.41E-19,
#          8.66E-30,
#          2.21E-14,
#          2.00E-17,
#          7.11E-64,
#          1.01E-40),
# gamma = -0.04396946,
# alpha = 0.001,
# b1 = 56.82505,
# b2 = -1.529048,
# c1 = 161.3135,
# c2 = 1.689175,
# e1 = 1628.088,
# e2 = 19.99451,
# mu = 0.8712209
# )
#-----------------------------------------------------------------------------#
# 


#-----------------------------------------------------------------------------#
# Run the model
#-----------------------------------------------------------------------------#
setwd(source_directory)
source("Sim_master.R")

setwd(source_directory)
source("Sim_graphics.R")
makeMaps(start_year, end_year, paste0(output_directory, "/", output_root))

