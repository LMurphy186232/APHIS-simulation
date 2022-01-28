###############################################################################
# Functions for preparing graphics from outputs of the APHIS simulation
###############################################################################

#-----------------------------------------------------------------------------#
# Make a set of maps from APHIS simulation output. For each year, a map is
# created showing all trees in black, removed trees in gray, and infested trees
# in red. Survey areas are also shown.
# 
# This assumes that the output is called "trees". Passing a dataset as an 
# argument to this function would make a copy of it, and these datasets are 
# often already extremely large.
# 
# Maps are saved as PNG files. This can easily be changed on the savePlot 
# line below.
# 
# Arguments:
# start_year, end_year: years to start and end mapping
# filename: what you want output files to be named. They will be in the 
# format filename[year].png
#-----------------------------------------------------------------------------#
makeMaps <- function(start_year, end_year, filename) {
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
    
    savePlot(paste0(filename, formatC(year, width=2, flag=0)), type="png")
  }
}

#-----------------------------------------------------------------------------#
# Make a GIF of the output maps produced with the makeMaps function. This 
# looks for all PNG files in a directory and makes them into a single GIF.
#
# Arguments:
# path: path to PNG map files
# filename: what you want your output GIF to be named
#-----------------------------------------------------------------------------#
makeGif <- function(path, filename) {
  library(dplyr)
  library(raster)
  library(magick)
  library(tidyverse)
  
  list.files(path, pattern = "png$", full.names = TRUE) %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps= 2,loop = 1) %>% # animates, can opt for number of loops
  image_write(paste0(filename,".gif")) # write to current dir
}

