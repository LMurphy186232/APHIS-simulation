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
    if (!is.null(surveys_done[[year]]) && mode(surveys_done[[year]]) ==  "S4") 
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
  library(magick)
  
  ff <- list.files(path, pattern = "png$", full.names = TRUE)
  ir <- image_read(ff)
  ij <- image_join(ir)
  im <- image_animate(ij, fps= 2,loop = 1)
  image_write(im, paste0(filename,".gif"))
}


#-----------------------------------------------------------------------------#
# Make line graphs for output.
#
# Arguments:
# path: path to PNG map files
# filename: what you want your output GIF to be named
#-----------------------------------------------------------------------------#
makeLineGraphs <- function(start_year, end_year, filename) {
  library(ggplot2)
  library(ggpubr)
  
  total_infested <- rep(0, length(start_year:end_year))
  newly_infested <- rep(0, length(start_year:end_year))
  total_removed  <- rep(0, length(start_year:end_year))
  newly_removed  <- rep(0, length(start_year:end_year))
  
  # Collect the statistics year-by-year
  for (year in start_year:end_year) {
    total_infested[year] <- length(which(trees$year_infested <= year))
    newly_infested[year] <- length(which(trees$year_infested == year))
    total_removed [year] <- length(which(trees$year_removed <= year))
    newly_removed [year] <- length(which(trees$year_removed == year))  
  }
  
  # Graph this nicely
  df <- data.frame(Year = start_year:end_year,
                   total_infested = total_infested,
                   newly_infested = newly_infested,
                   total_removed  = total_removed,
                   newly_removed  = newly_removed)
  
  p1 <- ggplot(data = df, aes(x = Year, y = total_infested)) +
    geom_line(lwd=2) +
    labs(y = "Cumulative total infested") +
    theme_bw()
  
  p2 <- ggplot(data = df, aes(x = Year, y = newly_infested)) +
    geom_line(lwd=2) +
    labs(y = "Newly infested") +
    theme_bw()
  
  p3 <- ggplot(data = df, aes(x = Year, y = total_removed)) +
    geom_line(lwd=2) +
    labs(y = "Cumulative total removed") +
    theme_bw()
  
  p4 <- ggplot(data = df, aes(x = Year, y = newly_removed)) +
    geom_line(lwd=2) +
    labs(y = "Newly removed") +
    theme_bw()
  
  ggpubr::ggarrange(p1, p2, p3, p4, nrow=2, ncol=2)
}
