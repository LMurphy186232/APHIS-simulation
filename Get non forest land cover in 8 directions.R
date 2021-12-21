###############################################################################
# For each tree, this script finds the nearest non-forest edge in each of 8 
# directions in a landcover raster, and identifies what landcover it is. If the
# tree is not in a forest cover cell, all will be 0s.
# 
# Inputs: 
# -- A SpatialPointsDataFrame of trees, or a dataframe with X and Y coordinates
#    that can be turned into a SpatialPointsDataFrame
# -- A landcover raster. This is assumed to be the NLCD raster. It is highly
#    recommended that you give it a rough clip to the area covered by your 
#    trees. This script will further clip it.
#
# You must have rgdal installed as a package. This has extra requirements 
# compared to most R packages so consult the installation instructions.
###############################################################################

library(rgdal)
library(raster)

# What values represent forest types? Any other value encountered is assumed 
# to be non-forest. 
forest_cover <- 41:43

#-----------------------------------------------------------------------------#
# Load the tree dataset.  
#-----------------------------------------------------------------------------#
# By chance, is this a dataframe that we need to turn into a 
# SpatialPointsDataFrame? Here's how: 
# Load the thing, assume it's called "trees"
# Know the projection the coordinates are in, and how to specify that in R. 
# An example is below. A full explanation on determining projection strings is a
# bit outside the scope of this script, but if you have a GIS object in the same
# projection, you can load it into R and query its string.
# Locate the columns with the X and Y coordinates - example:
xcol <- 4
ycol <- 5
dat <- SpatialPointsDataFrame(coords = trees[,c(xcol, ycol)], data = trees, 
             proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# So now assume we have a SpatialPointsDataFrame called "dat"...
#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
# First off, let's trim down the landcover raster so that it's only a bit 
# bigger than our tree dataset. To be very sure, we'll do a 5 km buffer
#-----------------------------------------------------------------------------#


# Load the landcover raster, replacing filename as appropriate. 
cov <- raster("myraster.tif")

# Make the projection for the trees match the raster
dat <- spTransform(dat, proj4string(cov))

# Get the minimum and maximum tree coordinates, and add a 5 km buffer to get
# the extent to which we will clip the landcover raster
minX <- min(dat@coords[,1]) - 5000
maxX <- max(dat@coords[,1]) + 5000
minY <- min(dat@coords[,2]) - 5000
maxY <- max(dat@coords[,2]) + 5000

# Clip the raster
cov <- crop(cov, extent(minX, maxX, minY, maxY))
rm(minX, maxX, minY, maxY)

# Get the cell size of the raster - assuming it's square
cellsize <- res(cov)[1]
half_cell <- cellsize/2

#-----------------------------------------------------------------------------#


#-----------------------------------------------------------------------------#
# Function for getting the distance along a diagonal line, used
# in the raytrace function below
#
# Arguments:
# x, y - coordinates of origin point
# cell_center - cell coordinates of cell to get distance to
# horiz - whether we last made a horizonal (TRUE) or vertical (FALSE)
#         move to get to this cell
# dir - one of NE, NW, SE, or SW
#-----------------------------------------------------------------------------#
getdistance <- function(x, y, cell_center, horiz, dir) {
  
  # We need to know the distance to the cell we reached. If we last made a
  # horizontal move, we know the X coordinate of the intersecting point.
  # Otherwise, we know the Y coordinate.
  if (!horiz) {
    # We made a vertical move. Find the Y distance we moved
    
    if (dir == "NE" || dir == "NW") {
      # Moving north - get the southern edge of the cell
      celly = cell_center[2] - half_cell
    } else {
      # Moving south - get the northern edge of the cell
      celly = cell_center[2] + half_cell
    }
    
    # Trig - we have a right triangle and we know angle (45 deg)
    # and delta y. Find hypotenuse
    dist = abs(y - celly) / sin(pi/4)
    
  } else {
    
    # We made a horizontal move. Find the X distance we moved
    if (dir == "NE" || dir == "SE") {
      # Moving east - get the western edge of the cell
      cellx = cell_center[1] - half_cell
    } else {
      # Moving west - get the eastern edge of the cell
      cellx = cell_center[1] + half_cell
    }
    
    # Trig - we have a right triangle and we know angle (45 deg)
    # and delta x. Find hypotenuse
    dist = abs(x - cellx) / cos(pi/4)
    
  }
  dist
}


#-----------------------------------------------------------------------------#
# Raytrace function for tracing lines along the diagonal.
# Since the slope of our line is always +- 1, there is a consistent
# pattern to which cells you intersect. You always alterate 
# advancing one cell horizontally and one cell vertically. Where
# you start depends on where in the origin cell the point is.
#
# We'll find the x and y distances to the cell edges - the shortest
# one is the direction we'll move first.
#
# Arguments: x, y = starting point; homecell - cell we're starting in;
# dir = one of NE, NW, SE, SW.
#-----------------------------------------------------------------------------#
raytrace <- function(x, y, homecell, dir) {
  
  # Home cell row and column
  col <- colFromCell(cov, homecell)
  row <- rowFromCell(cov, homecell)
  
  # Home cell center
  cell_center <- xyFromCell(cov, homecell)
  
  #---------------------------------------------------------------------------#
  # Figure out if our first move is horizontal or vertical.
  # To do this, find the X and Y distances between our point
  # and the cell edge in the appropriate direction. The
  # shorter distance determines which way we go.
  #---------------------------------------------------------------------------#
  
  # Find the x and y distances to the cell edge
  if (dir == "NE" || dir == "NW") {
    # Find the y distance to the top of the cell
    ydist = (cell_center[2] + half_cell) - y
  } else {
    # Find the y distance to the bottom of the cell
    ydist = y - (cell_center[2] - half_cell)
  }
  
  if (dir == "NE" || dir == "SE") {
    # Find the x distance to the right of the cell
    xdist = (cell_center[1] + half_cell) - x
  } else {
    # Find the x distance to the left of the cell
    xdist = x - (cell_center[1] - half_cell)
  }
  
  # This is a toggle for which direction to advance next
  horiz = xdist < ydist
  
  
  #---------------------------------------------------------------------------#
  # Set up the increment for rows and columns so we walk the
  # correct direction
  #---------------------------------------------------------------------------#
  if (dir == "NE" || dir == "NW") {
    y.incr = -1
  } else {
    y.incr = 1
  }
  
  if (dir == "NE" || dir == "SE") {
    x.incr = 1
  } else {
    x.incr = -1
  }
  
  if (horiz) {
    col = col + x.incr
  } else {
    row = row + y.incr
  }
  # Flip the toggle to move the other way next time
  horiz = !horiz
  
  #---------------------------------------------------------------------------#
  # Walking loop
  #---------------------------------------------------------------------------#
  firstdist <- NA
  firstcov <- NA
  
  while (validCol(cov, col) && validRow(cov, row) && is.na(firstdist)) {
    
    # Are we out of forest type?
    cell_cover <- extract(cov, cellFromRowCol(cov, row, col))
    if (!(cell_cover %in% forest_cover)) {
      
      firstdist <- getdistance(x, y, xyFromCell(cov, cellFromRowCol(cov, row, col)), !horiz, dir)
      firstcov <- cell_cover
    }
    
    if (horiz) {
      col = col + x.incr
    } else {
      row = row + y.incr
    }
    
    # Flip the toggle to move the other way next time
    horiz = !horiz
  }
  
  list(firstdist = firstdist, firstcov = firstcov)
}
#-----------------------------------------------------------------------------#









#-----------------------------------------------------------------------------#
# Add columns for the results 
#-----------------------------------------------------------------------------#
dat$Nnoforestdist <- NA
dat$Nnoforestcov <- NA

dat$Snoforestdist <- NA
dat$Snoforestcov <- NA

dat$Enoforestdist <- NA
dat$Enoforestcov <- NA

dat$Wnoforestdist <- NA
dat$Wnoforestcov <- NA

dat$NEnoforestdist <- NA
dat$NEnoforestcov <- NA

dat$NWnoforestdist <- NA
dat$NWnoforestcov <- NA

dat$SEnoforestdist <- NA
dat$SEnoforestcov <- NA

dat$SWnoforestdist <- NA
dat$SWnoforestcov <- NA

#-----------------------------------------------------------------------------#
# Find the cover type each tree is in and its home cell
#-----------------------------------------------------------------------------#
dat$tree_cov <- extract(cov, dat, method = "simple")

# Any trees that aren't in forest, set their distances to 0
p <- which(!dat$tree_cov %in% forest_cover)
dat$Nnoforestdist[p] <- 0
dat$Snoforestdist[p] <- 0
dat$Enoforestdist[p] <- 0
dat$Wnoforestdist[p] <- 0
dat$NEnoforestdist[p] <- 0
dat$NWnoforestdist[p]<- 0
dat$SEnoforestdist[p] <- 0
dat$SWnoforestdist[p] <- 0

dat$Nnoforestcov[p]  <- dat$tree_cov[p]
dat$Snoforestcov[p]  <- dat$tree_cov[p]
dat$Enoforestcov[p]  <- dat$tree_cov[p]
dat$Wnoforestcov[p]  <- dat$tree_cov[p]
dat$NEnoforestcov[p] <- dat$tree_cov[p]
dat$NWnoforestcov[p] <- dat$tree_cov[p]
dat$SEnoforestcov[p] <- dat$tree_cov[p]
dat$SWnoforestcov[p] <- dat$tree_cov[p]

# Get the raster cell of each tree
dat$cell <- cellFromXY(cov, dat@coords)
#-----------------------------------------------------------------------------#







#-----------------------------------------------------------------------------#
# Process the cardinal directions cell-by-cell
# Cell numbers start at 1 in the upper left corner, and increase from left to 
# right, and then from top to bottom. The last cell number equals the number of
# cells of the Raster* object.
# 
# Rows are the same. Row 1 is at the top (north). Columns go left-to-right.
#-----------------------------------------------------------------------------#
cells_to_do <- unique(dat$cell[dat$tree_cov %in% forest_cover])
for (cell in cells_to_do) {
  
  # Get the center of the cell
  cell_center <- xyFromCell(cov, cell)
  
  # Get its row and column
  homerow <- rowFromCell(cov, cell)
  homecol <- colFromCell(cov, cell)
  
  # Get all the trees in this cell
  cell_trees <- which(dat$cell == cell)
  
  #---------------------------------------------------------------------------#
  # Walk north until we find our edge
  #---------------------------------------------------------------------------#
  row       <- homerow - 1
  col       <- homecol
  edge_cell <- -1
  
  while(validRow(cov, row) && edge_cell < 0) {
    
    # Are we out of forest cover?
    cell_cover <- extract(cov, cellFromRowCol(cov, row, col))
    if (!(cell_cover %in% forest_cover)) {
      edge_cell <- cellFromRowCol(cov, row, col)
    }
    row <- row - 1
  }
  
  # If we didn't find a cell, throw an error!
  if (edge_cell < 0) stop("Didn't find a north cell!")
  
  # Get all the trees that are in this cell, and have them calculate their
  # distance to the southern edge of the edge cell
  dat$Nnoforestdist[cell_trees] <- xyFromCell(cov, edge_cell)[2] - half_cell - dat@coords[cell_trees,2]
  dat$Nnoforestcov [cell_trees] <- cell_cover
  #---------------------------------------------------------------------------#
  
  
  
  
  
  
  #---------------------------------------------------------------------------#
  # Walk south until we find our edge
  #---------------------------------------------------------------------------#
  row       <- homerow + 1
  col       <- homecol
  edge_cell <- -1
  
  while(validRow(cov, row) && edge_cell < 0) {
    
    # Are we out of forest cover?
    cell_cover <- extract(cov, cellFromRowCol(cov, row, col))
    if (!(cell_cover %in% forest_cover)) {
      edge_cell <- cellFromRowCol(cov, row, col)
    }
    row <- row + 1
  }
  
  # If we didn't find a cell, throw an error!
  if (edge_cell < 0) stop("Didn't find a south cell!")
  
  # Get all the trees that are in this cell, and have them calculate their
  # distance to the northern edge of the edge cell
  dat$Snoforestdist[cell_trees] <- dat@coords[cell_trees,2] - (xyFromCell(cov, edge_cell)[2] + half_cell)
  dat$Snoforestcov [cell_trees] <- cell_cover
  #---------------------------------------------------------------------------#  
    
    
    
  
  
  
  #---------------------------------------------------------------------------#
  # Walk east until we find our edge
  #---------------------------------------------------------------------------#
  row       <- homerow
  col       <- homecol + 1
  edge_cell <- -1
  
  while(validCol(cov, col) && edge_cell < 0) {
    
    # Are we out of forest cover?
    cell_cover <- extract(cov, cellFromRowCol(cov, row, col))
    if (!(cell_cover %in% forest_cover)) {
      edge_cell <- cellFromRowCol(cov, row, col)
    }
    col <- col + 1
  }
  
  # If we didn't find a cell, throw an error!
  if (edge_cell < 0) stop("Didn't find an east cell!")
  
  # Get all the trees that are in this cell, and have them calculate their
  # distance to the western edge of the edge cell
  dat$Enoforestdist[cell_trees] <- xyFromCell(cov, edge_cell)[1] - half_cell - dat@coords[cell_trees,1]
  dat$Enoforestcov [cell_trees] <- cell_cover
  #---------------------------------------------------------------------------#
  
  
  
  
  
  
  #---------------------------------------------------------------------------#
  # Walk west until we find our edge
  #---------------------------------------------------------------------------#
  row       <- homerow
  col       <- homecol - 1
  edge_cell <- -1
  
  while(validCol(cov, col) && edge_cell < 0) {
    
    # Are we out of forest cover?
    cell_cover <- extract(cov, cellFromRowCol(cov, row, col))
    if (!(cell_cover %in% forest_cover)) {
      edge_cell <- cellFromRowCol(cov, row, col)
    }
    col <- col - 1
  }
  
  # If we didn't find a cell, throw an error!
  if (edge_cell < 0) stop("Didn't find a west cell!")
  
  # Get all the trees that are in this cell, and have them calculate their
  # distance to the eastern edge of the edge cell
  dat$Wnoforestdist[cell_trees] <- dat@coords[cell_trees,1] - (xyFromCell(cov, edge_cell)[1] + half_cell)
  dat$Wnoforestcov [cell_trees] <- cell_cover
  #---------------------------------------------------------------------------#
}
    




#-----------------------------------------------------------------------------#
# The intercardinal directions must be done tree-by-tree
#-----------------------------------------------------------------------------#
for (i in 1:nrow(dat)) {
  
  if (dat$tree_cov[i] %in% forest_cover) {
  
    # Walk northeast
    results <- raytrace(dat@coords[i,1], dat@coords[i,2], dat$cell[i], "NE")
    dat$NEnoforestdist[i] <- results$firstdist
    dat$NEnoforestcov [i] <- results$firstcov
    
    # Walk northwest
    results <- raytrace(dat@coords[i,1], dat@coords[i,2], dat$cell[i], "NW")
    dat$NWnoforestdist[i] <- results$firstdist
    dat$NWnoforestcov [i] <- results$firstcov
    
    # Walk southeast
    results <- raytrace(dat@coords[i,1], dat@coords[i,2], dat$cell[i], "SE")
    dat$SEnoforestdist[i] <- results$firstdist
    dat$SEnoforestcov [i] <- results$firstcov
    
    # Walk southwest
    results <- raytrace(dat@coords[i,1], dat@coords[i,2], dat$cell[i], "SW")
    dat$SWnoforestdist[i] <- results$firstdist
    dat$SWnoforestcov [i] <- results$firstcov
  }
}

# Flag for whether the tree was in forest cover or not
dat$in_forest <- ifelse(dat$tree_cov %in% 41:43, 1, 0)
# Mean distance to nonforest cover
forestcols <- c(which(names(dat) == "Nnoforestdist"),
                which(names(dat) == "Snoforestdist"),
                which(names(dat) == "Enoforestdist"),
                which(names(dat) == "Wnoforestdist"),
                which(names(dat) == "NEnoforestdist"),
                which(names(dat) == "NWnoforestdist"),
                which(names(dat) == "SEnoforestdist"),
                which(names(dat) == "SWnoforestdist"))

# Adjust units if necessary! The value should be in m. Check the linear
# units of your projection.
dat$mean_noforestdist <- rowMeans(dat[,forestcols])
dat$mean_noforestdist[dat$in_forest==0] <- 0

# Clean up columns if desired - not required
dat$cell <- NULL
dat$in_forest <- NULL
dat$Nnoforestdist <- NULL
dat$Snoforestdist <- NULL
dat$Enoforestdist <- NULL
dat$Wnoforestdist <- NULL
dat$NEnoforestdist <- NULL
dat$NWnoforestdist <- NULL
dat$SEnoforestdist <- NULL
dat$SWnoforestdist <- NULL