#-----------------------------------------------------------------------------#
# This function will hex a 100 m cell grid on top of the trees, identify which
# cell each tree is in and add that info to the "trees" SpatialDataFrame, and
# then create the "cells" dataframe to use for distance calculations.
# 
# This uses a mix of the sp and sf packages, mostly because only sf has a good
# way to calculate polygon centroids.
# 
# Since this is being used in the auspices of our own back-end model functions,
# this just assumes that "trees" is in the workspace already.
#-----------------------------------------------------------------------------#
createHexGridOnTrees <- function() {
  #---------------------------------------------------------------------------#
  # We need to create the hexagonal grid that sp is based on.
  # Create a rectangle spatial polygon bounding all points
  #---------------------------------------------------------------------------#
  minX <- min(trees$x) - (cellsize*2)
  maxX <- max(trees$x) + (cellsize*2)
  minY <- min(trees$y) - (cellsize*2)
  maxY <- max(trees$y) + (cellsize*2)
  coords <- rbind(c(minX, minY), 
                  c(minX, maxY), 
                  c(maxX, maxY),
                  c(maxX, minY),
                  c(minX, minY))
  
  # Make a Polygon
  p <- Polygon(coords, hole = F)
  
  # Make a Polygons - separate class
  poly_list <- list(Polygons(list(p), ID = 1))
  
  # Make our SpatialPolygons object
  sp_poly <- SpatialPolygons(poly_list)                                                             
  
  #-----------------------------------------------------------------------------#
  # Hex a grid over the top. 
  #-----------------------------------------------------------------------------#
  hex_pts <-spsample(sp_poly,type="hexagonal",cellsize=cellsize)
  hex_polys <- HexPoints2SpatialPolygons(hex_pts)
  #id <- names(hex_polys)
  #pid <- sapply(slot(hex_polys, "polygons"), function(x) slot(x, "ID"))
  #p.df <- data.frame( ID=1:length(hex_polys), row.names = pid)
  #writeOGR(SpatialPolygonsDataFrame(hex_polys, data=p.df), dsn=dsn, "oldpolys", driver="ESRI Shapefile")
  
  
  # Get the centroids of the hex polys
  cell_centers <- st_centroid(st_as_sf(hex_polys))#, byid = T)
  #df <- SpatialPointsDataFrame(as_Spatial(cell_centers), data=data.frame(ID=1:nrow(cell_centers)))
  #writeOGR(df, dsn=dsn, "oldpoints", driver="ESRI Shapefile")
  
  #-----------------------------------------------------------------------------#
  # Extract the polygon for each tree
  #-----------------------------------------------------------------------------#
  # x will be a vector of indexes of hex_polys that each tree falls into
  x <- over(trees, hex_polys)
  cell_id_for_tree <- names(hex_polys)[x]
  trees$cell <<- cell_id_for_tree
  
  #-----------------------------------------------------------------------------#
  # Filter and package cells data
  #-----------------------------------------------------------------------------#
  # Remove hex polys with no trees
  cell_centers <- cell_centers[unique(x),]
  hex_polys <- hex_polys[unique(x)]
  
  #save(cell_centers, hex_polys, file = "GIS data.Rdata")
  
  cells <<- as.data.frame(st_coordinates(cell_centers))
  cells$ID <<- names(hex_polys)
  cells <<- cells[cells$ID %in% trees$cell,]
}