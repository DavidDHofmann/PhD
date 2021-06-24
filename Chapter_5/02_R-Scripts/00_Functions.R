################################################################################
#### Function to Get the PDF odf a Mixed von Mises Distribution
################################################################################
# Function to determine the pdf of a mixed von mises distribution
dvonMises <- function(x, k1, k2){
  exp(k1 * cos(x)) / (4 * pi * besselI(k1, nu = 0)) +
  exp(k2 * cos(x - pi)) / (4 * pi * besselI(k2, nu = 0))
}

# Function to randomly sample from a mixed von mises distribution
rvonMises <- function(size, k1, k2, by = 0.01){
  x <- seq(-pi, +pi, by = by)
  probs <- dvonMises(x, k1 = k1, k2 = k2)
  random <- sample(x, size = size, prob = probs, replace = T)
  return(random)
}

################################################################################
#### Function to Simulate Movement
################################################################################
# Function to simulate wolf walk
move <- function(
      n       = 10      # Number of coordinates to return
    , dist_sl = NULL    # Distribution of step lengths
    , dist_ta = NULL    # Distribution of turning angles
    , ext     = NULL    # Extent (bounding box)
    , start   = c(0, 0) # Initial coordinates
  ){

  # Make sure the extent is a spatial polygon
  if (!is.null(ext)){
    ext <- as(ext, "SpatialPolygons")
  }

  # Initiate vectors to store coordinates
  x <- rep(NA, n)
  y <- rep(NA, n)

  # Insert start coordinates
  x[1] <- start[1]
  y[1] <- start[2]

  # Generate random initial absolute turning angle
  absta <- runif(1, -pi, +pi)

  # Simulate movement
  for (i in 2:n){

    # Sample new steps until one is inside the study area
    inside <- F
    while (!inside){

      # Sample turning angle and step length
      sl <- rgamma(n = 1, scale = dist_sl$scale, shape = dist_sl$shape)
      ta <- rvonMises(size = 1, k1 = dist_ta$k1, k2 = dist_ta$k2)

      # Calculate new absolute turning angle
      absta_new <- absta + ta

      # We need to make sure that the absolute turning angle ranges from 0 to 2
      # * pi
      if (absta_new > 2 * pi){
        absta_new <- absta_new - 2 * pi
      } else if (absta_new < 0){
        absta_new <- absta_new + 2 * pi
      }

      # Calculate new endpoints
      x[i] <- x[i - 1] + sin(absta_new) * sl
      y[i] <- y[i - 1] + cos(absta_new) * sl

      # If an extent is provided, check if the new endpoints are within the
      # study extent
      if (!is.null(ext)){
        pts <- SpatialPoints(cbind(x[i], y[i]))
        crs(pts) <- crs(ext)
        inside <- gContains(ext, pts)
      } else {
        inside <- T
      }
    }
    absta <- absta_new
  }

  # Coerce the coordinates to a spatial line
  path <- spLines(SpatialPoints(cbind(x = x, y = y)))

  # Return the path
  return(path)
}

################################################################################
#### Function to Cut a Line to a Given Start and Endpoint
################################################################################
# Function to cut line to given start and endpoint (measured in distance from
# origin)
lineSubstring <- function(
      line  = NULL # Line to be cut
    , start = 0    # Distance at which the line should start
    , end   = 1    # Distance at which the line should end
  ){

  # Identify distances at which new lines begin
  dists <- gProject(line, as(line, "SpatialPoints"), normalized = T)

  # Identify lines that are within the desired distances
  index <- (dists >= start) & (dists <= end)

  # Retrieve coordinates of startpoint, midpoints, and endpoint
  res <- list(
      gInterpolate(line, start, normalized = T)
    , as(line, "SpatialPoints")[index, ]
    , gInterpolate(line, end, normalized = T)
  )

  # Bind them together and coerce them to "SpatialLines"
  res <- as(do.call(rbind, res), "SpatialLines")

  # Return the final line
  return(res)
}

################################################################################
#### Function to Cut Line into Segments at Intersections
################################################################################
# Function to cut line into segments at intersections
segmentLine <- function(
      line = NULL # Line to segment
    , poly = NULL # Polygon according to which the line is segmented
  ){

  # Identify intersection points
  ints <- gIntersection(as(poly, "SpatialLines"), line)

  # Can only segment if there is at least one intersection
  if (!is.null(ints)){

    # Calculate distance of intersections on line
    dists <- sort(gProject(line, ints, normalized = TRUE))
    dists <- c(0, dists, 1)

    # Segmentize line at intersections
    segments <- lapply(1:(length(dists) - 1), function(x){
      lineSubstring(line, dists[x], dists[x + 1])
    })
    segments <- do.call(rbind, segments)

    # Return the segments
    return(segments)

  # If there are no intersections...
  } else {

    # ... return the original line
    return(line)
  }
}

################################################################################
#### Function to Retain only Parts Between First and Last Intersection
################################################################################
# Function to segment a line to keep only the part between first and last
# intersection with another object
cropSequence <- function(line, polygon){

  # If the Line never intersects with the polygon, throw an error
  if (!gIntersects(line, polygon)){stop("Line does not intersect with polygon")}

  # Check if the line lies starts or ends within the polygon
  starts_inside <- gContains(polygon, gInterpolate(line, 0, normalized = T))
  ends_inside <- gContains(polygon, gInterpolate(line, 1, normalized = T))

  # Identify intersections (if there are any)
  poly <- as(polygon, "SpatialLines")
  line <- as(line, "SpatialLines")
  ints <- raster::intersect(poly, line)
  nints <- length(ints)

  # If the line lies fully inside the polygon, return everything
  if (length(ints) == 0){return(line)}

  # Project intersections onto line
  dists <- sort(gProject(line, ints, normalized = TRUE))

  # If the line starts inside the polygon and ends outside it, crop from 0 to
  # last intersection
  if (starts_inside & !ends_inside){return(lineSubstring(line, 0, dists[nints]))}

  # If the line starts inside the polygon and ends inside it, don't crop
  if (starts_inside & ends_inside){return(line)}

  # If the line starts outside the polygon and ends inside it, crop from first
  # intersection to 1
  if (!starts_inside & !ends_inside){return(lineSubstring(line, dists[1], 1))}

  # If the line starts outside the polygon and ends outside it, crop from first
  # intersection to the last one
  if (!starts_inside & !ends_inside){return(lineSubstring(line, dists[1], dists[ndists]))}
}

################################################################################
#### Function to Distribute Camera Traps in Study Area
################################################################################
# Function to distribute camera traps in study area
distributeCams <- function(
      area  = NULL  # Area in which cameras should be distributed
    , n     = 2     # Number of cameras per row/column
    , range = 2     # Range of the cameras
  ){

  # Create a raster with the desired number of traps
  grid <- raster(area, nrow = n, ncol = n, vals = 1:n**2)

  # Coerce the raster to SpatialPolygons and SpatialPoints
  polys   <- as(grid, "SpatialPolygons")
  points  <- SpatialPoints(coordinates(grid))

  # Buffer points to desired range
  range <- gBuffer(points, width = range, byid = T)

  # Return all spatial objects them
  return(list(Grid = polys, Points = points, Range = range))

}

################################################################################
#### Function to Calculate Distance to Walk between Points
################################################################################
# Function to calculate distance to walk between points
getDistance <- function(points = NULL){
  dist <- gDistance(points, byid = T)
  if (length(points) > 1){
    diag(dist) <- NA
    dist <- min(dist, na.rm = T)
    dist <- (length(points) - 1) * dist
    return(dist)
  } else {
    return(0)
  }
}
