
################################################################################
#### Function to calculate Distances on a Raster Efficiently
################################################################################
#' Calculate Raster Distance
#'
#' Function to calculate the distance of a raster cell to the nearest cell
#' containing a certain value
#' @export
#' @param x \code{RasterLayer} on which distances should be calculated
#' @param value numeric to which the distance should be calculated
#' @return \code{RasterLayer}
distanceTo <- function(x, value = 1) {

  # Convert to regular raster if terra raster is provided
  if (inherits(x, "SpatRaster")) {
    x <- raster(x)
    returnterra <- T
  }

  # Function requires two packages
  require(spatstat)
  require(maptools)

  # Convert raster to ppp object
  ppp <- x %>%
    rasterToPoints(., fun = function(x){x == value}, spatial = TRUE) %>%
    spTransform(., CRS("+init=epsg:32734")) %>%
    as(., "ppp")

  # Create empty raster onto which the distances are calculated
  distance <- raster(x)
  values(distance) <- distance %>%
    as(., "SpatialPoints") %>%
    spTransform(., CRS("+init=epsg:32734")) %>%
    as(., "ppp") %>%
    nncross(., ppp) %>%
    .[["dist"]]

  # Return the distance raster
  return(distance)
}

################################################################################
#### Function to Extend Raster and Fill with Observed Values
################################################################################
#' Extend Raster and Fill with Observed Values
#'
#' Function to extend a \code{SpatRaster} an fill the extended area with values
#' sampled from the original \code{SpatRaster}
#' @export
#' @param x \code{SpatRaster}
#' @param y extent to which the raster should be extended
#' @return \code{SpatRaster}
extendRaster <- function(x, y) {
  x_mask  <- is.na(x)
  x_mask  <- extend(x_mask, y, fill = F)
  x_large <- extend(x, y, fill = NA)
  for (i in 1:nlyr(x_large)) {
    indices <- which(is.na(x_large[[i]][]))
    n       <- length(indices)
    values  <- na.omit(x[[i]][])
    x_large[[i]][indices] <- sample(values, size = n, replace = T)
    x_large[[i]] <- mask(x_large[[i]], x_mask[[i]], maskvalues = 1, updatevalue = NA)
  }
  names(x_large) <- names(x)
  return(x_large)
}

################################################################################
#### Function to Prepare Movement Model for Dispersal Simulation
################################################################################
#' Prepare Movement Model for Dispersal Simulation
#'
#' Function to prepare movement model for dispersal simulation
#' @export
#' @param model Movement model
#' @return \code{list}
prepareModel <- function(model = NULL){

  # Extract parameter estimates from model
  coeffs <- fixef(model)$cond

  # Remove intercept
  coeffs <- coeffs[names(coeffs) != "(Intercept)"]

  # Obtain original model formula and coerce it to a terms object
  form <- model
  form <- formula(form)
  form <- terms(form)

  # Identify and remove unnecessary terms
  newform <- c("Intercept", "step_id_", "0 +")
  newform <- paste(newform, collapse = "|")
  newform <- grep(newform, attr(form, "term.labels"))
  newform <- drop.terms(form, newform, keep.response = F)
  newform <- formula(newform)

  # Return the coefficients and the formula
  return(list(
      coefficients  = coeffs
    , formula       = newform
  ))
}

################################################################################
#### Function to Prepare Raster Layers for Dispersal Simulation
################################################################################
#' Prepare Covariates for Simulation
#'
#' Function to prepare covariate layers for dispersal simulation
#' @export
#' @param layers \code{RasterStack} of covariates
#' @return \code{list}
prepareCovars <- function(layers){

  # Extract layer names
  names <- names(layers)

  # If it is a terra raster, convert it
  if (inherits(layers, "SpatRaster")) {
    layers <- as(layers, "Raster")
  }

  # Extract layer extent
  extent <- extent(layers)
  extent <- list(
      xmin = extent@xmin
    , xmax = extent@xmax
    , ymin = extent@ymin
    , ymax = extent@ymax
  )

  # Convert layers to velox raster
  layers <- velox(layers)

  # Return a named list
  return(list(
      covars      = layers
    , covar_names = names
    , extent      = extent
  ))
}

################################################################################
#### Function to Check if Points are within or outside an extent
################################################################################
#' Check if Coordinates are within Extent
#'
#' Function to Check if Points are within or outside an extent
#' @export
#' @param \code{data.frame} or \code{matrix} containing columns representing x
#' and y coordinates for which you want to know if they lie within an extent or
#' not.
#' @param \code{list} containing xmin, xmax, ymin, and ymax
#' @return Logical (\code{TRUE} = within extent, \code{FALSE} = outside extent)
pointsInside <- function(xy = NULL, extent = NULL){

  # Check if x value is within boundaries
  xy[, 1] > extent$xmin & xy[, 1] < extent$xmax &

    # Check if y value is within boundaries
    xy[, 2] > extent$ymin & xy[, 2] < extent$ymax
}

################################################################################
#### Function to Extract Covariates (Fast)
################################################################################
#' Extract covariates quickly
#'
#' Function to extract raster values and calculate their average coverage along
#' a spatial feature. This function requires the velox package which allows a
#' much quicker value extraction.
#' @export
#' @param raster \code{velox raster} containing covariates
#' @param feature Spatial feature (\code{sp} package) below which covariates
#' should be extracted
#' @return \code{data.frame}
extrCov <- function(raster = NULL, feature = NULL) {

  # Copy velox raster
  raster_copy <- raster$copy()

  # Crop the raster
  raster_copy$crop(feature)

  # Extract values
  extracted <- raster_copy$extract(feature, small = TRUE)

  # Calculate average values
  extracted <- vapply(extracted, colMeans, numeric(raster_copy$nbands))

  # Coerce to dataframe
  extracted <- as.data.frame(t(extracted))

  # Return the final dataframe
  return(extracted)
}

################################################################################
#### Function to Scale a Covariate Using a Scaling Table
################################################################################
#' Scale Covariates
#'
#' Function to scale covariates using a scaling table
#' @export
#' @param covars Either a \code{data.frame} or a \code{RasterStack} containing
#' continuous covariates.
#' @param scaling \code{data.frame} containing a "Center" column and a "Scale"
#' column. The \code{data.frame} must contain rownames to indicating the
#' covariates to which the scaling parameters belong to.
#' @return \code{data.frame} or \code{RasterStack}
scaleCovars <- function (covars, scaling) {
  scaled <- sapply(1:ncol(covars), function(x) {
    scale(covars[, x]
      , center  = scaling$center[names(covars)[x]]
      , scale   = scaling$scale[names(covars)[x]]
    )
  })
  scaled <- as.data.frame(scaled)
  names(scaled) <- names(covars)
  return(scaled)
}

################################################################################
#### Function to Reproject Coordinates
################################################################################
#' Reproject coordinates to another projection
#'
#' Function to transform coordinates into coordinates of another projection
#' @export
#' @param xy \code{data.frame} or matrix of coordinates that need to be
#' transformed.
#' @param from CRS of the coordinates
#' @param to CRS to which the coordinates should be reprojected
#' @return numeric value or vector
#' @examples
#' degreesToMeters(0.01)
reprojCoords <- function(xy, from = NULL, to = NULL, multicore = F) {
  xy <- as.matrix(xy)
  xy <- terra::project(xy, from, to)
  return(xy)
}

################################################################################
#### Function to Coerce Single Simulated Trajectories to Proper Lines
################################################################################
#' Coerce multiple Simulated Trajectories to SpatialLines
#'
#' Function to coerce simulated coordinates to lines after a desired number of
#' steps. Basically a wrapper around \code{sim2tracks()}
#' @export
#' @param simulation \code{data.frame} resulting from the function
#' \code{disperse()}
#' @param keep.data logical Should the input dataframe be conserved?
#' @return \code{SpatialLinesDataFrame}
sim2tracks <- function(simulation = NULL, crs = NULL, keep.data = F) {

  # Unless keep.data is desired, remove all unnecessary columns
  if (!keep.data){
    simulation <- simulation[, c("x", "y")]
  }

  # Calculate number of steps
  steps <- nrow(simulation) - 1

  # If data does not need to be kept
  if (!keep.data){
    coordinates(simulation) <- c("x", "y")
    lines <- spLines(simulation)
    if (!is.null(crs)){
      crs(lines) <- crs
    }
    return(lines)

  # If data needs to be kept
  } else {
    pts <- simulation
    coordinates(pts) <- c("x", "y")
    line <- spLines(pts)
    lines <- createSegments(line)
    lines <- as(lines, "SpatialLinesDataFrame")
    lines@data <- simulation[1:steps, ]
    if (!is.null(crs)){
      crs(lines) <- crs
    }
    return(lines)
  }
}

################################################################################
#### Function to Coerce Multiple Simulated Trajectories to Proper Lines
################################################################################
#' Coerce multiple Simulated Trajectories to SpatialLines
#'
#' Function to coerce simulated coordinates to lines after a desired number of
#' steps. Basically a wrapper around \code{sim2tracks()}
#' @export
#' @param simulations \code{data.frame} resulting from the function
#' \code{disperse}
#' @param id character column name of the column containing the track id
#' @param keep.data logical Should the input dataframe be conserved?
#' @param mc.cores integer Number of cores used for the coercion. Set to
#' 1 by default.
#' @param messages logical Should a progress bar and messages be printed?
#' @return \code{SpatialLinesDataFrame}
sims2tracks <- function(
    simulations = NULL
  , id          = "TrackID"
  , crs         = NULL
  , keep.data   = F
  , mc.cores    = 1
  , messages    = T
  ) {

  # Keep only relevant data
  if (!keep.data){
    simulations <- simulations[, c("x", "y", "TrackID")]
  }

  # Nest data by id
  simulations <- nest(simulations, data = -all_of(id))

  # Coerce each trajectory to a spatial lines object
  if (messages){
    lines <- pbmclapply(1:nrow(simulations)
    , mc.cores           = mc.cores
    , ignore.interactive = T
    , FUN                = function(x){
      l <- sim2tracks(
          simulation = simulations$data[[x]]
        , crs        = crs
        , keep.data  = keep.data
      )
      l$TrackID <- simulations$TrackID[x]
      return(l)
    })
  } else {
    lines <- mclapply(1:nrow(simulations)
    , mc.cores           = mc.cores
    , FUN                = function(x){
      l <- sim2tracks(
          simulation = simulations$data[[x]]
        , keep.data  = keep.data
      )
      l$TrackID <- simulations$TrackID[x]
      return(l)
    })
  }

  # Bind them
  lines <- do.call(rbind, lines)

  # Return the lines
  return(lines)

}

################################################################################
#### Function to Rasterize Using SpatStat
################################################################################
#' Rasterize Shapes Using spatstat
#'
#' Function to rasterize using spatstat
#' @export
#' @param l \code{SpatialPoints}, \code{SpatialLines}, or \code{SpatialPolygons}
#' to be rasterized
#' @param r \code{RasterLayer} onto which the objects should be rasterized
#' @param mc.cores numeric How many cores should be used?
#' @return \code{RasterLayer}
rasterizeSpatstat <- function(l, r, mc.cores = 1){

  # In case we run the rasterization on a single core, run a loop
  if (mc.cores == 1){

    # Create im layer
    values(r) <- 0
    im <- as.im.RasterLayer(r)
    summed <- im

    # Prepare progress bar
    pb <- txtProgressBar(
        min     = 0
      , max     = length(l)
      , initial = 0
      , style   = 3
      , width   = 55
    )

    # Go through all lines and rasterize them
    for (y in 1:length(l)){
      line    <- as.psp(l[y, ], window = im)
      line    <- as.mask.psp(line)
      line_r  <- as.im.owin(line, na.replace = 0)
      summed  <- Reduce("+", list(summed, line_r))
      setTxtProgressBar(pb, y)
    }

    # Return heatmap as a raster
    return(raster(summed))

  # If multicore is desired
  } else {

    # Split lines into a package for each core
    lines <- splitShape(l, n = mc.cores)

    # Run in parallel
    heatmap <- mclapply(lines, mc.cores = mc.cores, function(x){

      # Create im layer
      im <- as.im.RasterLayer(r)
      summed <- im

      # Rasterize each line
      for (y in 1:length(x)){
        line    <- as.psp(x[y, ], window = im)
        line    <- as.mask.psp(line)
        line_r  <- as.im.owin(line, na.replace = 0)
        summed  <- Reduce("+", list(summed, line_r))
      }

      # Return the sum
      return(summed)

    })

    # Combine heatmaps
    heatmap <- Reduce("+", heatmap)

    # Return as raster
    return(raster(heatmap))

  }
}

################################################################################
#### Function to Rasterize Simulations
################################################################################
#' Rasterize Simulations
#'
#' Function to rasterize using spatstat
#' @export
#' @param simulations \code{data.frame} containing the simulated data
#' @param raster \code{RasterLayer} onto which the objects should be rasterized
#' @param steps numeric, how many steps should be rasterized
#' @param area numeric, ID of the source areas that should be considered
#' @param flood character, one of "Min", "Mean", "Max"
#' @param messages boolean, should messages be printed during the rasterization
#' @param mc.cores numeric How many cores should be used?
#' @return \code{RasterLayer}
rasterizeSims <- function(
      simulations = NULL      # Simulated trajectories
    , raster      = NULL      # Raster onto which we rasterize
    , steps       = 500       # How many steps should be considered
    , area        = NULL      # Simulations from which areas?
    , flood       = "Min"     # Which flood level?
    , messages    = T         # Print update messages?
    , filename    = tempfile()
    , mc.cores    = detectCores() - 1
  ) {

  # Subset to corresponding data
  sub <- simulations[
    simulations$StepNumber <= steps &
    simulations$FloodLevel == flood, ]
  if (!is.null(area)) {
    sub <- sub[sub$SourceArea %in% area, ]
  }

  # Make sure raster values are all 0
  values(raster) <- 0

  # Create spatial lines
  sub_traj <- sims2tracks(
      simulations = sub
    , id          = "TrackID"
    , messages    = messages
    , mc.cores    = mc.cores
  )

  # Convert into spatial lines
  sub_traj <- as(sub_traj, "SpatialLines")
  crs(sub_traj) <- crs(raster)

  # Rasterize lines onto the cropped raster
  if (messages) {
    cat("Rasterizing spatial lines...\n")
  }
  heatmap <- rasterizeSpatstat(
      l        = sub_traj
    , r        = raster
    , mc.cores = 1
  )

  # Store the map to a temporary file and return it
  crs(heatmap) <- crs(raster)
  heatmap <- writeRaster(heatmap, filename)
  return(heatmap)
}

################################################################################
#### Function to Generate Visitation Histories
################################################################################
#' Generate Visitation Histories
#'
#' Function to extract the visitation history of a simulated trajectory
#' @export
#' @param x \code{vector} containing the visited raster cells
#' @param singlecount boolean, should double visits be considered?
#' @return \code{data.frame}
visitHist <- function(x, singlecount = F) {
  transitions <- data.frame(from = lag(x), to = x) %>%
    group_by(from, to) %>%
    subset(from != to) %>%
    na.omit() %>%
    summarize(TotalConnections = n(), .groups = "drop")
  if (singlecount){
    transitions$TotalConnections <- 1
  }
  return(transitions)
}


################################################################################
#### Function to Interpolate a Path
################################################################################
#' Interpolate a path
#'
#' Function that takes a sequence of xy coordinates and generates interpolated
#' points
#' @export
#' @param x \code{vector} x-coordinates
#' @param y \code{vector} y-coordinates
#' @param eps distance at which interpolated points should be placed
#' @return \code{matrix} of interpolated xy-coordinates
interpolatePath <- function(x, y, eps = 0.1) {
  inter <- lapply(1:(length(x) - 1), function(i) {
    xy_new <- interpolatePoints(
        x1 = x[i]
      , x2 = x[i + 1]
      , y1 = y[i]
      , y2 = y[i + 1]
      , by = eps
    )
    return(xy_new)
  }) %>% do.call(rbind, .)
  return(inter)
}


################################################################################
#### Function to Compute Betweenness
################################################################################
#' Calculate Betweenness for Simulations
#'
#' Function to compute betweenness on simulated data
#' @export
#' @param simulations \code{data.frame} containing the simulated data
#' @param raster \code{RasterLayer} based on which betweenness should be computed
#' @param steps numeric, how many steps should be considered?
#' @param area numeric, ID of the source areas that should be considered
#' @param flood character, one of "Min", "Mean", "Max"
#' @param messages boolean, should messages be printed during the rasterization
#' @param eps numeric, the interpolation distance (if desired) for point interpolation
#' @param filename character, the filename to which the raster should be stored (tempfile() by default)
#' @param mc.cores numeric, How many cores should be used?
#' @return \code{RasterLayer}
betweenSims <- function(
      simulations = NULL      # Simulated trajectories
    , raster      = NULL      # Raster onto which we rasterize
    , steps       = 500       # How many steps should be considered
    , area        = NULL      # Simulations from which areas?
    , flood       = "Min"     # Which flood level?
    , messages    = T         # Print update messages?
    , eps         = NULL      # Interpolation distance
    , filename    = tempfile()
    , mc.cores    = detectCores() - 1
  ) {

  # Subset to corresponding data
  sub <- simulations[
    simulations$StepNumber <= steps &
    simulations$FloodLevel == flood, ]
  if (!is.null(area)) {
    sub <- sub[sub$SourceArea %in% area, ]
  }

  # Prepare vertices and layout of the network for betweenness
  raster[] <- 1:ncell(raster)
  ver      <- values(raster)
  lay      <- as.matrix(as.data.frame(raster, xy = T)[, c(1, 2)])

  # Nest tracks by their IDs
  sub <- nest(sub, data = -TrackID)

  # Determine the visitation history of each path. Note that we will
  # "interpolate" each of the simulated steps. This will allow us to determine
  # cell-transitions at a much finer scale than if we would simply use the start
  # and endpoint of each step.
  if (messages) {
    cat("Computing visitation history...\n")
  }
  if (mc.cores > 1) {
    history <- pbmclapply(
        X                  = sub$data
      , ignore.interactive = T
      , mc.cores           = 1
      , FUN                = function(path) {
        if (!is.null(eps)) {
          path <- interpolatePath(path$x, path$y, eps = eps)
        }
        visits <- raster::extract(raster, path)
        visits <- visitHist(visits, singlecount = T)
        return(visits)
      })
  } else {
    if (messages) {
      pb <- txtProgressBar(min = 0, max = length(sub$data), style = 3)
    }
    history <- lapply(1:length(sub$data), function(x) {
      path <- sub$data[[x]]
      if (!is.null(eps)) {
        path <- interpolatePath(path$x, path$y, eps = eps)
      }
      visits <- raster::extract(raster, path)
      visits <- visitHist(visits, singlecount = T)
      if (messages) {
        setTxtProgressBar(pb, x)
      }
      return(visits)
    })
  }

  # Aggregate visitation histories across all paths
  history <- history %>%
    do.call(rbind, .) %>%
    group_by(from, to) %>%
    summarize(TotalConnections = sum(TotalConnections), .groups = "drop") %>%
    ungroup() %>%
    mutate(weight = mean(TotalConnections) / TotalConnections)

  # Create network, compute betweenness, and put values onto the raster
  if (messages) {
    cat("Computing betweenness...\n")
  }
  net <- graph_from_data_frame(history, vertices = ver)
  betweenness <- raster
  values(betweenness) <- betweenness(net)

  # Store the map to a temporary file and return it
  crs(betweenness) <- crs(raster)
  betweenness <- writeRaster(betweenness, filename)
  return(betweenness)
}


################################################################################
#### Function to Create Centroids that lie WITHIN Polygons
################################################################################
#' Create centroid within a polygon
#'
#' Function to create centroid within a polygon
#' @export
#' @param pol \code{SpatialPolygons} or \code{SpatialPolygonsDataFrame} in which
#' the centroids should be identified
#' @return \code{SpatialPointsDataFrame}
gCentroidWithin <- function(pol){

  # Load required packages
  require(rgeos)

  # Identify the number of polygons
  pol$.tmpID <- 1:length(pol)

  # Calculate centroids
  initialCents <- gCentroid(pol, byid = T)

  # Put data from the polygons to the centroids
  centsDF <- SpatialPointsDataFrame(initialCents, pol@data)

  # Indicate that these are true centroids
  centsDF$isCentroid <- TRUE

  # Check if the centroids are within the polygons
  centsInOwnPoly <- sapply(1:length(pol), function(x){
    gIntersects(pol[x,], centsDF[x, ])
  })

  # In case all centroids are within the polygons. We're done
  if(all(centsInOwnPoly) == TRUE){
        return(centsDF)
  } else {

    # We substitue outside centroids with points inside the polygon
    newPoints <- SpatialPointsDataFrame(
        gPointOnSurface(
            pol[!centsInOwnPoly, ]
          , byid = T
        )
      , pol@data[!centsInOwnPoly,]
    )

    # Indicate that these are not true centroids
    newPoints$isCentroid <- FALSE

    # Replace outside entrouds
    centsDF <- rbind(centsDF[centsInOwnPoly, ], newPoints)

    # Order points according to their polygon counterparts
    centsDF <- centsDF[order(centsDF$.tmpID), ]

    # Remove temporary ID column
    centsDF@data <- centsDF@data[, - which(names(centsDF@data) == ".tmpID")]

    # Return points
    return(centsDF)
  }
}

################################################################################
#### Function to darken a color
################################################################################
#' Darken a color
#'
#' Function to darken a color
#' @export
#' @param color A character string naming a color
#' @param factor Factor by which the color should be darkened. 1.4 by default
#' @return Hexadecimal code of the darkened color
#' @examples
#' darken("blue")
#' plot(1:2, 1:2, cex = 70, pch = 20, col = c("blue", darken("blue", 3)))
darken <- function(color, factor = 1.4){
    col <- col2rgb(color)
    col <- col / factor
    col <- rgb(t(col), maxColorValue = 255)
    col
}

################################################################################
#### Function to lighten a color
################################################################################
#' Lighten a color
#'
#' Function to lighten a color
#' @export
#' @param color A character string naming a color
#' @param factor Factor by which the color should be lightened. 1.4 by default
#' @return Hexadecimal code of the lightened color
#' @examples
#' lighten("blue")
#' plot(1:2, 1:2, cex = 70, pch = 20, col = c("blue", lighten("blue", 3)))
lighten <- function(color, factor = 1.4){
    col <- col2rgb(color)
    col <- col * factor
    col <- rgb(
        t(as.matrix(apply(col, 1, function(x) if (x > 255) 255 else x)))
      , maxColorValue = 255
    )
    col
}

################################################################################
#### Function to Rotate Coordinates
################################################################################
#' Rotate coordinates
#'
#' Function to rotate coordinates
#' @export
#' @param coords, matrix of coordinates to be rotated
#' @param degrees, degrees by which the coordinates should be rotated
#' @return matrix of rotated coordinates
#' @examples
#' rotate_coords(cbind(1, 1))
rotate_coords <- function(coords, degrees) {
  rad <- -1 * degrees * pi / 180
  x <- coords[, 1]
  y <- coords[, 2]
  x_new <- x * cos(rad) - y * sin(rad)
  y_new <- x * sin(rad) + y * cos(rad)
  coords_new <- cbind(x = x_new, y = y_new)
  return(coords_new)
}

################################################################################
#### Function to Generate an Ellipse
################################################################################
#' Generate an Ellipse
#'
#' Function to generate an ellipse
#' @export
#' @param xc numeric, x-coordiante of the centerpoint
#' @param yc numeric, y-coordinate of the centerpoint
#' @param a numeric, width
#' @param b numeric, height
#' @param phi numeric, rotation
#' @param spatial boolean, should a spatial object be returned?
#' @return either a matrix of ellipse coordiantes or an object of type \code{SpatVector}
#' @examples
#' elli <- ellipse(0, 0, 3, 2, 90)
#' plot(elli)
ellipse <- function(xc = 0, yc = 0, a = 2, b = 1, phi = 0, spatial = F) {
  x <- a * sin(seq(pi / 2, 5 * pi / 2, 0.1))
  y <- b * cos(seq(pi / 2, 5 * pi / 2, 0.1))
  x[length(x) + 1] <- x[1]
  y[length(y) + 1] <- y[1]
  # scale <- 1 / max(abs(c(x, y)))
  # x <- x * scale
  # y <- y * scale
  xy <- cbind(x = x, y = y)
  xy <- rotate_coords(xy, degrees = phi)
  xy[, 1] <- xy[, 1] + xc
  xy[, 2] <- xy[, 2] + yc
  if (spatial) {
    xy <- vect(xy, type = "polygons")
  }
  return(xy)
}

################################################################################
#### Function to Split Polygon Randomly
################################################################################
#' Split a Polygon randomly
#'
#' Function to split a polygon randomly
#' @export
#' @param x SpatVect
#' @return SpatVect Split polygon
splitPoly <- function(x) {

  # Get extent from the object
  lims <- ext(x)

  # Decide (randomly) whether to split vertically or horizontally
  split_vert <- rbernoulli(n = 1)

  # Coordinates of the splitting shape
  if (split_vert) {
    xy <- rbind(
        c(xmin(lims), ymin(lims))
      , c(xmin(lims), ymax(lims))
      , c(runif(1, min = xmin(lims), max = xmax(lims)), ymax(lims))
      , c(runif(1, min = xmin(lims), max = xmax(lims)), ymin(lims))
    )

  } else {
    xy <- rbind(
        c(xmin(lims), ymin(lims))
      , c(xmax(lims), ymin(lims))
      , c(xmax(lims), runif(1, min = ymin(lims), max = ymax(lims)))
      , c(xmin(lims), runif(1, min = ymin(lims), max = ymax(lims)))
    )
  }

  # Convert to spatial object and use it to clip the input shape
  shape <- vect(xy, type = "polygons")
  clip1 <- terra::intersect(x, shape)
  clip2 <- x - shape

  # Put all together and disaggregate into separate geometries
  clipped <- rbind(clip1, clip2)
  clipped <- disagg(clipped)
  return(clipped)
}

################################################################################
#### Function to Split Polygon Until its Area is < km
################################################################################
#' Split a Polygon randomly until its area falls below a threshold
#'
#' Split a Polygon randomly until its area falls below a threshold
#' @export
#' @param x SpatVect Polygon to be split
#' @param max_area numeric area above which polygon should be split
#' @return SpatVect Split polygon
splitPolyUntil <- function(x, max_area = NULL) {
  x$Area <- expanse(x, unit = "km")
  while (any(x$Area > max_area)) {

    # Split areas that are too small
    index <- which(x$Area > max_area)
    split <- list()
    for (i in 1:length(x)) {
      if (i %in% index) {
        split[[i]] <- splitPoly(x[i, ])
      } else {
        split[[i]] <- x[i, ]
      }
    }
    x      <- do.call(rbind, split)
    x$Area <- expanse(x, unit = "km")
  }
  return(x)
}

################################################################################
#### Function to Simulate a HR using Triangulation
################################################################################
#' Generate Home Ranges
#'
#' Function to simulate home ranges using triangulation algorithm
#' @export
#' @param n numeric. number of random points used for triangulation
#' @param m numeric. SpatRast of areas that should be removed
#' @param r SpatRast onto which the HRs should be returned/rasterized
#' @param buff numeric. Buffer width used around the extent of the input raster
#' @param max_area numeric area above which polygon should be split
#' @param min_area numeric area above which a polygon should be removed
#' @return SpatRast of simulated homeranges
simHR <- function(
      n        = NULL # Number of random points
    , m        = NULL # Mask of areas that shall not be suitable
    , r        = NULL # Reference raster onto which the hrs are rasterized
    , buff     = NULL # Should the input raster be buffered a little bit?
    , min_area = 100
    , max_area = 1000
  ) {

  # Get an extent from the input raster, and buffer if desired
  area <- ext(r)
  area <- as.polygons(area, crs = crs(r))
  if (!is.null(buff)) {
    area <- buffer(area, width = buff)
  }
  if (!is.null(m)) {
    area <- area - m
  }

  # Randomly place points in areas that shall not be removed
  pts <- spatSample(area, size = n)

  # Run triangulation to get polygons
  pols    <- voronoi(pts)
  if (!is.null(m)) {
    pols    <- pols - m
  }
  pols    <- crop(pols, ext(r))
  pols    <- disagg(pols)

  # Split polygons that are too big and remove polygons that are too small
  if (crs(pols) == "") { crs(pols) <- crs(r) }
  pols      <- splitPolyUntil(pols, max_area = max_area)
  pols      <- pols[pols$Area >= min_area, ]
  pols      <- disagg(pols)
  pols$ID   <- sample(1:length(pols))

  # Create a raster, rasterize polygons, and return the raster
  pols_r <- rasterize(pols, r, field = "ID")
  return(pols_r)
}

################################################################################
#### Function to Get Metrics of HRs
################################################################################
#' Get HR-Metrics from a raster of home ranges
#'
#' Function to calculate HR-Metrics from a raster of home ranges
#' @export
#' @param x SpatRaster of the home ranges
#' @return tibble containing the metrics
metricsHR <- function(x) {
  metrics <- x %>%
    expanse(unit = "km", byValue = T) %>%
    as.tibble() %>%
    summarize(
        minArea  = min(area)
      , meanArea = mean(area)
      , maxArea  = max(area)
      , Number   = length(area)
    )
  return(metrics)
}

################################################################################
#### Dispersal Simulation Function
################################################################################
#' Simulate Dispersal
#'
#' Function to simulate dispersal using a step selection function
#' @export
#' @param source matrix of start coordinates
#' @param covars spatial covariates prepared using the custom function
#' @param model issf model prepared using the custom function
#' @param sl_dist Step Length Distribution
#' @param sl_max What is the largest possible step?
#' @param date Start date
#' @param n_steps Number of steps simulated
#' @param n_rsteps Number of random steps proposed
#' @param scaling Dataframe to scale extracted covariates
#' @param stop Should the simulation stop at a boundary?
#' @param x SpatRaster of the home ranges
#' @return tibble containing the metrics
disperse <- function(
    source              = NULL    # Start Coordinates
  , covars              = NULL    # Spatial Covariates, prepared with our funct.
  , model               = NULL    # iSSF Model, prepared with our funct.
  , sl_dist             = NULL    # Step Length Distribution
  , sl_max              = Inf     # What is the largest possible step?
  , date                = as.POSIXct("2015-06-15 07:00:00", tz = "UTC")
  , n_steps             = 10      # Number of steps simulated
  , n_rsteps            = 25      # Number of random steps proposed
  , scaling             = NULL    # Dataframe to scale extracted covariates
  , stop                = F) {    # Should the simulation stop at a boundary?

  # # For testing
  # source <- cbind(x = c(23), y = c(-19))
  # covars <- cov$Min
  # model <- mod
  # date <- as.POSIXct("2015-06-15 07:00:00", tz = "UTC")
  # sl_max <- 35000
  # n_rsteps <- 25
  # stop <- F

  # Create a new dataframe indicating the first location. Note that we draw
  # random turning angles to start off
  track <- data.frame(
      x           = coordinates(source)[, 1]
    , y           = coordinates(source)[, 2]
    , absta_      = runif(1, min = 0, max = 2 * pi) # tentative
    , ta_         = NA
    , sl_         = NA
    , Timestamp   = date
    , BoundaryHit = FALSE
    , inactive    = NA
  )

  # Simulate random steps
  for (i in 1:n_steps) {

    # Check if the timestamp corresponds to low or high activity
    inactive <- strftime(date, tz = "UTC", format = "%H:%M:%S")
    inactive <- ifelse(inactive %in% c("07:00:00"), 1, 0)

    # Draw random turning angles
    ta_new <- runif(n_rsteps
      , min = - pi
      , max = + pi
    )

    # Draw random step lengths
    sl_new <- rgamma(n_rsteps
      , shape = sl_dist$params$shape
      , scale = sl_dist$params$scale
    )

    # In case the sl_ should be capped, do so
    if (sl_max != Inf) {
      sl_new <- pmin(sl_new, sl_max)
    }

    # Identify origin of track
    begincoords <- track[i, c("x", "y")]

    # Calculate new absolute turning angles
    absta_new <- getAbsNew(
        absta = track$absta_[i]
      , ta    = ta_new
    )

    # Calculate new endpoints
    endpoints_new <- calcEndpoints(
        xy    = as.matrix(track[i, c("x", "y")])
      , absta = absta_new
      , sl    = sl_new
    )

    # Check which endpoints leave the study extent
    inside <- pointsInside(
        xy     = endpoints_new
      , extent = covars$extent
    )

    # In case some steps are not inside the study area and we want the loop to
    # break
    if (sum(!inside) > 0 & stop) {

        # Break the loop
        break

      # In case some steps are not inside the study area and we DONT want the
      # loop to break
      } else if (sum(!inside) > 0 & !stop) {

        # Keep only steps inside the study area
        endpoints_new <- endpoints_new[inside, ]
        absta_new     <- absta_new[inside]
        ta_new        <- ta_new[inside]
        sl_new        <- sl_new[inside]

    }

    # Create spatial lines from origin to new coordinates
    l <- vector("list", nrow(endpoints_new))
    for (j in seq_along(l)){
        l[[j]] <- Lines(
          list(
            Line(
              rbind(
                  begincoords[1, ]
                , endpoints_new[j,]
              )
            )
          ), as.character(j)
        )
    }

    # Coerce to spatial lines
    steps <- SpatialLines(l)

    # Extract covariates along each step
    extracted <- extrCov(covars$covars, steps)

    # Put some nice column names
    names(extracted) <- covars$covar_names

    # Put everything into a dataframe
    rand <- data.frame(
        x           = endpoints_new[, 1]
      , y           = endpoints_new[, 2]
      , absta_      = absta_new
      , ta_         = ta_new
      , sl_         = sl_new
      , BoundaryHit = sum(!inside) > 0
    )

    # Put all covariates into a dataframe. We will use this to calculate
    # selection scores
    covariates <- data.frame(
        extracted
      , cos_ta_  = cos(ta_new)
      , log_sl_  = log(sl_new)
      , sl_      = sl_new
    )

    # Scale covariates
    covariates <- scaleCovars(covariates, scaling)

    # Put the activity phase into the covariate table as well
    covariates$inactive <- inactive

    # Prepare model matrix (and remove intercept)
    mat <- model.matrix(model$formula, covariates)
    mat <- mat[ , -1]

    # Calculate selection scores
    score <- exp(mat %*% model$coefficients)

    # Update date
    date <- date + hours(4)

    # Note that we assume that no fix exists at 11:00. In this case we add
    # another 4 hours
    if(strftime(date, tz = "UTC", format = "%H:%M:%S") == "11:00:00") {
      date <- date + hours(4)
    }

    # Coerce selection scores to probabilities
    Probs <- score / sum(score)

    # Sample a step according to the above predicted probabilities
    rand <- rand[sample(1:nrow(rand), size = 1, prob = Probs), ]

    # Add updated values to current step
    track$absta_[i]       <- rand$absta_[1]
    track$ta_[i]          <- rand$ta_[1]
    track$sl_[i]          <- rand$sl_[1]
    track$BoundaryHit[i]  <- rand$BoundaryHit[1]
    track$inactive[i]     <- inactive

    # Add new endpoints and new (tentative) absolute turning angle to the next
    # step
    track[i + 1, "x"]         <- rand$x[1]
    track[i + 1, "y"]         <- rand$y[1]
    track[i + 1, "absta_"]    <- rand$absta_[1]
    track[i + 1, "Timestamp"] <- date
  }
  return(track)
}
