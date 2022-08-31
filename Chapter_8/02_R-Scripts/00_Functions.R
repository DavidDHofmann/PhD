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

  # Create mask of NA values
  na_mask <- is.na(x)

  # Extract values from raster
  vals <- values(x, data.frame = F)

  # Remove na's
  vals <- na.omit(vals)

  # Extend raster to new extent
  r <- extend(x, y)

  # Also extend the mask
  na_mask <- extend(na_mask, y)

  # Identify NAs in the new extent
  indices <- which(is.na(values(r)))

  # Replace them with values sampled from old raster
  values(r)[indices] <- vals[runif(length(indices), min = 1, max = length(vals))]

  # Make sure original NA's are put back on the map
  r <- mask(r, na_mask, maskvalue = 1, updatevalue = NA)

  # Return the new raster
  return(r)
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
  xy <- suppressWarnings(vect(xy, crs = from))
  xy <- terra::project(xy, to)
  return(crds(xy))
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
    na.omit() %>%
    summarize(TotalConnections = n(), .groups = "drop")
  if (singlecount){
    transitions$TotalConnections <- 1
  }
  return(transitions)
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
#' @param mc.cores numeric How many cores should be used?
#' @return \code{RasterLayer}
betweenSims <- function(
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

  # Prepare vertices and layout of the network for betweenness
  raster[] <- 1:ncell(raster)
  ver <- values(raster)
  lay  <- as.matrix(as.data.frame(raster, xy = T)[, c(1, 2)])

  # Make coordinates of simulated trajectories spatial
  coordinates(sub) <- c("x", "y")
  crs(sub) <- crs(raster)

  # At each coordinate of the simulated trajectories we now extract the cell IDs
  # from the betweenness raster
  sub <- data.frame(
      TrackID    = sub$TrackID
    , StepNumber = sub$StepNumber
    , FloodLevel = sub$FloodLevel
    , x          = coordinates(sub)[, 1]
    , y          = coordinates(sub)[, 2]
    , r          = raster::extract(raster, sub)
  )

  # Get the visitation history
  sub <- nest(sub, data = -TrackID)
  if (messages) {
    cat("Computing visitation history...\n")
  }
  if (mc.cores > 1) {
    history <- pbmclapply(1:length(sub$data), mc.cores = mc.cores, ignore.interactive = T, function(x) {
      visitHist(sub$data[[x]]$r, singlecount = T)
    })
  } else {
    if (messages) {
      pb <- txtProgressBar(min = 0, max = length(sub$data), style = 3)
    }
    history <- lapply(1:length(sub$data), function(x) {
      hist <- visitHist(sub$data[[x]]$r, singlecount = T)
      if (messages) {
        setTxtProgressBar(pb, x)
      }
      return(hist)
    })
  }
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
