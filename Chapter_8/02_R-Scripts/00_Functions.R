################################################################################
#### Function to calculate Distances on a Raster Efficiently
################################################################################
#' Calculate Raster Distance
#'
#' Function to calculate the distance of a raster cell to the nearest cell
#' containing a certain value
#' @export
#' @param x \code{RasterLayer} on which distances should be calculated
#' @param value value to which the distance should be calculated
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
