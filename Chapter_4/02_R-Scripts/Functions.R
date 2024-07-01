################################################################################
#### Functions to Simulate Landscape
################################################################################
# Function to simulate center of attraction (using a distance)
simDistance <- function(n, x, y) {
  r    <- raster(nrows = n, ncols = n, xmn = 0, xmx = n, ymn = 0, ymx = n)
  cent <- SpatialPoints(t(c(x, y)))
  dist <- distanceFromPoints(r, cent)
  dist <- (dist - cellStats(dist, min)) / (cellStats(dist, max) - cellStats(dist, min))
  return(dist)
}

# Function to simulate a CONTINUOUS (e.g. elevation) layer
simElevation <- function(n, autocorr_range, seed = NULL) {
  r <- raster(res = 1, xmn = 0, xmx = n, ymn = 0, ymx = n)
  w <- focalWeight(r, d = autocorr_range, type = "circle")
  add_rows <- (nrow(w) - 1) / 2
  add_cols <- (ncol(w) - 1) / 2
  r <- raster(res = 1
    , xmn = 0 - add_cols
    , xmx = n + add_cols
    , ymn = 0 - add_rows
    , ymx = n + add_rows
  )
  if (is.null(seed)) {
      set.seed(round(runif(1, 0, 1e8)))
    } else {
      set.seed(seed)
  }
  r[] <- rnorm(n = ncell(r))
  r   <- focal(r, w = w, fum = "sum")
  r   <- trim(r)
  r <- (r - cellStats(r, min)) / (cellStats(r, max) - cellStats(r, min))
  return(r)
}

# Function to simulate a BINARY (e.g. water) layer
simForest <- function(n, autocorr_range, prop = 0.5, seed = NULL) {
  if (is.null(seed)) {
    seed <- round(runif(1, 0, 1e8))
  }
  r      <- simElevation(n, autocorr_range, seed)
  cutoff <- quantile(r, 1 - prop)
  r      <- r > cutoff
  return(r)
}

################################################################################
#### Functions to Compute Step Metrics
################################################################################
# Function to compute the absolute turning angle
absTA <- function(dx, dy) {
  absta <- atan2(dy, dx)
  absta <- (absta - pi / 2) * (-1)
  absta <- ifelse(absta < 0, 2 * pi + absta, absta)
  return(absta)
}

# Function to compute the relative turning angle
relTA <- function(absta) {
  relta <- (absta[-1] - absta[-length(absta)])
  relta <- c(NA, relta)
  relta <- ifelse(relta > +pi, relta - 2 * pi, relta)
  relta <- ifelse(relta < -pi, 2 * pi + relta, relta)
  return(relta)
}

# Function to compute step metrics
stepMet <- function(x, y) {

    # Compute distances moved in x and y direction
    dx <- c(x[-1], NA) - x
    dy <- c(y[-1], NA) - y

    # Calculate step length
    sl <- sqrt(dx ** 2 + dy ** 2)

    # Compute absolute turn angle
    absta <- absTA(dx, dy)

    # Compute relative turn angle
    relta <- relTA(absta)

    # Put metrics into data.frame
    metrics <- data.frame(sl = sl, absta = absta, relta = relta)

    # Return the metrics
    return(metrics)
}

################################################################################
#### Function to Interpolate Between Two Coordinates
################################################################################
# Extracting covariate values along spatial lines can be really really slow. To
# improve efficiency, you can therefore extract covariates along points that are
# distributed on the line instead. Note that you may want to modify this
# function if you work with non-projected data!

# Function to interpolate coordinates between two points
interpolatePoints <- function(x1, x2, y1, y2, by = 1) {

  # Calculate length of line between points
  length <- sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

  # Calculate how many segments we need
  nsegs <- max(ceiling(length / by), 1)

  # Interpolate between points
  x <- seq(x1, x2, length.out = nsegs + 1)
  y <- seq(y1, y2, length.out = nsegs + 1)
  return(cbind(x, y))
}

# # Try it
# interpolatePoints(x1 = 0, x2 = 10, y1 = 0, y2 = 10, by = 2)
#
# # Let's generate a random line
# start <- runif(n = 2, min = 0, max = 100)
# end <- runif(n = 2, min = 0, max = 100)
# line <- spLines(SpatialPoints(rbind(start, end)))
#
# # Now distribute points on that line using our custom function
# # The smaller the value for "by", the more points are put on the line
# line_inter <- interpolatePoints(
#     x1 = start[1]
#   , x2 = end[1]
#   , y1 = start[2]
#   , y2 = end[2]
#   , by = 0.1
# )
#
# # Plot the line and interpolated coordinates
# plot(line)
# points(line_inter, col = "red", pch = 20, cex = 0.5)
#
# # Are results the same?
# a <- raster::extract(cov, line, fun = "mean")
# b <- colMeans(raster::extract(cov, line_inter))
# rbind(along_line = as.vector(a), at_points = b)
#
# # Let's see how extraction speeds compare
# benchmark <- microbenchmark(
#     AlongLine = raster::extract(cov, line, fun = "mean")
#   , AtPoints  = colMeans(raster::extract(cov, line_inter, fun = "mean"))
#   , times     = 20
# )
# summary(benchmark)

################################################################################
#### Von Mises Distribution Functions
################################################################################
# Function to determine the pdf of a mixed von mises distribution
dvonmises <- function(x, kappa, mu) {
  exp(kappa * cos(x - mu)) / (2 * pi * besselI(kappa, nu = 0))
}

# Function to randomly sample from a mixed von mises distribution
rvonmises <- function(n, kappa, mu, by = 0.01) {
  x <- seq(-pi, +pi, by = by)
  probs <- dvonmises(x, kappa = kappa, mu = mu)
  random <- sample(x, size = n, prob = probs, replace = T)
  return(random)
}

# # Let's ensure that the function works as expected
# x <- seq(-pi, +pi, by = 0.01)
# y1 <- dvonmises(x, mu = 0, kappa = 0)
# y2 <- dvonmises(x, mu = 0, kappa = 1)
# y3 <- dvonmises(x, mu = 0, kappa = 2)
#
# # Visualize the density for different kappas
# plot(NA, xlim = c(-pi, +pi), ylim = c(0, 0.6), xlab = "Turning Angle", ylab = "Prob. Density")
# abline(v = 0, lty = 2, col = "gray80")
# lines(y1 ~ x, col = "blue")
# lines(y2 ~ x, col = "purple")
# lines(y3 ~ x, col = "red")
# legend("topright"
#   , lty    = 1
#   , col    = c("blue", "purple", "red")
#   , legend = c("kappa = 0", "kappa = 1", "kappa = 2")
# )

################################################################################
#### Function to Simulate Movement from a fitted SSF
################################################################################
# Function to simulate movement
move <- function(
      xy       = NULL    # Source point (in matrix form -> n * 2)
    , covars   = NULL    # Stack of covariate layers (names need to match formula!)
    , formula  = NULL    # Model formula used to predict selection score
    , prefs    = NULL    # Preferences used to predict selection score
    , sl_dist  = NULL    # Parameters describing the step length distribution
    , ta_dist  = NULL    # Parameters describing the turning angle distribution
    , ext      = NULL    # Extent on which animals are allowed to move
    , n_steps  = 10      # Number of steps simulated
    , n_rsteps = 25      # Number of random steps proposed at each step
    , stop     = TRUE    # Should the simulation stop at boundaries?
    , messages = F       # Would you like to print update messages?
    , extract  = "end"   # Where to extract covariates (at the end or along steps)
  ) {

  #     # For testing only
  #     xy       <- matrix(runif(2, xmin(ext2), xmax(ext2)), ncol = 2)
  #     covars    <- cov
  #     formula   <- ~ elev + dist
  #     prefs     <- c(0.5, -0.2) # They need to match the formula!!!
  #     sl_dist   <- list(name = "gamma", params = list(shape = 3, scale = 0.5))
  #     ta_dist   <- list(name = "vonmises", params = list(kappa = 0.2, mu = 0))
  #     n_rsteps  <- 1000
  #     n_steps   <- 200
  #     stop      <- F
  #     extract   <- "end"

  # Create a new dataframe based on the source point. Note that we draw random
  # turning angles to start off
  track <- data.frame(
      x     = c(NA, xy[, 1])
    , y     = c(NA, xy[, 2])
    , absta = c(runif(1, min = 0, max = 2 * pi), NA)
    , ta    = NA
    , sl    = NA
  )

  # Simulate random steps
  if (messages) {
    pb <- txtProgressBar(0, n_steps, style = 3)
  }
  for (i in 2:n_steps) {

    # # For testing only
    # i <- 2

    # Draw random turning angles
    ta_new <- rvonmises(n_rsteps
      , mu    = ta_dist$params$mu
      , kappa = ta_dist$params$kappa
    )

    # Draw random step lengths
    sl_new <- rgamma(n_rsteps
      , shape = sl_dist$params$shape
      , scale = sl_dist$params$scale
    )

    # Make sure that the steps cover at least a minimal distance (this is
    # relevant if we need to compute the log of it)
    sl_new[sl_new < 0.0001] <- 0.0001

    # Put the step lengths and turning angles into a new dataframe. These are
    # our proposed random steps.
    rand <- data.frame(
        absta  = track$absta[i - 1] + ta_new
      , ta     = ta_new
      , sl     = sl_new
    )

    # We need to make sure that the absolute turning angle ranges from 0 to 2 *
    # pi
    rand$absta[rand$absta > 2 * pi] <-
      rand$absta[rand$absta > 2 * pi] - 2 * pi
    rand$absta[rand$absta < 0] <-
      rand$absta[rand$absta < 0] + 2 * pi

    # Calculate new endpoints
    rand$x <- track$x[i] + sin(rand$absta) * rand$sl
    rand$y <- track$y[i] + cos(rand$absta) * rand$sl

    # Create spatial points from endpoints
    coordinates(rand) <- c("x", "y")

    # Depending on the answer in the beginning, the loop breaks if one of the
    # new coordinates is outside the map boundaries
    if (stop){
      if (nrow(rand[ext, ]) != n_rsteps){
        break
      }
    } else {
      rand <- rand[ext, ]
    }

    # Coerce back to regular dataframe
    rand <- as.data.frame(rand, xy = T)
    rand$xy <- NULL

    # Prepare a "line" for each random step. We first need the coordinates of
    # the steps for this
    begincoords <- track[i, c("x", "y")]
    endcoords   <- rand[, c("x", "y")]

    # Interpolate coordinates and extract covariates
    if (extract == "end") {
      extracted <- raster::extract(covars, endcoords)
      rand <- cbind(rand, extracted)
    } else {
      extracted <- sapply(1:nrow(endcoords), function(x){
        line <- interpolatePoints(
            x1 = begincoords[1, 1]
          , x2 = endcoords[x, 1]
          , y1 = begincoords[1, 2]
          , y2 = endcoords[x, 2]
          , by = 0.1
        )
        extr <- raster::extract(covars, line)
        extr <- colMeans(extr)
        return(extr)
      })
      rand <- cbind(rand, t(extracted))
    }

    # Calculate cos_ta and log_sl
    rand$cos_ta <- cos(rand$ta)
    rand$log_sl <- log(rand$sl)

    # Prepare model matrix (and remove intercept)
    mat <- model.matrix(formula, rand)
    mat <- mat[ , -1]

    # Calculate selection scores
    score <- exp(mat %*% prefs)

    # Convert scores to probabilities
    probs <- score / sum(score)

    # Sample one of the steps based on predicted probabilities
    rand <- rand[sample(nrow(rand), 1, prob = probs), ]

    # Add the step to our track
    track$absta[i] <- rand$absta
    track$ta[i] <- rand$ta
    track$sl[i] <- rand$sl
    track[i + 1, "x"] <- rand$x
    track[i + 1, "y"] <- rand$y

    # Print update
    if (messages) {
      setTxtProgressBar(pb, i)
    }
  }

  # Assign step numbers
  track$step_number <- 0:(nrow(track) - 1)

  # Return track, yet remove initial pseudo-fix
  return(track[-1, ])
}

################################################################################
#### Function to Extend a Covariate Layer
################################################################################
# Function that extends a covariate layer and fills the added border with values
# sampled from the layer
extendRaster <- function(x, y) {
  x_mask  <- is.na(x)
  x_mask  <- extend(x_mask, y, value = 0)
  x_large <- extend(x, y, value = NA)
  for (i in 1:nlayers(x_large)) {
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
#### Main Analysis Functions
################################################################################
# Function that takes gps data and resamples fixes to a desired duration
resFix <- function(data, hours, start = 1, tol = 0.5) {

  # Nest by id and resample fixes
  data_resampled <- data %>%
    nest(Data = -id) %>%
    mutate(Data = map(Data, function(x) {
      first <- (min(x$timestamp) - days(1)) %>%
        update(., hour = start, min = 0, sec = 0)
      last <- (max(x$timestamp) + days(1)) %>%
        update(., hour = 24, min = 0, sec = 0)
      dates <- seq(first, last, by = paste0(hours, " hours")) %>%
        as.data.frame() %>%
        set_names("timestamp")
      closest <- sapply(1:nrow(dates), function(y) {
        index <- which.min(abs(dates$timestamp[y] - x$timestamp))[1]
        close <- as.numeric(abs(dates$timestamp[y] - x$timestamp[index]), units = "hours") <= tol
        if (close){
          return(index)
        } else {
          return(NA)
        }
      })
      closest <- na.omit(closest)
      return(x[closest, ])
    })) %>% unnest(Data)

  # Return
  return(data_resampled)

}

# Function to compute durations between fixes (by id)
computeDurations <- function(data, units = "hours") {
  data_durations <- data %>%
      nest(Data = -id) %>%
      mutate(Data = map(Data, function(x) {
        x$duration <- difftime(lead(x$timestamp), x$timestamp, units = units)
        return(x)
      })) %>% unnest(Data)
    return(data_durations)
}

# Function to create a dataset with rarified observations (the same function is
# also used to subsample x individuals from all individuals)
rarifyData <- function(data, missingness) {
  data_rarified <- data %>%
    nest(Data = -id) %>%
    mutate(Data = map(Data, function(x) {
      keep <- sample(1:nrow(x), size = nrow(x) * (1 - missingness))
      rarified <- x[sort(keep), ]
      return(rarified)
    })) %>% unnest(Data)
  return(data_rarified)
}

# Function to impute missing fixes to get a regularized dataframe. For
# simplicity, only conduct a single imputation.
imputeFixes <- function(data) {
  inds <- unique(data$id)
  pred <- lapply(inds, function(x) {

    # Fit model but suppress the unneccessary messages
    sink(tempfile())
    mod <- suppressMessages(suppressWarnings(crwMLE(
        mov.model   = ~ 1
      , data        = data[data$id == x, ]
      , Time.name   = "timestamp"
      , coord       = c("x", "y")
      , timeStep    = "hour"
      , method      = "L-BFGS-B"
      , control     = list(maxit = 100, trace = 0, REPORT = 1)
      , initialSANN = list(maxit = 100, trace = 1, REPORT = 1, temp = 0.5, tmax = 10)
      , attemps     = 1000
    )))
    sink()

    # Make prediction
    pre <- suppressMessages(suppressWarnings(crwPredict(
        object.crwFit = mod
      , predTime      = "1 hour"
      , flat          = TRUE
    )))

    # Clean the output
    pre <- data.frame(
        id          = pre$id
      , x           = pre$mu.x
      , y           = pre$mu.y
      , step_number = 1:nrow(pre)
      , step_id     = NA
      , timestamp   = pre$timestamp
    )

    # Indicate if a fix was imputed or not
    pre$imputed <- !(pre$timestamp %in% data$timestamp[data$id == x])

    # Return the imputed dataframe
    return(pre)
  })

  # Bind all individuals together, assign unique step id, return all
  pred <- do.call(rbind, pred)
  pred$step_id <- 1:nrow(pred)
  return(pred)
}

# Function to compute bursts (per id, depending on the fogriveness). A new burst
# always starts if the step-duration exceeds the maximally accepted duration
# (forgiveness * regular duration).
computeBursts <- function(data, max_duration) {

  # Nest data by id (a burst cannot expand across multiple ids)
  data_bursted <- data %>%
    nest(Data = -id) %>%

    # Compute bursts by id. A new burst is defined if the duration is >
    # max_duration
    mutate(Data = map(Data, function(x) {
      x$duration <- difftime(lead(x$timestamp), x$timestamp, units = "hours")
      x$irregular <- x$duration > max_duration
      x$burst <- NA
      x$burst[1] <- 1
      x$burst[2:nrow(x)] <- lag(cumsum(x$irregular) + 1)[-1]
      return(x)

    # Unnest the data again
    })) %>% unnest(Data) %>% ungroup()

  # Return the "bursted" data
  return(data_bursted)
}

# Function to compute step metrics (step length, relative turning angle,
# absolute turning angle). Step metrics are calculated on bursts.
computeMetrics <- function(data) {

  # Nest by id and burst
  data_metrics <- data %>%
    nest(Data = -c(id, burst)) %>%

    # Compute step metrics
    mutate(Data = map(Data, function(z) {
      metrics <- stepMet(x = z$x, y = z$y)
      metrics <- cbind(z, metrics)
      return(metrics)
    })) %>% unnest(Data) %>%

    # Tidy up
    dplyr::select(burst, step_number, step_id, everything()) %>%
    ungroup()

  # Return the data containing step metrics
  return(data_metrics)
}

# Function to fit step length and turning angle distributions. Statically, as
# well as dynamically to different step durations
fitDists <- function(data
    , durations   = 1
    , regular     = NULL
    , dynamic     = T
    , multicore   = F
    , resample    = T
    , start       = 0
    , rarify      = F
    , missingness = NULL
  ) {

  # Use first duration as regular duration, unless specified is provided
  if (is.null(regular)) {
    regular <- min(durations)
  }

  # Ensure that all parameters are provided
  if (rarify & is.null(missingness)) {
      stop("You need to provide a missingness value")
    } else if (!rarify) {
      missingness <- 0
  }

  # Make sure only unique durations are checked for
  durations <- sort(unique(durations))

  # Fit distribution to the specified duration
  metrics <- data %>%
    rarifyData(missingness = 0) %>%
    # resFix(hours = regular, start = start) %>%
    computeBursts(max_duration = regular) %>%
    computeMetrics() %>%
    subset(duration == regular) %>%
    dplyr::select(sl, relta) %>%
    tibble()
  sl_dist <- suppressWarnings(fit_distr(metrics$sl, dist_name = "gamma")$params)
  ta_dist <- suppressWarnings(fit_distr(metrics$relta, dist_name = "vonmises")$params)
  dists <- list(static = list(
      sl = list(shape = sl_dist$shape, scale = sl_dist$scale)
    , ta = list(kappa = ta_dist$kappa, mu = as.numeric(ta_dist$mu))
  ))

  # If desired, also get dynamic distributional parameters
  if (dynamic) {

    # Setup apply function to be used (depending on multicore use or not)
    if (multicore) {
        applier <- function(...) {
          pbmclapply(ignore.interactive = T, mc.cores = detectCores() - 1, ...)
        }
      } else {
        applier <- function(...) {
          lapply(...)
        }
    }

    # Fit distributions to different step-durations
    dists_dynamic <- tibble(duration = durations)
    dists_dynamic$Parameters <- applier(dists_dynamic$duration, function(i) {

      # Resample and rarify track if so desired
      metrics <- data %>%
        rarifyData(missingness = missingness) %>%
        {
          if (resample) {
            resFix(., hours = i, start = start)
          } else {
            .
          }
        } %>%
        computeBursts(max_duration = i) %>%
        computeMetrics() %>%
        subset(duration == i) %>%
        dplyr::select(sl, relta) %>%
        tibble()
      sl_dist <- suppressWarnings(fit_distr(metrics$sl, dist_name = "gamma")$params)
      ta_dist <- suppressWarnings(fit_distr(metrics$relta, dist_name = "vonmises")$params)
      results <- tibble(
          shape     = sl_dist$shape
        , scale     = sl_dist$scale
        , kappa     = ta_dist$kappa
        , mu        = as.numeric(ta_dist$mu)
        , nsteps_sl = length(na.omit(metrics$sl))
        , nsteps_ta = length(na.omit(metrics$relta))
      )
    })
    dists_dynamic <- unnest(dists_dynamic, Parameters)
    dists$dynamic <- list(
        sl = dplyr::select(dists_dynamic, duration, shape, scale, nsteps_sl)
      , ta = dplyr::select(dists_dynamic, duration, kappa, mu, nsteps_ta)
    )
  }

  # Return parameters
  return(dists)
}

# Function to generate random steps. The approach parameter is used to specify
# the approach that should be used to generate random steps
computeSSF <- function(data, n_rsteps, approach, dists) {

  # # TESTING
  # data <- testing[2, ]
  # data$x <- 0
  # data$y <- 0
  # # # data$duration <- c(1, 3)
  # data$duration <- 3
  # data$sl <- 1
  # data$absta <- 0
  # data$relta <- 0
  # n_rsteps <- 1
  # approach <- "multistep"

  # In case of the "model" approach, we want to avoid issues with convergence.
  # Hence, we only keep steps where the duration is represented a minimum number
  # of times. This avoids that we have to estimate a parameter for a duration of
  # 5, when there is only a single step with duration = 5
  keep <- data %>%
    count(duration) %>%
    subset(n > 5) %>%
    pull(duration)
  data <- subset(data, duration %in% keep)

  # Generate a new column that indicates that the steps are "observed" steps
  data$case <- 1

  # Cannot work with steps that have no turning angle, so remove them
  data <- subset(data, !is.na(relta))

  # Create a new dataframe into which we can put alternative/random steps
  rand <- data[rep(1:nrow(data), each = n_rsteps), ]

  # Indicate that these steps are random steps (case = 0)
  rand$case <- 0

  ##############################################################################
  #### Approach 1 - Uncorrected
  ##############################################################################
  # Step lengths sampled from "minimal" step-duration distributions
  if (approach == "uncorrected" | approach == "imputed") {
    rand$sl <- rgamma(n = nrow(rand)
      , scale = dists$static$sl$scale
      , shape = dists$static$sl$shape
    )
    rand$relta_new <- rvonmises(n = nrow(rand)
      , kappa = dists$static$ta$kappa
      , mu    = dists$static$ta$mu
      , by    = 0.01
    )

  ##############################################################################
  #### Approach 2 - Naive
  ##############################################################################
  # Step lengths adjusted merely by multiplying with the step duration
  } else if (approach == "naive") {
    rand$sl <- rgamma(n = nrow(rand)
      , scale = dists$static$sl$scale * rand$duration
      , shape = dists$static$sl$shape
    )
    rand$relta_new <- rvonmises(n = nrow(rand)
      , kappa = dists$static$ta$kappa
      , mu    = dists$static$ta$mu
      , by    = 0.01
    )

  ##############################################################################
  #### Approach 3 - Dynamic or Modelled
  ##############################################################################
  # Step lengths and turning angles sampled from distributions fitted to
  # different step-durations
  } else if (approach == "dynamic+model") {
    rand$sl <- rgamma(n = nrow(rand)
      , scale = dists$dynamic$sl$scale[match(rand$duration, dists$dynamic$sl$duration)]
      , shape = dists$dynamic$sl$shape[match(rand$duration, dists$dynamic$sl$duration)]
    )
    rand$relta_new <- sapply(1:nrow(rand), function(z) {
      rvonmises(n = 1
        , kappa = dists$dynamic$ta$kappa[dists$dynamic$ta$duration == rand$duration[z]]
        , mu    = dists$dynamic$ta$mu[dists$dynamic$ta$duration == rand$duration[z]]
        , by    = 0.01
      )
    })

  ##############################################################################
  #### Approach 4 - Multistep
  ##############################################################################
  # Step lengths and turning angles resulting from multiple random steps, where
  # the number of generated steps matches the step duration of the observed step
  } else if (approach == "multistep") {
    rand <- rand %>% group_by(duration) %>% nest()
    rand$data <- lapply(1:nrow(rand), function(z) {

      # Extract important info of the current iteration (i.e. the "z")
      duration <- rand$duration[z]
      x        <- rand$data[[z]]$x
      y        <- rand$data[[z]]$y
      relta    <- rand$data[[z]]$relta
      absta    <- rand$data[[z]]$absta

      # Generate sequence of random steps (if duration = 1, generate only one
      # step, if duration = 2, generate two, etc.)
      for (i in 1:duration) {
        sl <- rgamma(n = nrow(rand$data[[z]])
          , scale = dists$static$sl$scale
          , shape = dists$static$sl$shape
        )
        relta_new <- rvonmises(n = nrow(rand$data[[z]])
          , kappa = dists$static$ta$kappa
          , mu    = dists$static$ta$mu
          , by    = 0.01
        )

        # # TESTING
        # sl <- c(1, 1, 1)[i]
        # relta_new <- c(pi/2, 0, -pi/2)[i]
        # sl <- 1
        # relta_new <- 0

        # In case we are looking at the first step, we can only calculate the
        # new absolute turning angle by first deriving the difference in
        # relative turning angles
        if (i == 1) {
          relta_diff <- relta_new - relta
          absta <- absta + relta_diff
        } else {
          absta <- absta + relta_new
        }

        # Compute new endpoints
        absta <- ifelse(absta < 0, 2 * pi + absta, absta)
        absta <- ifelse(absta > 2 * pi, absta - 2 * pi, absta)
        x <- x + sin(absta) * sl
        y <- y + cos(absta) * sl
        relta <- relta_new
      }

      # Compute step length and relative turning angle from the origins to the
      # end coordinate of the last multi-random-steps
      dx <- x - rand$data[[z]]$x
      dy <- y - rand$data[[z]]$y
      sl <- sqrt(dx ** 2 + dy ** 2)
      absta <- atan2(dy, dx)
      absta <- (absta - pi / 2) * (-1)
      absta <- ifelse(absta < 0, 2 * pi + absta, absta)
      absta <- ifelse(absta > 2 * pi, absta - 2 * pi, absta)
      # relta <- absta - rand$data[[z]]$absta           # THIS IS CAUSING ISSUES!!!!
      relta <- rand$data[[z]]$relta - (rand$data[[z]]$absta - absta)
      relta <- ifelse(relta > +pi, relta - 2 * pi, relta)
      relta <- ifelse(relta < -pi, 2 * pi + relta, relta)
      rand$data[[z]]$sl <- sl
      rand$data[[z]]$relta_new <- relta
      return(rand$data[[z]])
    })
    rand <- unnest(rand, data)
    rand <- ungroup(rand)
    rand <- arrange(rand, id, step_id, desc(case))

  } else {
    stop("Provide valid input for the desired approach")
  }

  # Calculate new "absolute" turning angle
  rand$relta_diff <- rand$relta_new - rand$relta
  rand$absta <- rand$absta + rand$relta_diff
  rand$absta <- ifelse(rand$absta < 0, 2 * pi + rand$absta, rand$absta)
  rand$absta <- ifelse(rand$absta > 2 * pi, rand$absta - 2 * pi, rand$absta)
  rand$relta <- rand$relta_new

  # Remove undesired stuff
  rand$relta_new <- NULL
  rand$relta_diff <- NULL

  # Put steps together
  all <- rbind(data, rand)
  all <- arrange(all, id, step_id, desc(case))

  # Calculate new endpoints
  all$x_to <- all$x + sin(all$absta) * all$sl
  all$y_to <- all$y + cos(all$absta) * all$sl

  # Sort and ungroup
  all <- dplyr::select(
      all
    , id
    , step_number
    , step_id
    , x
    , y
    , x_to
    , y_to
    , everything()
  )
  all <- ungroup(all)
  return(all)
}

# Function to extract covariates along steps and compute step covariates
computeCovars <- function(data, covars, extract = "end", multicore = F) {

  # Run the covariate extraction on multiple cores
  if (multicore & extract == "along") {
    extracted <- pbmclapply(1:nrow(data), ignore.interactive = T, mc.cores = detectCores() - 1, function(x) {
      ints <- interpolatePoints(
          x1 = data$x[x]
        , x2 = data$x_to[x]
        , y1 = data$y[x]
        , y2 = data$y_to[x]
        , by = 1
      )
      extr <- raster::extract(covars, ints)
      extr <- colMeans(extr)
      return(extr)
    })
    extracted <- as.data.frame(do.call(rbind, extracted))

  # Run the covariate extraction on a single core
  } else if (extract == "along") {
    extracted <- lapply(1:nrow(data), function(x) {
      ints <- interpolatePoints(
          x1 = data$x[x]
        , x2 = data$x_to[x]
        , y1 = data$y[x]
        , y2 = data$y_to[x]
        , by = 1
      )
      extr <- raster::extract(covars, ints)
      extr <- colMeans(extr)
      return(extr)
    })
    extracted <- as.data.frame(do.call(rbind, extracted))
  } else {
    extracted <- raster::extract(covars, cbind(data$x_to, data$y_to))
  }

  # Put extracted covariates into a dataframe and bind them to the original data
  data <- cbind(data, extracted)

  # Ensure that step lengths cover a minimal distance
  data$sl[data$sl == 0] <- min(data$sl[data$sl != 0])

  # Calculate derived movement metrics
  data$log_sl <- log(data$sl)
  data$cos_ta <- cos(data$relta)

  # Return the results
  return(data)
}

# Function to run the step selection model using two different approaches
runModel <- function(data, approach, covariates) {

  # Run the step selection model
  if (approach == "dynamic+model" & length(unique(data$duration)) > 1) {
    data$duration <- as.factor(data$duration)
    form <- update(covariates, . ~ .
      + sl
      + log_sl
      + cos_ta
      + sl:duration
      + log_sl:duration
      + cos_ta:duration
      + strata(step_id)
    )
    mod <- clogit(form, data = data)
  } else {
    form <- update(covariates, . ~ .
      + sl
      + log_sl
      + cos_ta
      + strata(step_id)
    )
    mod <- clogit(form, data = data)
  }

  # Extract model coefficients and put them into a dataframe. Also, compute
  # confidence intervals.
  ci <- confint(mod, level = 0.95)
  coefs <- summary(mod)$coefficients
  coefs <- data.frame(
      Coefficient = rownames(coefs)
    , Estimate    = coefs[, "coef"]
    , SE          = coefs[, "se(coef)"]
    , Z           = coefs[, "z"]
    , LCI         = ci[, 1]
    , UCI         = ci[, 2]
  )
  rownames(coefs) <- NULL

  # Return them
  return(coefs)
}

################################################################################
#### Function to update Step Distirbutions given Parameters
################################################################################
# Function to udpate distributions
updateDists <- function(shape, scale, kappa, mu, beta_log_sl, beta_sl, beta_cos_ta) {

  # Update distribution parameters
  shape_upd <- shape + beta_log_sl
  scale_upd <- 1 / (1 / scale - beta_sl)
  kappa_upd <- kappa + beta_cos_ta
  mu_upd    <- mu

  # Return the updated distribution parameters
  updated <- list(
      sl = list(shape = shape_upd, scale = scale_upd)
    , ta = list(kappa = kappa_upd, mu = mu_upd)
  )
  return(updated)
}

################################################################################
#### Function to Correct Tentative Parameters given a Model Result
################################################################################
# Function to obtain the corrected distributional parameters from our results
correctDists <- function(result) {

  # Extract distribution and coefficients
  dists <- result$Dists
  coefs <- result$Coefs

  # Check if the distribution has been fitted dynamically
  dynamic <- !is.null(dists$dynamic)

  # Check if the step duration has been included as covariate
  model <- any(grepl(coefs$Coefficient, pattern = "duration"))

  # Extract the tentative parameters (we'll only do this for a step-duration of
  # one, i.e. the first row of the dynamic step parameters)
  if (dynamic) {
      shape_tent <- dists$dynamic$sl$shape[1]
      scale_tent <- dists$dynamic$sl$scale[1]
      kappa_tent <- dists$dynamic$ta$kappa[1]
      mu_tent    <- dists$dynamic$ta$mu[1]
    } else {
      shape_tent <- dists$static$sl$shape
      scale_tent <- dists$static$sl$scale
      kappa_tent <- dists$static$ta$kappa
      mu_tent    <- dists$static$ta$mu
  }

  # Extract model estimates (again, only for the regular step-duration)
  beta_cos_ta <- coefs$Estimate[coefs$Coefficient == "cos_ta"]
  beta_log_sl <- coefs$Estimate[coefs$Coefficient == "log_sl"]
  beta_sl     <- coefs$Estimate[coefs$Coefficient == "sl"]

  # Update the distributions
  updated <- updateDists(
      shape       = shape_tent
    , scale       = scale_tent
    , kappa       = kappa_tent
    , mu          = mu_tent
    , beta_log_sl = beta_log_sl
    , beta_sl     = beta_sl
    , beta_cos_ta = beta_cos_ta
  )

  # Return the updated distributions
  result$Dists$updated <- updated
  return(result)

}

################################################################################
#### Function to convert a piece of text to TeX for latex2exp
################################################################################
# Function to convert text to TeX
totex <- function(x) {
  TeX(paste0(x))
}
