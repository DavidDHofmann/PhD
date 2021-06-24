################################################################################
#### Simulating Movement Data with Known Preferences
################################################################################
# Description: Use ISSF analysis to simulate movement trajectories with known
# preferences

# Clear R's brain
rm(list = ls())

# Load required packages
library(RandomFields)  # To simulate covariates
library(raster)        # To handle spatial data
library(Rcpp)          # For faster point interpolation
library(pbmcapply)     # To simulate in parallel
library(tidyverse)     # For data wrangling
library(survival)      # To run conditional logistic regression

# Function to interpolate between spatial points
sourceCpp("/home/david/ownCloud/Dokumente/Bibliothek/Wissen/R-Scripts/interpolatepoints.cpp")

################################################################################
#### Simulate Covariates
################################################################################
# Generate random covariate (elevation)
elev <- RMexp(var = 5, scale = 10) +
  RMnugget(var = 1) +
  RMtrend(mean = 0)
elev <- RFsimulate(elev, x = 1:200, y = 1:200)
elev <- raster(elev)

# Add a point of attraction
center <- coordinates(elev)
center <- colMeans(center)
center <- SpatialPoints(t(center))

# Calculate distance to center
dist <- distanceFromPoints(elev, center)

# Normalize both covariates
elev <- (elev - cellStats(elev, mean)) / cellStats(elev, sd)
dist <- (dist - cellStats(dist, mean)) / cellStats(dist, sd)

# Put covariate layers into a stack
covars <- stack(elev, dist)
names(covars) <- c("elev", "dist")

# Visualize all covariates
rasterVis::levelplot(covars, at = seq(-5, 5, length = 50))

################################################################################
#### Function To Simulate Movement
################################################################################
# Function to simulate movement
move <- function(
      xy       = NULL
    , covars   = NULL
    , prefs    = NULL
    , sl_dist  = NULL
    , n_steps  = 10
    , n_rsteps = 25
    , stop     = TRUE
  ){

  # create a new dataframe based on the source point. Note that we draw random
  # turning angles to start off
  track <- data.frame(
      x     = c(NA, xy[, 1])
    , y     = c(NA, xy[, 2])
    , absta = c(runif(1, min = 0, max = 2 * pi), NA)
    , ta    = c(runif(1, min = 0, max = 2 * pi), NA)
    , sl    = NA
  )

  # Get the extent of the covariates
  extent <- as(extent(covars), "SpatialPolygons")

  # Simulate random steps
  for (i in 2:n_steps){

    # Prepare an empty list in which we can store the random steps
    rand <- list()

    # Draw random turning angles
    ta_new <- runif(n_rsteps
      , min = -pi
      , max = +pi
    )

    # Draw random step lengths
    sl_new <- rgamma(n_rsteps
      , shape = sl_dist["shape"]
      , scale = sl_dist["scale"]
    )

    # Make sure that the steps cover at least a minimal distance
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
      if (nrow(rand[extent, ]) != n_rsteps){
        break
      }
    } else {
      rand <- rand[extent, ]
    }

    # Coerce back to regular dataframe
    rand <- as.data.frame(rand, xy = T)
    rand$xy <- NULL

    # Prepare a "line" for each random step. We first need the coordinates of
    # the steps for this
    begincoords <- track[i, c("x", "y")]
    endcoords   <- rand[, c("x", "y")]

    # Interpolate coordinates and extract covariates
    extracted <- sapply(1:nrow(endcoords), function(x){
      line <- interpolatePointsC(
          x1 = begincoords[1, 1]
        , x2 = endcoords[x, 1]
        , y1 = begincoords[1, 2]
        , y2 = endcoords[x, 2]
        , by = 1
      )
      extr <- raster::extract(covars, line)
      extr <- colMeans(extr)
      return(extr)
    })

    # Bind with other data
    rand <- cbind(rand, t(extracted))

    # Calculate cos_ta and log_sl
    rand$cos_ta <- cos(rand$ta)
    rand$log_sl <- log(rand$sl)

    # Prepare model matrix
    mat <- model.matrix(~ elev + dist + sl + log_sl + cos_ta, rand)
    mat <- mat[ , 2:ncol(mat)]

    # Calculate selection scores
    score <- exp(mat %*% prefs)

    # Convert scores to probabilities
    probs <- score / sum(score)

    # Keep only the step with the highest score
    rand <- rand[sample(nrow(rand), 1, prob = probs), ]

    # Add the step to our track
    track$absta[i] <- rand$absta
    track$ta[i] <- rand$ta
    track$sl[i] <- rand$sl
    track[i + 1, "x"] <- rand$x
    track[i + 1, "y"] <- rand$y
  }

  # Assign step numbers
  track$step_number <- 0:(nrow(track) - 1)

  # Return track, yet remove initial pseudo-fix
  return(track[-1, ])
}

################################################################################
#### Single Trajectory
################################################################################
# Simulation Parameters
stop      <- F
n_rsteps  <- 25
n_steps   <- 200
sl_dist   <- c(shape = 3, scale = 1)
prefs     <- c(
    elev   = 0.5
  , dist   = -3
  , sl     = 0.1
  , log_sl = 0.1
  , cos_ta = 1
)

# Simulate a test trajectory
sim <- move(
    xy       = matrix(c(100, 100), ncol = 2)
  , covars   = covars
  , stop     = stop
  , n_rsteps = n_rsteps
  , n_steps  = n_steps
  , sl_dist  = sl_dist
  , prefs    = prefs
)

# Visualize the simulation
plot(covars[[1]])
points(sim$y ~ sim$x, type = "o", pch = 16, cex = 0.5)
plot(center, add = T, col = "red", pch = 20)

################################################################################
#### Simulate Multiple Trajectories
################################################################################
# Number of simulated individuals
ndisp <- 10

# Simulate track for each individual
sims <- pbmclapply(
    X                  = 1:ndisp
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x){

  # Simulate trajectory
  sim <- move(
      xy       = matrix(c(100, 100), ncol = 2)
    , covars   = covars
    , stop     = stop
    , n_rsteps = n_rsteps
    , n_steps  = n_steps
    , sl_dist  = sl_dist
    , prefs    = prefs
  )

  # Assign unique simulation ID
  sim$ID <- x

  # Return the simulation
  return(sim)
})

# Bind simulations
sims <- do.call(rbind, sims)
rownames(sims) <- NULL

# Assign unique step ids
sims$step_id <- 1:nrow(sims)

# Visualize them
ids <- unique(sims$ID)
col <- rainbow(length(ids))
plot(covars[[1]])
for (i in 1:length(unique(sims$ID))){
  sub <- subset(sims, ID == ids[i])
  points(sub$y ~ sub$x, col = col[i], pch = 16, cex = 0.4)
  lines(sub$y ~ sub$x, col = col[i], lwd = 0.3)
}

################################################################################
#### Generate "Observed" Data
################################################################################
# Assume that we only observed xy data + ID
obs <- dplyr::select(sims, x, y, ID, step_number, step_id)

# Function to calculate the step length
stepLength <- function(x, y){
  length <- sqrt((x - lead(x)) ** 2 + (y - lead(y)) ** 2)
  return(length)
}

# Function to calculate the absolute turning angle
absAngle <- function(x, y, rad = T){
  xx <- lead(x) - x
  yy <- lead(y) - y
  xx <- na.omit(xx)
  yy <- na.omit(yy)
  b <- sign(xx)
  b[b == 0] <- 1
  tempangle <- b * (yy < 0) * pi + atan(xx / yy)
  tempangle[tempangle < 0] <- tempangle[tempangle < 0] + 2 * pi
  if (!rad){
    tempangle <- tempangle * 180 / pi
  }
  tempangle <- c(tempangle, NA)
  return(tempangle)
}

# Calculate step length, absolute and relative turning angles for along each
# trajectory
obs <- obs %>%
  group_by(ID) %>%
  nest() %>%
  mutate(data = map(data, function(x){
    x$sl <- stepLength(x$x, x$y)
    x$absta <- absAngle(x$x, x$y)
    x$ta <- x$absta - lag(x$absta)
    return(x)
  })) %>%
  unnest(data) %>%
  mutate(ta = case_when(
      ta < -pi ~ ta + 2 * pi
    , ta > +pi ~ ta - 2 * pi
    , TRUE ~ ta
  )) %>%
  dplyr::select(ID, step_number, step_id, everything())

# Let's ensure that we calculated metrics correctly
cbind(sims$sl, obs$sl)
cbind(sims$absta, obs$absta)
cbind(sims$ta, obs$ta)

################################################################################
#### Step Selection Analysis
################################################################################
# Number of random steps?
n_rsteps <- 25

# Indicate case steps
obs$case <- 1

# Cannot work with steps that have no turning angle
obs <- subset(obs, !is.na(ta))

# Create a new dataframe for alternative steps
rand <- obs[rep(1:nrow(obs), each = n_rsteps), ]

# Indicate that they are control steps (case = 0)
rand$case <- 0

# Sample new step lengths and turning angles
rand$sl <- rgamma(n = nrow(rand)
  , scale = sl_dist["scale"]
  , shape = sl_dist["shape"]
)
rand$ta_new <- runif(n = nrow(rand), min = -pi, max = +pi)

# Calculate new "absolute" turning angle
rand$ta_diff <- rand$ta - rand$ta_new
rand$absta <- rand$absta - rand$ta_diff
rand$absta[rand$absta > 2 * pi] <-
  rand$absta[rand$absta > 2 * pi] - 2 * pi
rand$absta[rand$absta < 0] <-
  rand$absta[rand$absta < 0] + 2 * pi
rand$ta <- rand$ta_new

# Remove undesired stuff
rand$ta_new <- NULL
rand$ta_diff <- NULL

# Put steps together
all <- rbind(obs, rand)
all <- arrange(all, ID, step_number, desc(case))

# Calculate new endpoints
all$x_to <- all$x + sin(all$absta) * all$sl
all$y_to <- all$y + cos(all$absta) * all$sl

# Sort
all <- dplyr::select(all, ID, step_number, step_id, x, y, x_to, y_to, everything())
all <- ungroup(all)

# Create interpolated coordinates for each step
extracted <- sapply(1:nrow(all), function(i){
  ints <- interpolatePointsC(
      x1 = all$x[i]
    , x2 = all$x_to[i]
    , y1 = all$y[i]
    , y2 = all$y_to[i]
    , by = 1
  )
  extr <- raster::extract(covars, ints)
  extr <- colMeans(extr)
  return(extr)
})

# Bind with other data
all <- cbind(all, t(extracted))

# Calculate movement metrics
all$log_sl <- log(all$sl)
all$cos_ta <- cos(all$ta)

# # Create lines from coordinates
# begincoords <- all[, c("x", "y")]
# endcoords <- all[, c("x_to", "y_to")]
# l <- vector("list", nrow(begincoords))
# for (i in seq_along(l)){
#   l[[i]] <- Lines(list(Line(rbind(
#       as.matrix(begincoords[i, ])
#     , as.matrix(endcoords[i,])
#   ))), as.character(i))
# }
#
# # Make the lines spatial
# lines <- SpatialLines(l)
# lines <- as(lines, "SpatialLinesDataFrame")
# lines@data <- all
#
# # Visualize
# head(lines)
# plot(lines)
#
# # Store to file
# library(rgdal)
# writeOGR(lines, ".", "test", driver = "ESRI Shapefile", overwrite = T)

################################################################################
#### Run Model
################################################################################
# Run model (here, we're not going to fit a mixed model)
mod <- clogit(case ~
  + elev
  + dist
  + sl
  + log_sl
  + cos_ta
  + strata(step_id)
  , data = all
)
summary(mod)

# Plot
mod_plot <- broom::tidy(mod)

# Compare estimates to true values
cbind(coef(mod), prefs)

# Visualize
preferences <- as.data.frame(prefs)
preferences$term <- rownames(preferences)
ggplot(mod_plot, aes(x = estimate, y = term)) +
  geom_vline(xintercept = 0, lty = 2, col = "gray80") +
  geom_errorbarh(aes(xmin = estimate - 1.96 * std.error, xmax = estimate + 1.96 * std.error), height = 0.2) +
  geom_point() +
  geom_point(data = preferences, aes(x = prefs, y = term), col = "red", pch = 1, cex = 2) +
  theme_classic()
