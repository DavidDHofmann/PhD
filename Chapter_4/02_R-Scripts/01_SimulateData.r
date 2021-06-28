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
library(pbmcapply)     # For parallel computing
library(tidyverse)     # For data wrangling

# Function to interpolate between spatial points
sourceCpp("/home/david/ownCloud/Dokumente/Bibliothek/Wissen/R-Scripts/interpolatepoints.cpp")

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

################################################################################
#### Simulate Covariates
################################################################################
# Generate random covariate (elevation)
elev <- RMexp(var = 5, scale = 10) +
  RMnugget(var = 1) +
  RMtrend(mean = 0)
elev <- RFsimulate(elev, x = 1:300, y = 1:300)
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

# Let's also specify an extent on which the animals are allowed to move
ext <- extent(c(50, 250, 50, 250))
ext <- as(ext, "SpatialPolygons")

# Visualize all covariates + extent
rasterVis::levelplot(covars, at = seq(-5, 5, length = 50))

################################################################################
#### Function To Simulate Movement
################################################################################
# Function to simulate movement
move <- function(
      xy       = NULL
    , covars   = NULL
    , ext      = NULL
    , prefs    = NULL
    , sl_dist  = NULL
    , n_steps  = 10
    , n_rsteps = 25
    , stop     = TRUE
  ){

  # If no extent is provided, use the extent of the covariates
  if (is.null(ext)){
    ext <- extent(covars)
  }
  ext <- as(ext, "SpatialPolygons")

  # create a new dataframe based on the source point. Note that we draw random
  # turning angles to start off
  track <- data.frame(
      x     = c(NA, xy[, 1])
    , y     = c(NA, xy[, 2])
    , absta = c(runif(1, min = 0, max = 2 * pi), NA)
    , ta    = c(runif(1, min = 0, max = 2 * pi), NA)
    , sl    = NA
  )

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

    # Prepare model matrix (and remove intercept)
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
    elev   = 0.5    # Preference towards elevation
  , dist   = -10.0   # Preference for distance to point of attraction
  , sl     = 0      # Preference for step length
  , log_sl = 0      # Preference for step length
  , cos_ta = 0      # Preference for turning angle
)

# # Prepare a function that resamples preferences from given distributions
# randomPrefs <- function(){
#   prefs <- c(
#       elev   = rnorm(n = 1, mean = 0.5, sd = 0.2)
#     , dist   = rnorm(n = 1, mean = -0.5, sd = 0.5)
#     , sl     = 0
#     , log_sl = 0
#     , cos_ta = 0
#   )
#   return(prefs)
# }
#
# # Try it
# randomPrefs()

# Simulate a test trajectory
sim <- move(
    xy       = matrix(runif(2, 50, 150), ncol = 2)
  , covars   = covars
  , ext      = ext
  , stop     = stop
  , n_rsteps = n_rsteps
  , n_steps  = n_steps
  , sl_dist  = sl_dist
  , prefs    = prefs
)

# Add an ID for the individual and an ID for each step
sim$ID <- 1
sim$step_id <- 1:nrow(sim)

# Visualize the simulation
plot(covars[[1]])
plot(ext, add = T, border = "red", lwd = 2)
points(sim$y ~ sim$x, type = "o", pch = 16, cex = 0.5)
plot(center, add = T, col = "red", pch = 20)

################################################################################
#### Multiple Trajectories
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
      xy       = matrix(runif(2, 50, 150), ncol = 2)
    , covars   = covars
    , ext      = ext
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

# Assign unique ID to each step
sims$step_id <- 1:nrow(sims)

# Visualize them
ids <- unique(sims$ID)
col <- rainbow(length(ids))
plot(covars[[1]])
plot(ext, add = T, border = "red", lwd = 2)
for (i in 1:length(unique(sims$ID))){
  sub <- subset(sims, ID == ids[i])
  points(sub$y ~ sub$x, col = col[i], pch = 16, cex = 0.4)
  lines(sub$y ~ sub$x, col = col[i], lwd = 0.3)
}

# In reality, we don't observe step lengths, turning angles etc., but we derive
# them from xy coordinates. Hence, let's that we only observed xy data + ID
obs <- dplyr::select(sims, x, y, ID, step_number, step_id)

# This is the type of data that we would observe in reality. Let's store it to
# file. Also store the simulated spatial layers to file
write_csv(obs, "03_Data/01_RawData/ObservedMovements.csv")
writeRaster(covars, "03_Data/01_RawData/CovariateLayers.tif")

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

################################################################################
#### Estimate Beta: USING CLOGIT
################################################################################
# Run model (only one variable)
mod <- clogit(case ~
  + elev
  + dist
  + strata(step_id)
  , data = all
)
summary(mod)
beta_cl <- coef(mod)

################################################################################
#### Estimate Beta: MAXIMUM LIKELIHOOD BY HAND
################################################################################
# Function to calculate the log-likelihoood (in parallel)
loglik <- function(beta){

  # Identify number of strata
  strata <- unique(all$step_id)

  # Go through each strata and calculate log-likelihood for proposed beta (here
  # we use parallel computing)
  loglik <- mclapply(strata, mc.cores = detectCores() - 1, function(x){
    sub <- all[all$step_id == x, ]
    mat <- model.matrix(case ~ elev + dist, sub)
    score <- exp(mat %*% c(0, beta))
    lik <- score[1] / sum(score)
    log(lik)
  })
  loglik <- do.call(rbind, loglik)

  # Return the summed log-likelihood
  sum <- sum(loglik)
  return(sum)

}

# Calculate log-likelihood for different betas
betas <- seq(0, 1, by = 0.01)
betas <- expand.grid(beta1 = seq(0, 1, by = 0.1), beta2 = seq(-1, 0, by = 0.1))
betas <- as.matrix(betas)
logliks <- sapply(betas, function(x){loglik(beta = x)})
logliks <- sapply(1:nrow(betas), function(x){loglik(beta = betas[x, ])})

# Find maximum likelihood estimate
beta_ml <- betas[which.max(logliks), ]

# Compare to model
cbind(ML_clogit = beta_cl, ML_hand = beta_ml)

# Plot
plot(logliks ~ betas[, 1], type = "o", xlab = "Beta", ylab = "Log-Likelihood", axes = F, pch = 16)
abline(v = beta_ml, lty = 2, col = "red")
abline(v = coef(mod), lty = 2, col = "blue")
legend("bottomright", col = c("red", "blue"), legend = c("Brute ML", "Model ML"), lty = c(2, 2))
axis(1)
axis(2)

################################################################################
#### Estimate Beta: MAXIMUM LIKELIHOOD USING MLE
################################################################################
# To use the mle function, we need to be able to calculate the negative
# log-likelihood instead of the log-likelihood
negloglik <- function(beta){
  loglik <- loglik(beta)
  negloglik <- -loglik
  return(negloglik)
}
beta_mle <- coef(mle(negloglik, start = list(beta = 0), method = "BFGS"))

# Compare to model
cbind(ML_clogit = beta_cl, ML_hand = beta_ml, ML_mle = beta_mle)

################################################################################
#### Estimate Beta: BAYESIAN
################################################################################
dat <- all[, c("step_id", "case", "elev", "dist")]
y <- dat %>%
  dplyr::select(step_id, case) %>%
  group_by(step_id) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  spread(key = row, value = case) %>%
  dplyr::select(-step_id) %>%
  as.matrix()
elev <- dat %>%
  dplyr::select(step_id, elev) %>%
  group_by(step_id) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  spread(key = row, value = elev) %>%
  dplyr::select(-step_id) %>%
  as.matrix()
dist <- dat %>%
  dplyr::select(step_id, dist) %>%
  group_by(step_id) %>%
  mutate(row = row_number()) %>%
  ungroup() %>%
  spread(key = row, value = dist) %>%
  dplyr::select(-step_id) %>%
  as.matrix()

# Bundle data
dat_jags <- list(
    y           = y
  , elev        = elev
  , dist        = dist
  , strata      = nrow(y)
  , strata_size = ncol(y)
)

# Write JAGS model file
cat(file = "model.txt", "model {

  # Prior
  beta_1 ~ dnorm(0, 0.001)
  beta_2 ~ dnorm(0, 0.001)

  # Likelihood
  for (i in 1:strata){
    for (j in 1:strata_size){
      phi[i, j] <- exp(beta_1 * elev[i, j] + beta_2 * dist[i, j])
    }
    for (j in 1:strata_size){
      lambda[i, j] <- phi[i, j] / sum(phi[i, 1:strata_size])
    }
    y[i, 1:strata_size] ~ dmulti(lambda[i, 1:strata_size], 1)
  }

}")

# Function to sample initial values
inits <- function(){
  list(
      beta_1 = rnorm(1)
    , beta_2 = rnorm(1)
  )
}

# Define parameters to be monitored (i.e. estimated)
params <- c("beta_1", "beta_2")

# MCMC settings (usually defined using trial and error)
na <- 2000        # Number of iterations in the adaptive phase
ni <- 2500        # Number of draws from the posterior (in each chain)
nb <- 1000        # Number of draws to discard as burn-in
nc <- 5           # Number of chains
nt <- 5           # Thinning rate (nt = 1 means we do not thin)

# Run the model
mod_jags <- jags(
    data               = dat_jags
  , inits              = inits
  , parameters.to.save = params
  , model.file         = "model.txt"
  , n.iter             = ni
  , n.burnin           = nb
  , n.chains           = nc
  , n.thin             = nt
  , n.adapt            = na
  , parallel           = T
)

# Show traceplots
par(mfrow = c(2, 2))
jagsUI::traceplot(mod_jags)

# Visualize distribution of beta_1
hist(mod_jags$sims.list$beta_1
  , breaks = 20
  , col    = "cornflowerblue"
  , main   = "Distribution of Beta"
  , xlab   = expression(beta)
  , xlim   = c(-1.5, 1.5)
)
abline(v = 0, lty = 2, col = "red")

# Visualize distribution of beta_2
hist(mod_jags$sims.list$beta_2
  , breaks = 20
  , col    = "cornflowerblue"
  , main   = "Distribution of Beta"
  , xlab   = expression(beta)
  , xlim   = c(-3, 3)
)
abline(v = 0, lty = 2, col = "red")

# Summary of output, rounded to 3 digits
print(mod_jags, 3)

# Let's look at the estimates
beta_jags <- c(mod_jags$mean$beta_1, mod_jags$mean$beta_2)

# Compare estimates to true parameters
cbind(
    ML_clogit = beta_cl
  , ML_hand   = beta_ml
  , ML_mle    = beta_mle
  , JAGS      = beta_jags
)

################################################################################
#### TESTING
################################################################################
# For testing only
# rmultinom(n = 1, prob = c(1, 1, 5), size = 1)
# dmultinom(x = c(0, 0, 1), prob = c(1, 1, 5))

# strata <- nrow(y)
# strata_size <- ncol(y)
# beta <- 0.5
# phi <- matrix(rep(NA, strata * strata_size), ncol = strata_size)
# lambda <- matrix(rep(NA, strata * strata_size), ncol = strata_size)

# for (i in 1:strata){
#   for (j in 1:strata_size){
#     phi[i, j] <- exp(beta * elev[i, j])
#   }
#   for (j in 1:strata_size){
#     lambda[i, j] <- phi[i, j] / sum(phi[i, 1:strata_size])
#   }
#   y[i, 1:strata_size] ~ dmultinom(x = y[i, ], prob = lambda[i, 1:strata_size])
# }
