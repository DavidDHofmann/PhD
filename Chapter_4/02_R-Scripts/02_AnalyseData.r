################################################################################
#### Simulating Movement Data with Known Preferences
################################################################################
# Description: Use ISSF analysis to simulate movement trajectories with known
# preferences

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)        # To handle spatial data
library(Rcpp)          # For faster point interpolation
library(pbmcapply)     # For parallel computing
library(tidyverse)     # For data wrangling
library(survival)      # To run conditional logistic regression

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

################################################################################
#### Helpful Functions
################################################################################
# Function to interpolate between spatial points
sourceCpp("/home/david/ownCloud/Dokumente/Bibliothek/Wissen/R-Scripts/interpolatepoints.cpp")

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

################################################################################
#### Prepare Data
################################################################################
# Load observed movement data
obs <- read_csv("03_Data/01_RawData/ObservedMovements.csv")

# Specify the different design combinations
dat <- expand_grid(
    Missingness = seq(0, 0.5, by = 0.1)
  , Futility    = 1:5
  , Replicate   = 1:10
)

# Create datasets with rarified observations. That is, randomly remove fixes,
# from 0% to 50%.
dat <- mutate(dat, Observations = map(Missingness, function(x){
  obs[sort(sample(1:nrow(obs), size = nrow(obs) * (1 - x))), ]
}))

# Prepare data for regular step selection analysis
dat <- mutate(dat, SSF = map2(Observations, Futility, function(x, y){

  # Calculate temporal difference between steps (in steps)
  x$duration <- x$step_number - lag(x$step_number)
  x$duration[1] <- 1

  # Define bursts. A new burst starts when the temporal lag is more than the
  # futility
  x$irregular <- x$duration > y
  x$burst <- cumsum(x$irregular) + 1
  x$irregular <- NULL

  # Calculate step length, absolute and relative turning angles for along each
  # burst
  x <- x %>%
    group_by(ID, burst) %>%
    nest() %>%
    mutate(data = map(data, function(z){
      z$sl <- stepLength(z$x, z$y)
      z$absta <- absAngle(z$x, z$y)
      z$ta <- z$absta - lag(z$absta)
      return(z)
    })) %>%
    unnest(data) %>%
    mutate(ta = case_when(
        ta < -pi ~ ta + 2 * pi
      , ta > +pi ~ ta - 2 * pi
      , TRUE ~ ta
    )) %>%
    dplyr::select(ID, burst, step_number, step_id, everything()) %>%
    ungroup()

}))

# Unnest all
test <- dat %>%
  dplyr::select(-Observations) %>%
  unnest(SSF)
nrow(test) / 1e6 * 25

# Check some of the datasets
dat %>%
  dplyr::select(-Observations) %>%
  subset(Missingness == 0.5 & Futility == 2 & Replicate == 1) %>%
  unnest(SSF)

################################################################################
#### CONTINUE HERE
################################################################################


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
