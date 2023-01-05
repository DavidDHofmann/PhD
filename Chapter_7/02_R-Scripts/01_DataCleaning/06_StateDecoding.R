################################################################################
#### Analysis of Average Activity using HMMs
################################################################################
# Description: Analysis of activity data using HMMs

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_7"
setwd(wd)

# Load required packages
library(tidyverse)   # For data wrangling
library(lubridate)   # To handle dates
library(hms)         # To handle times
library(pbmcapply)   # For multicore use with progress bar

# Load cleaned activity data (note that the timestamps are all in UTC)
dat <- read_csv("03_Data/02_CleanData/ActivityDataCovariates.csv")
dat <- subset(dat, State == "Resident")

# # For now, keep only the first 2000 entries per individual
# dat <- dat %>%
#   group_by(DogID) %>%
#   slice(1:2000) %>%
#   ungroup()

# Function to transform model parameters
transformTheta  <- function(theta, N, direction) {
  cut <- (N - 1) * N
  if (direction == "forward") {
    theta_trans <- c(
        qlogis(theta[1:cut])
      , log(theta[(cut + 1):length(theta)])
    )
  } else if (direction == "backward") {
    theta_trans <- c(
      plogis(theta[1:cut])
    , exp(theta[(cut + 1):length(theta)])
    )
  } else {
    stop("Can only transform 'forward' or 'backward'")
  }
  return(theta_trans)
}

# Compute likelihood for a given theta
lik <- function(theta_star, x, N) {

  # Backtransform theta
  # N <- 3
  theta <- transformTheta(theta = theta_star, N = N, direction = "backward")

  # Transition matrix and steady state
  cut     <- (N - 1) * N
  t       <- diag(N)
  t[!t]   <- theta[1:cut]
  t       <- t(t)
  diag(t) <- 2 - rowSums(t)
  d <- solve(t(diag(N) - t + 1), rep(1, N))

  # Distributional parameters for the two states
  shape <- theta[(cut + 1):(cut + N)]
  scale <- theta[(cut + 1 + N):length(theta)]

  # Compute probabilities outside the loop
  ind <- which(!is.na(x))
  probs <- matrix(1, length(x), N)
  for (i in 1:N) {
    probs[ind, i] <- dgamma(x[ind], shape = shape[i], scale = scale[i])
  }

  # Forward algorithm
  foo <- d %*% diag(probs[1, ])
  l <- log(sum(foo))       # to avoid numerical issues
  phi <- foo / sum(foo)    # to avoid numerical issues
  for (i in 2:length(x)) {
    foo <- phi %*% t %*% diag(probs[i, ])
    l <- l + log(sum(foo)) # to avoid numerical issues
    phi <- foo / sum(foo)  # to avoid numerical issues
  }

  # We want to minimize the log-likelihood
  return(-l)
}

# Function to decode the states
viterbi <- function(x, theta, N) {

  # Transition matrix and steady state
  cut     <- (N - 1) * N
  t       <- diag(N)
  t[!t]   <- theta[1:cut]
  t       <- t(t)
  diag(t) <- 2 - rowSums(t)
  d <- solve(t(diag(N) - t + 1), rep(1, N))

  # Distributional parameters for the two states
  shape <- theta[(cut + 1):(cut + N)]
  scale <- theta[(cut + 1 + N):length(theta)]

  # Compute probabilities outside the loop
  ind <- which(!is.na(x))
  probs <- matrix(1, length(x), N)
  for (i in 1:N) {
    probs[ind, i] <- dgamma(x[ind], shape = shape[i], scale = scale[i])
  }

  # Backwards algorithm
  n <- length(x)
  xi <- matrix(0, n, N)
  foo <- d * probs[1, ]
  xi[1, ] <- foo / sum(foo)
  for (i in 2:n) {
    foo <- apply(xi[i - 1, ] * t, 2, max) * probs[i, ]
    xi[i, ] <- foo / sum(foo)
  }
  iv <- numeric(n)
  iv[n] <- which.max(xi[n, ])
  for (i in (n - 1):1) {
    iv[i] <- which.max(t[, iv[i + 1]] * xi[i, ])
  }
  return(iv)
}

################################################################################
#### Single Individual
################################################################################
# Subset to a desired individual
print(unique(dat$DogID))
sub <- subset(dat, DogID == "Sishen")
sub <- sub[1:2000, ]

# Initial values
theta_star <- c(rep(0.2, 6), 1, 8, 15, 0.5, 2, 10)
theta_star <- transformTheta(theta_star, N = 3, direction = "forward")

# Using optim
ml_optim <- optim(par = theta_star, fn = lik, x = sub$Act30, N = 3)
ml_optim <- transformTheta(ml_optim$par, N = 3, direction = "backward")

# Get states
sub$State <- viterbi(sub$Act30, theta = ml_optim, N = 3)

# Plot
ggplot(sub, aes(x = as_hms(Timestamp), y = Act30, col = as.factor(State))) +
  geom_point(alpha = 0.5, size = 0.75) +
  theme_minimal() +
  scale_color_manual(values = c("cornflowerblue", "orange", "purple")) +
  facet_wrap(~Date, ncol = 1)

# Let's check if we get nice unimodal distributions within the states
# hist(sub$Act30[sub$State == 3])
# hist(sub$Act30[sub$State == 2])
# hist(sub$Act30[sub$State == 1])
hist(sub$Act[sub$State == 3])
hist(sub$Act[sub$State == 2])
hist(sub$Act[sub$State == 1])

################################################################################
#### Multiple Individuals
################################################################################
# Nest data by dog
dat_nested <- nest(dat, Data = -DogID)

# Let's run the HMM model for each individual
dat_nested$ModelParams <- pbmclapply(
    X                  = 1:nrow(dat_nested)
  , mc.cores           = detectCores() / 2
  , ignore.interactive = T
  , FUN                = function(i) {
    sub <- dat_nested$Data[[i]]
    ml_optim <- optim(par = theta_star, fn = lik, x = sub$Act30, N = 3)
    ml_optim <- transformTheta(ml_optim$par, N = 3, direction = "backward")
    return(ml_optim)
  })

# Let's decode the states for all individuals
dat_nested$Data <- pbmclapply(
    X                  = 1:nrow(dat_nested)
  , mc.cores           = detectCores() - 1
  , ignore.interactive = T
  , FUN                = function(i) {
    sub <- dat_nested$Data[[i]]
    par <- dat_nested$ModelParams[[i]]
    sub$State <- viterbi(sub$Act30, theta = par, N = 3)
    return(sub)
  })

# Store the results
write_rds(dat_nested, "03_Data/02_CleanData/ActivityDataCovariatesStates.rds")
