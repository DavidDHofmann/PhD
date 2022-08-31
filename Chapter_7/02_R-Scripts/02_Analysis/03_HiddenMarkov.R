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
library(runner)      # To compute moving windows

# Load cleaned activity data (note that the timestamps are all in UTC)
dat <- read_csv("03_Data/02_CleanData/ActivityDataCovariates.csv")

################################################################################
#### Functions
################################################################################
# Write the likelihood function
lik <- function(theta_star, x) {

  # Backtransform theta
  theta <- c(
      plogis(theta_star[1])
    , plogis(theta_star[2])
    , exp(theta_star[3])
    , exp(theta_star[4])
    , exp(theta_star[5])
    , exp(theta_star[6])
  )

  # Transition matrix and steady state
  t <- matrix(c(theta[1], 1 - theta[2], 1 - theta[1], theta[2]), nrow = 2)
  d <- solve(t(diag(2) - t + 1), c(1, 1))

  # Distributional parameters for the two states
  shape <- theta[3:4]
  scale <- theta[5:6]

  # Compute probabilities outside the loop -> Maybe this needs to be inside the
  # loop for the bayesian approach
  ind <- which(!is.na(x))
  probs <- matrix(1, length(x), 2)
  probs[ind, ] <- cbind(
      dgamma(x[ind], shape = shape[1], scale = scale[1])
    , dgamma(x[ind], shape = shape[2], scale = scale[2])
  )

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

# Write a function to backtransform model estimates
getTheta <- function(theta_est) {
  c(
    plogis(theta_est[1])
  , plogis(theta_est[2])
  , exp(theta_est[3])
  , exp(theta_est[4])
  , exp(theta_est[5])
  , exp(theta_est[6])
  )
}

# Function to decode the states
viterbi <- function(x, theta) {

  # Transition matrix and steady state
  t <- matrix(c(theta[1], 1 - theta[2], 1 - theta[1], theta[2]), nrow = 2)
  d <- solve(t(diag(2) - t + 1), c(1, 1))

  # Distributional parameters for the two states
  shape <- theta[3:4]
  scale <- theta[5:6]

  # Compute probabilities outside the loop -> Maybe this needs to be inside the
  # loop for the bayesian approach
  ind <- which(!is.na(x))
  probs <- matrix(1, length(x), 2)
  probs[ind, ] <- cbind(
      dgamma(x[ind], shape = shape[1], scale = scale[1])
    , dgamma(x[ind], shape = shape[2], scale = scale[2])
  )

  # Backwards algorithm
  n <- length(x)
  xi <- matrix(0, n, 2)
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
#### WHAT ABOUT MISSING DATA?
################################################################################

################################################################################
#### THINK ABOUT DELTA-T
################################################################################

# Subset to first individual
dat <- mutate(dat, Time = as_hms(Timestamp))
dat <- subset(dat, DogID == "Abel")
dat <- dat[1:1000, ]

# Function to apply a moving window to the data
windowActivity <- function(x, timestamps, window_size = "60_mins") {
  x_coarse <- runner(x
    , k = "60 mins"
    , idx = timestamps
    , f = function(x) {mean(x)}
  )
  return(x_coarse)
  return(x)
}

# Apply the function
dat$ActX60 <- windowActivity(dat$ActX, dat$Timestamp, window_size = "60 mins")

# Visualize
dat %>% ggplot(aes(x = Time, y = ActX60, group = as.factor(Date))) +
    geom_line(alpha = 0.5, lwd = 0.5) +
    theme_minimal() +
    scale_color_viridis_d() +
    theme(legend.position = "none")

# Initial values for the maximization
theta_star <- c(
    qlogis(0.5)
  , qlogis(0.5)
  , log(1)
  , log(5)
  , log(1)
  , log(5)
)

# Run maximization
ml_nlm <- nlm(lik, theta_star, x = dat$ActX60)
ml_nlm <- getTheta(ml_nlm$estimate)

# Decode the states
states <- viterbi(dat$ActX60, ml_nlm)

# Plot
all <- cbind(dat, StatesEstimated = states)
ggplot(all, aes(x = Timestamp, y = ActX60, col = StatesEstimated)) +
  geom_line() +
  theme_minimal()


################################################################################
#### Testing
################################################################################
# Run estimation for the different individuals
dat_nested <- nest(dat, Data = -DogID)
dat_nested <- mutate(dat_nested, Estimates = map(Data, function(x) {

  # Return the results
  return(results)

}))
print(dat_nested)

par(mfrow = c(2, 1))
plot(dat_nested$Data[[1]]$ActX60[1:200], type = "l")
plot(dat_nested$Estimates[[1]]$States[1:200], type = "l")

#### NEED TO WEIGH THE DISTRIBUTIONS!!!
# Visualize the distributions
x <- seq(0, 255, by = 0.1)
y1 <- dgamma(x, shape = ml_nlm[3], scale = ml_nlm[5])
y2 <- dgamma(x, shape = ml_nlm[4], scale = ml_nlm[6])
hist(sub$ActX60, freq = F)
lines(y1 ~ x)
lines(y2 ~ x)

# Apply it
states <- viterbi(sub$ActX60, ml_nlm)

# Visualize
table(states)
par(mfrow = c(3, 1))
plot(sub$ActX[501:1000], type = "l")
plot(sub$ActX60[501:1000], type = "l")
plot(states[501:1000], type = "o", pch = 20)





# # Write a custom gamma distribution
# dgammacust <- function(x, shape, scale, invert = F) {
#   if (invert) {
#     x <- x * (-1) + 255
#   }
#   probs <- dgamma(x, shape = shape, scale = scale)
#   probs[x < 0 | x > 255] <- 0
#   return(probs)
# }
#
# # Should be the same
# dgammacust(255, scale = 2, shape = 2, invert = T)
# dgammacust(0, scale = 2, shape = 2)
#
# # Should be the same
# dgammacust(254, scale = 2, shape = 2, invert = T)
# dgammacust(1, scale = 2, shape = 2)
#
# # Write the likelihood function
# lik <- function(theta_star, x) {
#
#   # Backtransform theta
#   theta <- c(
#       plogis(theta_star[1])
#     , plogis(theta_star[2])
#     , exp(theta_star[3])
#     , exp(theta_star[4])
#     , exp(theta_star[5])
#     , exp(theta_star[6])
#   )
#
#   # Transition matrix and steady state
#   t <- matrix(c(theta[1], 1 - theta[2], 1 - theta[1], theta[2]), nrow = 2)
#   d <- solve(t(diag(2) - t + 1), c(1, 1))
#
#   # Distributional parameters for the two states
#   shape <- theta[3:4]
#   scale <- theta[5:6]
#   lambd <- theta[7:8]
#
#   # Compute probabilities outside the loop -> Maybe this needs to be inside the
#   # loop for the bayesian approach
#   probs <- cbind(
#       dgammacust(x, shape = shape[1], scale = scale[1], invert = F)
#     , dgammacust(x, shape = shape[2], scale = scale[2], invert = T)
#   )
#
#   # Forward algorithm
#   foo <- d %*% diag(probs[1, ])
#   l <- log(sum(foo))       # to avoid numerical issues
#   phi <- foo / sum(foo)    # to avoid numerical issues
#   for (i in 2:length(x)) {
#     foo <- phi %*% t %*% diag(probs[i, ])
#     l <- l + log(sum(foo)) # to avoid numerical issues
#     phi <- foo / sum(foo)  # to avoid numerical issues
#   }
#
#   # We want to minimize the log-likelihood
#   return(-l)
# }
#
# # Try it!
# theta_star <- c(
#     qlogis(0.7)
#   , qlogis(0.7)
#   , log(1.5)
#   , log(10)
#   , log(1)
#   , log(10)
# )
#
# # Write a function to backtransform model estimates
# getTheta <- function(theta_est) {
#   c(
#     plogis(theta_est[1])
#   , plogis(theta_est[2])
#   , exp(theta_est[3])
#   , exp(theta_est[4])
#   , exp(theta_est[5])
#   , exp(theta_est[6])
#   )
# }
#
# # Using nlm
# ml_nlm <- nlm(lik, theta_star, x = sub$ActX)
# ml_nlm <- getTheta(ml_nlm$estimate)
#
# # Or using optim
# ml_optim <- optim(par = theta_star, fn = lik, x = sub$ActX)
# ml_optim <- getTheta(ml_optim$par)
#
# # Compare to truth (note that states may be switched!)
# cbind(ml_nlm, ml_optim)
#
# # Function to decode the states
# viterbi <- function(x, theta) {
#
#   # Transition matrix and steady state
#   t <- matrix(c(theta[1], 1 - theta[2], 1 - theta[1], theta[2]), nrow = 2)
#   d <- solve(t(diag(2) - t + 1), c(1, 1))
#
#   # Distributional parameters for the two states
#   shape <- theta[3:4]
#   scale <- theta[5:6]
#
#   # Compute probabilities outside the loop -> Maybe this needs to be inside the
#   # loop for the bayesian approach
#   probs <- cbind(
#       dgammacust(x, shape = shape[1], scale = scale[1], invert = F)
#     , dgammacust(x, shape = shape[2], scale = scale[2], invert = T)
#   )
#
#   # Backwards algorithm
#   n <- length(x)
#   xi <- matrix(0, n, 2)
#   foo <- d * probs[1, ]
#   xi[1, ] <- foo / sum(foo)
#   for (i in 2:n) {
#     foo <- apply(xi[i - 1, ] * t, 2, max) * probs[i, ]
#     xi[i, ] <- foo / sum(foo)
#   }
#   iv <- numeric(n)
#   iv[n] <- which.max(xi[n, ])
#   for (i in (n - 1):1) {
#     iv[i] <- which.max(t[, iv[i + 1]] * xi[i, ])
#   }
#   return(iv)
# }
#
# # Apply it
# states <- viterbi(sub$x, ml_optim)
#
# # Visualize
# table(states)
# plot(states)
