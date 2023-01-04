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

# Unnest
unique(dat$DogID)
sub <- subset(dat, DogID = "Aspen")
par(mfrow = c(2, 2))
hist(sub$Act15, breaks = 25)
hist(sub$Act30, breaks = 25)
hist(sub$Act45, breaks = 25)
hist(sub$Act60, breaks = 25)

# Let's go with the 30 minutes cutoff
par(mfrow = c(2, 2))
hist(dat$Act15, breaks = 25)
hist(dat$Act30, breaks = 25)
hist(dat$Act45, breaks = 25)
hist(dat$Act60, breaks = 25)

# To be able to work with a gamma distribution, we need to get rid of 0s
dat$Act15[dat$Act15 == 0] <- min(dat$Act15[dat$Act15 > 0])
dat$Act30[dat$Act30 == 0] <- min(dat$Act30[dat$Act30 > 0])
dat$Act45[dat$Act45 == 0] <- min(dat$Act45[dat$Act45 > 0])
dat$Act60[dat$Act60 == 0] <- min(dat$Act60[dat$Act60 > 0])

# Likelihood for a custom gamma function
lik_gamma <- function(theta_star, x) {

  # Backtransform theta
  theta <- c(
      plogis(theta_star[1]) # Transition
    , plogis(theta_star[2]) # Transition
    , exp(theta_star[3])    # shape
    , exp(theta_star[4])    # shape
    , exp(theta_star[5])    # scale
    , exp(theta_star[6])    # scale
  )

  # Transition matrix and steady state
  t <- matrix(c(theta[1], 1 - theta[2], 1 - theta[1], theta[2]), nrow = 2)
  d <- solve(t(diag(2) - t + 1), c(1, 1))

  # Distributional parameters for the two states
  shape1 <- theta[3]
  shape2 <- theta[4]
  scale1 <- theta[5]
  scale2 <- theta[6]

  # Compute probabilities
  ind <- which(!is.na(x))
  probs <- matrix(1, length(x), 2)
  probs[ind, ] <- cbind(
      dgamma(x[ind], shape = shape1, scale = scale1)
    , dgamma(x[ind], shape = shape2, scale = scale2)
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
getTheta_gamma <- function(theta_est) {
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
viterbi_gamma <- function(x, theta) {

  # Transition matrix and steady state
  t <- matrix(c(theta[1], 1 - theta[2], 1 - theta[1], theta[2]), nrow = 2)
  d <- solve(t(diag(2) - t + 1), c(1, 1))

  # Distributional parameters for the two states
  shape1 <- theta[3]
  shape2 <- theta[4]
  scale1 <- theta[5]
  scale2 <- theta[6]

  # Compute probabilities outside the loop -> Maybe this needs to be inside the
  # loop for the bayesian approach
  ind <- which(!is.na(x))
  probs <- matrix(1, length(x), 2)
  probs[ind, ] <- cbind(
      dgamma(x[ind], shape = shape1, scale = scale1)
    , dgamma(x[ind], shape = shape2, scale = scale2)
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

# Initial values for the maximization
theta_star_gamma <- c(
    qlogis(0.5)
  , qlogis(0.5)
  , log(10)
  , log(5)
  , log(20)
  , log(15)
)

# Write a custom likelihood function for a fancy poisson distribution that
# allows us to handle the heavy tail.
lik_poisson <- function(theta_star, x) {

  # Original x (0 to 255) and inverted x (255 to 0). Allows us to account for
  # the fact that we have a regular and one mirrored poisson distribution.
  # Instead of mirroring the distribution, we mirror the x-values.
  x1 <- x
  x2 <- x1 * (-1) + 255

  # Backtransform theta
  theta <- c(
      plogis(theta_star[1]) # Transition
    , plogis(theta_star[2]) # Transition
    , exp(theta_star[3])    # Lambda
    , exp(theta_star[4])    # Lambda
  )

  # Transition matrix and steady state
  t <- matrix(c(theta[1], 1 - theta[2], 1 - theta[1], theta[2]), nrow = 2)
  d <- solve(t(diag(2) - t + 1), c(1, 1))

  # Distributional parameters for the two states
  lambda1 <- theta[3]
  lambda2 <- theta[4]

  # Compute probabilities
  ind <- which(!is.na(x))
  probs <- matrix(1, length(x), 2)
  probs[ind, ] <- cbind(
      dpois(x1[ind], lambda = lambda1)
    , dpois(x2[ind], lambda = lambda2)
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
getTheta_poisson <- function(theta_est) {
  c(
    plogis(theta_est[1])
  , plogis(theta_est[2])
  , exp(theta_est[3])
  , exp(theta_est[4])
  )
}

# Function to decode the states
viterbi_poisson <- function(x, theta) {

  # Convert x
  x1 <- x
  x2 <- x1 * (-1) + 255

  # Transition matrix and steady state
  t <- matrix(c(theta[1], 1 - theta[2], 1 - theta[1], theta[2]), nrow = 2)
  d <- solve(t(diag(2) - t + 1), c(1, 1))

  # Distributional parameters for the two states
  lambda1 <- theta[3]
  lambda2 <- theta[4]

  # Compute probabilities outside the loop -> Maybe this needs to be inside the
  # loop for the bayesian approach
  ind <- which(!is.na(x))
  probs <- matrix(1, length(x), 2)
  probs[ind, ] <- cbind(
      dpois(x1[ind], lambda = lambda1)
    , dpois(x2[ind], lambda = lambda2)
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

# Initial values for the maximization
theta_star_poisson <- c(
    qlogis(0.5)
  , qlogis(0.5)
  , log(5)
  , log(10)
)

################################################################################
#### Single Individual
################################################################################
# Let's test this with one individual
print(unique(dat$DogID))
sub <- subset(dat, DogID == "Odzala")
sub <- sub[1:2000, ]

# Initial values for the maximization
theta_star_gamma <- c(
    qlogis(0.5)
  , qlogis(0.5)
  , log(1)
  , log(5)
  , log(10)
  , log(50)
)

# Run maximization
ml_optim <- optim(par = theta_star_gamma, fn = lik_gamma, x = sub$Act30)
ml_optim <- getTheta_gamma(ml_optim$par)

# Decode the states
sub$State <- viterbi_gamma(sub$Act30, ml_optim)

# Plot
ggplot(sub[1:2000, ], aes(x = as_hms(Timestamp), y = Act30, col = as.factor(State))) +
  geom_point(alpha = 0.5, size = 0.75) +
  theme_minimal() +
  scale_color_manual(values = c("cornflowerblue", "orange")) +
  facet_wrap(~Date, ncol = 1)
hist(sub$Act30[sub$State == 2])
hist(sub$Act30[sub$State == 1])

################################################################################
#### Multiple Individuals
################################################################################
# Nest data by dog
dat_nested <- nest(dat, Data = -DogID)

# Let's run the HMM model for each individual
dat_nested$ModelParams <- pbmclapply(
    X                  = 1: nrow(dat_nested)
  , mc.cores           = detectCores() - 1
  , ignore.interactive = T
  , FUN                = function(i) {
    sub <- dat_nested$Data[[i]]
    ml_nlm <- suppressWarnings(nlm(lik, theta_star, x = sub$ActX))
    ml_nlm <- getTheta(ml_nlm$estimate)
    return(ml_nlm)
  })
print(dat_nested)

# Let's decode the states
dat_nested$Data <- pbmclapply(
    X                  = 1: nrow(dat_nested)
  , mc.cores           = detectCores() - 1
  , ignore.interactive = T
  , FUN                = function(i) {
    sub       <- dat_nested$Data[[i]]
    ml_nlm    <- dat_nested$ModelParams[[i]]
    sub$State <- viterbi(sub$ActX, ml_nlm)
    return(sub)
  })

# Store the results
write_rds(dat_nested, "03_Data/03_Results/99_DecodedHiddenMarkov.rds")
dat_nested <- read_rds("03_Data/03_Results/99_DecodedHiddenMarkov.rds")
sub <- dat_nested$Data[[2]]
sub$ActX60 <- windowActivity(sub$ActX, sub$Timestamp, window_size = "30 mins")
par(mfrow = c(2, 1))
hist(sub$ActX, breaks = 50)
hist(sub$ActX60, breaks = 50)
ggplot(sub[1000:3000, ], aes(x = Timestamp, y = ActX, col = as.factor(State))) +
  geom_point()
