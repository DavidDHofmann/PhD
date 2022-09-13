################################################################################
#### Hidden Markov Analysis of Dispersal Data
################################################################################
# Clear R's Brain
rm(list = ls())

# Load required packages
library(wilddogr)
library(tidyverse)
library(lubridate)
library(amt)
library(moveHMM)

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_9"
setwd(wd)

################################################################################
#### Data Download
################################################################################
# Identify files on dropbox (including the rvc data)
files <- dog_files(rvc = F)

# Let's take a look at the files we could download
head(files)

# Subset to the files of interest
todownload <- subset(files, DogName %in% c("Abel", "Everest", "Odzala", "Mirage"))

# Download the data
downloaded <- dog_download(
    x         = todownload
  , clean     = T
  , overwrite = T
  , outdir    = file.path(getwd(), "03_Data/01_RawData")
  , printpath = T
)

# Compute step metrics
dat <- downloaded %>%
  read_csv() %>%
  subset(!is.na(x) & !is.na(y)) %>%
  arrange(DogName, CollarID, Timestamp) %>%
  nest(Data = - c(DogName, CollarID)) %>%
  mutate(Data = map(Data, function(x) {
    x <- make_track(x, .x = "x", .y = "y", .t = "Timestamp", crs = 4326, all_cols = T)
    x <- transform_coords(x, crs_to = 32734)
    index <- which(hour(round(x$t_, "hour")) == 7)[1]
    x <- track_resample(x, rate = hours(24), tolerance = hours(1), start = index)
    x <- steps_by_burst(x, zero_dir = "N", clockwise = T, keep_cols = "start")
    return(x)
  })) %>%
  unnest(Data)

# Fit global step length and turning angles distributions to the two states
dists <- list(
    sl = list(Resident = NA, Disperser = NA)
  , ta = list(Resident = NA, Disperser = NA)
)
dist_sl_resi <- fit_distr(dat$sl_[dat$State == "Resident"], dist_name = "gamma")$params
dist_sl_disp <- fit_distr(dat$sl_[dat$State == "Disperser"], dist_name = "gamma")$params
dist_ta_resi <- fit_distr(dat$ta_[dat$State == "Resident"], dist_name = "vonmises")$params
dist_ta_disp <- fit_distr(dat$ta_[dat$State == "Disperser"], dist_name = "vonmises")$params

# Compute mean and sd for the gamma distribution
dist_sl_resi$mean <- dist_sl_resi$shape * dist_sl_resi$scale
dist_sl_disp$mean <- dist_sl_disp$shape * dist_sl_disp$scale
dist_sl_resi$sd <- sqrt(dist_sl_resi$shape * dist_sl_resi$scale ** 2)
dist_sl_disp$sd <- sqrt(dist_sl_disp$shape * dist_sl_disp$scale ** 2)

# Define initial values for step length distribution(s)
mu0 <- c(dist_sl_resi$mean, dist_sl_disp$mean)    # state 1 and state 2 parameters
sd0 <- c(dist_sl_resi$sd, dist_sl_disp$sd)        # state 1 and state 2 parameters
stepPar <- c(mu0, sd0)

# Define intial values for turning angle distribution(s)
anglemean0 <- c(0, 0)                                  # state 1 and state 2 parameters
anglecon0 <- c(dist_ta_resi$kappa, dist_ta_disp$kappa) # state 1 and state 2 parameters
anglePar <- c(anglemean0, anglecon0)

# Function to generate a "moveData" object
moveData <- function(data) {
  if(is.null(data$ID) | is.null(data$step) | is.null(data$x) | is.null(data$y)) {
    stop("Can't construct moveData object: fields are missing")
  }
  obj <- data
  class(obj) <- append("moveData",class(obj))
  return(obj)
}

# Function to fit an HMM for a single individual
runHMM <- function(dat) {
  dat_move <- dplyr::select(dat, ID = DogName, step = sl_, angle = ta_, x = x1_, y = y1_)
  dat_move <- moveData(dat_move)
  mod <- fitHMM(
      data      = dat_move
    , nbStates  = 2
    , stepPar0  = stepPar
    , anglePar0 = anglePar
  )
  return(mod)
}

# Run through the individuals


# Try it
test <- runHMM(subset(dat, DogName == "Abel"))


# Fit the model
dat_move <- dplyr::select(dat, ID = DogName, step = sl_, angle = ta_, x = x1_, y = y1_)
dat_move <- moveData(dat_move)
mod <- fitHMM(
    data      = dat_move
  , nbStates  = 2
  , stepPar0  = stepPar
  , anglePar0 = anglePar
)
dat$StateModel <- viterbi(mod)
table(dat$State, dat$StateModel, dat$DogName)
