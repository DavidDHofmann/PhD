################################################################################
#### Case Study: Employing the Dynamic + Model Approach to GPS of Hyenas
################################################################################
# Description: Application of the dynamic + model approach to hyena data. We'll
# run through three scenarios:
# - Forgiveness = 1
# - Forgiveness = 3, Habitat-Duration-Interactions = F
# - Forgiveness = 3, Habitat-Duration-Interactions = T

# Clear R's brain
rm(list = ls())

# Load required packages
library(amt)               # To do step-selection analysis
library(tidyverse)         # To wrangle data
library(pbmcapply)         # To apply stuff in parallel
library(survival)          # To fit conditional logistic regression model
library(raster)            # To handle spatial data

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_4"
setwd(wd)

# Load custom functions
source("02_R-Scripts/Functions.R")

# Parameters
n_rsteps <- 200

# Seed for reproducability
set.seed(12345)

# Load hyena gps data and covariates
dat <- read_rds("03_Data/Hyena_GPS.rds")
cov <- stack("03_Data/Hyena_Covariates.tif")

# Missingness is really low in this dataset, so we'll ever so slightly increase
# it.
dat <- dat %>% subset(DOP < 10) %>% rarifyData(missingness = 0.25)

# Visualize the tracks
ggplot() +
  geom_raster(data = na.omit(as.data.frame(cov[["Water"]], xy = T)), aes(x = x, y = y, fill = as.factor(Water))) +
  geom_path(data = dat, aes(x = x, y = y, col = timestamp), linewidth = 0.1) +
  geom_point(data = dat, aes(x = x, y = y, col = timestamp), pch = 20, size = 0.5) +
  scale_fill_manual(values = c("white", "cornflowerblue"), name = "Water") +
  scale_color_datetime(low = "red", high = "green", name = "Timestamp") +
  coord_sf() +
  theme_minimal()

# Resample the fixes to ensure they are all spread at least two hours apart
dat <- resFix(dat, hours = 2, start = 0, tol = 0.5)

# Round the timestamps to the nearest hour and compute durations between fixes
dat$timestamp <- round_date(dat$timestamp, "hour")
dat <- computeDurations(dat, units = "hour")

# Assess the median sampling rate for each individual (should be 2 hours for
# every individual)
summary(as.numeric(dat$duration))
count(dat, duration)
dat$duration <- NULL

# Define the regular step-duration
regular <- 2

# Fit static and dynamic distributions
dis <- fitDists(dat
  , dynamic     = T
  , durations   = c(2, 4, 6)
  , multicore   = T
  , resample    = F
  , rarify      = T
  , missingness = 0.1
)

# Compute bursts and associated metrics for a forgiveness of 1 and 3
dat_forg1 <- dat %>%
  computeBursts(max_duration = regular * 1) %>%
  computeMetrics() %>%
  computeSSF(n_rsteps = n_rsteps, dists = dis, approach = "dynamic+model") %>%
  computeCovars(covars = cov) %>%
  mutate(duration = as.factor(duration))
dat_forg3 <- dat %>%
  computeBursts(max_duration = regular * 3) %>%
  computeMetrics() %>%
  computeSSF(n_rsteps = n_rsteps, dists = dis, approach = "dynamic+model") %>%
  computeCovars(covars = cov) %>%
  mutate(duration = as.factor(duration))

# Model with Forgiveness = 1
mod_forg1 <- clogit(case ~
  + sl
  + log_sl
  + cos_ta
  + Water
  + DistanceToWater
  + Trees
  + strata(step_id)
  , data = dat_forg1
)

# Model with Forgiveness = 3, Interactions = F
mod_forg3_nohabitatints <- clogit(case ~
  + sl
  + log_sl
  + cos_ta
  + sl:duration
  + log_sl:duration
  + cos_ta:duration
  + Water
  + DistanceToWater
  + Trees
  + strata(step_id)
  , data = dat_forg3
)

# Model with Forgiveness = 3, Interactions = T
mod_forg3_habitatints <- clogit(case ~
  + sl
  + log_sl
  + cos_ta
  + sl:duration
  + log_sl:duration
  + cos_ta:duration
  + Water
  + DistanceToWater
  + Trees
  + Water:duration
  + DistanceToWater:duration
  + Trees:duration
  + strata(step_id)
  , data = dat_forg3
)

# Put them together into a single tibble and compute AIC
mod <- tibble(
    Forgiveness                 = c(1, 3, 3)
  , HabitatDurationInteractions = c(F, F, T)
  , Model                       = list(mod_forg1, mod_forg3_nohabitatints, mod_forg3_habitatints)
  , AIC                         = sapply(Model, AIC)
)

# Derivative of the correction function for the scale parameter. This is
# required to apply the delta method and retrieve the SE for the corrected
# scale parameter.
df_dsl <- function(scale, beta_sl) {
  ((1 / scale) - beta_sl) ** (-2)
}

# Go through the models and use results to compute updated distribution
# parameters (for the regular step duration), as well as the associated
# confidence intervals. Then, compile the habitat-selection estimates and
# updated movement parameters for the regular step duration.
mod$Results <- lapply(mod$Model, function(x) {

  # Extract model coefficients
  coefs <- summary(x)$coefficients
  coefs <- tibble(
      Coefficient = rownames(coefs)
    , Estimate    = coefs[, "coef"]
    , SE          = coefs[, "se(coef)"]
  )

  # Extract tentative movement parameters, correctors and their associated SE
  shape_0      <- dis$dynamic$sl$shape[1]
  scale_0      <- dis$dynamic$sl$scale[1]
  kappa_0      <- dis$dynamic$ta$kappa[1]
  mu_0         <- dis$dynamic$ta$mu[1]
  beta_sl      <- coefs$Estimate[coefs$Coefficient == "sl"]
  beta_log_sl  <- coefs$Estimate[coefs$Coefficient == "log_sl"]
  beta_cos_ta  <- coefs$Estimate[coefs$Coefficient == "cos_ta"]
  se_sl        <- coefs$SE[coefs$Coefficient == "sl"]
  se_log_sl    <- coefs$SE[coefs$Coefficient == "log_sl"]
  se_cos_ta    <- coefs$SE[coefs$Coefficient == "cos_ta"]

  # Compute corrected distribution parameters
  dis_updated <- updateDists(
      shape       = shape_0
    , scale       = scale_0
    , kappa       = kappa_0
    , mu          = mu_0
    , beta_sl     = beta_sl
    , beta_log_sl = beta_log_sl
    , beta_cos_ta = beta_cos_ta
  )

  # Compute the associated standard errors. The standard error for the shape and
  # kappa parameter are simply the standard errors associated with log_sl and
  # cos_ta. To compute the standard error for the scale, however, we need to
  # employ the delta method.
  se_scale <- sqrt(df_dsl(scale = scale_0, beta_sl = beta_sl) ** 2 * se_sl ** 2)
  se_shape <- se_log_sl
  se_kappa <- se_cos_ta

  # Put all results together
  coefs_habitat  <- coefs %>%
    subset(Coefficient %in% c("Water", "DistanceToWater", "Trees")) %>%
    mutate(Kernel = "Habitat")
  coefs_movement <- tibble(
      Coefficient = c("scale", "shape", "kappa")
    , Estimate    = c(dis_updated$sl$scale, dis_updated$sl$shape, dis_updated$ta$kappa)
    , SE          = c(se_scale, se_shape, se_kappa)
  ) %>% mutate(Kernel = "Movement")
  coefs <- rbind(coefs_habitat, coefs_movement) %>%
    dplyr::select(Variable = Coefficient, Value = Estimate, everything())

  # Compute 90%, 95%, and 99% confidence intervals
  coefs <- coefs %>%
    mutate(
        LCI_90 = Value - qnorm(1 - (1 - 0.90) / 2) * SE
      , UCI_90 = Value + qnorm(1 - (1 - 0.90) / 2) * SE
      , LCI_95 = Value - qnorm(1 - (1 - 0.95) / 2) * SE
      , UCI_95 = Value + qnorm(1 - (1 - 0.95) / 2) * SE
      , LCI_99 = Value - qnorm(1 - (1 - 0.99) / 2) * SE
      , UCI_99 = Value + qnorm(1 - (1 - 0.99) / 2) * SE
    )

  # Return them
  return(coefs)

})

# Store the results to file
write_rds(mod, "03_Data/CaseStudyResults.rds")

# Store session information
writeLines(capture.output(sessionInfo()), "02_R-Scripts/02_CaseStudy/SessionInfo/01_CaseStudy.txt")
