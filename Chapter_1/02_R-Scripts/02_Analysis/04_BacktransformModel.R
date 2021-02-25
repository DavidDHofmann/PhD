################################################################################
#### Backtransform Scaled Model
################################################################################
# Description: Here, I want to backtransform the scaled model to retrieve
# coefficients for the movement metrics on the untransformed scale. However,
# this script is only supposed to make sure the function that I use to
# backtansform the data actually works. I can verify this by fitting an unscaled
# model and then compare the coefficients to a model with scaled covariates that
# I backtransform.

# Clear R's brain
rm(list = ls())

# Surpress scientific notation
options(scipen = 999)

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)    # For data wrangling
library(davidoff)     # Custom functions
library(glmmTMB)      # For modelling

################################################################################
#### Loading Data
################################################################################
# Load step selection data
dat <- read_csv("03_Data/02_CleanData/00_General_Dispersers_POPECOL(SSF_Extracted).csv")

# Keep only desired columns
dat <- dat %>% select(c(
      id
    , step_id_
    , case_
    , sl_
    , ta_
    , inactive
    , Water
    , DistanceToWater
    , Trees
    , Shrubs
    , HumansBuff5000 = HumanInfluenceBuffer_5000
))

# This is a crucial step. In the original model, step lengths were included in
# meters. This resulted in convergence issues, which I overcame by standardizing
# also the movement metrics. Alternatively, I could have also converted step
# lengths to kilometers, yet I realized that the model fit was better using
# standardized step lengths. Here, I will convert the step lengths to kilometers
# which will allow me to refit the full model and see if my
# "back-transformation" of the full model is correct.
dat$sl_ <- dat$sl_ / 1000

# We want to add the log of the step speed and the cosine of the turning angle
# and calculate the sqrt of the DistanceToWater
dat <- dat %>% mutate(
    log_sl_             = log(sl_)
  , cos_ta_             = cos(ta_)
  , SqrtDistanceToWater = sqrt(DistanceToWater)
)

# Let's also move all movement metrics to the front
dat <- dat %>% select(c(
  id, step_id_, case_, sl_, log_sl_, ta_, cos_ta_, inactive, everything()
))

# Create scaled covariates
dat <- transform(dat
  , sl_s                  = scale(sl_)
  , cos_ta_s              = scale(cos_ta_)
  , log_sl_s              = scale(log_sl_)
  , DistanceToWater_s     = scale(DistanceToWater)
  , SqrtDistanceToWater_s = scale(SqrtDistanceToWater)
  , Water_s               = scale(Water)
  , Trees_s               = scale(Trees)
  , Shrubs_s              = scale(Shrubs)
  , HumansBuff5000_s      = scale(HumansBuff5000)
)

################################################################################
#### Function to Backtransform a Model with Scaled Covariates
################################################################################
# Function to "back-transform" a model
backTransform <- function(
    model     = NULL
  , variables = NULL
  , means     = NULL
  , sds       = NULL){

  # Extract model coefficients
  coefs <- summary(model)$coefficients$cond[, 1]
  coefs_frame <- enframe(coefs)
  names(coefs_frame) <- c("Term", "Coefficient")

  # Split interaction terms and retrieve interacting variables
  coefs_frame <- separate(coefs_frame
    , col    = 1
    , into   = c("Var1", "Var2")
    , sep    = ":"
    , fill   = "right"
    , remove = F
  )

  # Generate column indicating interaction effects
  coefs_frame$Interaction <- ifelse(is.na(coefs_frame$Var2), F, T)

  # Find interactions where both variables need to be unscaled
  coefs_frame$BothScaled <-
    coefs_frame$Var1 %in% variables &
    coefs_frame$Var2 %in% variables

  # Prepare named vector for the backtransformed coefficients
  coefs_back <- coefs

  # Calculate new intercept
  # Part_A = old intercept
  # Part_B = coefficient * mean / sd for all variables that need to be unscaled
  # Part_C = (coefficient * mean_1 * mean_2) / (sd_1 * sd_2) for all
  # interactions where both variables need to be unscaled
  part_A <- coefs[1]
  part_B <- sapply(variables, function(x){
    coefs[x] * means[x] / sds[x]
  })
  part_C <- sapply(coefs_frame$Term[coefs_frame$BothScaled], function(x){
    covar1 <- coefs_frame$Var1[coefs_frame$Term == x]
    covar2 <- coefs_frame$Var2[coefs_frame$Term == x]
    (coefs[x] * means[covar1] * means[covar2]) / (sds[covar1] * sds[covar2])
  })
  part_C <- ifelse(length(part_C) == 0, 0, sum(part_C))
  coefs_back[1] <- part_A - sum(part_B) + sum(part_C)

  # Remove the intercept for simplicity
  coefs_frame <- coefs_frame[coefs_frame$Term != "(Intercept)", ]

  # Calculate new main effects. Formulas are different depending on whether the
  # respective main effect should be unscaled as well or not
  main <- coefs_frame$Term[!coefs_frame$Interaction]
  for (x in main){
    if (x %in% variables){
      part_A <- coefs[x] / sds[x]
      ints <- coefs_frame[
        (x == coefs_frame$Var1 | x == coefs_frame$Var2) &
        coefs_frame$BothScaled, ]
      if (nrow(ints) != 0){
        part_B <- sapply(1:nrow(ints), function(z){
          othervar <- c(ints$Var1[z], ints$Var2[z])
          othervar <- othervar[othervar != x]
          ints$Coefficient[z] * means[othervar] /
            (sds[ints$Var1[z]] * sds[ints$Var2[z]])
        })
      } else {
        part_B <- 0
      }
      coefs_back[x] <- part_A - sum(part_B)
    } else {
      part_A <- coefs[x]
      ints <- coefs_frame[
        (x == coefs_frame$Var1 | x == coefs_frame$Var2) &
        coefs_frame$Interaction, ]
      if (nrow(ints) != 0){
        part_B <- sapply(1:nrow(ints), function(z){
          othervar <- c(ints$Var1[z], ints$Var2[z])
          othervar <- othervar[othervar != x]
          ints$Coefficient[z] * means[othervar] / sds[othervar]
        })
      } else {
        part_B <- 0
      }
      coefs_back[x] <- part_A - sum(part_B)
    }
  }

  # Calculate new interaction effects. This is always the same
  ints <- coefs_frame$Term[coefs_frame$Interaction]
  for (x in ints){
    covar1 <- coefs_frame$Var1[coefs_frame$Term == x]
    covar2 <- coefs_frame$Var2[coefs_frame$Term == x]
    covars <- c(covar1, covar2)
    covars <- covars[covars %in% variables]
    coefs_back[x] <- coefs[x] / prod(sds[covars])
  }

  # Return the back transformed covariates
  return(coefs_back)
}

################################################################################
#### Show that Back-Transformation Works in a Simple Model
################################################################################
# Run an unscaled model
mod1 <- glmm_clogit(case_ ~
  + cos_ta_
  + sl_
  + log_sl_
  + Water
  + log_sl_:Water
  + log_sl_:inactive
  + (1 | step_id_)
  + (0 + cos_ta_ | id)
  + (0 + sl_ | id)
  + (0 + log_sl_ | id)
  + (0 + Water | id)
  , data = dat
)

# Run a model where all covariates are scaled
mod2 <- glmm_clogit(case_ ~
  + cos_ta_s
  + sl_s
  + log_sl_s
  + Water_s
  + log_sl_s:Water_s
  + log_sl_s:inactive
  + (1 | step_id_)
  + (0 + cos_ta_s | id)
  + (0 + sl_s | id)
  + (0 + log_sl_s | id)
  + (0 + Water_s | id)
  , data = dat
)

# Compare model summaries
summary(mod1)
summary(mod2)

# Do the backtransformation
variables <- c("cos_ta_s", "sl_s", "log_sl_s", "Water_s")
means <- c(
    "cos_ta_s" = mean(dat$cos_ta_)
  , "sl_s"     = mean(dat$sl_)
  , "log_sl_s" = mean(dat$log_sl_)
  , "Water_s"  = mean(dat$Water)
)
sds <- c(
    "cos_ta_s" = sd(dat$cos_ta_)
  , "sl_s"     = sd(dat$sl_)
  , "log_sl_s" = sd(dat$log_sl_)
  , "Water_s"  = sd(dat$Water)
)
back <- backTransform(mod2, variables = variables, means = means, sds = sds)

# Compare coefficients
cbind(
    summary(mod1)$coefficients$cond[, 1]
  , back
)

# Note: The intercept will not be correct. Although I don't know the reason for
# this, it doesn't matter as the intercept will not be used for anything.

################################################################################
#### Show that Back-Transformation Works in a Simple Model
################################################################################
# Run a model where only environmental covariates are scaled. This strucutre is
# identical to the proposed workflow in Fieberg et al. 2020
mod1 <- glmm_clogit(case_ ~
  + cos_ta_
  + log_sl_
  + Water_s
  + log_sl_:Water_s
  + log_sl_:inactive
  + (1 | step_id_)
  + (0 + cos_ta_ | id)
  + (0 + sl_ | id)
  + (0 + log_sl_ | id)
  + (0 + Water_s | id)
  , data = dat
)

# Run a model where all covariates are scaled
mod2 <- glmm_clogit(case_ ~
  + cos_ta_s
  + log_sl_s
  + Water_s
  + log_sl_s:Water_s
  + log_sl_s:inactive
  + (1 | step_id_)
  + (0 + cos_ta_s | id)
  + (0 + sl_s | id)
  + (0 + log_sl_s | id)
  + (0 + Water_s | id)
  , data = dat
)

# Compare model summaries
summary(mod1)
summary(mod2)

# Do the backtransformation of the log of the step length
variables <- c("log_sl_s", "cos_ta_s")
means <- c("log_sl_s" = mean(dat$log_sl_), "cos_ta_s" = mean(dat$cos_ta_))
sds <- c("log_sl_s" = sd(dat$log_sl_), "cos_ta_s" = sd(dat$cos_ta_))
back <- backTransform(mod2, variables = variables, means = means, sds = sds)

# Compare coefficients
cbind(
    summary(mod1)$coefficients$cond[, 1]
  , back
)

################################################################################
#### Show that Back-Transformation Works in the Full Model (sl_ in km)
################################################################################
# Only environmental predictors scaled (similar to Fieberg et al. 2020). This is
# basically the model we would have liked to fit with step lengths in meters.
mod1 <- glmm_clogit(case_ ~
  + cos_ta_
  + sl_
  + log_sl_
  + (1 | step_id_)
  + (0 + cos_ta_ | id)
  + (0 + sl_ | id)
  + (0 + log_sl_ | id)
  + Shrubs_s
  + SqrtDistanceToWater_s
  + Trees_s
  + Water_s
  + HumansBuff5000_s
  + (0 + Shrubs_s | id)
  + (0 + SqrtDistanceToWater_s | id)
  + (0 + Trees_s | id)
  + (0 + Water_s | id)
  + (0 + HumansBuff5000_s | id)
  + sl_:inactive
  + sl_:Water_s
  + cos_ta_:log_sl_
  + cos_ta_:SqrtDistanceToWater_s
  + sl_:Trees_s
  + cos_ta_:HumansBuff5000_s
  + sl_:Shrubs_s
  + sl_:SqrtDistanceToWater_s
  + cos_ta_:sl_
  , data = dat
)

# All scaled
mod2 <- glmm_clogit(case_ ~
  + cos_ta_s
  + sl_s
  + log_sl_s
  + (1 | step_id_)
  + (0 + cos_ta_s | id)
  + (0 + sl_s | id)
  + (0 + log_sl_s | id)
  + Shrubs_s
  + SqrtDistanceToWater_s
  + Trees_s
  + Water_s
  + HumansBuff5000_s
  + (0 + Shrubs_s | id)
  + (0 + SqrtDistanceToWater_s | id)
  + (0 + Trees_s | id)
  + (0 + Water_s | id)
  + (0 + HumansBuff5000_s | id)
  + sl_s:inactive
  + sl_s:Water_s
  + cos_ta_s:log_sl_s
  + cos_ta_s:SqrtDistanceToWater_s
  + sl_s:Trees_s
  + cos_ta_s:HumansBuff5000_s
  + sl_s:Shrubs_s
  + sl_s:SqrtDistanceToWater_s
  + cos_ta_s:sl_s
  , data = dat
)

# Compare model summaries
summary(mod1)
summary(mod2)

# Backtransform the movement metrics
variables <- c(
    "sl_s"
  , "log_sl_s"
  , "cos_ta_s"
)
means <- c(
    "sl_s"     = mean(dat$sl_)
  , "log_sl_s" = mean(dat$log_sl_)
  , "cos_ta_s" = mean(dat$cos_ta_)
)
sds <- c(
    "sl_s"     = sd(dat$sl_)
  , "log_sl_s" = sd(dat$log_sl_)
  , "cos_ta_s" = sd(dat$cos_ta_)
)
back <- backTransform(mod2, variables = variables, means = means, sds = sds)

# Compare covariates
cbind(
    summary(mod1)$coefficients$cond[, 1]
  , back
  , summary(mod2)$coefficients$cond[, 1]
)

################################################################################
#### Backtransform Full Model
################################################################################
# Convert step lengths back to meters
dat$sl_ <- dat$sl_ * 1000
dat$log_sl_ <- log(dat$sl_)
dat$sl_s <- scale(dat$sl_)
dat$log_sl_s <- scale(dat$log_sl_)

# Reload the movement model
models <- read_rds("03_Data/03_Results/99_MovementModel.rds")
best <- models$Model[[1]]

# Backtransform the movement metrics
variables <- c(
    "sl_"
  , "log_sl_"
  , "cos_ta_"
)
means <- c(
    "sl_"     = mean(dat$sl_)
  , "log_sl_" = mean(dat$log_sl_)
  , "cos_ta_" = mean(dat$cos_ta_)
)
sds <- c(
    "sl_"     = sd(dat$sl_)
  , "log_sl_" = sd(dat$log_sl_)
  , "cos_ta_" = sd(dat$cos_ta_)
)
back <- backTransform(best, variables = variables, means = means, sds = sds)

# As we will see, the back-transformed model corresponds to the model we fitted
# with step lengths in km, simply with coefficients of "sl" adjusted by 1/1000
cbind(
    BackTransformed             = back
  , EnvironmentScaledKilometers = summary(mod1)$coefficients$cond[, 1]
  , AllScaledMeters             = summary(best)$coefficients$cond[, 1]
)

# Let's store the back-transformed model coefficients
write_rds(back, "03_Data/03_Results/99_MovementModelBacktransformed.rds")
