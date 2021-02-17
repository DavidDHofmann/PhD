################################################################################
#### Step Selection Function - Model Selection
################################################################################
# Description: In this script I run forward model selection using the 4-hourly
# fixes of our dispersers

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
library(parallel)     # To run stuff in parallel
library(corrplot)     # To plot correlations
library(cowplot)      # For nice plots
library(ggpubr)       # For nice plots
library(glmmTMB)      # For modelling

################################################################################
#### Loading Data
################################################################################
# Load data (and put iSSF and TiSSF data together)
dat <- read_csv("03_Data/02_CleanData/00_General_Dispersers_POPECOL(SSF_Extracted).csv")

# Remove columns that we don't need
dat <- dat %>% select(-c(
  X1, State, x1_, x2_, y1_, y2_, t1_, t2_, dt_, absta_
))

# We want to add the log of the step speed and the cosine of the turning angle
# and the square root of "DistanceToWater"
dat <- dat %>% mutate(
    log_sl_             = log(sl_)
  , cos_ta_             = cos(ta_)
  , SqrtDistanceToWater = sqrt(DistanceToWater)
)

# Let's also move all movement metrics to the front
dat <- dat %>% select(c(
  id, step_id_, case_, sl_, log_sl_, ta_, cos_ta_, inactive, everything()
))

################################################################################
#### Scaling Data
################################################################################
# Scale the covariates. Depending on the method (iSSF, TiSSF) we'll scale the
# covariates differrently.
dat <- transform(dat
  , log_sl_                   = scale(log_sl_)
  , cos_ta_                   = scale(cos_ta_)
  , sl_                       = scale(sl_)
  , DistanceToWater           = scale(DistanceToWater)
  , SqrtDistanceToWater       = scale(SqrtDistanceToWater)
  , Water                     = scale(Water)
  , Trees                     = scale(Trees)
  , Shrubs                    = scale(Shrubs)
  , HumanInfluenceBuffer_0000 = scale(HumanInfluenceBuffer_0000)
  , HumanInfluenceBuffer_1000 = scale(HumanInfluenceBuffer_1000)
  , HumanInfluenceBuffer_2000 = scale(HumanInfluenceBuffer_2000)
  , HumanInfluenceBuffer_3000 = scale(HumanInfluenceBuffer_3000)
  , HumanInfluenceBuffer_4000 = scale(HumanInfluenceBuffer_4000)
  , HumanInfluenceBuffer_5000 = scale(HumanInfluenceBuffer_5000)
)

################################################################################
#### Function to Contrast Covariates
################################################################################
# Function to compare SSF models when using different covariates
compCovars <- function(covars = NULL, data){

  # Create comparison table
  comp <- tibble(Covars  = covars)

  # Run models
  comp$Model <- mclapply(1:nrow(comp), mc.cores = detectCores() - 1, function(x){
    glmm_clogit(
        writeForm(comp$Covars[x])
      , data = data
    )
  })

  # Extract model AICs
  comp$AIC <- sapply(comp$Model, AIC)

  # Return the tibble
  return(comp)

}

################################################################################
#### Contrast Covariates
################################################################################
# Distance to Water
comp1 <- compCovars(c("DistanceToWater", "SqrtDistanceToWater"), data = dat)
print(comp1)

# Compare different buffers for human influence
comp2 <- compCovars(c(
    "HumanInfluenceBuffer_0000"
  , "HumanInfluenceBuffer_1000"
  , "HumanInfluenceBuffer_2000"
  , "HumanInfluenceBuffer_3000"
  , "HumanInfluenceBuffer_4000"
  , "HumanInfluenceBuffer_5000"
), data = dat)

# Add some helper columns so that we can easily plot the results
comp2$Buffer <- as.numeric(substr(comp2$Covars, 22, nchar(comp2$Covars)))

# Visualize
ggplot(comp2, aes(x = Buffer, y = AIC)) +
  geom_point() +
  geom_line()
