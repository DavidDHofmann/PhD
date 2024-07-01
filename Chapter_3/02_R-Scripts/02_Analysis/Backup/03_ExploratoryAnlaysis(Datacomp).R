################################################################################
#### Comparing New to Old Data
################################################################################
# Description: In this script, I wanted to compare the data used in chapter 1 to
# the data from chapter 2 to better understand the differences in resulting
# model estimates and to pinpoint potential errors. The different data types are
# dat_osoc: old steps (os) and old covariates (oc)
# dat_osnc: old steps (os) and new covariates (nc)
# dat_nsnc: new steps (ns) and new covariates (nc)
# -> dat_nsncind: regular iSSF (i), new dispersers
# -> dat_nsnctnd: time varying iSSF (t), new dispersers
# -> dat_nsnciod: regular iSSF (i), old dispersers (od)
# -> dat_nsnctod: time varying iSSF (t), old dispersers (od)


# Clear R's brain
rm(list = ls())

# Surpress scientific notation
options(scipen = 999)

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)      # For data wrangling
library(davidoff)       # Custom functions
library(corrplot)       # To plot correlations
library(cowplot)        # For nice plots
library(ggpubr)         # For nice plots
library(glmmTMB)        # For modelling
library(gridExtra)      # For laying out multiple plots
library(survival)       # For conditional logistic regression
library(parallel)       # To allow parallelization
library(broom)          # To clean model results
library(TwoStepCLogit)  # Fieberg's two-step approach
library(ggdark)         # Dark ggplot themes

# Identify number of cores we can use
cores <- detectCores()

################################################################################
#### Loading Data
################################################################################
# Load the data
dat_osoc <- read_csv("/home/david/ownCloud/University/15. PhD/Chapter_0/03_Data/02_CleanData/00_General_Dispersers_Popecol(SSF4Hours).csv")
dat_osnc <- read_csv("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/00_General_Dispersers_POPECOL(Extracted).csv")
dat_nsncind <- read_csv("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/00_General_Dispersers_POPECOL(iSSF_Extracted).csv")
dat_nsnctnd <- read_csv("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/00_General_Dispersers_POPECOL(TiSSF_Extracted).csv")

# Replace the step length with the step speed in the last dataset
dat_nsnctnd$sl_ <- dat_nsnctnd$speed_

# There are step lengths == 0, make them positive
dat_osoc$sl_[dat_osoc$sl_ == 0] <- 1
dat_osnc$sl_[dat_osnc$sl_ == 0] <- 1
dat_nsncind$sl_[dat_nsncind$sl_ == 0] <- 1
dat_nsnctnd$sl_[dat_nsnctnd$sl_ == 0] <- 1

# Let's create a subset of the new data where we only keep those steps that are
# also present in the old data
index <- do.call(paste0, dat_nsncind[, c("id", "t1_")]) %in%
  do.call(paste0, dat_osoc[, c("id", "t1_")])
dat_nsnciod <- dat_nsncind[index, ]

index <- do.call(paste0, dat_nsnctnd[, c("id", "t1_")]) %in%
  do.call(paste0, dat_osoc[, c("id", "t1_")])
dat_nsnctod <- dat_nsnctnd[index, ]

# Keep only the columns we are actually interested in
dat_osoc <- select(dat_osoc, c(id, step_id_, case_, sl_, ta_, absta_, Water,
  DistanceToWater, Trees, Shrubs, HumansBuff5000))
dat_osnc <- select(dat_osnc, c(id, step_id_, case_, sl_, ta_, absta_, Water,
  DistanceToWater, Trees, Shrubs, HumansBuff5000 =
  Facebook_HumanInfluenceBuffer_5000))
dat_nsncind <- select(dat_nsncind, c(id, step_id_, case_, sl_, ta_, absta_, Water,
  DistanceToWater, Trees, Shrubs, HumansBuff5000 =
  Facebook_HumanInfluenceBuffer_5000))
dat_nsnctnd <- select(dat_nsnctnd, c(id, step_id_, case_, sl_, ta_, absta_, Water,
  DistanceToWater, Trees, Shrubs, HumansBuff5000 =
  Facebook_HumanInfluenceBuffer_5000))
dat_nsnciod <- select(dat_nsnciod, c(id, step_id_, case_, sl_, ta_, absta_, Water,
  DistanceToWater, Trees, Shrubs, HumansBuff5000 =
  Facebook_HumanInfluenceBuffer_5000))
dat_nsnctod <- select(dat_nsnctod, c(id, step_id_, case_, sl_, ta_, absta_, Water,
  DistanceToWater, Trees, Shrubs, HumansBuff5000 =
  Facebook_HumanInfluenceBuffer_5000))

# Make case_ logical
dat_osoc$case_ <- as.logical(dat_osoc$case_)
dat_osnc$case_ <- as.logical(dat_osnc$case_)
dat_nsncind$case_ <- as.logical(dat_nsncind$case_)
dat_nsnctnd$case_ <- as.logical(dat_nsnctnd$case_)
dat_nsnciod$case_ <- as.logical(dat_nsnciod$case_)
dat_nsnctod$case_ <- as.logical(dat_nsnctod$case_)

# Create an indicator column for the datasource
dat_osoc$datacat <- "osoc"
dat_osnc$datacat <- "osnc"
dat_nsncind$datacat <- "nsnc_issf_nd"
dat_nsnctnd$datacat <- "nsnc_tissf_nd"
dat_nsnciod$datacat <- "nsnc_issf_od"
dat_nsnctod$datacat <- "nsnc_tissf_od"

# Put them together
dat <- rbind(dat_osoc, dat_osnc, dat_nsncind, dat_nsnctnd, dat_nsnciod, dat_nsnctod)

################################################################################
#### Check variation in covariates
################################################################################
# Calculate variance in covariate values for each datacategory
variation <- dat %>%
  gather(key = Covariate, value = Value, 4:11) %>%
  group_by(id, Covariate, datacat) %>%
  summarize(Variance = var(Value))

# Let's see for which indivs which covariates don't vary
variation[variation$Variance == 0, ]

# Store indivs for whom covariates don't vary
novary <- variation[variation$Variance == 0, ]

################################################################################
#### Compare Data Descriptively
################################################################################
# Compare number of observations (case steps)
dat %>%
  subset(case_) %>%
  group_by(datacat) %>%
  summarize(n = n())

# Compare number of observations by individual
dat %>%
  subset(case_) %>%
  group_by(datacat, id) %>%
  summarize(n = n()) %>%
  spread(key = datacat, value = n)

# Compare step lengths: We are quite close to the true mean, yet the sampled
# data generates a median that is quite a bit lower than the observed median.
dat %>%
  group_by(datacat, case_) %>%
  summarize(
      sl_mean   = mean(sl_)
    , sl_sd     = sd(sl_)
    , sl_median = median(sl_)
  )

# Compare cosine of the turning angles: True steps are much more directional
# than the sampled steps (i.e. they have a higher cos(ta) value)
dat %>%
  group_by(datacat, case_) %>%
  summarize(
      ta_mean   = mean(cos(ta_))
    , ta_sd     = sd(cos(ta_))
    , ta_median = median(cos(ta_))
  )

# Compare water cover: It looks like in the new data the difference in water
# cover between case & control steps has become much smaller. This could suggest
# that more steps are further away from water.
dat %>%
  group_by(datacat, case_) %>%
  summarize(
      water_mean   = mean(Water)
    , water_sd     = sd(Water)
    , water_median = median(Water)
  )

# Compare distance to water: We see that steps are much further away from water
# compared to earlier. This may be caused by the pronounced flood in 2019.
dat %>%
  group_by(datacat, case_) %>%
  summarize(
      dist2water_mean   = mean(DistanceToWater)
    , dist2water_sd     = sd(DistanceToWater)
    , dist2water_median = median(DistanceToWater)
  )

# Compare tree cover: On average, tree cover is slightly lower in comparison to
# old data
dat %>%
  group_by(datacat, case_) %>%
  summarize(
      trees_mean   = mean(Trees)
    , trees_sd     = sd(Trees)
    , trees_median = median(Trees)
  )

# Compare shrub cover: Almost identical between new and old data
dat %>%
  group_by(datacat, case_) %>%
  summarize(
      shrubs_mean   = mean(Shrubs)
    , shrubs_sd     = sd(Shrubs)
    , shrubs_median = median(Shrubs)
  )

# Compare HumansBuff: Here's where it gets weird. It looks like case_ steps tend
# to be in higher human influence than control steps. Whuuuut?
dat %>%
  group_by(datacat, case_) %>%
  summarize(
      humans_mean   = mean(HumansBuff5000)
    , humans_sd     = sd(HumansBuff5000)
    , humans_median = median(HumansBuff5000)
  )

################################################################################
#### Scale Covariates for Modelling
################################################################################
# Scale covariates
dat_osoc <- dat_osoc %>% mutate(
    cos_ta_         = scale(cos(ta_))
  , log_sl_         = scale(log(sl_))
  , sl_             = scale(sl_)
  , Water           = scale(Water)
  , DistanceToWater = scale(sqrt(DistanceToWater))
  , HumansBuff5000  = scale(HumansBuff5000)
  , Shrubs          = scale(Shrubs)
  , Trees           = scale(Trees)
)
dat_osnc <- dat_osnc %>% mutate(
    cos_ta_         = scale(cos(ta_))
  , log_sl_         = scale(log(sl_))
  , sl_             = scale(sl_)
  , Water           = scale(Water)
  , DistanceToWater = scale(sqrt(DistanceToWater))
  , HumansBuff5000  = scale(HumansBuff5000)
  , Shrubs          = scale(Shrubs)
  , Trees           = scale(Trees)
)
dat_nsncind <- dat_nsncind %>% mutate(
    cos_ta_         = scale(cos(ta_))
  , log_sl_         = scale(log(sl_))
  , sl_             = scale(sl_)
  , Water           = scale(Water)
  , DistanceToWater = scale(sqrt(DistanceToWater))
  , HumansBuff5000  = scale(HumansBuff5000)
  , Shrubs          = scale(Shrubs)
  , Trees           = scale(Trees)
)
dat_nsnctnd <- dat_nsnctnd %>% mutate(
    cos_ta_         = scale(cos(ta_))
  , log_sl_         = scale(log(sl_))
  , sl_             = scale(sl_)
  , Water           = scale(Water)
  , DistanceToWater = scale(sqrt(DistanceToWater))
  , HumansBuff5000  = scale(HumansBuff5000)
  , Shrubs          = scale(Shrubs)
  , Trees           = scale(Trees)
)
dat_nsnciod <- dat_nsnciod %>% mutate(
    cos_ta_         = scale(cos(ta_))
  , log_sl_         = scale(log(sl_))
  , sl_             = scale(sl_)
  , Water           = scale(Water)
  , DistanceToWater = scale(sqrt(DistanceToWater))
  , HumansBuff5000  = scale(HumansBuff5000)
  , Shrubs          = scale(Shrubs)
  , Trees           = scale(Trees)
)
dat_nsnctod <- dat_nsnctod %>% mutate(
    cos_ta_         = scale(cos(ta_))
  , log_sl_         = scale(log(sl_))
  , sl_             = scale(sl_)
  , Water           = scale(Water)
  , DistanceToWater = scale(sqrt(DistanceToWater))
  , HumansBuff5000  = scale(HumansBuff5000)
  , Shrubs          = scale(Shrubs)
  , Trees           = scale(Trees)
)

# Combine the scaled data
dat <- rbind(dat_osoc, dat_osnc, dat_nsncind, dat_nsnctnd, dat_nsnciod, dat_nsnctod)

# Make the datacategory factorial
dat$datacat <- factor(dat$datacat, levels = c("osoc", "osnc", "nsnc_issf_od", "nsnc_tissf_od", "nsnc_issf_nd", "nsnc_tissf_nd"))

# Remove Pula as there are only two fixes for her
dat <- subset(dat, id != "Pula")

################################################################################
#### Mixed Effects Models
################################################################################
# Nest data for modelling
dat_nest <- dat %>%
  group_by(datacat) %>%
  nest()

# Run a basic model
covars <- c("Water")
dat_nest$Models <- mclapply(dat_nest$data, mc.cores = cores - 1, function(x){
  glmm_clogit(writeForm(covars), data = x)
})

# Extract coefficients
dat_nest$Coeffs <- mclapply(dat_nest$Models, mc.cores = cores - 1, function(x){
  getCoeffs(x)
})

# Plot them
dat_nest %>%
  select(datacat, Coeffs) %>%
  unnest() %>%
  subset(Covariate != "(Intercept)") %>%
  ggplot(aes(x = Coefficient, y = Covariate, col = datacat)) +
    geom_point(position = position_dodge(-0.4)) +
    geom_errorbarh(
      aes(
          xmin = Coefficient - 1.96 * SE
        , xmax = Coefficient + 1.96 * SE
        , height = 0.2), position = position_dodge(-0.4)
      ) +
    geom_vline(xintercept = 0, linetype = "dashed", col = "darkgray")

# Run model from chapter 1
covars <- c("Water", "DistanceToWater", "HumansBuff5000", "Shrubs", "Trees")
dat_nest$Models <- mclapply(dat_nest$data, mc.cores = cores - 1, function(x){
  glmm_clogit(writeForm(covars), data = x)
})

# Extract coefficients
dat_nest$Coeffs <- mclapply(dat_nest$Models, mc.cores = cores - 1, function(x){
  getCoeffs(x)
})

# Plot them
dat_nest %>%
  select(datacat, Coeffs) %>%
  unnest() %>%
  subset(Covariate != "(Intercept)") %>%
  ggplot(aes(x = Coefficient, y = Covariate, col = datacat)) +
    geom_point(position = position_dodge(-0.6)) +
    geom_errorbarh(
      aes(
          xmin = Coefficient - 1.96 * SE
        , xmax = Coefficient + 1.96 * SE
        , height = 0.2), position = position_dodge(-0.6)
      ) +
    geom_vline(xintercept = 0, linetype = "dashed", col = "darkgray")

# We can see that adding the new dispersers (nd) to our data results in
# coefficients that are closer to 0. Hence, I'll try to identify the individuals
# that cause this.
################################################################################
#### Individual Models
################################################################################
# We cant fit the same model to all individuals because for some there are
# non-varying covariates.
omit <- novary %>%
  group_by(id, datacat) %>%
  mutate(omit = paste0(Covariate, collapse = ", ")) %>%
  select(-c(Covariate, Variance)) %>%
  distinct()

# Nest data by data category and individual
dat_nest <- dat %>%
  group_by(id, datacat) %>%
  nest()

# Join the variables we need to omit
dat_nest <- left_join(dat_nest, omit, by = c("id", "datacat"))

# Write the full model
covars <- c("Water", "DistanceToWater", "HumansBuff5000", "Shrubs", "Trees")

# For each row in the nested table, keep only the covariates that should not be
# removed
dat_nest$Covariates <- lapply(dat_nest$omit, function(x){
  covars[!covars %in% x]
})

# Function to write the model formula
writeForm2 <- function(covars){
  form <- as.formula(
    paste0(
        "case_ ~ + sl_ + log_sl_ + cos_ta_ + "
      , paste0(covars, collapse = " + ")
      , " + strata(step_id_)"
    )
  )
  return(form)
}

# Run some models for each individual
dat_nest$Models <- lapply(1:nrow(dat_nest), function(x){
  clogit(writeForm2(dat_nest$Covariates[[x]]), data = dat_nest$data[[x]])
})

# Extract coefficients
dat_nest$Coeffs <- lapply(dat_nest$Models, function(x){
  tidy(x)
})

# Plot results for each individual
p <- dat_nest %>%
  select(id, datacat, Coeffs) %>%
  unnest(Coeffs) %>%
  subset(abs(estimate) < 5) %>%
  ggplot(aes(x = estimate, y = term, col = datacat)) +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbarh(aes(xmin = estimate - 1.96 * std.error, xmax = estimate + 1.96 * std.error, height = 0.2), position = position_dodge(width = 0.5)) +
    facet_wrap("id", scales = "free")
# ggsave("test.png", scale = 2.5)

# Summarize the selection coefficients
se <- function(x){
  sd(x) / sqrt(length(x))
}

# Unnest and visualize
dat_nest %>%
  select(id, datacat, Coeffs) %>%
  unnest(Coeffs) %>%
  subset(abs(estimate) < 5) %>%
  filter(term != "(Intercept)") %>%
  group_by(term, datacat) %>%
  mutate(weight = (1 / std.error ** 2) / sum((1 / std.error ** 2))) %>%
  # summarize(test = sum(weight))
  summarize(mean = sum(weight * estimate), se = sqrt((sum(weight * (estimate - mean) ** 2)) / (n() - 1))) %>%
  # summarize(mean = mean(estimate), se = se(estimate)) %>%
  ggplot(aes(y = term, x = mean, col = datacat)) +
    geom_errorbarh(aes(xmin = mean - 1.96 * se, xmax = mean + 1.96 * se, height = 0.2), position = position_dodge(width = 0.5)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgray") +
    geom_point(position = position_dodge(width = 0.5))
