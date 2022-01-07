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
dat1 <- read_csv("03_Data/02_CleanData/00_General_Dispersers_POPECOL(iSSF_Extracted).csv")
dat2 <- read_csv("03_Data/02_CleanData/00_General_Dispersers_POPECOL(TiSSF_Extracted).csv")
dat <- rbind(dat1, dat2)

# IMPORTANT: For simplicity I'll replace the step length by the step speed
dat$sl_ <- dat$speed_
dat$speed_ <- NULL

# Remove columns that we don't need
dat <- dat %>% select(-c(
  X1, burst, State, x1_, x2_, y1_, y2_, t1_, t2_, dt_, absta_
))

# Make sure there are no 0 step lengths
dat$sl_[dat$sl_ == 0] <- 1

# We want to add the log of the step speed and the cosine of the turning angle
dat <- dat %>% mutate(
    log_sl_ = log(sl_)
  , cos_ta_ = cos(ta_)
)

# Let's also move all movement metrics to the front
dat <- dat %>% select(c(
  id, step_id_, case_, sl_, log_sl_, ta_, cos_ta_, everything()
))

# Merge protection categories
dat$Protected <- dat$ForestReserve + dat$Protected + dat$NationalPark
dat$ForestReserve <- NULL
dat$NationalPark  <- NULL
dat$Unprotected   <- NULL

# Calculate square rooted distances
dat <- dat %>% mutate(
    SqrtDistanceToWater             = sqrt(DistanceToWater)
  , SqrtFacebook_DistanceToVillages = sqrt(Facebook_DistanceToVillages)
  , SqrtWorldpop_DistanceToVillages = sqrt(Worldpop_DistanceToVillages)
  , SqrtDistanceToRoads             = sqrt(DistanceToRoads)
  , SqrtFacebook_DistanceToHumans   = sqrt(Facebook_DistanceToHumans)
  , SqrtWorldpop_DistanceToHumans   = sqrt(Worldpop_DistanceToHumans)
)

# Compare the number of data for the two methods
dat <- dat %>%
  group_by(method) %>%
  nest()

# Compare the number of individuals in the two datasets
dat$NoDogs <- sapply(dat$data, function(x){length(unique(x$id))})
dat$NoSteps <- sapply(dat$data, function(x){sum(x$case_)})
diff(dat$NoDogs)
diff(dat$NoSteps)
print(dat)

# Let's remove the columns again
dat$NoDogs  <- NULL
dat$NoSteps <- NULL

################################################################################
#### Scaling Data
################################################################################
# Scale the covariates. Depending on the method (iSSF, TiSSF) we'll scale the
# covariates differrently.
dat <- dat %>% mutate(data = map(data, function(x){
  x %>%
    transform(
        log_sl_                             = scale(log_sl_)
      , cos_ta_                             = scale(cos_ta_)
      , sl_                                 = scale(sl_)
      , Globeland_Water                     = scale(Globeland_Water)
      , Globeland_Urban                     = scale(Globeland_Urban)
      , Globeland_Cropland                  = scale(Globeland_Cropland)
      , Globeland_Forest                    = scale(Globeland_Forest)
      , Globeland_Shrubs                    = scale(Globeland_Shrubs)
      , Globeland_Grassland                 = scale(Globeland_Grassland)
      , Copernicus_Water                    = scale(Copernicus_Water)
      , Copernicus_Urban                    = scale(Copernicus_Urban)
      , Copernicus_Cropland                 = scale(Copernicus_Cropland)
      , Copernicus_Forest                   = scale(Copernicus_Forest)
      , Copernicus_Shrubs                   = scale(Copernicus_Shrubs)
      , Copernicus_Grassland                = scale(Copernicus_Grassland)
      , DistanceToWater                     = scale(DistanceToWater)
      , SqrtDistanceToWater                 = scale(SqrtDistanceToWater)
      , Water                               = scale(Water)
      , Trees                               = scale(Trees)
      , Shrubs                              = scale(Shrubs)
      , Protected                           = scale(Protected)
      , Facebook_Villages                   = scale(Facebook_Villages)
      , Worldpop_Villages                   = scale(Worldpop_Villages)
      , Facebook_DistanceToVillages         = scale(Facebook_DistanceToVillages)
      , SqrtFacebook_DistanceToVillages     = scale(SqrtFacebook_DistanceToVillages)
      , Worldpop_DistanceToVillages         = scale(Worldpop_DistanceToVillages)
      , SqrtWorldpop_DistanceToVillages     = scale(SqrtWorldpop_DistanceToVillages)
      , DistanceToRoads                     = scale(DistanceToRoads)
      , SqrtDistanceToRoads                 = scale(SqrtDistanceToRoads)
      , Facebook_DistanceToHumans           = scale(Facebook_DistanceToHumans)
      , SqrtFacebook_DistanceToHumans       = scale(SqrtFacebook_DistanceToHumans)
      , Worldpop_DistanceToHumans           = scale(Worldpop_DistanceToHumans)
      , SqrtWorldpop_DistanceToHumans       = scale(SqrtWorldpop_DistanceToHumans)
      , Facebook_HumanDensity               = scale(Facebook_HumanDensity)
      , Worldpop_HumanDensity               = scale(Worldpop_HumanDensity)
      , Facebook_HumanInfluenceBuffer_0000  = scale(Facebook_HumanInfluenceBuffer_0000)
      , Facebook_HumanInfluenceBuffer_1000  = scale(Facebook_HumanInfluenceBuffer_1000)
      , Facebook_HumanInfluenceBuffer_2000  = scale(Facebook_HumanInfluenceBuffer_2000)
      , Facebook_HumanInfluenceBuffer_3000  = scale(Facebook_HumanInfluenceBuffer_3000)
      , Facebook_HumanInfluenceBuffer_4000  = scale(Facebook_HumanInfluenceBuffer_4000)
      , Facebook_HumanInfluenceBuffer_5000  = scale(Facebook_HumanInfluenceBuffer_5000)
      , Worldpop_HumanInfluenceBuffer_0000  = scale(Worldpop_HumanInfluenceBuffer_0000)
      , Worldpop_HumanInfluenceBuffer_1000  = scale(Worldpop_HumanInfluenceBuffer_1000)
      , Worldpop_HumanInfluenceBuffer_2000  = scale(Worldpop_HumanInfluenceBuffer_2000)
      , Worldpop_HumanInfluenceBuffer_3000  = scale(Worldpop_HumanInfluenceBuffer_3000)
      , Worldpop_HumanInfluenceBuffer_4000  = scale(Worldpop_HumanInfluenceBuffer_4000)
      , Worldpop_HumanInfluenceBuffer_5000  = scale(Worldpop_HumanInfluenceBuffer_5000)
    )
}))

################################################################################
#### Investigate Correlation Among Covariates (iSSF)
################################################################################
# Investigate correlations among covariates
correlations <- dat$data[[1]] %>%
  select(., c(Globeland_Water:ncol(.)), - RoadCrossing) %>%
  cor(., use = "pairwise.complete.obs")

# Visualize
corrplot(correlations, type = "upper")
corrplot(correlations, method = "number", type = "upper")

# We only need to look at the upper triangle of the derived matrix
correlations[upper.tri(correlations, diag = TRUE)] <- NA

# Now look at the correlation matrix
correlations

# Identify the covariates with a correlation above a specific value (some say
# that correlation should not lie above 0.6, others say that 0.7 is fine)
mat <- abs(correlations) > 0.6

# Prepare a dataframe that shows the covariates with a correlation above 0.6
df <- mat %>%
  as.table() %>%
  as.data.frame(., stringsAsFactors = FALSE)

# Remove the NAs
correlated <- na.omit(df[df$Freq, ])

# Add the correlation values to the table
correlated$Corr <- na.omit(correlations[mat])

# Look at the result
correlated

################################################################################
#### Investigate Correlation Among Covariates (TiSSF)
################################################################################
# Investigate correlations among covariates
correlations <- dat$data[[2]] %>%
  select(., c(Globeland_Water:ncol(.)), - RoadCrossing) %>%
  cor(., use = "pairwise.complete.obs")

# Visualize
corrplot(correlations, type = "upper")
corrplot(correlations, method = "number", type = "upper")

# We only need to look at the upper triangle of the derived matrix
correlations[upper.tri(correlations, diag = TRUE)] <- NA

# Now look at the correlation matrix
correlations

# Identify the covariates with a correlation above a specific value (some say
# that correlation should not lie above 0.6, others say that 0.7 is fine)
mat <- abs(correlations) > 0.6

# Prepare a dataframe that shows the covariates with a correlation above 0.6
df <- mat %>%
  as.table() %>%
  as.data.frame(., stringsAsFactors = FALSE)

# Remove the NAs
correlated <- na.omit(df[df$Freq, ])

# Add the correlation values to the table
correlated$Corr <- na.omit(correlations[mat])

# Look at the result
correlated

################################################################################
#### Function to Contrast Covariates
################################################################################
# Function to compare SSF models when using different covariates
compCovars <- function(covars = NULL){

  # Create comparison table
  comp <- tibble(
    expand_grid(
        Method  = c("iSSF", "TiSSF")
      , Covars  = covars
    )
  )

  # Run models
  comp$Model <- mclapply(1:nrow(comp), mc.cores = detectCores() - 1, function(x){
    glmm_clogit(
        writeForm(comp$Covars[x])
      , data = dat$data[[which(dat$method == comp$Method[x])]]
    )
  })

  # Extract model AICs
  comp$AIC <- sapply(comp$Model, AIC)

  # Return the tibble
  return(comp)

}

################################################################################
#### Contrast Covariates: DistanceTo vs. SqrtDistanceTo
################################################################################
# Distance to Roads
comp1 <- compCovars(c("DistanceToRoads", "SqrtDistanceToRoads"))
print(comp1)

# Distance to Water
comp2 <- compCovars(c("DistanceToWater", "SqrtDistanceToWater"))
print(comp2)

# Distance to Humans (Facebook)
comp3 <- compCovars(c("Facebook_DistanceToHumans", "SqrtFacebook_DistanceToHumans"))
print(comp3)

# Distance to Humans (Worldpop)
comp4 <- compCovars(c("Worldpop_DistanceToHumans", "SqrtWorldpop_DistanceToHumans"))
print(comp4)

# Distance to Villages (Facebook)
comp5 <- compCovars(c("Facebook_DistanceToVillages", "SqrtFacebook_DistanceToVillages"))
print(comp5)

# Distance to Villages (Worldpop)
comp6 <- compCovars(c("Worldpop_DistanceToVillages", "SqrtWorldpop_DistanceToVillages"))
print(comp6)

################################################################################
#### Compare Covariates: LandCover (Globeland, Copernicus, MODIS)
################################################################################
# Compare different land cover datasets (using iSSF data)
mod1 <- glmm_clogit(writeForm(c("Globeland_Water", "Globeland_Shrubs",
  "Globeland_Grassland", "Globeland_Forest")), data = dat$data[[1]])
mod2 <- glmm_clogit(writeForm(c("Copernicus_Water", "Copernicus_Shrubs",
  "Copernicus_Grassland", "Copernicus_Forest")), data = dat$data[[1]])
mod3 <- glmm_clogit(writeForm(c("Water", "Shrubs", "Trees")), data = dat$data[[1]])
summary(mod1)
summary(mod2)
summary(mod3)
AIC(mod1, mod2, mod3)

# Compare different land cover datasets (using TiSSF data)
mod1 <- glmm_clogit(writeForm(c("Globeland_Water", "Globeland_Shrubs",
  "Globeland_Grassland", "Globeland_Forest")), data = dat$data[[2]])
mod2 <- glmm_clogit(writeForm(c("Copernicus_Water", "Copernicus_Shrubs",
  "Copernicus_Grassland", "Copernicus_Forest")), data = dat$data[[2]])
mod3 <- glmm_clogit(writeForm(c("Water", "Shrubs", "Trees")), data = dat$data[[2]])
summary(mod1)
summary(mod2)
summary(mod3)
AIC(mod1, mod2, mod3)

# Let's compare the water layers in unimodal models
comp7 <- compCovars(c("Globeland_Water", "Copernicus_Water", "Water"))
print(comp7)

# Let's compare the shrubs layers in unimodal models
comp8 <- compCovars(c("Globeland_Shrubs", "Copernicus_Shrubs", "Shrubs"))
print(comp8)

# Let's compare the forest layers in unimodal models
comp9 <- compCovars(c("Globeland_Forest", "Copernicus_Forest", "Trees"))
print(comp9)

################################################################################
#### Contrast Covariates: Humans Buffered
################################################################################
# Compare different buffers for human influence
comp10 <- compCovars(c(
    "Facebook_HumanInfluenceBuffer_0000"
  , "Facebook_HumanInfluenceBuffer_1000"
  , "Facebook_HumanInfluenceBuffer_2000"
  , "Facebook_HumanInfluenceBuffer_3000"
  , "Facebook_HumanInfluenceBuffer_4000"
  , "Facebook_HumanInfluenceBuffer_5000"
  , "Worldpop_HumanInfluenceBuffer_0000"
  , "Worldpop_HumanInfluenceBuffer_1000"
  , "Worldpop_HumanInfluenceBuffer_2000"
  , "Worldpop_HumanInfluenceBuffer_3000"
  , "Worldpop_HumanInfluenceBuffer_4000"
  , "Worldpop_HumanInfluenceBuffer_5000"
))

# Add some helper columns so that we can easily plot the results
comp10$Source <- substr(comp10$Covars, 1, 8)
comp10$Covariates <- substr(comp10$Covars, 10, nchar(comp10$Covars))
comp10$Buffer <- as.numeric(substr(comp10$Covars, 31, nchar(comp10$Covars)))

# Visualize
ggplot(comp10, aes(x = Buffer, y = AIC, col = Source)) +
  geom_point() +
  geom_line() +
  facet_wrap("Method", scales = "free")

# Prepare a more "cheesy" version of this plot
ggplot(comp10, aes(x = Buffer, y = AIC, col = Source, lty = Source)) +
  geom_point() +
  geom_line() +
  facet_wrap("Method", scales = "free", ncol = 1) +
  theme_cowplot() +
  scale_color_manual(values = c("black", "gray"))

################################################################################
#### Contrast Covariates: Distance To Humans, Human Density, Villages
################################################################################
# Distance to Humans
comp11 <- compCovars(c(
    "SqrtFacebook_DistanceToHumans"
  , "SqrtWorldpop_DistanceToHumans"
))
print(comp11)

# Human Density
comp12 <- compCovars(c(
    "Facebook_HumanDensity"
  , "Worldpop_HumanDensity"
))
print(comp12)

# Villages
comp13 <- compCovars(c(
    "Facebook_Villages"
  , "Worldpop_Villages"
))
print(comp13)

################################################################################
#### Contrast Distance To Humans, Human Density in a more complex model
################################################################################
# Run different models similar to Hofmann et al. 2021
mod1 <- glmm_clogit(writeForm(c("Water", "SqrtDistanceToWater", "Trees",
  "Shrubs", "Facebook_HumanInfluenceBuffer_5000")), data = dat$data[[2]])
mod2 <- glmm_clogit(writeForm(c("Water", "SqrtDistanceToWater", "Trees",
  "Shrubs", "Worldpop_HumanInfluenceBuffer_5000")), data = dat$data[[2]])
mod3 <- glmm_clogit(writeForm(c("Water", "SqrtDistanceToWater", "Trees",
  "Shrubs", "SqrtFacebook_DistanceToHumans")), data = dat$data[[2]])
mod4 <- glmm_clogit(writeForm(c("Water", "SqrtDistanceToWater", "Trees",
  "Shrubs", "SqrtWorldpop_DistanceToHumans")), data = dat$data[[2]])

# Check results
summary(mod1)
summary(mod2)
summary(mod3)
summary(mod4)

# Compare AIC
AIC(mod1, mod2, mod3, mod4)

# Visualize models
showCoeffs(getCoeffs(mod1)[-1, ])
showCoeffs(getCoeffs(mod2)[-1, ])
showCoeffs(getCoeffs(mod3)[-1, ])
showCoeffs(getCoeffs(mod4)[-1, ])

################################################################################
#### Store Comparisons
################################################################################
# Put all comparisons together
comp10 <- select(comp10, Method, Covars, Model, AIC)
comp <- rbind(
    comp1, comp2, comp3, comp4, comp5, comp6, comp7
  , comp8, comp9, comp10, comp11, comp12, comp13
)

# Check object size
format(object.size(comp), units = "Mb")

# Remove the models
comp$Model <- NULL
format(object.size(comp), units = "Mb")

# Store the comparisons
write_csv(comp, "03_Data/03_Results/99_CovariateComparisons.csv")
