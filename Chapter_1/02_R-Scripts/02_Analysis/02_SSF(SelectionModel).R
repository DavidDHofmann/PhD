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
#### Comparison to Old Data
################################################################################
# I want to compare the current dataset to the dataset used in Hofmann et al
# 2020. Let's thus load the two
old <- read_csv("/home/david/ownCloud/University/15. PhD/Chapter_0/03_Data/02_CleanData/00_General_Dispersers_Popecol(SSF4Hours).csv")
new <- read_csv("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/00_General_Dispersers_POPECOL(iSSF_Extracted).csv")

# Subset to case_ steps only
old <- subset(old, case_)
new <- subset(new, case_ == 1)

# Remove undesired columns
old <- select(old, c(id, x1_, x2_, y1_, y2_, t1_, t2_))
new <- select(new, c(id, x1_, x2_, y1_, y2_, t1_, t2_))

# Let's only compare the data of the common individuals
new <- subset(new, id %in% old$id)
old <- subset(old, id %in% new$id)

# Create new column indicating dataset
old$set <- "old"
new$set <- "new"

# Combine them
comb <- rbind(old, new)

# Compare number of observations
table(comb$set)
diff(as.data.frame(table(comb$set))$Freq)

# Compare fixes per individual
table(comb$id, comb$set)

# Join the dataframes
joined <- full_join(old, new, by = c("id", "t1_", "t2_"))

# Sort the columns
joined <- select(joined, c(
    id
  , t1_
  , t2_
  , x1_old = x1_.x
  , x1_new = x1_.y
  , x2_old = x2_.x
  , x2_new = x2_.y
  , y1_old = y1_.x
  , y1_new = y1_.y
  , y2_old = y2_.x
  , y2_new = y2_.y
  , everything()
))

# Let's also do an antijoin
anti_join(old, new, by = c("id", "t1_", "t2_"))
anti_join(new, old, by = c("id", "t1_", "t2_"))

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
  X1, burst, State, x1_, x2_, y1_, y2_, t1_, t2_, dt_, drctn_p, absta_
))

# We want to add the log of the step speed and the cosine of the turning angle
dat <- dat %>% mutate(
    log_sl_ = log(sl_)
  , cos_ta_ = cos(ta_)
)

# Let's also move all movement metrics to the front
dat <- dat %>% select(c(
  id, step_id_, case_, sl_, log_sl_, ta_, cos_ta_, everything()
))

# Simplify protection categories
dat$Protected <- dat$ForestReserve + dat$Protected + dat$NationalPark
dat$ForestReserve <- NULL
dat$NationalPark  <- NULL
dat$Unprotected   <- NULL

# We're not going to consider human influences buffered beyond 5000 meters as
# this would be simply too much
dat <- dat %>% select(-c(
    contains(c("7500", "10000", "12500", "15000", "17500", "20000"))
))

# Calculate square rooted distances
dat <- dat %>% mutate(
    SqrtDistanceToWater             = sqrt(DistanceToWater)
  , SqrtFacebook_DistanceToVillages = sqrt(Facebook_DistanceToVillages)
  , SqrtWorldpop_DistanceToVillages = sqrt(Worldpop_DistanceToVillages)
  , SqrtDistanceToRoads             = sqrt(DistanceToRoads)
  , SqrtFacebook_DistanceToHumans   = sqrt(Facebook_DistanceToHumans)
  , SqrtWorldpop_DistanceToHumans   = sqrt(Worldpop_DistanceToHumans)
)

# Nest the data
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
      , Facebook_HumanInfluenceBuffer_2500  = scale(Facebook_HumanInfluenceBuffer_2500)
      , Facebook_HumanInfluenceBuffer_5000  = scale(Facebook_HumanInfluenceBuffer_5000)
      , Worldpop_HumanInfluenceBuffer_0000  = scale(Worldpop_HumanInfluenceBuffer_0000)
      , Worldpop_HumanInfluenceBuffer_2500  = scale(Worldpop_HumanInfluenceBuffer_2500)
      , Worldpop_HumanInfluenceBuffer_5000  = scale(Worldpop_HumanInfluenceBuffer_5000)
    )
}))

# Extract the scaling parameters from the two datasets
scales <- lapply(dat$data, function(x){

  # Prepare dataframe into which we store the values
  scaling <- data.frame(
      ColumnName  = names(x)
    , Center      = NA
    , Scale       = NA
  )

  # Extract the values from the data
  for (i in 1:ncol(x)){
    center <- attr(x[, i], "scaled:center")
    scale <- attr(x[, i], "scaled:scale")
    scaling$Center[i] <- ifelse(is.null(center), NA, center)
    scaling$Scale[i] <- ifelse(is.null(scale), NA, scale)
  }

  # Return the dataframe
  return(scaling)
})

# Look at the resulting table
names(scales) <- c("iSSF", "TiSSF")
scales

# Store the object to an RDS
write_rds(scales, "03_Data/03_Results/99_Scaling.rds")

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
  , "Facebook_HumanInfluenceBuffer_2500"
  , "Facebook_HumanInfluenceBuffer_5000"
  , "Worldpop_HumanInfluenceBuffer_0000"
  , "Worldpop_HumanInfluenceBuffer_2500"
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
#### Contrast Covariates: Villages, Human Density, Human Buffer, Dist To Humans
################################################################################
# We still have some covariates that are correlated. Let's compare them
comp14 <- compCovars(c(
    "Facebook_Villages"
  , "Facebook_HumanDensity"
  , "Facebook_HumanInfluenceBuffer_5000"
  , "SqrtFacebook_DistanceToVillages"
))
print(comp14)

################################################################################
#### Keep Only Desired Covariates
################################################################################
# Overall, the comparisons show that the covariates used in Hofmann et al. 2020
# were mostly adquate and that only minor improvements may be gained by using
# alternative data sources. For continuity, I'll thus keep only the covariates
# used in Hofmann et al. 2020.
dat <- dat %>% unnest()

# Keep only "SqrtDistanceTo", remove "DistanceTo"
dat <- dat %>%
  select(-c(
      DistanceToRoads
    , DistanceToWater
    , Facebook_DistanceToHumans
    , Worldpop_DistanceToHumans
    , Facebook_DistanceToVillages
    , Worldpop_DistanceToVillages
  ))

# Keep only MODISD data, remove Globeland and Copernicus land cover
dat <- dat %>% select(-contains(c("Copernicus", "Globeland")))

# Keep only 5km human influence buffer
dat <- dat %>% select(-contains(c("Worldpop_HumanInfluence", "0000", "2500")))

# Keep only Facebook human density data, remove worldpop
dat <- dat %>% select(-contains(c("Worldpop")))

# Keep only buffered human influence, remove "DistanceToHumans" etc.
dat <- dat %>% select(-c(
    Facebook_Villages
  , Facebook_HumanDensity
  , SqrtFacebook_DistanceToHumans
  , SqrtFacebook_DistanceToVillages
))

# Rename the remaining covariates nicely
dat$HumansBuff5000 <- dat$Facebook_HumanInfluenceBuffer_5000
dat$Facebook_HumanInfluenceBuffer_5000 <- NULL

# Nest data again
dat <- dat %>% group_by(method) %>% nest()

################################################################################
#### TESTING
################################################################################
test <- read_csv("/home/david/ownCloud/University/15. PhD/Chapter_0/03_Data/02_CleanData/00_General_Dispersers_Popecol(SSF4Hours).csv")
disps <- unique(test$id)
sub <- dat$data[[1]]
subset(sub, t1_ %in% test$t1_)
mod1 <- glmm_clogit(
  writeForm(
    c("Water", "SqrtDistanceToWater", "Trees", "Shrubs", "HumansBuff5000")
  )
  , data = sub
)
summary(mod1)
mod2 <- glmm_clogit(
  writeForm(
    c("Water", "SqrtDistanceToWater", "Trees", "Shrubs", "HumansBuff5000")
  )
  , data = subset(sub, id != "Chiounard")
)
mod3 <- glmm_clogit(
  writeForm(
    c("Water", "SqrtDistanceToWater", "Trees", "Shrubs", "HumansBuff5000")
  )
  , data = subset(sub, id %in% disps)
)
summary(mod1)
summary(mod2)
summary(mod3)

################################################################################
#### Forward Model Selection for Main Effects
################################################################################
# We now want to run forward model selection. That is, we iteratively increase
# the complexity of our models and keep track of the changes in the AIC.
# Ultimately, we will use the model with the lowest AIC as our most parsimonious
# model. To start, we need to know the set of covariates that we can choose
# from. So let's first identify all possible covariates.
covars <- c(
    "Water"
  , "SqrtDistanceToWater"
  , "Trees"
  , "Shrubs"
  , "Protected"
  , "HumansBuff5000"
  , "SqrtDistanceToRoads"
  , "RoadCrossing"
)

# Initate vector of selected covars
selected <- c()

# Initiate vector into which we store the corresponding models
model_sel <- list()

# Prepare loop counter
i <- 1

# Run loop until no more covariates are left for selection
while(length(covars) > 0){

  # Prepare model formula
  forformulas <- lapply(covars, function(x){
    c(selected, x)
  })

  # Prepare a vector that contains all covariates. This will allow us to easily
  # spot the covariates of a specific model
  covariates <- forformulas %>%
    do.call(rbind, .) %>%
    as.data.frame(.) %>%
    unite(., col = "Covariates", sep = ", ")

  # Prepare a tibble that we want to fill
  models <- tibble(
      ModelID     = 1:length(covars)
    , Covariates  = covariates
    , Formula     = lapply(forformulas, writeForm)
    , Model       = NA # Result of the model
    , AIC         = NA # AIC of the model
  )

  # Run all models in parallel
  models$Model <- mclapply(models$Formula, mc.cores = detectCores() - 1, function(x){
    glmm_clogit(x, data = dat$data[[2]])
  })

  # Extract the AIC value of each model
  models$AIC <- sapply(models$Model, AIC)

  # Identify the model with the lowest AIC and put the respective covariate
  # into the "selected" vector
  best_covar <- na.omit(covars[models$AIC == min(models$AIC, na.rm = TRUE)])
  selected <- c(selected, best_covar)

  # Update which covariates remain for selection
  covars <- covars[!(covars %in% selected)]

  # Store the model results into the list
  model_sel[[i]] <- models

  # Update the index
  i <- i + 1

  # Print how many covariates are left
  cat(length(covars), "covariates left...\n")

}

############################################################
#### Select Best Model
############################################################
# Put all relevant information of the model selection into a single tibble
models <- tibble(
    Covariates  = model_sel %>%
    lapply(., function(x){x$Covariates}) %>% do.call(rbind, .)
  , Formula     = model_sel %>%
    lapply(., function(x){x$Formula}) %>% do.call(c, .)
  , Model       = model_sel %>%
    lapply(., function(x){x$Model}) %>% do.call(c, .)
  , AIC         = model_sel %>%
    lapply(., function(x){x$AIC}) %>% do.call(c, .)
  , DeltaAIC    = NA # Difference in AIC to the "best" model
  , Weight      = NA # Weight assigned according to AIC
  , LogLik      = NA # Loglikelihood of the model
)

# Extract the LogLiks
models$LogLik <- lapply(models$Model, logLik) %>%
  do.call(rbind, .) %>%
  as.vector()

# Add a unique model id
models$ModelID <- 1:nrow(models)

# Calculate the DeltaAIC to the most parsimonious model
models$DeltaAIC <- models$AIC - min(models$AIC, na.rm = TRUE)

# Sort the tibble by the DeltaAIC
models <- arrange(models, DeltaAIC)

# Look at the best model
summary(models$Model[[1]])

# Plot the results
showCoeffs(getCoeffs(models$Model[[1]])[-1, ])

# ############################################################
# #### Forward Model Selection for Interaction Terms
# ############################################################
# # Let's now take the best model and use the selected covariates and check for
# # possible interactions that are important
# base <- models$Covariates[1, ] %>% strsplit(., split = ", ") %>% .[[1]]
#
# # Now clean the "selected" vector as we will use it to store the selected
# # interaction terms
# selected <- c()
#
# # Identify all possible 2-way interaction terms for the baseline covariates
# covars <- base %>%
#
#   # Find all possible 2-way combinations
#   combn(., 2) %>%
#
#   # Transpose resulting matrix
#   t(.) %>%
#
#   # Convert matrix to dataframe
#   as.data.frame(.) %>%
#
#   # Merge columns to create interaction terms
#   unite(., col, sep = ":") %>%
#
#   # Get rid of the dataframe format
#   .[["col"]]
#
# # Look at all interactions
# covars
#
# # We now remove all interaction terms which we don't want to consider since we
# # believe it does not make any sense to look at them
# covars <- covars[-1]
#
# # Run a loop that iteratively adds interaction terms until none are left for
# # selection
# while (length(covars) > 0){
#
#   # Write down all of the covariates for the different models
#   forformulas <- lapply(covars, function(x){
#     c(base, selected, x)
#   })
#
#   # Prepare a vector that makes the covariates of a model slightly better
#   # visible
#   covariates <- forformulas %>%
#
#     # Collapse the list
#     do.call(rbind, .) %>%
#
#     # Convert the matrix to a dataframe
#     as.data.frame(.) %>%
#
#     # Combine the columns into a vector that shows all covariates of a specific
#     # model
#     unite(., col = "Covariates", sep = ", ")
#
#   # Prepare an empty tibble to fill
#   models <- tibble(
#       ModelID     = 1:length(covars)
#     , Covariates  = covariates
#     , Formula     = lapply(covars, function(x){
#
#       # Write the model for the baseline covariates
#       formula <- writeForm(base)
#
#       # Check if we selected some interaction terms in the previous iteration
#       if (length(selected) > 0){
#
#         # If so, add them to the models as seperate terms (without random slope)
#         toadd <- selected %>%
#
#           # Transpose the matrix
#           t(.) %>%
#
#           # Convert the matrix to a dataframe
#           as.data.frame(.) %>%
#
#           # Combine the interactions using + as a seperator
#           unite(., col, sep = " + ") %>%
#
#           # Get rid of the dataframe structure
#           .[["col"]]
#
#         # Update the baseline model with the selected interaction terms
#         formula <- update(formula, paste("~ . +", toadd))
#       }
#
#       # Add the interaction terms we still need to check
#       formula <- update(formula, paste("~ . +", x))
#
#       # Return the model formula that we just created
#       return(formula)
#     })
#     , Model       = NA # Result of the model
#     , AIC         = NA # AIC of the model
#   )
#
#   # Run all models
#   models <- mutate(models, Model = map(Formula, function(x){
#     glmm_clogit(x, data = ssf)
#   }))
#
#   # Extract the AIC value of each model
#   models$AIC <- lapply(models$Model, AIC) %>%
#     do.call(rbind, .) %>%
#     as.vector()
#
#   # Identify the model with the lowest AIC and put the respective covariate
#   # into the "selected" vector
#   best_covar <- na.omit(covars[models$AIC == min(models$AIC, na.rm = TRUE)])
#   selected <- c(selected, best_covar)
#
#   # Update which covariates remain for selection
#   covars <- covars[!(covars %in% selected)]
#
#   # Store the model results into the list
#   model_sel[[i]] <- models
#
#   # Update the index
#   i <- i + 1
# }
#
# # Put all relevant information into a single tibble
# models <- tibble(
#     Covariates  = model_sel %>%
#     lapply(., function(x){x$Covariates}) %>% do.call(rbind, .)
#   , Formula     = model_sel %>%
#     lapply(., function(x){x$Formula}) %>% do.call(c, .)
#   , Model       = model_sel %>%
#     lapply(., function(x){x$Model}) %>% do.call(c, .)
#   , AIC         = model_sel %>%
#     lapply(., function(x){x$AIC}) %>% do.call(c, .)
#   , DeltaAIC    = NA # Difference in AIC to the "best" model
#   , Weight      = NA # Weight assigned according to AIC
#   , LogLik      = NA # Loglikelihood of the model
# )
#
# # Extract the LogLiks
# models$LogLik <- lapply(models$Model, logLik) %>%
#   do.call(rbind, .) %>%
#   as.vector()
#
# # Add a unique model id
# models$ModelID <- 1:nrow(models)
#
# # Calculate the DeltaAIC to the most parsimonious model
# models$DeltaAIC <- models$AIC - min(models$AIC, na.rm = TRUE)
#
# # Sort the tibble by the DeltaAIC
# models <- arrange(models, DeltaAIC)

############################################################
#### Calculating Model Weights
############################################################
# Calculate the AIC weights for the models with deltaAIC <= 2. Lets first write
# a function to get the AIC weight of the best models
AICweight <- function(x){
  (exp(-0.5 * x)) / sum(exp(-0.5 * x))
}

# Split our tibble into models with DeltaAIC <= 2 as the threshold
threshold <- 2
models1 <- subset(models, DeltaAIC <= threshold)
models2 <- subset(models, !DeltaAIC <= threshold | is.na(DeltaAIC))

# Calculate the model weights. Models of the second group will all get weight 0
models1$Weight <- AICweight(models1$DeltaAIC)
models2$Weight <- 0

# Put the models together in a tibble again
rownames(models1$Covariates) <- "0"
models <- rbind(models1, models2)

# Look at the result
print(models, n = 100)

# Look at the results of the best model
summary(models$Model[[1]])

# Plot the model coefficients of the best model
showCoeffs(getCoeffs(models$Model[[1]])[-1, ])

# Let's prepare a nice table to report in our results
table <- cbind(
    models$Covariates
  , select(models, -c(Covariates, Formula, Model))
)

# Store the table
write_rds(table, "03_Data/03_Results/99_PermeabilityModelAICs.rds")

# To save space when saving the model results we subset to only those models
# that get positive weight according to AIC
models <- subset(models, Weight > 0)

# Store the result for later
write_rds(models, "03_Data/03_Results/99_ModelSelection.rds")

# ############################################################
# #### Interpretation of Interactions
# ############################################################
# # We now want to prepare some plots that help us to understand how interactions
# # impact our dependent variable. Let's first isolate the best model
# best <- models$Model[[1]]
#
# # Look at the model summary
# summary(best)
#
# # Run the function on the interactions of interest
# visInt(best, "Water", "DistanceToHumans")
# visInt(best, "DistanceToWater", "HumansBuff5000")
# visInt(best, "HumansBuff5000", "DistanceToHumans")
# visInt(best, "Trees", "DistanceToHumans")
# visInt(best, "DistanceToWater", "Shrubs")
# visInt(best, "Water", "Shrubs")
#
# ############################################################
# #### Model Averaging
# ############################################################
# # Check the paper by Symonds and Moussalli (2011) for details on the formulas
# # In our case model averaging isn't really necessary since the first model gets
# # a weight of 1 anyways. However, this might still be a useful exercise. First
# # we subset to the models with a postive weight (only one in our case)
# predis <- subset(models, Weight > 0)
#
# # From each model we need to extract the coefficients
# predis <- predis %>% mutate(Coeffs = map(Model, function(x){
#   getCoeffs(x)
# }))
#
# # Add the weights to the coefficient tables. Let's also use this step to add the
# # model ID into each dataframe
# for (i in 1:nrow(predis)){
#   predis$Coeffs[[i]]$Weight <- predis$Weight[i]
#   predis$Coeffs[[i]]$ModelID <- predis$ModelID[i]
# }
#
# # Bind all dataframes of coefficients together
# predis <- do.call(rbind, predis$Coeffs)
#
# # Make sure that the weights of the covariates from the base model add up to 1
# predis %>%
#   group_by(Covariate) %>%
#   summarize(TotalWeight = sum(Weight))
#
# # Calculate the weighted coefficients
# predis <- predis %>% mutate(WeightedCoefficient = Coefficient * Weight)
#
# # Sum up weighted coefficients
# avg_coeffs <- predis %>%
#   group_by(Covariate) %>%
#   summarize(AveragedCoefficient = sum(WeightedCoefficient))
#
# # We still need to calculate the averaged variance. To do so we firsts put the
# # averaged coefficients back into the predis table
# predis <- left_join(predis, avg_coeffs, by = "Covariate")
#
# # We can now calculate the weighted SEs
# predis <- predis %>%
#   mutate(WeightedSE = Weight * sqrt(
#     SE ** 2 + (Coefficient - AveragedCoefficient) ** 2
#   ))
#
# # Sum up the weighted SEs
# avg_ses <- predis %>%
#   group_by(Covariate) %>%
#   summarize(AveragedSEs = sum(WeightedSE))
#
# # Put everything into one table
# coeffs <- full_join(avg_coeffs, avg_ses) %>%
#   rename(Coefficient = AveragedCoefficient, SE = AveragedSEs)
#
# # Look at the results
# showCoeffs(coeffs[-1, ])

############################################################
#### Model Validation: Necessary Function
############################################################
# Finally I validate the best model using the procedure described in Fortin et
# al. 2009. Let's reload the model results
models <- read_rds("03_Data/03_Results/99_ModelSelection.rds")

# Set a seed for reproducibility
set.seed(1234)

# Write a function that allows us to run the validation process
crossVal <- function(formula, data, ratio, random = FALSE){

  # Identify all unique step_ids
  step_ids <- unique(data$step_id_)

  # Split the step_ids into two groups with x% and (100-x)% of all data
  index <- !(sample.split(step_ids, SplitRatio = ratio))

  # Identify the step_ids of the training steps and those of the validation
  # steps
  steps_train <- step_ids[index]
  steps_valid <- step_ids[!index]

  # Split the data into training and validation datasets
  ssf_train <- subset(ssf, step_id_ %in% steps_train)
  ssf_valid <- subset(ssf, step_id_ %in% steps_valid)

  # Train the model using the training data
  model <- glmm_clogit(formula, data = ssf_train)

  # Identify the coefficients from the model
  coeffs <- getCoeffs(model)

  # We want to use the derived model to predict the selection scores in the
  # validation data. Let's create a simplified dataframe of the model coefficients
  pred <- as.data.frame(select(coeffs, Coefficient))

  # For easier indexing we assign the covariates as row-names
  rownames(pred) <- coeffs$Covariate

  # Predict the selection scores for the validation data
  ssf_valid$Scores <- exp(
    pred["cos(ta_)", ]        * cos(ssf_valid[, "ta_"]) +
    pred["log(sl_)", ]        * log(ssf_valid[, "sl_"]) +
    pred["Water", ]           * ssf_valid[, "Water"] +
    pred["DistanceToWater", ] * ssf_valid[, "DistanceToWater"] +
    pred["Shrubs", ]          * ssf_valid[, "Shrubs"] +
    pred["HumansBuff5000", ]  * ssf_valid[, "HumansBuff5000"] +
    pred["Trees", ]           * ssf_valid[, "Trees"]
  )

  # Depending on whether we want to randomize preferences, we include or exclude
  # realized steps
  if (random){

    # In case we want to randomize, get rid of case_ steps
    validation <- subset(ssf_valid, !case_) %>%

      # Group by cluster
      group_by(., step_id_) %>%

      # Identify the rank of each step
      mutate(., Rank = order(order(Scores, decreasing = TRUE))) %>%

      # Randomly select one of the steps per cluster
      sample_n(1)

  } else {

    # In case we dont want to randomize, we keep all steps
    validation <- ssf_valid %>%

      # Group by cluster
      group_by(., step_id_) %>%

      # Identify the rank of each step
      mutate(., Rank = order(order(Scores, decreasing = TRUE))) %>%

      # Subset to realized steps only
      subset(., case_)

  }

  # Now do some transformations
  validation <- validation %>%

    # Ungroup
    ungroup(.) %>%

    # Select the rank column
    select(., Rank) %>%

    # Prepare a table that gives the frequency of each rank
    table(.) %>%

    # Coerce table to dataframe
    as.data.frame(.) %>%

    # Rename columns nicer
    set_names(., c("Rank", "Frequency")) %>%

    # Make the rank a numeric variable
    mutate(., Rank = as.numeric(Rank))

  # Run a pearson rank corrleation test
  spearman <- cor.test(
      x       = validation$Rank
    , y       = validation$Frequency
    , method  = "spearman"
    , exact   = FALSE
  )

  # Return the spearman test result as well as the data
  return(
    list(
        "Data" = validation
      , "SpearmanCorrelation" = as.vector(spearman$estimate)
    )
  )
}

############################################################
#### Validation with observed preferences
############################################################
# Run the function 100 times without randomized preferences, but skip an
# iteration if there is an erorr (due to convergence)
valid_pref <- list()

# We run the loop until there are 100 models
i <- 1
while (i < 101){
   output <- tryCatch(crossVal(
      formula = models$Formula[[1]]
    , data    = ssf
    , ratio   = 0.2
    , random  = FALSE
  ), warning = function(w){return("skip")})

  # In case there was a warning, we skip the iteration, otherwise we store the
  # output and increase the counter
  if (class(output) != "list"){
    next
  } else {
    valid_pref[[i]] <- output
    i <- i + 1
  }

  # Print status
  cat(i, "of", 100, "done... \n")
}

# Extract all the rhos
rho_pref <- lapply(valid_pref, function(x){
  x[[2]]
}) %>% unlist()

# Extract all the data
dat_pref <- lapply(valid_pref, function(x){
  x[[1]]
}) %>% do.call(rbind, .)


############################################################
#### Validation with random preferences
############################################################
# Run the function 100 times with randomized preferences, but skip an iteration
# if there is an erorr (due to convergence)
valid_rand <- list()

# We run the loop until there are 100 models
i <- 1
while (i < 101){
   output <- tryCatch(crossVal(
      formula = models$Formula[[1]]
    , data    = ssf
    , ratio   = 0.2
    , random  = TRUE
  ), warning = function(w){return("skip")})

  # In case there was a warning, we skip the iteration, otherwise we store the
  # output and increase the counter
  if (class(output) != "list"){
    next
  } else {
    valid_rand[[i]] <- output
    i <- i + 1
  }
}

# Extract all the rhos
rho_rand <- lapply(valid_rand, function(x){
  x[[2]]
}) %>% unlist()

# Extract all the data
dat_rand <- lapply(valid_rand, function(x){
  x[[1]]
}) %>% do.call(rbind, .)

############################################################
#### Validation: Testing
############################################################
# Find the means and confidence intervals of the rho's from the two datasets
test_pref <- t.test(rho_pref)
test_rand <- t.test(rho_rand)

# Put everything together into a single dataframe
validation <- data.frame(
    Group = c("rho_pref", "rho_rand")
  , Mean  = c(test_pref$estimate, test_rand$estimate)
  , LCL   = c(test_pref$conf.int[1], test_rand$conf.int[1])
  , UCL   = c(test_pref$conf.int[2], test_rand$conf.int[2])
)

# Put the results to file
write_rds(
    validation
  , "03_Data/03_Results/99_ModelValidation.rds"
)
write_rds(
    list(dat_pref, dat_rand)
  , "03_Data/03_Results/99_ModelValidation(Data).rds"
)
