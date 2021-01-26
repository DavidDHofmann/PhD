################################################################################
#### Descriptive Stats of Extracted Covariates
################################################################################
# Description: Get an understanding of the different covariates

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
library(corrplot)     # To plot correlations
library(cowplot)      # For nice plots
library(ggpubr)       # For nice plots
library(glmmTMB)      # For modelling
library(gridExtra)    # For laying out multiple plots

################################################################################
#### Comparison to Old Data
################################################################################
# I want to compare the current dataset to the dataset used in Hofmann et al
# 2021. Let's thus load the two
old <- read_csv("/home/david/ownCloud/University/15. PhD/Chapter_0/03_Data/02_CleanData/00_General_Dispersers_Popecol(SSF4Hours).csv")
# new <- read_csv("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/00_General_Dispersers_POPECOL(iSSF_Extracted).csv")
new <- read_csv("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/00_General_Dispersers_POPECOL(Extracted).csv")

# Subset to case_ steps
old <- subset(old, case_)
new <- subset(new, case_ == 1)

# Let's join the data by timestamp
joined <- inner_join(old, new, by = c("id", "t1_"))

# Compare covariates: Step Length
cbind(
    summary(joined$sl_.x)
  , summary(joined$sl_.y)
)
cor(joined$sl_.x, joined$sl_.y)

# Compare covariates: Turning Angle
cbind(
    summary(joined$ta_.x)
  , summary(joined$ta_.y)
)
cor(joined$ta_.x, joined$ta_.y)

# Compare covariates: Humans
cbind(
    summary(joined$HumansBuff5000)
  , summary(joined$Facebook_HumanInfluenceBuffer_5000)
)
cor(joined$HumansBuff5000, joined$Facebook_HumanInfluenceBuffer_5000)

# Compare covariates: Distance to Humans
cbind(
    summary(joined$DistanceToHumans)
  , summary(joined$Facebook_DistanceToHumans)
)
cor(joined$DistanceToHumans, joined$Facebook_DistanceToHumans)

# Compare covariates: Water
cbind(
    summary(joined$Water.x)
  , summary(joined$Water.y)
)
cor(joined$Water.x, joined$Water.y)

# Compare covariates: Distance To Water
cbind(
    summary(joined$DistanceToWater.x)
  , summary(joined$DistanceToWater.y)
)
cor(joined$DistanceToWater.x, joined$DistanceToWater.y)

# Compare covariates: Shrubs
cbind(
    summary(joined$Shrubs.x)
  , summary(joined$Shrubs.y)
)
cor(joined$Shrubs.x, joined$Shrubs.y)

# Compare covariates: Roads
cbind(
    summary(as.numeric(joined$RoadCrossing.x))
  , summary(as.numeric(joined$RoadCrossing.y))
)
cor(joined$RoadCrossing.x, joined$RoadCrossing.y)

# Compare covariates: Distance to Roads
cbind(
    summary(as.numeric(joined$DistanceToRoads.x))
  , summary(as.numeric(joined$DistanceToRoads.y))
)
cor(joined$DistanceToRoads.x, joined$DistanceToRoads.y)

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
#### Compare Models
################################################################################
# Reload data
old <- read_csv("/home/david/ownCloud/University/15. PhD/Chapter_0/03_Data/02_CleanData/00_General_Dispersers_Popecol(SSF4Hours).csv")
# new <- read_csv("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/00_General_Dispersers_POPECOL(iSSF_Extracted).csv")
new <- read_csv("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/00_General_Dispersers_POPECOL(Extracted).csv")
old <- subset(old, t1_ %in% new$t1_ & id %in% new$id)
new <- subset(new, t1_ %in% old$t1_ & id %in% old$id)

# Compare fixes by individual (almost identical)
tab1 <- table(old$id[old$case_]) %>% as.data.frame()
tab2 <- table(new$id[new$case_ == 1]) %>% as.data.frame()
tab <- full_join(tab1, tab2, by = "Var1")
tab$Freq.x[is.na(tab$Freq.x)] <- 0
tab$Diff <- abs(tab$Freq.x - tab$Freq.y)
tab

# Remove some individuals
# old <- subset(old, !(id %in% c("Kalahari", "Dell", "Aspen", "Chiounard", "Encinitas")))
# new <- subset(new, !(id %in% c("Kalahari", "Dell", "Aspen", "Chiounard", "Sishen")))

# Let's run the model from Hofmann et al. 2020 again and compare the results
# using the two different datasets
old <- old %>% mutate(
    cos_ta_         = scale(cos(ta_))
  , log_sl_         = scale(log(sl_))
  , sl_             = scale(sl_)
  , Water           = scale(Water)
  , DistanceToWater = scale(sqrt(DistanceToWater))
  , HumansBuff5000  = scale(HumansBuff5000)
  , Shrubs          = scale(Shrubs)
  , Trees           = scale(Trees)
)
new <- new %>% mutate(
    cos_ta_         = scale(cos(ta_))
  , log_sl_         = scale(log(sl_))
  , sl_             = scale(sl_)
  , Water           = scale(Water)
  , DistanceToWater = scale(sqrt(DistanceToWater))
  , HumansBuff5000  = scale(Facebook_HumanInfluenceBuffer_5000)
  , Shrubs          = scale(Shrubs)
  , Trees           = scale(Trees)
)

# Run models
covars <- c("HumansBuff5000", "Water", "DistanceToWater", "Shrubs", "Trees")
mod1 <- glmm_clogit(writeForm(covars), data = old)
mod2 <- glmm_clogit(writeForm(covars), data = new)

# Compare results
summary(mod1)
summary(mod2)
coefs1 <- getCoeffs(mod1)[-1, ]
coefs2 <- getCoeffs(mod2)[-1, ]
list(coefs1, coefs2)

# Visualize coefficients
p1 <- showCoeffs(coefs1[-1, ])
p2 <- showCoeffs(coefs2[-1, ])
grid.arrange(p1, p2)

# Load step length distributions
summary(old$sl_[old$case_])
summary(new$sl_[new$case_ == 1])
