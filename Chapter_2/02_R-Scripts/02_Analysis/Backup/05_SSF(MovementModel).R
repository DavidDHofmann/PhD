############################################################
#### Step Selection Function - Movement Model
############################################################
# Description: We will now use the most parsimonious model and add interaction
# terms with movement metrics. This will allow us to identify how individuals
# move differently depending on the prevailing habitat. We will use the model
# results to simulate dispersers

# Clear R's brain
rm(list = ls())

# Surpress scientific notation
options(scipen = 999)

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/00_WildDogs"
setwd(wd)

# Load required packages
library(raster)
library(rgeos)
library(data.table)
library(tidyverse)
library(glmmTMB)
library(lubridate)
library(viridis)
library(scales)
library(ggpubr)

# Load custom functions
source("Functions.r")

############################################################
#### Data Preparation and Exploratory Analyses
############################################################
# Load the 4 hourly fixes (we will load the shapefile first)
ssf <- "03_Data/02_CleanData/00_General_Dispersers_Popecol(SSF4Hours).csv" %>%

  # We need to manually set some of the column types
  read_csv(., col_types = cols(
        Homerange     = "d"
      , OtherPack     = "d"
      , id            = "f"
      , tod_          = "f"
      , RoadCrossing  = "f"
      , State         = "f")) %>%

  # Remove undesired columns
  select(-c("X1")) %>%

  # Make some adjustments to covariates and scale them. Remember that one cannot
  # scale categorical variables.
  mutate(
      DistanceToWater   = sqrt(DistanceToWater)
    , DistanceToHumans  = sqrt(DistanceToHumans)
    , DistanceToRoads   = sqrt(DistanceToRoads)
    , log_sl_           = log(sl_)
    , cos_ta_           = cos(ta_)
  ) %>%
  transform(
      Water                 = scale(Water)
    , DistanceToWater       = scale(DistanceToWater)
    , Trees                 = scale(Trees)
    , Shrubs                = scale(Shrubs)
    , Protected             = scale(Protected)
    , HumansBase            = scale(HumansBase)
    , HumansAverage         = scale(HumansAverage)
    , HumansBuff5000        = scale(HumansBuff5000)
    , DistanceToHumans      = scale(DistanceToHumans)
    , DistanceToRoads       = scale(DistanceToRoads)
    , OtherPack             = scale(OtherPack)
    , Homerange             = scale(Homerange)
    , NoHomeranges          = scale(NoHomeranges)
)

# We know that wild dogs typically follow a diurnal activity pattern. That is,
# they are typically most active in the early morning and evening. To account
# for this, we may want to derive an indicator of the time of the day.
ssf$tod2_ <- ssf$t1_ %>%
  round_date("30 minutes") %>%
  strftime(tz = "UTC", format = "%H:%M:%S") %>%
  factor(., levels = c("03:00:00", "07:00:00", "15:00:00", "19:00:00", "23:00:00"))

# Let's see if the step length is somewhat correlated with this metric
ssf %>%
  subset(case_) %>%
  ggplot(aes(y = sl_, col = tod2_)) +
    geom_boxplot() +
    scale_y_log10()

# Maybe there are differences among individuals
ssf %>%
  subset(case_) %>%
  ggplot(aes(y = sl_, col = tod2_)) +
   geom_boxplot() +
   scale_y_log10() +
   facet_wrap("id")

# We may want to create two classes that separate these two distinct periods
ssf$Activity <- ifelse(ssf$tod2_ %in% c("03:00:00", "15:00:00")
  , yes = "MainActivity"
  , no  = "LowActivity"
)

# Visualize again by activity phase
ssf %>%
  subset(case_) %>%
  ggplot(aes(y = sl_, col = Activity)) +
    geom_boxplot() +
    scale_y_log10()

# We can also easily check for significant differences among means
ggboxplot(subset(ssf, case_), x = "Activity", y = "sl_") +
  stat_compare_means(
      comparison = list(c("LowActivity", "MainActivity"))
    , label = "p.signif"
  ) +
  scale_y_log10()

# Maybe there is a comparable effect for the turning angle?
ssf %>%
  subset(case_) %>%
  ggplot(aes(y = cos(ta_), col = Activity)) +
    geom_boxplot()

# Maybe there are differences among individuals
ssf %>%
  subset(case_) %>%
  ggplot(aes(y = cos(ta_), col = Activity)) +
   geom_boxplot() +
   facet_wrap("id")

# I would also expect that seasonality eventually comes into play. Let'
# therefore also get an indicator of the current season.
ssf$Season <- getSeason(ssf$t1_)

# Visualize step length by season
ssf %>%
  subset(case_) %>%
  ggplot(aes(y = sl_, col = Season)) +
    geom_boxplot() +
    scale_y_log10()

# Let's see how many fixes we have in each season
table(subset(ssf, case_)$Season)

# Maybe dispersers have to rest from time to time. To examine this, we need to
# identify the distance traveled during the past few fixes and see whether large
# distances are followed by smaller steps. To achieve this, I first identify
# bursts within which there is coherent GPS data. Let's assign a BurstID to each
# step
sub <- ssf %>%
  subset(case_) %>%
  group_by(id) %>%
  nest() %>%
  mutate(data = map(data, function(x){
    x$BurstID <- stepBursts(x)
    return(x)
  })) %>%
  unnest(data)

# Within each burst we can now calculate the distance moved during the past 6
# and 12 fixes (approx. 24 and 48 hours). Note that there is a slight bias as we
# assume that steps from 0700-1500 are four hours, rather than eight. This is
# also why we have to specify fixes rather than a time period
sub <- sub %>%
  group_by(id, BurstID) %>%
  nest() %>%
  mutate(data = map(data, function(x){
    x$CumDistancePast1Fixes <- distanceTraveled(x, fixes = 1)
    x$CumDistancePast6Fixes <- distanceTraveled(x, fixes = 6)
    x$CumDistancePast12Fixes <- distanceTraveled(x, fixes = 12)
    x$CumDistancePast18Fixes <- distanceTraveled(x, fixes = 18)
    x$CumDistancePast50Fixes <- distanceTraveled(x, fixes = 50)
    return(x)
  })) %>%
  unnest(data)

# Split the timestamp into a date and into time
sub$Date <- as.Date(sub$t1_)
sub$Time <- strftime(sub$t1_, format = "%H:%M:%S", tz = "UTC")
sub$Time <- as.POSIXct(sub$Time, format = "%H:%M:%S", tz = "UTC")

# Visualize distribution of step lengths depending on previous steps
ggplot(sub, aes(x = CumDistancePast1Fixes, y = sl_)) +
  geom_point(size = 0.1)

# It is difficult to see what is going on. It may be easier to visualize the log
# transformed data
ggplot(sub, aes(x = log(CumDistancePast1Fixes), y = log(sl_), color = Time)) +
  geom_point() +
  scale_color_gradientn(colors = viridis(50), trans = "time") +
  ggtitle("Past 1 Fix")

# Is there a similar pattern if we look at more previous fixes?
ggplot(sub, aes(x = log(CumDistancePast50Fixes), y = log(sl_), color = Time)) +
  geom_point() +
  scale_color_gradientn(colors = viridis(50), trans = "time") +
  ggtitle("Past 50 Fixes")

# Let's also calculate the distance to the first fix
################################################################################
#### I THINK THIS IS WRONG -> Also need to take into account y2_ - y1_
################################################################################
sub <- sub %>%
  group_by(id) %>%
  nest() %>%
  mutate(data = map(data, function(x){
    x$DistanceToFirst <- abs(x$x2_ - x$x1_[1])
    return(x)
  })) %>%
  unnest(data)

# Check how distance to first fix evolves over time
ggplot(sub, aes(x = t1_, y = DistanceToFirst)) +
  geom_line() +
  facet_wrap("id", scales = "free")

# We may want to check if the realized step length decreases over time
ggplot(sub, aes(x = t1_, y = sl_)) +
  geom_line() +
  facet_wrap("id", scales = "free")

# Let's see if there is a daily pattern (for single individuals)
ggplot(subset(sub, id == "Mirage"), aes(x = Time, y = sl_, col = factor(Date))) +
  geom_line() +
  facet_wrap("id", scales = "free") +
  theme(legend.position = "none") +
  scale_x_datetime(labels = date_format("%H"))

# Let's see if there is a daily pattern (on a population level)
ggplot(sub, aes(x = Time, y = sl_, col = factor(Date))) +
  geom_line() +
  facet_wrap("id", scales = "free") +
  theme(legend.position = "none") +
  scale_x_datetime(labels = date_format("%H"))

# Make the tod_ and season factorial
ssf$Activity <- as.factor(ssf$Activity)
ssf$Season <- as.factor(ssf$Season)

# Also scale the movement metrics now
ssf <- ssf %>% transform(
    sl_                   = scale(sl_)
  , log_sl_               = scale(log_sl_)
  , cos_ta_               = scale(cos_ta_)
)

############################################################
#### Forward Selection of Movement Interaction Terms
############################################################
# Reload the results from the model selection script
models <- readRDS("03_Data/03_Results/99_ModelSelection.rds")
summary(models$Model[[1]])

# # Run a model
# mod <- glmm_clogit(case_ ~
#   + cos_ta_
#   + log_sl_
#   + Water
#   + Trees
#   + DistanceToWater
#   + HumansBuff5000
#   + Shrubs
#   + log_sl_:Activity
#   + log_sl_:Water
#   + log_sl_:Trees
#   + cos_ta_:HumansBuff5000
#   + cos_ta_:DistanceToWater
#   + (1|step_id_)
#   + (0 + cos_ta_|id)
#   + (0 + log_sl_|id)
#   + (0 + Water|id)
#   + (0 + Trees|id)
#   + (0 + DistanceToWater|id)
#   + (0 + HumansBuff5000|id)
#   + (0 + Shrubs|id)
#   , data = ssf
# )
#
# # Check model results
# summary(mod)
# showCoeffs(getCoeffs(mod)[-1, ], xlim = c(-3, 3))


# Let's retrieve all covariates from the best model
base <- models$Covariates[1, ] %>% strsplit(., ", ") %>% .[[1]]

# We will only test for two-way movement interactions. Let us therefore get rid
# of any interaction terms in the covariates
remove <- grep(":", base)
if (length(remove) > 0){
  covars <- base[-grep(":", base)]
} else {
  covars <- base
}

# Now combine all of these covariates with the cos(ta_) and the log(sl_)
covars <- c(paste0(covars, ":", "log_sl_"), paste0(covars, ":", "cos_ta_"))

# We also want to test for an interaction between the step length and the
# activity phase
covars <- c(covars, "log_sl_:Activity")

# Look at the covariates that we are going to test for
covars

# Prepare a vector into which we put the covariates that are selected during our
# model selection loop
selected <- c()

# Also prepare a list into which we store all model results
model_sel <- list()

# Initiate a counter
i <- 1

# Run a loop that iteratively adds interaction terms until none are left for
# selection
while (length(covars) > 0){

  # Write down all of the covariates for the different models
  forformulas <- lapply(covars, function(x){
    c(base, selected, x)
  })

  # Prepare a vector that makes the covariates of a model slightly better
  # visible
  covariates <- forformulas %>%

    # Collapse the list
    do.call(rbind, .) %>%

    # Convert the matrix to a dataframe
    as.data.frame(.) %>%

    # Combine the columns into a vector that shows all covariates of a specific
    # model
    unite(., col = "Covariates", sep = ", ")

  # Prepare an empty tibble to fill
  models <- tibble(
      ModelID     = 1:length(covars)
    , Covariates  = covariates
    , Formula     = lapply(covars, function(x){

      # Write the model for the baseline covariates
      formula <- writeForm(base)

      # Check if we selected some interaction terms in the previous iteration
      if (length(selected) > 0){

        # If so, add them to the models as seperate terms (without random slope)
        toadd <- selected %>%

          # Transpose the matrix
          t(.) %>%

          # Convert the matrix to a dataframe
          as.data.frame(.) %>%

          # Combine the interactions using + as a seperator
          unite(., col, sep = " + ") %>%

          # Get rid of the dataframe structure
          .[["col"]]

        # Update the baseline model with the selected interaction terms
        formula <- update(formula, paste("~ . +", toadd))
      }

      # Add the interaction terms we still need to check
      formula <- update(formula, paste("~ . +", x))

      # Return the model formula that we just created
      return(formula)
    })
    , Model       = NA # Result of the model
    , AIC         = NA # AIC of the model
  )

  # Run all models
  models <- mutate(models, Model = map(Formula, function(x){
    glmm_clogit(x, data = ssf)
  }))

  # Extract the AIC value of each model
  models$AIC <- lapply(models$Model, AIC) %>%
    do.call(rbind, .) %>%
    as.vector()

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

# Put all relevant information into a single tibble
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
rownames(models1$Covariates) <- 1:nrow(models1)
rownames(models2$Covariates) <- (1 + nrow(models1)):(nrow(models2) + nrow(models1))
models <- rbind(models1, models2)

# Look at the result
print(models, n = 100)

# Look at the results of the best model
summary(models$Model[[1]])

# Plot the model coefficients of the best model
showCoeffs(getCoeffs(models$Model[[1]])[-1, ], xlim = c(-2, 2))

# Let's prepare a nice table to report in our results
table <- cbind(
    models$Covariates
  , select(models, -c(Covariates, Formula, Model))
)

# Store the table
write_rds(table, "03_Data/03_Results/99_MovementModelAICs.rds")

# To save space when saving the model results we subset to only those models
# that get positive weight according to AIC
models <- subset(models, Weight > 0)

# Store the result for later
write_rds(models, "03_Data/03_Results/99_MovementModel.rds")
