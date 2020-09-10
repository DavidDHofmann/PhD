############################################################
#### Step Selection Function - Model Selection
############################################################
# Description: In this script I run forward model selection using the 4-hourly
# fixes of our dispersers

# Clear R's brain
rm(list = ls())

# Surpress scientific notation
options(scipen = 999)

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_0"
setwd(wd)

# Load required packages
library(raster)
library(data.table)
library(rgeos)
library(tidyverse)
library(glmmTMB)
library(caTools)
library(lemon)
library(corrplot)
library(davidoff)

############################################################
#### Loading and Examining Data
############################################################
# Load the 4 hourly fixes
ssf <- "03_Data/02_CleanData/00_General_Dispersers_Popecol(SSF4Hours).csv" %>%

  # We need to manually set some of the column types
  read_csv(., col_types = cols(
        id            = "f"
      , tod_          = "f"
      , RoadCrossing  = "f"
      , State         = "f")) %>%

  # Remove undesired columns
  select(-c("X1", "x1_", "x2_", "y1_", "y2_", "t1_", "t2_", "dt_"))

# Let's check whether there is variation in all covariates for all dispersers
# (obviously not the case)
summary <- ssf %>%

  # Select variables for which we want to know the variance
  select(., id, Water:RoadCrossing) %>%

  # Turn dataframe to long format
  gather(., "Covariate", "Value", 2:ncol(.)) %>%

  # Group by id and covariate
  group_by(., id, Covariate) %>%

  # Calculate variances
  summarize(., Var = var(Value)) %>%

  # Turn back to wide format
  spread(., Covariate, Var)

# Show the summary
summary

# Make some adjustments to covariates and scale them. Remember that one cannot
# scale categorical variables.
ssf <- ssf %>%
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
    , log_sl_               = scale(log_sl_)
    , cos_ta_               = scale(cos_ta_)
    , sl_                   = scale(sl_)
)

# Make sure that the structure of the data is fine
str(ssf)

# Note that the scaling parameters of a specific covariate can be accessed using
# the following code. We will need the parameters to rescale the covariates when
# we run predictions. So lets store these parameters. First we prepare an empty
# dataframe
scaling <- data.frame(
    ColumnName  = names(ssf)
  , Center      = NA
  , Scale       = NA
)

# Extract the scaling and centering parameters from the dataframe
for (i in 1:ncol(ssf)){
  center <- attr(ssf[, i], "scaled:center")
  scale <- attr(ssf[, i], "scaled:scale")
  scaling$Center[i] <- ifelse(is.null(center), NA, center)
  scaling$Scale[i] <- ifelse(is.null(scale), NA, scale)
}

# Look at the resulting table
scaling

# Store the object to an RDS
write_rds(scaling, "03_Data/03_Results/99_ScalingSSF.rds")

############################################################
#### Correlation Analysis
############################################################
# We need to exclude covariates that are overly correlated. Let's thus check for
# correlation among them (not for factorial variables)
correlations <- ssf %>%
  select(., c(Water:RoadCrossing), - RoadCrossing) %>%
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

############################################################
#### Single Model
############################################################
# Write model forumla
form <- writeForm(
  c("Water", "DistanceToWater", "Trees", "Shrubs", "HumansBuff5000")
)

# Fit model
mod <- glmm_clogit(form, data = ssf)

# Check and plot results
summary(mod)
showCoeffs(getCoeffs(mod)[-1, ])

############################################################
#### Forward Model Selection for Main Effects
############################################################
# We now want to run forward model selection. That is, we iteratively increase
# the complexity of our models and keep track of the changes in the AIC.
# Ultimately, we will use the model with the lowest AIC as our most parsimonious
# model. If there are multiple models within a range of 2 AICs, we can average
# the corresponding models' coefficients if we feel so. To start, we need to
# know the set of covariates that we can choose from. So let's first identify
# all possible covariates.
covars <- c(
    "Water"
  , "DistanceToWater"
  , "Trees"
  , "Shrubs"
  , "Protected"
  # , "HumansBase"
  # , "HumansAverage"
  , "HumansBuff5000"
  , "DistanceToRoads"
  , "RoadCrossing"
)

# Prepare a vector into which we put the covariates that are selected during our
# model selection loop
selected <- c()

# Also prepare a list into which we store all model results
model_sel <- list()

# Prepare a counter that increases in every iteration of the loop(s)
i <- 1

# Run a while loop that iteratively adds covariates until there are no more
# covariates left for selection
while (length(covars) > 0){

  # Combine the selected covariates with each of the remaining covariates.
  # We can use the resulting list to prepare our model formulas. This step only
  # starts working in the second iteration of the loop
  forformulas <- lapply(covars, function(x){
    c(selected, x)
  })

  # Just for easier visibility, we also prepare a vector that contains all
  # covariates. This will allow us to easily spot the covariates of a specific
  # model
  covariates <- forformulas %>%
    do.call(rbind, .) %>%
    as.data.frame(.) %>%
    unite(., col = "Covariates", sep = ", ")

  # Prepare an empty tibble to fill
  models <- tibble(
      ModelID     = 1:length(covars)
    , Covariates  = covariates
    , Formula     = lapply(forformulas, writeForm)
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

##############################################################
#### Examine Random Effects
##############################################################
# Reload model
mod <- read_rds("03_Data/03_Results/99_ModelSelection.rds")
mod <- models$Model[[1]]

# Check out coefficients per individual
coeffs <- coef(mod)
ranefs <- ranef(mod, condVar = T)
coeffs$cond$id
ranefs$cond$id

# Note: ranef yields the difference between the individual specific effect and
# the mean level effect. coef, on the other hand, yields the individual specific
# effect. Thus, the followin two lines yield (approximately) the same
mean(coeffs$cond$id$cos_ta_) + ranefs$cond$id$cos_ta_[3]

# We now want to visualize the individual variation. There are two possibilities
# for this: lme4::dotplot() or a ggplot. The dotplot is easier, yet not
# customizable. Let's first do the dotplot, then recreate it in ggplot.
lme4:::dotplot.ranef.mer(ranef(mod)$cond)

# Maybe scalefree?
lme4:::dotplot.ranef.mer(ranef(mod)$cond, scales = list(x = list(relation = "free")))

# Prepare dataframe that we need to plot the same in ggplot
rfs <- ranefs$cond$id %>%
  rownames_to_column() %>%
  gather(key = Covariate, value = Mean, 2:9)
names(rfs)[1] <- "id"

# We need to add the conditional variance
condVar <- attributes(ranefs$cond$id)$condVar
names(condVar) <- attributes(ranefs$cond$id)$names
condVar <- as.data.frame(do.call(rbind, condVar))
names(condVar) <- attributes(ranefs$cond$id)$row.names
condVar <- rownames_to_column(condVar)
names(condVar)[1] <- "Covariate"
condVar <- gather(condVar, key = id, value = Variance, 2:17)

# Join data to rfs dataframe
rfs <- left_join(rfs, condVar)

# Rename stuff nicely
rfs$Covariate <- gsub(rfs$Covariate, pattern = "cos_ta_", replacement = "cos(ta)")
rfs$Covariate <- gsub(rfs$Covariate, pattern = "log_sl_", replacement = "log(sl)")
rfs$Covariate <- gsub(rfs$Covariate, pattern = "sl_", replacement = "sl")
rfs$Covariate <- gsub(rfs$Covariate, pattern = "HumansBuff5000", replacement = "HumanInfluence")

# Make covariates a factor
rfs$Covariate <- factor(rfs$Covariate, levels = c(
  "cos(ta)", "sl", "log(sl)", "Water", "DistanceToWater", "Shrubs"
  , "Trees", "HumanInfluence"
))

# Visualize. Note that I am transforming the variance using mean - 2 *
# sqrt(Variance). This was taken from here: https://stackoverflow.com/questions
# /13847936/plot-random-effects-from-lmer-lme4-package-using-qqmath-or-dotplot-
# how-to-mak
ggplot(rfs, aes(x = Mean, y = id)) +
  geom_point() +
  facet_wrap("Covariate", nrow = 2) +
  geom_errorbarh(aes(
      xmin = Mean - 2 * sqrt(Variance)
    , xmax = Mean + 2 * sqrt(Variance)
  ), colour = "black", height = 0) +
  xlim(-2.1, 2.1)

# Check out the variation in effects
as.data.frame(apply(coeffs$cond$id, 2, mean))
as.data.frame(apply(coeffs$cond$id, 2, sd))
as.data.frame(apply(coeffs$cond$id, 2, var))

# Note: Variance of some random terms is mereley 0. This shouldn't be an issue
# though -> Check this post: # https://stats.stackexchange.com/questions/115090/
# why-do-i-get-zero-variance-of-a-random-effect-in-my-mixed-model-despite-some-va

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

  # Extract coefficients
  coeffs <- fixef(model)$cond

  # Prepare model matrix
  modeldat <- model.matrix(terms(model), data = ssf_valid)

  # Calculate selection score
  ssf_valid$Scores <- as.vector(exp(modeldat %*% coeffs))

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
