################################################################################
#### Step Selection Model
################################################################################
# Description: Fit the movement model and run the validation procedure. I
# thought about splitting the two tasks into separate scripts, yet both require
# the covariates to be scaled and utilize the same framework. The script
# therefore appears a bit lengthy and complicated.

# Clear R's brain
rm(list = ls())

# Surpress scientific notation
options(scipen = 999)

# Load required packages
library(tidyverse)    # For data wrangling
library(glmmTMB)      # For modelling
library(pbmcapply)    # For progress bar parallel
library(ggh4x)        # For nested ggplot-facettes

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Custom functions
source("02_R-Scripts/00_Functions.R")

# Specify file names
file_models <- "03_Data/03_Results/MovementModels.rds"
file_valid  <- "03_Data/03_Results/Validation.rds"
file_form   <- "03_Data/03_Results/Formula.rds"
file_scale  <- "03_Data/03_Results/Scaling.rds"
dir_valid   <- "03_Data/03_Results/Validation"

# Specify if movement models and validation should be rerun
rerun_models <- F
rerun_valid  <- F

# Define the maximum number of random steps that should be considered
# n_rsteps_consider <- 25

# Remove files if necessary
if (rerun_models & file.exists(file_models)) {
  file.remove(file_models)
}
if (rerun_models & file.exists(file_form)) {
  file.remove(file_form)
}
if (rerun_valid & file.exists(file_valid)) {
  file.remove(file_valid)
}
if (rerun_valid & dir.exists(dir_valid)) {
  unlink(dir_valid, recursive = T)
}

################################################################################
#### Loading Data
################################################################################
# Load data
dat <- read_rds("03_Data/02_CleanData/SSFExtracted.rds")

# Keep only desired columns
dat <- dat %>% select(
    id              = ID
  , Timestamp
  , TimestampRounded
  , step_id
  , step_id_within
  , case
  , sl
  , ta              = relta
  , dt
  , inactive
  , Sex
  , SeasonClimate
  , SeasonHerbivores
  , Night
  , Moonlight       = meanMoonlight
  , LightType

  # Static Covariates
  , HumansStatic
  , ForestStatic
  , TreesStatic
  , ShrubsStatic
  , WaterStatic
  , NDVIStatic
  , DistanceToWaterStatic
  , DistanceToPansStatic
  , TemperatureStatic
  , PrecipitationStatic

  # Dynamic Covariates
  , HumansDynamic
  , ForestDynamic
  , TreesDynamic
  , ShrubsDynamic
  , WaterDynamic
  , NDVIDynamic
  , DistanceToWaterDynamic
  , DistanceToPansDynamic
  , TemperatureDynamic
  , PrecipitationDynamic
)

# Add log of the step length, cos of the turning angle, and sqrt of the
# distance to water
dat <- dat %>% mutate(
    log_sl                     = log(sl)
  , cos_ta                     = cos(ta)
  , SqrtDistanceToWaterStatic  = sqrt(DistanceToWaterStatic)
  , SqrtDistanceToPansStatic   = sqrt(DistanceToPansStatic)
  , SqrtDistanceToWaterDynamic = sqrt(DistanceToWaterDynamic)
  , SqrtDistanceToPansDynamic  = sqrt(DistanceToPansDynamic)
)

# Let's also move all movement metrics to the front
dat <- dat %>% select(c(
  id, step_id, step_id_within, case, sl, log_sl, ta, cos_ta, inactive, everything()
))

################################################################################
#### Scale Data
################################################################################
# Note: Scaling is rather tricky in the situation considered. Since each
# covariate is represented by a "static" and a "dynamic" version, we need to
# think a bit harder how to best scale, while preserving comparability. I see
# two options to achieve this. First, to scale across all covariate values (i.e.
# ignoring if static or dynamic), or to first compute scaling parameters for the
# two datasets separately, and then take the mean. To avoid that one dataset has
# more power than the other, I'll go with the second approach.

# Get a summary over the (numeric) covariates
# summarizeData(dat)

# Define the columns that we want to standardize
vars <- c(
    "sl"
  , "log_sl"
  , "cos_ta"
  , "Night"
  , "Moonlight"

  # Dynamic Variables
  , "HumansDynamic"
  , "ForestDynamic"
  , "TreesDynamic"
  , "ShrubsDynamic"
  , "WaterDynamic"
  , "NDVIDynamic"
  , "DistanceToWaterDynamic"
  , "DistanceToPansDynamic"
  , "SqrtDistanceToWaterDynamic"
  , "SqrtDistanceToPansDynamic"
  , "TemperatureDynamic"
  , "PrecipitationDynamic"

  # Static Variables
  , "HumansStatic"
  , "ForestStatic"
  , "TreesStatic"
  , "ShrubsStatic"
  , "WaterStatic"
  , "NDVIStatic"
  , "DistanceToWaterStatic"
  , "DistanceToPansStatic"
  , "SqrtDistanceToWaterStatic"
  , "SqrtDistanceToPansStatic"
  , "TemperatureStatic"
  , "PrecipitationStatic"
)

# Find scaling parameters (mean and sd)
scaling <- dat %>%
  dplyr::select(vars) %>%
  pivot_longer(sl:PrecipitationStatic, names_to = "Covariate", values_to = "Value") %>%
  group_by(Covariate) %>%
  summarize(Mean = mean(Value, na.rm = T), SD = sd(Value, na.rm = T), .groups = "drop") %>%
  mutate(Covariate2 = gsub(Covariate, pattern = "Static", replacement = "")) %>%
  mutate(Covariate2 = gsub(Covariate2, pattern = "Dynamic", replacement = "")) %>%
  group_by(Covariate2) %>%
  mutate(center = mean(Mean), scale = mean(SD)) %>%
  ungroup() %>%
  dplyr::select(-c(Covariate2, Mean, SD)) %>%
  as.data.frame()

# Apply scaling
dat[, vars] <- scaleCovars(dat[, vars], scaling)

# # Run the standardization on them
# dat <- mutate(dat, across(all_of(vars), ~standardize(., operation = "scale")))

# Check if it worked
# summarizeData(subset(dat, case == 1))
#
# # Let's write the scaling parameters to a separate file
# standardized <- lapply(1:ncol(dat), function(x) {
#   center <- attr(dat[, x, drop = T], "standardize:center")
#   scale  <- attr(dat[, x, drop = T], "standardize:scale")
#   if (!is.null(center) & !is.null(scale)) {
#     return(data.frame(Variable = names(dat)[x], Center = center, Scale = scale))
#   } else {
#     return(NA)
#   }
# }) %>% do.call(rbind, .) %>% subset(!is.na(Variable))
# standardized <- list(
#       center = setNames(standardized[, 2], standardized$Variable)
#     , scale  = setNames(standardized[, 3], standardized$Variable)
# )
#
# # Write scaling factors and scaled data to file
write_rds(scaling, file_scale)
# write_rds(dat, "03_Data/02_CleanData/SSFExtractedScaled.rds")

################################################################################
#### Main Movement Models
################################################################################
# If model outputs already exist, load them. Otherwise, prepare a tibble into
# which we'll store all model results
if (file.exists(file_models)) {
    cat(file_models, "already exists and will be reloaded\n")
    design <- read_rds(file_models)
  } else {
    cat(file_models, "does not exist and will be created\n")
    design <- expand_grid(
          FittingCovariates    = c("Static", "Dynamic")
        , ModelSeasons         = c("Single", "Multi")
        , Formula              = c("Simple", "Full")
        , NumberRandomSteps    = c(10, 25, 50, 75, 100)
      ) %>%
      left_join(., tibble(
          ModelSeasons = c("Single", "Multi")
        , Season       = list(c("All"), c("Wet", "Dry")))
        , by           = "ModelSeasons"
      ) %>%
      unnest(Season)
    write_rds(design, file_models)
}

# Model formulas (simple formula, and then a more complicated one)
form_simple <- formula(case ~
  + (1 | step_id)
  + cos_ta
  + (0 + cos_ta | id)
  + sl
  + (0 + sl | id)
  + log_sl
  + (0 + log_sl | id)
  + Humans
  + (0 + Humans | id)
  + Shrubs
  + (0 + Shrubs | id)
  + Trees
  + (0 + Trees | id)
  + Water
  + (0 + Water | id)
  + SqrtDistanceToWater
  + (0 + SqrtDistanceToWater | id)
)
form_full <- update(form_simple, . ~ .
  + Temperature:sl
  + Temperature:log_sl
  + LightType:sl
  + LightType:log_sl
  + sl:Water
  + sl:Trees
  + sl:Shrubs
  + cos_ta:SqrtDistanceToWater
)

# Put them together
forms <- list(Simple = form_simple, Full = form_full)

# Write the formula to file
write_rds(forms, file_form)

# Loop through the design, fit the respective models, and extract the
# coefficients
design <- slice_sample(design, n = nrow(design), replace = F)
done <- "Models" %in% names(design)
if (!done) {

  # Run through the design in parallel
  cat("Fitting movement models\n")
  design$Results <- pbmclapply(
      X                  = seq_len(nrow(design))
    , ignore.interactive = T
    , mc.cores           = detectCores() - 1
    , FUN                = function(i) {

    # Create duplicate dataset and subset to desired number of random steps
    dat_dup <- subset(dat, step_id_within <= design$NumberRandomSteps[i])

    # Keep required covariates (either the static or the dynamic ones)
    if (design$FittingCovariates[i] == "Static") {
        dat_dup <- select(dat_dup, -contains("Dynamic"))
        names(dat_dup) <- gsub(names(dat_dup), pattern = "Static", replacement = "")
      } else {
        dat_dup <- select(dat_dup, -contains("Static"))
        names(dat_dup) <- gsub(names(dat_dup), pattern = "Dynamic", replacement = "")
    }

    # Split the data by season (this can either be two seasons or a single one)
    if (design$Season[i] == "Wet") {
        dat_dup <- subset(dat_dup, SeasonClimate == "Wet")
      } else if (design$Season[i] == "Dry") {
        dat_dup <- subset(dat_dup, SeasonClimate == "Dry")
    }

    # Select appropriate model formula
    form <- forms[[design$Formula[i]]]

    # Drop entries where the distance to pans is NA (also make sure that for the
    # remaining steps we have at least x random steps)
    # dat_dup <- subset(dat_dup, !is.na(SqrtDistanceToPans))
    # dat_dup <- subset(dat_dup, step_id %in% subset(count(dat_dup, step_id), n >= 26)$step_id)
    # dat_dup <- subset(dat_dup, !is.na(SqrtDistanceToWater))
    # dat_dup <- subset(dat_dup, step_id %in% subset(count(dat_dup, step_id), n >= min(design$NumberRandomSteps[i], 26))$step_id)

    # Fit model. Let's also stop the time as it will give us an indication of
    # how time-hungry the validation might be.
    mod <- glmm_clogit(form, dat_dup)

    # Extract coefficients tables
    mod <- getCoeffs(mod, zvalue = T, pvalue = T, ranefs = T)

    # Return the results
    return(mod)

  })

  # Unnest the results
  design <- unnest(design, Results)

  # Write results to file
  write_rds(design, file_models)
  design <- read_rds(file_models)

}

# # Plot of estimates depending on number of random steps
# design %>%
#   unnest(Models) %>%
#   ggplot(aes(x = Covariate, y = Coefficient, col = factor(NumberRandomSteps), group = factor(NumberRandomSteps))) +
#     geom_point(
#         shape    = 3
#       , size     = 1
#       , position = position_dodge(width = 0.8)
#     ) +
#     facet_nested(FittingCovariates + factor(Season) ~ "Nothing", scales = "free_x", space = "free_x") +
#     theme_awesome() +
#     labs(y = expression(beta*"-Coefficient")) +
#     scale_color_viridis_d(begin = 0.3, name = "Number of Random Steps", direction = -1) +
#     scale_linewidth_manual(
#         name   = "Confidence Level"
#       , values = c(1.5, 0.75, 0.3)
#     ) +
#     theme(
#       , legend.position   = "bottom"
#       , legend.margin     = margin(0, 50, 0, -20)
#       , legend.box.margin = margin(-5, -10, -5, -10)
#       , legend.text       = element_text(face  = 3)
#       , legend.title      = element_text(face  = 3)
#       , panel.grid.minor  = element_blank()
#       , axis.text.x       = element_text(angle = 45, hjust = 1)
#       , axis.title.y      = element_text(angle = 90)
#     ) +
#     guides(
#         colour    = guide_legend(title.position = "top", title.hjust = 0.5)
#       , linewidth = guide_legend(title.position = "top", title.hjust = 0.5, override.aes = list(colour = "#3CBB75FF"))
#     )

################################################################################
#### Validation
################################################################################
# Given the plot above, it is obvious that fewer than 100 random steps suffice
# to achieve stable model estimates. To reduce the runtime of the validation, we
# can subset to fewer random steps.

# Validation parameters
ratio <- 0.8
steps <- 25
reps  <- 100

# Validation design to run through
design <- expand_grid(
      FittingCovariates    = c("Static", "Dynamic")
    , ModelSeasons         = c("Single", "Multi")
    , PredictionCovariates = c("Static", "Dynamic")
    , Formula              = c("Simple", "Full")
    , NumberRandomSteps    = steps
    , Replicate            = 1:reps
  ) %>%
  left_join(., tibble(
      ModelSeasons = c("Single", "Multi")
    , Season       = list(c("All"), c("Wet", "Dry")))
    , by           = "ModelSeasons"
  ) %>%
  unnest(Season) %>%
  subset(!(FittingCovariates == "Static" & PredictionCovariates == "Dynamic")) %>%
  mutate(ModelCode = paste0(
      substr(FittingCovariates, 1, 1)
    , substr(ModelSeasons, 1, 1)
    , substr(PredictionCovariates, 1, 1)
    , "_"
    , substr(Season, 1, 1)
    , substr(Formula, 1, 1)
  )) %>%
  mutate(Filename = paste0("03_Data/03_Results/Validation/", ModelCode, "_", sprintf("%04d", Replicate), ".rds")) %>%
  mutate(Done = file.exists(Filename))

# Store this design to file
write_rds(design, file_valid)

# Create directory if needed
if (!dir.exists(dir_valid)) {
    cat("Validation directory does not exist and will be created\n")
    dir.create(dir_valid)
  } else {
    cat("Validation directory exists and will not be created\n")
}

# Subset to entries that still need to be processed and shuffle the entries to
# get more accurate prediction of processing time
design <- subset(design, !Done)
design <- slice_sample(design, n = nrow(design), replace = F)
cat("There are", nrow(design), "validation scenarios to iterate through\n")

# Set seed for multicore processes
set.seed(1234)
RNGkind("L'Ecuyer-CMRG")
mc.reset.stream()

# Run through the design and apply the validation procedure
cat("Running validation procedure\n")
design$Validation <- pbmclapply(
    X                  = seq_len(nrow(design))
  , mc.cores           = detectCores() - 1
  , ignore.interactive = T
  , FUN                = function(i) {

  # Split the data by season (this can either be two seasons or a single one)
  dat_dup <- subset(dat, step_id_within <= design$NumberRandomSteps[i])
  if (design$Season[i] == "Wet") {
      dat_dup <- subset(dat_dup, SeasonClimate == "Wet")
    } else if (design$Season[i] == "Dry") {
      dat_dup <- subset(dat_dup, SeasonClimate == "Dry")
  }

  # Split data into training and prediction
  step_ids  <- unique(dat_dup$step_id)
  sampled   <- sample(step_ids, size = length(step_ids) * ratio, replace = F)
  dat_train <- subset(dat_dup, step_id %in% sampled)
  dat_valid <- subset(dat_dup, !(step_id %in% sampled))

  # Keep required covariates for training (either the static or the dynamic
  # ones)
  if (design$FittingCovariates[i] == "Static") {
      dat_train <- select(dat_train, -contains("Dynamic"))
      names(dat_train) <- gsub(names(dat_train), pattern = "Static", replacement = "")
    } else {
      dat_train <- select(dat_train, -contains("Static"))
      names(dat_train) <- gsub(names(dat_train), pattern = "Dynamic", replacement = "")
  }

  # Keep required covariates for prediction (either the static or the dynamic
  # ones)
  if (design$PredictionCovariates[i] == "Static") {
      dat_valid <- select(dat_valid, -contains("Dynamic"))
      names(dat_valid) <- gsub(names(dat_valid), pattern = "Static", replacement = "")
    } else {
      dat_valid <- select(dat_valid, -contains("Static"))
      names(dat_valid) <- gsub(names(dat_valid), pattern = "Dynamic", replacement = "")
  }

  # Select appropriate model formula
  form <- forms[[design$Formula[i]]]

  # Drop entries where the distance to pans is NA (also make sure that for the
  # remaining steps we have at least 25 random steps)
  # dat_train <- subset(dat_train, !is.na(SqrtDistanceToPans))
  # dat_train <- subset(dat_train, step_id %in% subset(count(dat_train, step_id), n >= 26)$step_id)

  # Run cross validation
  success <- F
  attempt <- 0
  while (!success & attempt <= 10) {
    valid <- tryCatch(crossVal(
        formula   = form
      , dat_train = dat_train
      , dat_valid = dat_valid
    ), warning = function(w){ return("skip") })
    success <- is.data.frame(valid)
    attempt <- attempt + 1
  }

  # Write it to file
  write_rds(valid, design$Filename[i])

})

# Now reload the validation data and ensure that all are done
valid <- file_valid %>%
  read_rds() %>%
  mutate(Done = file.exists(Filename))
write_rds(valid, file_valid)

# Once all is finished, load the data. Check how many validations were
# successful
valid <- file_valid %>%
  read_rds() %>%
  mutate(Data = map(Filename, read_rds)) %>%
  dplyr::select(-c(Filename, Done)) %>%
  mutate(Keep = sapply(Data, is.data.frame))

# Unnest the ones that were successful
valid <- valid %>%
  subset(Keep) %>%
  dplyr::select(-Keep) %>%
  unnest(Data) %>%
  nest(Data = -c(FittingCovariates, ModelSeasons, PredictionCovariates, Formula, Replicate, Preferences)) %>%
  ungroup()

# Run through the validation data and compute the rank-frequency
valid <- valid %>%
  mutate(Data = map(Data, function(x) {
    x <- x %>%
      ungroup(.) %>%
      select(., Rank) %>%
      table(.) %>%
      as.data.frame(.) %>%
      set_names(., c("Rank", "Frequency")) %>%
      mutate(., Rank = as.numeric(Rank))
    return(x)
  }))

# Compute spearman's rank correlation
valid$Spearman <- sapply(valid$Data, function(x) {

  # Run a pearson rank corrleation test
  spearman <- cor.test(
      x       = x$Rank
    , y       = x$Frequency
    , method  = "spearman"
    , exact   = FALSE
  )
  spearman <- as.vector(spearman$estimate)
  return(spearman)
})

# Create Model Codes
valid <- valid %>%
  select(-Data) %>%
  mutate(ModelCode = paste0(
      substr(FittingCovariates, 1, 1)
    , substr(ModelSeasons, 1, 1)
    , substr(PredictionCovariates, 1, 1)
  )) %>%
  mutate(ModelCode = factor(ModelCode, levels = c("SSS", "SMS", "DSS", "DSD", "DMS", "DMD")))

# # Put the results to file
write_rds(valid, "03_Data/03_Results/RankFrequency.rds")

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/02_Analysis/02_MovementModel.rds")
cat("Done :)\n")
