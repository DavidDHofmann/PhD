################################################################################
#### Analysis of Simulated Data using Different Approaches
################################################################################
# Description: Analysis of the generated data using different approaches

# Clear R's brain
rm(list = ls())

# Load required packages
library(pbmcapply)     # For parallel computing
library(tidyverse)     # For data wrangling
library(raster)        # To handle spatial data
library(sf)            # For plotting spatial features
library(amt)           # To fit distributions
library(survival)      # To run conditional logistic regression
library(crawl)         # To impute missing fixes

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_4"
setwd(wd)

# Load custom functions
source("02_R-Scripts/Functions.R")

# Function to fit step length and turning angle distributions. Statically, as
# well as dynamically to different step durations

# Reload simulated data
dat <- read_rds("03_Data/Simulation.rds")

# Load all simulated data into memory
dat$Covariates <- lapply(dat$Covariates, readAll)
dat$Movement   <- lapply(dat$MovementFilename, read_rds)

# Remove any columns that are not needed anymore
dat <- select(dat, -c(CovariateFilename, MovementFilename))

################################################################################
#### Global Model Parameters
################################################################################
n_rsteps         <- 200                    # Number of random steps to be generated
extract_where    <- "end"                  # Extract covariates at the end or along steps?
formula          <- ~ forest + elev + dist # Formula used to predict step-selection scores
maxattempts_mod  <- 20                     # Maximum number of times to repeat an analysis (in case of convergence issues)
maxattempts_imp  <- 100                    # Maximum number of times to try to impute fixes (note that this will be multiplied by the value above)
regular          <- 1                      # Regular step duration

################################################################################
#### Ensure All Functions Work as Intended
################################################################################
# Extract some example data
i <- sample(1:nrow(dat), size = 1)
cov <- dat$Covariates[[i]]
obs <- dat$Movement[[i]]

# Rarify data
obs_missing <- rarifyData(obs, missingness = 0.5)

# Try to impute missing fixes
obs_imputed <- imputeFixes(obs_missing)
ggplot(obs_imputed, aes(x = x, y = y, group = id, col = imputed)) +
  geom_path(lwd = 0.2) +
  geom_point(size = 0.2) +
  coord_sf() +
  theme_minimal() +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5))

# Compute bursts from the (non-imputed) data -> Try with different forgiveness
# values!
obs_bursted <- computeBursts(obs_missing, max_duration = 5 * regular)
ggplot(obs_bursted, aes(x = x, y = y, group = id, col = as.factor(burst))) +
  geom_path(lwd = 0.2) +
  geom_point(size = 0.2) +
  coord_sf() +
  theme_minimal() +
  theme(axis.title.y = element_text(angle = 0, vjust = 0.5), legend.position = "none")

# Compute step metrics
obs_metrics <- computeMetrics(obs_bursted)

# Fit step- and turning-angle distributions to original and rarified data
fit_distris <- fitDists(obs_missing, durations = 1, regular = 1, dynamic = F)
fit_distris <- fitDists(obs_missing, durations = 1:5, regular = 1, dynamic = T, multicore = T)

# Generate random steps
obs_stepsel <- computeSSF(obs_metrics, n_rsteps = 10, approach = "dynamic+model", dists = fit_distris)

# Extract covariates
obs_covaris <- computeCovars(obs_stepsel, cov, multicore = T, extract = "end")

# Run the model
mod <- runModel(obs_covaris, approach = "dynamic+model", covariates = case ~ dist + elev + forest)
mod

# Try to updated tentative model parameters
result <- list(Dists = fit_distris, Coefs = mod)
correctDists(result)

# Clear up space and remov unused objects
rm(
    cov
  , obs
  , mod
  , obs_missing
  , obs_imputed
  , obs_bursted
  , obs_metrics
  , obs_stepsel
  , obs_covaris
  , fit_distris
)
gc()

################################################################################
#### Proper Analysis
################################################################################
# Let's indicate the filenames into which we will store the analysis results
dat$ResultsFilename <- paste0(
    "03_Data/SimulationResults/"
  , "A", sprintf("%02d", dat$AutocorrRange), "_"
  , "R", sprintf("%03d", dat$Replicate)
  , ".rds"
)

# Create the respective folder
dir.create("03_Data/SimulationResults", showWarnings = F, recursive = T)

# For each set of covariates and simulations we will need to loop through
# different treatments. Let's create a dataframe to keep track of those
design <- expand_grid(
    Missingness    = seq(0, 0.5, by = 0.1)                                                 # Fraction of the fixes that is removed
  , Forgiveness    = 1:5                                                                   # Allowed step-duration
  , Approach       = c("uncorrected", "naive", "dynamic+model", "multistep", "imputed")    # Which approach to use to analyse the data
)

# Write the design to file (don't save the covariates and simualated movements
# again though)
write_rds(select(dat, -c(Covariates, Movement)), "03_Data/Analysis.rds")

# Let's randomize the order in which we will go through the different
# combinations
dat <- dat[sample(nrow(dat), replace = F), ]

# Subset to rows that haven't been run yet
dat_sub <- subset(dat, !file.exists(ResultsFilename))

# Loop through the simulations and apply the different analyses
lapply(seq_len(nrow(dat_sub)), function(x) {

  # ################################################################################
  # #### TESTING
  # ################################################################################
  # dat_sub <- subset(dat
  #   , AutocorrRange == 10
  #   & Replicate     == 45
  #   # & Missingness   == 0.4
  #   # & Forgiveness   == 3
  #   # & Approach      == "imputed"
  # )
  # x <- 1
  # ################################################################################
  # #### TESTING
  # ################################################################################

  # Print progress
  cat("Iteration", x, "out of", nrow(dat_sub), "\n")

  # Extract relevant information
  cov <- dat_sub$Covariates[[x]]
  obs <- dat_sub$Movement[[x]]
  fil <- dat_sub$ResultsFilename[x]

  # Analyse the data using the various approaches
  desi <- design

  # ################################################################################
  # #### TESTING
  # ################################################################################
  # desi <- subset(design
  #   # , AutocorrRange == 10
  #   # & Replicate     == 29
  #   , Missingness   == 0.4
  #   & Forgiveness   == 5
  #   & Approach      == "imputed"
  # )
  # y <- 1
  # desi$Results <- NULL
  # ################################################################################
  # #### TESTING
  # ################################################################################
  desi$Results <- pbmclapply(
      X                  = 1:nrow(desi)
    , ignore.interactive = T
    , mc.cores           = detectCores() - 1
    , FUN                = function(y) {

    # Extract relevant design parameters
    miss <- desi$Missingness[[y]]
    forg <- desi$Forgiveness[[y]]
    appr <- desi$Approach[[y]]

    # Reset objects
    coefs   <- NULL
    dis     <- NULL
    obs_sub <- NULL
    obs_ssf <- NULL
    imputed <- NULL

    # The model may not converge. Repeat the analysis if this is the case
    attempts_mod   <- 0
    while (!is.data.frame(coefs) & attempts_mod < maxattempts_mod) {

      # Reset loop objects
      coefs   <- NULL
      dis     <- NULL
      obs_sub <- NULL
      imputed <- NULL

      # Rarify the data
      obs_sub <- rarifyData(obs, missingness = miss)

      # Fit step length distributions
      if (appr == "dynamic+model") {
          dis <- fitDists(obs_sub, dynamic = T, durations = 1:forg, multicore = F, resample = T)
        } else {
          dis <- fitDists(obs_sub, dynamic = F, durations = 1:forg, multicore = F, resample = T)
      }

      # If required, impute fixes
      attempts_imp <- 0
      if (appr == "imputed") {
        while (!is.data.frame(imputed) & attempts_imp < maxattempts_imp) {
          imputed <- tryCatch(imputeFixes(obs_sub)
            , warning = function(w) {return (NA)}
            , error   = function(e) {return (NA)}
          )
          attempts_imp <- attempts_imp + 1
        }
        obs_sub <- imputed
      }

      # If the imputation failed, skip to the next iteration
      if (!is.data.frame(obs_sub)) {
        attempts_mod <- attempts_mod + 1
        next
      }

      # Compute bursts and on those the relevant step metrics
      obs_sub <- computeBursts(obs_sub, max_duration = forg * regular)
      obs_sub <- computeMetrics(obs_sub)

      # Generate random steps
      obs_ssf <- computeSSF(obs_sub
        , n_rsteps = n_rsteps
        , approach = appr
        , dists    = dis
      )

      # Extract covariates
      obs_ssf <- computeCovars(obs_ssf
        , covars    = cov
        , multicore = F
      )

      # Run the model
      coefs <- tryCatch(runModel(obs_ssf, approach = appr, covariates = case ~ dist + elev + forest)
        , warning = function(w) {return (NA)}
        , error   = function(e) {return (NA)}
      )

      # Update number of attempts
      attempts_mod <- attempts_mod + 1
    }

    # Also keep track of the distribution used to draw random steps
    results <- list(
        Dists              = dis
      , Coefs              = coefs
      , ModelAttempts      = attempts_mod
      , ImputationAttempts = attempts_imp
    )

    # Return the results
    return(results)
  })

  # Write results to file
  write_rds(desi, fil)

  # Store results to file
  return(NULL)
})

################################################################################
#### Updating Tentative Parameters
################################################################################
# Reload data
dat <- read_rds("03_Data/Analysis.rds") %>%
  mutate(Data = map(ResultsFilename, read_rds)) %>%
  unnest(Data)

# Check if all is valid
dat <- mutate(dat, Valid = map_lgl(Results, function(x) {
  !any(is.na(x$Coefs$LCI)) & !is.null(x$Coefs)
}))
count(dat, Valid)
# subset(dat, !Valid) %>% pull(ResultsFilename) %>% unique() %>% file.remove()

# Correct tentative parameters
dat <- dat %>%
  mutate(Results = map(Results, function(x) {
    correctDists(x)
  }))

# Create a column that contains "cleaned" results
dat$ResultsCleaned <- pbmclapply(
    X                  = dat$Results
  , ignore.interactive = T
  , mc.cores           = detectCores() - 1
  , FUN                = function(x) {

  # Extract the results (distributions and coefficients)
  dists <- x$Dists
  coefs <- x$Coefs

  # Compile the updated distribution parameters (for step duration of one)
  dist_params <- tibble(
      shape = as.numeric(dists$updated$sl$shape)
    , scale = as.numeric(dists$updated$sl$scale)
    , kappa = as.numeric(dists$updated$ta$kappa)
    , mu    = as.numeric(dists$updated$ta$mu)
    ) %>%
    pivot_longer(1:4, names_to = "Variable", values_to = "Value") %>%
    mutate(Kernel = "Movement")

  # Compile the habitat selection coefficients
  coef_params <- coefs %>%
    select(Variable = Coefficient, Value = Estimate) %>%
    subset(Variable %in% c("forest", "dist", "elev")) %>%
    mutate(Kernel = "Habitat")
  results <- rbind(dist_params, coef_params)
  return(results)
})

# Store everything into a file
write_rds(dat, "03_Data/SimulationResultsConsolidated.rds")

# Store session information
writeLines(capture.output(sessionInfo()), "02_R-Scripts/01_SimulationAnalysis/SessionInfo/02_Analysis.txt")
