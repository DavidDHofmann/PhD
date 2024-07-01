################################################################################
#### Movement Model
################################################################################
# Plots to help with the interpretation of the movement model

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Surpress scientific notation
options(scipen = 999)

# Load required packages
library(tidyverse)    # For data wrangling
library(viridis)      # For nice colors
library(ggpubr)       # To arrange plots
library(metR)         # To be able to get nice colorscales for factorial data
library(ggh4x)        # For nested facettes
library(latex2exp)    # To use latex

# Custom functions
source("02_R-Scripts/00_Functions.R")

# Let's write a function that allows us to predict the relative selection score
# df1 is the reference, df2 is the one where the covariate is varied
predictRSS <- function(model, df1, df2, ci = NULL, return_data = F) {

  # Testing
  # model <- model_prep
  # ci <- 0.95

  # Span model matrices
  df1_mm <- model.matrix(model$formula, df1)[, -1, drop = F]
  df2_mm <- model.matrix(model$formula, df2)[, -1, drop = F]

  # Predict scores
  s1 <- as.vector(exp(df1_mm %*% model$betas))
  s2 <- as.vector(exp(df2_mm %*% model$betas))

  # Compute rss
  rss <- s2 / s1

  # Calculate confidence interval if desired
  if (!is.null(ci)) {

    # Get difference matrix between df1 and df2
    diff <- sweep(as.matrix(df2_mm), 2, as.matrix(df1_mm))
    logrss_se <- as.vector(sqrt((diff %*% model$ses) ** 2))


    # Prepare table to store confidence intervals
    intervals <- lapply(1:length(ci), function(x) {

      # Get critical value
      p <- 1 - ((1 - ci[x]) / 2)
      zstar <- qnorm(p)

      # Compute intervals
      lwr <- exp(log(rss) - zstar * logrss_se)
      upr <- exp(log(rss) + zstar * logrss_se)

      # Prepare dataframe
      cis <- as.data.frame(cbind(lwr, upr))
      names(cis) <- paste0(c("Lower_", "Upper_"), 100 * ci[x])

      # Return them
      return(cis)

    }) %>% do.call(cbind, .)

    rss <- data.frame(RSS = rss, intervals)

  }

  # Return data if desired
  if (return_data){
    rss <- cbind(df1, rss)
  }

  # Return the final data
  return(rss)
}

# Function to generate prediction dataset
predictionSets <- function(covariate, length = 100) {

  # Create reference data
  df1 <- tibble(
      sl                  = dat_ranges$Mean[dat_ranges$Covariate == "sl"]
    , log_sl              = log(sl)
    , cos_ta              = cos(0)
    , Shrubs              = dat_ranges$Mean[dat_ranges$Covariate == "Shrubs"]
    , Water               = dat_ranges$Mean[dat_ranges$Covariate == "Water"]
    , SqrtDistanceToWater = sqrt(dat_ranges$Mean[dat_ranges$Covariate == "DistanceToWater"])
    , Trees               = dat_ranges$Mean[dat_ranges$Covariate == "Trees"]
    , Humans              = dat_ranges$Mean[dat_ranges$Covariate == "Humans"]
    , Temperature         = dat_ranges$Mean[dat_ranges$Covariate == "Temperature"]
  )

  # Create data with varying covariate
  df2 <- df1[rep(1, length), ]
  df2[, covariate] <- seq(
      from       = dat_ranges$Min[dat_ranges$Covariate == covariate]
    , to         = dat_ranges$Max[dat_ranges$Covariate == covariate]
    , length.out = length
  )

  # Scale both with same scaling as the models
  df1 <- scaleCovars(df1, scal)
  df2 <- scaleCovars(df2, scal)

  # Add light types
  df1$LightType <- factor("Bright", levels = c("Bright", "Dark"))
  df2$LightType <- factor("Bright", levels = c("Bright", "Dark"))

  # Return
  return(list(df1 = df1, df2 = df2))
}

# Load our movement models, scaling variables, and tentative gamma
gamma   <- read_rds("03_Data/03_Results/StepLengthDistribution.rds")
forms   <- read_rds("03_Data/03_Results/Formula.rds")
scal    <- read_rds("03_Data/03_Results/Scaling.rds")
model   <- read_rds("03_Data/03_Results/MovementModels.rds") %>% subset(NumberRandomSteps == 25)
sl_dist <- read_rds("03_Data/03_Results/StepLengthDistribution.rds")
ta_dist <- list(
    name   = "vonmises"
  , params = list(kappa = 0)
)

# We can simplify the scaling dataframe to either static or dynamic (they are
# the same parameters anyway)
scal <- scal %>%
  mutate(Covariate = gsub(Covariate, pattern = "Dynamic|Static", replacement = "")) %>%
  distinct()

# Load original (unscaled) data used to fit the model
dat_orig <- "03_Data/02_CleanData/SSFExtracted.rds" %>%
  read_rds() %>%
  subset(case == 1) %>%
  mutate(
      SqrtDistanceToWaterStatic = sqrt(DistanceToWaterStatic)
    , SqrtDistanceToWaterDynamic = sqrt(DistanceToWaterDynamic)
  )

# Find ranges for observed values
dat_ranges <- dat_orig %>%
  dplyr::select(where(is.numeric)) %>%
  pivot_longer(everything(), names_to = "Covariate", values_to = "Value") %>%
  mutate(Covariate = gsub(Covariate, pattern = "Static", replacement = "")) %>%
  mutate(Covariate = gsub(Covariate, pattern = "Dynamic", replacement = "")) %>%
  group_by(Covariate) %>%
  summarize(
      Min  = min(Value, na.rm = T)
    , Max  = max(Value, na.rm = T)
    , Mean = mean(Value, na.rm = T)
  )

# Design through which to loop
design <- expand_grid(
    FittingCovariates = factor(c("Static", "Dynamic"), levels = c("Static", "Dynamic"))
  , Season            = c("All", "Dry", "Wet")
  , Formula           = factor(c("Simple", "Complex"), levels = c("Simple", "Complex"))
  , Covariate         = c("Water", "Shrubs", "SqrtDistanceToWater", "Trees", "Humans")
)

# Loop through design
design$Predictions <- lapply(1:nrow(design), function(i) {

  # Get appropriate formula and extract name of covariate
  frm  <- forms[[design$Formula[i]]]
  covi <- design$Covariate[i]

  # Get appropriate model
  mod <- subset(model
    , FittingCovariates == design$FittingCovariates[i]
    & Season            == design$Season[i]
    & Formula           == ifelse(design$Formula[i] == "Complex", "Full", "Simple")
  ) %>% prepareModel(frm, .)

  # Prepare datasets, run predictions, and return everything
  dfs <- predictionSets(covariate = covi, length = 100)
  rss <- predictRSS(model = mod, df1 = dfs$df1, df2 = dfs$df2, ci = 0.95)

  # Backtransform the covariate
  covi_back <- dfs$df2[, covi] * scal$scale[scal$Covariate == covi] +
    scal$center[scal$Covariate == covi]
  result <- cbind(Value = covi_back, rss)
})

# Unnest
design <- design %>%
  unnest(Predictions) %>%
  mutate(
    Covariate = factor(Covariate
      , levels = c("Water", "SqrtDistanceToWater", "Shrubs", "Trees", "Humans")
      , labels = TeX(c("Water", "DistanceToWater$^{0.5}$", "Shrubs", "Trees", "Humans"))
    )
    , Formula = factor(Formula
        , levels = c("Simple", "Complex")
        , labels = c("Simple~Formula", "Complex~Formula")
    )
  )

# Prepare plot
p <- ggplot(design, aes(x = Value, y = RSS, color = Season, fill = Season)) +
  geom_line() +
  facet_nested(Formula + FittingCovariates ~ Covariate, scales = "free_x", labeller = label_parsed) +
  scale_color_viridis_d(begin = 0.3, name = "Season", direction = -1) +
  scale_fill_viridis_d(begin = 0.3, name = "Season", direction = -1) +
  theme_awesome() +
  coord_cartesian(ylim = c(0, 2)) +
  theme(
      axis.title.y     = element_text(angle = 90)
    , axis.text        = element_text(size = 6)
    , panel.grid.minor = element_blank()
  )

# Store the plot to file
ggsave("04_Manuscript/Figures/MovementModelInterpretation.png"
  , plot   = p
  , device = png
  , bg     = "white"
  , width  = 8
  , height = 5
  , scale  = 1
)
