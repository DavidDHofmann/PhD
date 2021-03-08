################################################################################
#### Movement Model
################################################################################
# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Surpress scientific notation
options(scipen = 999)

# Load required packages
library(tidyverse)    # For data wrangling
library(glmmTMB)      # To handle glmm models
library(davidoff)     # Custom functions
library(amt)          # For fitting distribution
library(lemon)        # For nice capped coords

# Load our two models (scaled and unscaled)
model_unscaled <- read_rds("03_Data/03_Results/99_MovementModelUnscaled.rds")
model_scaled   <- read_rds("03_Data/03_Results/99_MovementModel.rds")
model_scaled   <- model_scaled$Model[[1]]

################################################################################
#### Function to Predict the RSS
################################################################################
# Let's write a function that allows us to predict the relative selection score
predictRSS <- function(model, df1, df2, ci = NULL, return_data = F){

  # Use custom function to prepare model
  mod <- prepareModel(model)

  # Span model matrices
  df1_mm <- as.data.frame(model.matrix(mod$formula, df1))
  df2_mm <- as.data.frame(model.matrix(mod$formula, df2))

  # Remove intercept
  df1_mm <- df1_mm[, names(df1_mm) != "(Intercept)"]
  df2_mm <- df2_mm[, names(df2_mm) != "(Intercept)"]

  # Predict scores
  s1 <- predictScore(
      coefficients = mod$coefficients
    , formula      = mod$formula
    , data         = df1
  )
  s2 <- predictScore(
      coefficients = mod$coefficients
    , formula      = mod$formula
    , data         = df2
  )

  # Calculate relative selection scores
  rss <- s1 / s2

  # Calculate confidence interval if desired
  if (!is.null(ci)){

    # Get model variance matrix
    vars <- vcov(model)[[1]]

    # Remove intercept
    vars <- vars[rownames(vars) != "(Intercept)", colnames(vars) != "(Intercept)"]

    # Get difference matrix between df1 and df2
    diff <- sweep(as.matrix(df1_mm), 2, as.matrix(df2_mm))

    # Get variance of log-rss prediction
    varpred <- diag(diff %*% vars %*% t(diff))

    # Get the standard error
    logrss_se <- unname(sqrt(varpred))

    # Prepare table to store confidence intervals
    intervals <- lapply(1:length(ci), function(x){

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

################################################################################
#### Turning Angle
################################################################################
df1 <- data.frame(
    sl_                 = 1
  , cos_ta_             = seq(-2, 2, length.out = 1000)
  , log_sl_             = 0
  , Shrubs              = 0
  , Water               = 0
  , SqrtDistanceToWater = 0
  , Trees               = 0
  , HumansBuff5000      = 0
  , inactive            = F
)
df2 <- data.frame(
    sl_                 = 1
  , cos_ta_             = 0
  , log_sl_             = 0
  , Shrubs              = 0
  , Water               = 0
  , SqrtDistanceToWater = 0
  , Trees               = 0
  , HumansBuff5000      = 0
  , inactive            = F
)

# Predict scores
pred <- predictRSS(
    model       = model_unscaled
  , df1         = df1
  , df2         = df2
  , ci          = c(0.99, 0.95, 0.9)
  , return_data = T
)

# Make tidy
pred <- pred %>%
  gather(key = Interval, value = Boundary, Lower_99:Upper_90) %>%
  separate(Interval, into = c("Type", "Level"), sep = "_") %>%
  spread(key = Type, value = Boundary)

# Visualize
ggplot(pred, aes(x = cos_ta_, y = RSS)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, group = Level)
    , linetype  = "solid"
    , alpha     = 0.33
    , color     = "orange"
    , fill      = "orange"
    , lwd = 0.1
  ) +
  geom_line(size = 1) +
  xlab("cos(ta) (SD)") +
  ylab("RSS vs cos(ta)") +
  theme_classic()

################################################################################

################################################################################
#### Water
################################################################################
df1 <- data.frame(
    sl_                 = 0.001
  , cos_ta_             = 0
  , log_sl_             = log(0.001)
  , Shrubs              = 0
  , Water               = seq(-2, 2, length.out = 1000)
  , SqrtDistanceToWater = 0
  , Trees               = 0
  , HumansBuff5000      = 0
  , inactive            = F
)
df2 <- data.frame(
    sl_                 = 0.001
  , cos_ta_             = 0
  , log_sl_             = log(0.001)
  , Shrubs              = 0
  , Water               = 0
  , SqrtDistanceToWater = 0
  , Trees               = 0
  , HumansBuff5000      = 0
  , inactive            = F
)

# Predict scores
pred <- predictRSS(
    model       = model_unscaled
  , df1         = df1
  , df2         = df2
  , ci          = c(0.99, 0.95, 0.9)
  , return_data = T
)

# Make tidy
pred <- pred %>%
  gather(key = Interval, value = Boundary, Lower_99:Upper_90) %>%
  separate(Interval, into = c("Type", "Level"), sep = "_") %>%
  spread(key = Type, value = Boundary)

# Visualize
ggplot(pred, aes(x = Water, y = RSS)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, group = Level)
    , linetype  = "solid"
    , alpha     = 0.33
    , color     = "orange"
    , fill      = "orange"
    , lwd = 0.1
  ) +
  geom_line(size = 1) +
  xlab("Water (SD)") +
  ylab("RSS vs Water") +
  theme_classic()

################################################################################
#### Distance To Water
################################################################################
df1 <- data.frame(
    sl_                 = 0.001
  , cos_ta_             = 0
  , log_sl_             = log(0.001)
  , Shrubs              = 0
  , Water               = 0
  , SqrtDistanceToWater = seq(-2, 2, length.out = 1000)
  , Trees               = 0
  , HumansBuff5000      = 0
  , inactive            = F
)
df2 <- data.frame(
    sl_                 = 0.001
  , cos_ta_             = 0
  , log_sl_             = log(0.001)
  , Shrubs              = 0
  , Water               = 0
  , SqrtDistanceToWater = 0
  , Trees               = 0
  , HumansBuff5000      = 0
  , inactive            = F
)

# Predict scores
pred <- predictRSS(
  model         = model_scaled
  , df1         = df1
  , df2         = df2
  , ci          = c(0.99, 0.95, 0.9)
  , return_data = T
)

# Make tidy
pred <- pred %>%
  gather(key = Interval, value = Boundary, Lower_99:Upper_90) %>%
  separate(Interval, into = c("Type", "Level"), sep = "_") %>%
  spread(key = Type, value = Boundary)

# Visualize
ggplot(pred, aes(x = SqrtDistanceToWater, y = RSS)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, group = Level)
    , linetype  = "solid"
    , alpha     = 0.33
    , color     = "orange"
    , fill      = "orange"
    , lwd = 0.1
  ) +
  geom_line(size = 1) +
  xlab("SqrtDistanceToWater (SD)") +
  ylab("RSS vs SqrtDistanceToWater") +
  theme_classic()

################################################################################
#### Shrubs
################################################################################
df1 <- data.frame(
    sl_                 = 0.001
  , cos_ta_             = 0
  , log_sl_             = log(0.001)
  , Shrubs              = seq(-2, 2, length.out = 1000)
  , Water               = 0
  , SqrtDistanceToWater = 0
  , Trees               = 0
  , HumansBuff5000      = 0
  , inactive            = F
)
df2 <- data.frame(
    sl_                 = 0.001
  , cos_ta_             = 0
  , log_sl_             = log(0.001)
  , Shrubs              = 0
  , Water               = 0
  , SqrtDistanceToWater = 0
  , Trees               = 0
  , HumansBuff5000      = 0
  , inactive            = F
)

# Predict scores
pred <- predictRSS(
  model         = model_scaled
  , df1         = df1
  , df2         = df2
  , ci          = c(0.99, 0.95, 0.9)
  , return_data = T
)

# Make tidy
pred <- pred %>%
  gather(key = Interval, value = Boundary, Lower_99:Upper_90) %>%
  separate(Interval, into = c("Type", "Level"), sep = "_") %>%
  spread(key = Type, value = Boundary)

# Visualize
ggplot(pred, aes(x = Shrubs, y = RSS)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, group = Level)
    , linetype  = "solid"
    , alpha     = 0.33
    , color     = "orange"
    , fill      = "orange"
    , lwd = 0.1
  ) +
  geom_line(size = 1) +
  xlab("Shrubs (SD)") +
  ylab("RSS vs Shrubs") +
  theme_classic()

################################################################################
#### Trees
################################################################################
df1 <- data.frame(
    sl_                 = 0.001
  , cos_ta_             = 0
  , log_sl_             = log(0.001)
  , Shrubs              = 0
  , Water               = 0
  , SqrtDistanceToWater = 0
  , Trees               = seq(-2, 2, length.out = 1000)
  , HumansBuff5000      = 0
  , inactive            = F
)
df2 <- data.frame(
    sl_                 = 0.001
  , cos_ta_             = 0
  , log_sl_             = log(0.001)
  , Shrubs              = 0
  , Water               = 0
  , SqrtDistanceToWater = 0
  , Trees               = 0
  , HumansBuff5000      = 0
  , inactive            = F
)

# Predict scores
pred <- predictRSS(
  model         = model_scaled
  , df1         = df1
  , df2         = df2
  , ci          = c(0.99, 0.95, 0.9)
  , return_data = T
)

# Make tidy
pred <- pred %>%
  gather(key = Interval, value = Boundary, Lower_99:Upper_90) %>%
  separate(Interval, into = c("Type", "Level"), sep = "_") %>%
  spread(key = Type, value = Boundary)

# Visualize
ggplot(pred, aes(x = Trees, y = RSS)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, group = Level)
    , linetype  = "solid"
    , alpha     = 0.33
    , color     = "orange"
    , fill      = "orange"
    , lwd = 0.1
  ) +
  geom_line(size = 1) +
  xlab("Trees (SD)") +
  ylab("RSS vs Trees") +
  theme_classic()

################################################################################
#### Humans
################################################################################
df1 <- data.frame(
    sl_                 = 0.001
  , cos_ta_             = 0
  , log_sl_             = log(0.001)
  , Shrubs              = 0
  , Water               = 0
  , SqrtDistanceToWater = 0
  , Trees               = 0
  , HumansBuff5000      = seq(-2, 2, length.out = 1000)
  , inactive            = F
)
df2 <- data.frame(
    sl_                 = 0.001
  , cos_ta_             = 0
  , log_sl_             = log(0.001)
  , Shrubs              = 0
  , Water               = 0
  , SqrtDistanceToWater = 0
  , Trees               = 0
  , HumansBuff5000      = 0
  , inactive            = F
)

# Predict scores
pred <- predictRSS(
  model         = model_scaled
  , df1         = df1
  , df2         = df2
  , ci          = c(0.99, 0.95, 0.9)
  , return_data = T
)

# Make tidy
pred <- pred %>%
  gather(key = Interval, value = Boundary, Lower_99:Upper_90) %>%
  separate(Interval, into = c("Type", "Level"), sep = "_") %>%
  spread(key = Type, value = Boundary)

# Visualize
ggplot(pred, aes(x = HumansBuff5000, y = RSS)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, group = Level)
    , linetype  = "solid"
    , alpha     = 0.33
    , color     = "orange"
    , fill      = "orange"
    , lwd = 0.1
  ) +
  geom_line(size = 1) +
  xlab("HumansBuff5000 (SD)") +
  ylab("RSS vs HumansBuff5000") +
  theme_classic()

################################################################################
#### Step Length: Active vs Inactive
################################################################################
# Load tentative gamma distribution and adjust it to kilometers
sl_dist <- read_rds("03_Data/03_Results/99_GammaDistribution.rds")
sl_dist$params$scale <- sl_dist$params$scale / 1000

# Extract model parameters
coefs <- fixef(model_unscaled)$cond

# Step length distribution during activity
sl_active <- update_gamma(sl_dist
  , beta_sl     = coefs["sl_"]
  , beta_log_sl = coefs["log_sl_"]
)

# Step length distribution during inactivity
sl_inactive <- update_gamma(sl_dist
  , beta_sl     = coefs["sl_"] + coefs["sl_:inactiveTRUE"]
  , beta_log_sl = coefs["log_sl_"]
)

# Prepare dataframe for plot
plot_sl <- data.frame(sl_ = seq(from = 0.0, to = 35, length.out = 1000))
plot_sl$active <- dgamma(plot_sl$sl_
  , shape = sl_active$params$shape
  , scale = sl_active$params$scale
)
plot_sl$inactive <- dgamma(plot_sl$sl_
  , shape = sl_inactive$params$shape
  , scale = sl_inactive$params$scale
)
plot_sl <- pivot_longer(plot_sl, cols = -sl_)

# Visualize
ggplot(plot_sl, aes(x = sl_, y = value, color = factor(name))) +
  geom_line(size = 1) +
  scale_color_manual(values = c("green", "cornflowerblue")) +
  theme_classic()

################################################################################
#### Step Length: Water vs Dry
################################################################################
# Step length distribution in water
sl_dry <- update_gamma(sl_dist
  , beta_sl     = coefs["sl_"] + coefs["sl_:Water"] * -2
  , beta_log_sl = coefs["log_sl_"]
)

# Step length distribution in medium water
sl_medium <- update_gamma(sl_dist
  , beta_sl     = coefs["sl_"]
  , beta_log_sl = coefs["log_sl_"]
)

# Step length distribution in water
sl_water <- update_gamma(sl_dist
  , beta_sl     = coefs["sl_"] + coefs["sl_:Water"] * 2
  , beta_log_sl = coefs["log_sl_"]
)

# Prepare dataframe for plot
plot_sl <- data.frame(sl_ = seq(from = 0.0, to = 35, length.out = 1000))
plot_sl$dry <- dgamma(plot_sl$sl_
  , shape = sl_dry$params$shape
  , scale = sl_dry$scale
)
plot_sl$medium <- dgamma(plot_sl$sl_
  , shape = sl_medium$params$shape
  , scale = sl_medium$params$scale
)
plot_sl$water <- dgamma(plot_sl$sl_
  , shape = sl_water$params$shape
  , scale = sl_water$params$scale
)
plot_sl <- pivot_longer(plot_sl, cols = -sl_)

# Visualize
ggplot(plot_sl, aes(x = sl_, y = value, color = factor(name))) +
  geom_line(size = 1) +
  scale_color_manual(values = c("green", "cornflowerblue")) +
  theme_classic()
