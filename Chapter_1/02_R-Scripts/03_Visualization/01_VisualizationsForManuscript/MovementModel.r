################################################################################
#### Movement Model
################################################################################
# Interpreting the movement model. I have fitted two different model. One were
# all covariates (habitat AND movement covariates) have been scaled and one were
# only habitat covariates have been scaled, yet the step length was converted to
# kilometers. We'll use the first model to interpret the habitat kernel, yet the
# second model to investigate the movement kernel.

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
library(skimr)        # To skim through a dataframe
library(viridis)      # For nice colors

# Load our movement models
model1 <- read_rds("03_Data/03_Results/99_MovementModel.rds")$Model[[1]]
model2 <- read_rds("03_Data/03_Results/99_MovementModelUnscaled.rds")

# Load the scaling parameters
scaling <- read_rds("03_Data/03_Results/99_Scaling.rds")
sc_cent <- scaling$center
sc_scal <- scaling$scale

# Load the tentative gamma distribution and adjust it to kilometers
sl_dist <- read_rds("03_Data/03_Results/99_GammaDistribution.rds")
sl_dist$params$scale <- sl_dist$params$scale / 1000

# Prepare a turning angle distribution
ta_dist <- list(
    name   = "vonmises"
  , params = list(kappa = 0)
)

# Extract data from the movement model and unscale (un) the data
modeldat <- model1 %>%
  model.frame() %>%
  subset(case_) %>%
  dplyr::select(-c(case_, step_id_, id, inactive)) %>%
  mutate(
      sl_un = sl_ * sc_scal["sl_"] + sc_cent["sl_"]
    , log_sl_un = log_sl_ * sc_scal["log_sl_"] + sc_cent["log_sl_"]
    , cos_ta_un = cos_ta_ * sc_scal["cos_ta_"] + sc_cent["cos_ta_"]
    , Water_un = Water * sc_scal["Water"] + sc_cent["Water"]
    , SqrtDistanceToWater_un = SqrtDistanceToWater * sc_scal["SqrtDistanceToWater"] + sc_cent["SqrtDistanceToWater"]
    , Shrubs_un = Shrubs * sc_scal["Shrubs"] + sc_cent["Shrubs"]
    , Trees_un = Trees * sc_scal["Trees"] + sc_cent["Trees"]
    , HumansBuff5000_un = Trees * sc_scal["HumansBuff5000"] + sc_cent["HumansBuff5000"]
  ) %>%
  mutate_all(as.numeric) %>%
  dplyr::select(sort(names(.)))

# Identify the range on which we observed each covariate (scaled and unscaled)
ranges <- modeldat %>%
  gather(key = Covariate, value = Value) %>%
  group_by(Covariate) %>%
  summarize(Min = min(Value), Max = max(Value), Center = (Min + Max) / 2)

################################################################################
#### Turning Angle vs Step Length
################################################################################
# Check out model to see on which variables the effect of turning angle depends
summary(model2)

# Extract model Coefficients
coefs <- fixef(model2)$cond

# Sequence for different step lengths
seq_sl_ <- seq(
    ranges$Min[ranges$Covariate == "sl_un"]
  , ranges$Max[ranges$Covariate == "sl_un"]
  , length.out = 100
) / 1000

# Show turning angle for different values of sl_
dat <- lapply(seq_sl_, function(x){

  # Calculate updated vonmises distribution
  updated <- update_vonmises(dist = ta_dist
    , beta_cos_ta = coefs["cos_ta_"] +
      coefs["cos_ta_:sl_"] *
        x +
      coefs["cos_ta_:SqrtDistanceToWater"] *
        ranges$Center[ranges$Covariate == "SqrtDistanceToWater"] +
      coefs["cos_ta_:HumansBuff5000"] *
        ranges$Center[ranges$Covariate == "HumansBuff5000"]
  )

  # Prepare dataframe for plot
  plot_ta <- data.frame(ta_ = seq(from = -pi, to = +pi, length.out = 1000))

  # Insert the step length
  plot_ta$sl_ <- x

  # Get probabilities from updated distribution
  plot_ta$prob <- circular::dvonmises(plot_ta$ta_
    , kappa = updated$params$kappa
    , mu    = 0
  )

  # Return the final data
  return(plot_ta)

}) %>% do.call(rbind, .)

# Visualize using contour
ggplot(dat, aes(x = ta_, y = sl_, z = prob)) +
  geom_contour_filled() +
  geom_contour(col = "black") +
  scale_fill_viridis_d(option = "magma") +
  theme_classic() +
  coord_capped_cart(left = "both", bottom = "both") +
  labs(fill = "Probability")

# Visualize using ggplot
ggplot(dat, aes(x = ta_, y = prob, color = factor(sl_))) +
  geom_line(size = 1) +
  theme_classic() +
  scale_color_viridis_d(begin = 0.2, end = 0.8, option = "magma") +
  xlab("Turning Angle") +
  ylab("Probability Density") +
  labs(color = "Step Length") +
  theme(legend.position = "none")

################################################################################
#### Turning Angle vs SqrtDistanceToWater
################################################################################
# Sequence for different distances to water
seq_wat_ <- seq(
    ranges$Min[ranges$Covariate == "SqrtDistanceToWater"]
  , ranges$Max[ranges$Covariate == "SqrtDistanceToWater"]
  , length.out = 100
)

# Show turning angle for different values of sl_
dat <- lapply(seq_wat_, function(x){

  # Calculate updated vonmises distribution
  updated <- update_vonmises(dist = ta_dist
    , beta_cos_ta = coefs["cos_ta_"] +
      coefs["cos_ta_:sl_"] *
        ranges$Center[ranges$Covariate == "sl_un"] / 1000 +
      coefs["cos_ta_:SqrtDistanceToWater"] *
        x +
      coefs["cos_ta_:HumansBuff5000"] *
        ranges$Center[ranges$Covariate == "HumansBuff5000"]
  )

  # Prepare dataframe for plot
  plot_ta <- data.frame(ta_ = seq(from = -pi, to = +pi, length.out = 1000))

  # Insert the water cover
  plot_ta$SqrtDistanceToWater <- x

  # Get probabilities from updated distribution
  plot_ta$prob <- circular::dvonmises(plot_ta$ta_
    , kappa = updated$params$kappa
    , mu    = 0
  )

  # Return the final data
  return(plot_ta)

}) %>% do.call(rbind, .)

# Visualize using contour
ggplot(dat, aes(x = ta_, y = SqrtDistanceToWater, z = prob)) +
  geom_contour_filled() +
  geom_contour(col = "black") +
  scale_fill_viridis_d(option = "magma") +
  theme_classic() +
  coord_capped_cart(left = "both", bottom = "both") +
  labs(fill = "Probability")

# Visualize using ggplot
ggplot(dat, aes(x = ta_, y = prob, color = factor(SqrtDistanceToWater))) +
  geom_line(size = 1) +
  theme_classic() +
  scale_color_viridis_d(begin = 0.2, end = 0.8, option = "magma") +
  xlab("Turning Angle") +
  ylab("Probability Density") +
  labs(color = "SqrtDistanceToWater") +
  theme(legend.position = "none")

################################################################################
#### Turning Angle vs Humans
################################################################################
# Sequence for different human influences
seq_hum_ <- seq(
    ranges$Min[ranges$Covariate == "HumansBuff5000"]
  , ranges$Max[ranges$Covariate == "HumansBuff5000"]
  , length.out = 100
)

# Show turning angle for different values of sl_
dat <- lapply(seq_hum_, function(x){

  # Calculate updated vonmises distribution
  updated <- update_vonmises(dist = ta_dist
    , beta_cos_ta = coefs["cos_ta_"] +
      coefs["cos_ta_:sl_"] *
        ranges$Center[ranges$Covariate == "sl_un"] / 1000 +
      coefs["cos_ta_:SqrtDistanceToWater"] *
        ranges$Center[ranges$Covariate == "SqrtDistanceToWater"] +
      coefs["cos_ta_:HumansBuff5000"] *
        x
  )

  # Prepare dataframe for plot
  plot_ta <- data.frame(ta_ = seq(from = -pi, to = +pi, length.out = 1000))

  # Insert the water cover
  plot_ta$HumansBuff5000 <- x

  # Get probabilities from updated distribution
  plot_ta$prob <- circular::dvonmises(plot_ta$ta_
    , kappa = updated$params$kappa
    , mu    = 0
  )

  # Return the final data
  return(plot_ta)

}) %>% do.call(rbind, .)

# Visualize using contour
ggplot(dat, aes(x = ta_, y = HumansBuff5000, z = prob)) +
  geom_contour_filled() +
  geom_contour(col = "black") +
  scale_fill_viridis_d(option = "magma") +
  theme_classic() +
  coord_capped_cart(left = "both", bottom = "both") +
  labs(fill = "Probability")

# Visualize using lines
ggplot(dat, aes(x = ta_, y = prob, color = factor(SqrtDistanceToWater))) +
  geom_line(size = 1) +
  theme_classic() +
  scale_color_viridis_d(begin = 0.2, end = 0.8, option = "magma") +
  xlab("Turning Angle") +
  ylab("Probability Density") +
  labs(color = "SqrtDistanceToWater") +
  theme(legend.position = "none")

################################################################################
#### Step Length vs Water
################################################################################
# Check out model to see on which variables the effect of step lengths depends
summary(model2)

# Sequence for different distances to water
seq_wat_ <- seq(
    ranges$Min[ranges$Covariate == "SqrtDistanceToWater"]
  , ranges$Max[ranges$Covariate == "SqrtDistanceToWater"]
  , length.out = 100
)

# Show turning angle for different values of sl_
dat <- lapply(seq_wat_, function(x){

  # Calculate updated vonmises distribution
  updated <- update_gamma(dist = sl_dist
    , beta_sl = coefs["sl_"] +
      coefs["sl_:inactiveTRUE"] *
        0 +
      coefs["sl_:Water"] *
        x +
      coefs["sl_:Trees"] *
        ranges$Center[ranges$Covariate == "Trees"] +
      coefs["sl_:Shrubs"] *
        ranges$Center[ranges$Covariate == "Shrubs"] +
      coefs["sl_:SqrtDistanceToWater"] *
        ranges$Center[ranges$Covariate == "SqrtDistanceToWater"]
    , beta_log_sl = coefs["log_sl_"] +
      coefs["cos_ta_:log_sl_"] *
        ranges$Center[ranges$Covariate == "cos_ta_un"]
  )

  # Prepare dataframe for plot
  plot_ta <- data.frame(ta_ = seq(from = -pi, to = +pi, length.out = 1000))

  # Insert the water cover
  plot_ta$SqrtDistanceToWater <- x

  # Get probabilities from updated distribution
  plot_ta$prob <- circular::dvonmises(plot_ta$ta_
    , kappa = updated$params$kappa
    , mu    = 0
  )

  # Return the final data
  return(plot_ta)

}) %>% do.call(rbind, .)

# Visualize using contour
ggplot(dat, aes(x = ta_, y = SqrtDistanceToWater, z = prob)) +
  geom_contour_filled() +
  geom_contour(col = "black") +
  scale_fill_viridis_d(option = "magma") +
  theme_classic() +
  coord_capped_cart(left = "both", bottom = "both") +
  labs(fill = "Probability")

# Visualize using ggplot
ggplot(dat, aes(x = ta_, y = prob, color = factor(SqrtDistanceToWater))) +
  geom_line(size = 1) +
  theme_classic() +
  scale_color_viridis_d(begin = 0.2, end = 0.8, option = "magma") +
  xlab("Turning Angle") +
  ylab("Probability Density") +
  labs(color = "SqrtDistanceToWater") +
  theme(legend.position = "none")






















# Check out model
summary(model2)

# Load tentative gamma distribution and adjust it to kilometers
sl_dist <- read_rds("03_Data/03_Results/99_GammaDistribution.rds")
sl_dist$params$scale <- sl_dist$params$scale / 1000
coefs <- fixef(model2)$cond

# Step length distribution on dryland
sl_low <- update_gamma(sl_dist
  , beta_sl = coefs["sl_"] +
      coefs["sl_:Trees"] * ranges$Min[ranges$Covariate == "Trees"] +
      coefs["cos_ta_:sl_"] * ranges$Center[ranges$Covariate == "cos_ta_"]
  , beta_log_sl = coefs["log_sl_"]
)

# Step length distribution in medium water
sl_medium <- update_gamma(sl_dist
  , beta_sl = coefs["sl_"] +
      coefs["sl_:Trees"] * ranges$Center[ranges$Covariate == "Trees"] +
      coefs["cos_ta_:sl_"] * ranges$Center[ranges$Covariate == "cos_ta_"]
  , beta_log_sl = coefs["log_sl_"]
)

# Step length distribution in water
sl_high <- update_gamma(sl_dist
  , beta_sl = coefs["sl_"] +
      coefs["sl_:Trees"] * ranges$Max[ranges$Covariate == "Trees"] +
      coefs["cos_ta_:sl_"] * ranges$Center[ranges$Covariate == "cos_ta_"]
  , beta_log_sl = coefs["log_sl_"]
)

# Prepare dataframe for plot
plot_sl <- data.frame(sl_ = seq(from = 0.0, to = 35, length.out = 1000))
plot_sl$low <- dgamma(plot_sl$sl_
  , shape = sl_low$params$shape
  , scale = sl_low$params$scale
)
plot_sl$medium <- dgamma(plot_sl$sl_
  , shape = sl_medium$params$shape
  , scale = sl_medium$params$scale
)
plot_sl$high <- dgamma(plot_sl$sl_
  , shape = sl_high$params$shape
  , scale = sl_high$params$scale
)
plot_sl <- pivot_longer(plot_sl, cols = -sl_)
plot_sl$name <- factor(plot_sl$name, levels = c("high", "medium", "low"))

# Visualize
ggplot(plot_sl, aes(x = 1000 * sl_, y = value, color = name)) +
  geom_line(size = 1) +
  theme_classic() +
  scale_color_viridis_d(begin = 0.2, end = 0.8) +
  scale_y_sqrt(limits = c(0, 0.25)) +
  xlab("Step Length (m)") +
  ylab("Probability Density") +
  labs(color = "Tree Cover")
library(ggpubr)
ggarrange(p1, p2, nrow = 2)
# # Load observed steps
# steps <- read_csv("03_Data/02_CleanData/00_General_Dispersers_POPECOL(SSF_Extracted).csv")
# steps <- subset(steps, case_)
#
# # Categorize steps into small, medium, and large steps
# steps <- mutate(steps, step_size = case_when(
#     sl_ <= quantile(steps$sl_, 0.33) ~ "small"
#   , sl_ > quantile(steps$sl_, 0.33) & sl_ <= quantile(steps$sl_, 0.66) ~ "medium"
#   , sl_ > quantile(steps$sl_, 0.66) ~ "large"
# ))
#
# # Calculate means per group
# mean_steps <- steps %>%
#   group_by(step_size) %>%
#   summarize(sl_ = round(mean(sl_))) %>%
#   mutate(log_sl_ = log(sl_))
#
# # Scale them
# mean_steps <- mutate(mean_steps
#   , sl_s = scale(sl_
#     , center = scaling$center[["sl_"]]
#     , scale  = scaling$scale[["sl_"]]
#   )
#   , log_sl_s = scale(log_sl_
#     , center = scaling$center[["log_sl_"]]
#     , scale  = scaling$scale[["log_sl_"]]
#   )
# )

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

# df1 <- expand_grid(
#     sl_                 = mean_steps$sl_s[mean_steps$step_size == "medium"]
#   , cos_ta_             = seq(-2, 2, length.out = 100)
#   , log_sl_             = mean_steps$log_sl_s[mean_steps$step_size == "medium"]
#   , Shrubs              = 0
#   , Water               = 0
#   , SqrtDistanceToWater = 0
#   , Trees               = 0
#   , HumansBuff5000      = seq(-2, 2, length.out = 100)
#   , inactive            = F
# )
# df2 <- expand_grid(
#     sl_                 = mean_steps$sl_s[mean_steps$step_size == "medium"]
#   , cos_ta_             = 0
#   , log_sl_             = mean_steps$log_sl_s[mean_steps$step_size == "medium"]
#   , Shrubs              = 0
#   , Water               = 0
#   , SqrtDistanceToWater = 0
#   , Trees               = 0
#   , HumansBuff5000      = 0
#   , inactive            = F
# )
# library(raster)
# library(viridis)
# library(rasterVis)
# test <- predictRSS(model, df1, df2, ci = NULL, return_data = T)
# test <- test[, c("HumansBuff5000", "cos_ta_", "rss")]
# test$HumansBuff5000 <- test$HumansBuff5000 * scaling$scale["HumansBuff5000"] + scaling$center["HumansBuff5000"]
# test$cos_ta_ <- test$cos_ta_ * scaling$scale["cos_ta_"] + scaling$center["cos_ta_"]
# test <- rasterFromXYZ(test)
# plot(test, col = viridis(50), xlab = "HumanBuff5000", ylab = "cos_ta_")
# levelplot(test, xlab = "HumanBuff5000", ylab = "cos_ta_")
# contour(test, add = T)
# as.matrix(spread(test, cos_ta_, rss))

################################################################################
#### Water
################################################################################
df1 <- data.frame(
    sl_                 = mean_steps$sl_s[mean_steps$step_size == "medium"]
  , cos_ta_             = 0
  , log_sl_             = mean_steps$log_sl_s[mean_steps$step_size == "medium"]
  , Shrubs              = 0
  , Water               = seq(-2, 2, length.out = 1000)
  , SqrtDistanceToWater = 0
  , Trees               = 0
  , HumansBuff5000      = 0
  , inactive            = F
)
df2 <- data.frame(
    sl_                 = mean_steps$sl_s[mean_steps$step_size == "medium"]
  , cos_ta_             = 0
  , log_sl_             = mean_steps$log_sl_s[mean_steps$step_size == "medium"]
  , Shrubs              = 0
  , Water               = 0
  , SqrtDistanceToWater = 0
  , Trees               = 0
  , HumansBuff5000      = 0
  , inactive            = F
)

# Predict scores
pred <- predictRSS(
    model       = model
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
#### Turning Angle
################################################################################
# Check the model to see on what the effect of the turning angle depends
summary(model)

# Its influence depends on the values of "sl_", "HumansBuff5000",
# and "SqrtDistanceToWater". Let's prepare a dataframe for this
grid <- expand_grid(
    step_size           = mean_steps$step_size[2]
  , HumansBuff5000      = c(-2, 0, 2)
  , SqrtDistanceToWater = c(-2, 0, 2)[2]
)

# Run prediction for three different step sizes
preds <- lapply(1:nrow(grid), function(x){

  # Prepare data frames
  df1 <- data.frame(
      sl_                 = mean_steps$sl_s[mean_steps$step_size == grid$step_size[x]]
    , cos_ta_             = seq(-2, 2, length.out = 1000)
    , log_sl_             = mean_steps$sl_s[mean_steps$step_size == grid$step_size[x]]
    , Shrubs              = 0
    , Water               = 0
    , SqrtDistanceToWater = grid$SqrtDistanceToWater[x]
    , Trees               = 0
    , HumansBuff5000      = grid$HumansBuff5000[x]
    , inactive            = F
  )
  df2 <- data.frame(
      sl_                 = mean_steps$sl_s[mean_steps$step_size == grid$step_size[x]]
    , cos_ta_             = 0
    , log_sl_             = mean_steps$sl_s[mean_steps$step_size == grid$step_size[x]]
    , Shrubs              = 0
    , Water               = 0
    , SqrtDistanceToWater = grid$SqrtDistanceToWater[x]
    , Trees               = 0
    , HumansBuff5000      = grid$HumansBuff5000[x]
    , inactive            = F
  )

  # Predict scores
  pred <- predictRSS(
      model       = model
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

  # Indicate step size, human influence, SqrtDistanceToWater
  pred$step_size <- grid$step_size[x]
  pred$HumansBuff5000 <- grid$HumansBuff5000[x]
  pred$SqrtDistanceToWater <- grid$SqrtDistanceToWater[x]

  # Return the predictions
  return(pred)

}) %>% do.call(rbind, .)

# Visualize
ggplot(preds, aes(x = cos_ta_, y = RSS)) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "gray30") +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, group = Level)
    , linetype  = "solid"
    , alpha     = 0.33
    , color     = "orange"
    , fill      = "orange"
    , lwd = 0.1
  ) +
  geom_line(size = 1) +
  xlab("cos_ta_ (SD)") +
  ylab("RSS vs cos_ta_") +
  theme_classic() +
  facet_wrap( ~ step_size + HumansBuff5000)
