################################################################################
#### Plot Movement Model Results
################################################################################
# Description: A simple plot of the Dispersal Durations

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)  # For data wrangling
library(lubridate)  # To handle dates nicely
library(ggpubr)     # For nice plots
library(davidoff)   # Custom functions
library(glmmTMB)    # To handle model results
library(lemon)      # For nice axis labels
library(viridis)    # For nice colors
library(raster)     # To create rasters (visInt)
library(ggdark)     # Dark ggplot theme
library(latex2exp)  # For latex expressions

################################################################################
#### Prepare Data
################################################################################
# Old directory
# wd <- "/media/david/My Passport/Backups/WildDogs/15. PhD/00_WildDogs"
# setwd(wd)
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load the movement model
move.mod <- readRDS("03_Data/03_Results/99_MovementModel.rds")

# Select best model
best <- move.mod$Model[[1]]

# Summary of the best model
summary(best)

# Calculate Confidence intervals
coeffs <- getCoeffs(best, pvalue = TRUE)[-1, ] %>%
  mutate(
      LCI_90 = Coefficient - qnorm(1 - (1 - 0.90) / 2) * SE
    , UCI_90 = Coefficient + qnorm(1 - (1 - 0.90) / 2) * SE
    , LCI_95 = Coefficient - qnorm(1 - (1 - 0.95) / 2) * SE
    , UCI_95 = Coefficient + qnorm(1 - (1 - 0.95) / 2) * SE
    , LCI_99 = Coefficient - qnorm(1 - (1 - 0.99) / 2) * SE
    , UCI_99 = Coefficient + qnorm(1 - (1 - 0.99) / 2) * SE
  )

# Add stars indicating the significance
coeffs$Significance <- sapply(1:nrow(coeffs), function(x){
  if (coeffs$pvalue[x] <= 0.01){
    return("***")
  } else if (coeffs$pvalue[x] <= 0.05){
    return("**")
  } else if (coeffs$pvalue[x] <= 0.1){
    return("*")
  }
})

# Rename covariates
coeffs$Covariate <- gsub(coeffs$Covariate
  , pattern     = "cos_ta_"
  , replacement = "cos(ta)"
)
coeffs$Covariate <- gsub(coeffs$Covariate
  , pattern     = "log_sl_"
  , replacement = "log(sl)"
)
coeffs$Covariate <- gsub(coeffs$Covariate
  , pattern     = "sl_"
  , replacement = "sl"
)
coeffs$Covariate <- gsub(coeffs$Covariate
  , pattern     = "Shrubs"
  , replacement = "Shrubs/Grassland"
)
coeffs$Covariate <- gsub(coeffs$Covariate
  , pattern     = "Trees"
  , replacement = "Woodland"
)
coeffs$Covariate <- gsub(coeffs$Covariate
  , pattern     = "HumansBuff5000"
  , replacement = "HumanInfluence"
)
coeffs$Covariate <- gsub(coeffs$Covariate
  , pattern     = "inactiveTRUE"
  , replacement = "LowActivity"
)
coeffs$Covariate <- gsub(coeffs$Covariate
  , pattern     = "SqrtDistanceToWater"
  , replacement = "DistanceToWater '*m^0.5'"
)
coeffs$Preference <- ifelse(coeffs$Coefficient > 0, "Preferred", "Avoided")
coeffs$Preference <- factor(coeffs$Preference, levels = c("Preferred", "Avoided"))

# Specify the order in which the coefficients should be plotted
order <- c(
      "Water"
    , "DistanceToWater '*m^0.5'"
    , "Woodland"
    , "Shrubs/Grassland"
    , "HumanInfluence"
    , "sl"
    , "cos(ta)"
    , "log(sl)"
    , "cos(ta):sl"
    , "cos(ta):log(sl)"
    , "sl:LowActivity"
    , "sl:Water"
    , "sl:Woodland"
    , "sl:Shrubs/Grassland"
    , "sl:DistanceToWater '*m^0.5'"
    , "cos(ta):HumanInfluence"
    , "cos(ta):DistanceToWater '*m^0.5'"
)

# Prepare groups with "expressions"
groups <- c(
      "Water"
    , expression(DistanceToWater^0.5)
    , "Woodland"
    , "Shrubs/Grassland"
    , "HumanInfluence"
    , "sl"
    , "cos(ta)"
    , "log(sl)"
    , "cos(ta):sl"
    , "cos(ta):log(sl)"
    , "sl:LowActivity"
    , "sl:Water"
    , "sl:Woodland"
    , "sl:Shrubs/Grassland"
    , expression(sl:DistanceToWater^0.5)
    , "cos(ta):HumanInfluence"
    , expression(cos(ta):DistanceToWater^0.5)
)

# Prepare plot with Covariates on the y-axis and the corresponding
# coefficients on the x-axis
ggplot(data = coeffs, aes(y = Covariate, x = Coefficient, col = factor(Preference))) +
  geom_point(shape = 3, size = 2.5) +
  geom_errorbarh(aes(
      xmin = Coefficient - 1.645 * SE
    , xmax = Coefficient + 1.645 * SE)
    , height = 0, size = 1.5, alpha = 0.5
  ) +
  geom_errorbarh(aes(
      xmin = Coefficient - 1.96 * SE
    , xmax = Coefficient + 1.96 * SE)
    , height = 0, size = 0.8, alpha = 0.75
  ) +
  geom_errorbarh(aes(
      xmin = Coefficient - 2.575 * SE
    , xmax = Coefficient + 2.575 * SE)
    , height = 0, size = 0.2, alpha = 1
  ) +
  geom_text(aes(label = Significance, hjust = 0.5, vjust = -0.2), show.legend = F) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
  scale_y_discrete(labels = rev(groups), limits = rev(order)) +
  theme_classic() +
  xlim(c(-1.3, 1.3)) +
  coord_capped_cart(left = "both", bottom = "both") +
  labs(x = expression(beta*"-Coefficient")) +
  scale_color_manual(values = c("#5B9BD5", "orange")) +
  # dark_theme_classic() +
  theme(
    #   panel.grid.minor = element_line(size = 0.0)
    # , panel.grid.major = element_line(size = 0.0)
    # , panel.border     = element_blank()
    # , axis.line        = element_line()
    # , axis.ticks       = element_line(colour = "black")
    , legend.title     = element_blank()
    # , axis.text.y      = element_text(color = rev(labcols))
  )

# Store the plot
ggsave("test.png", device = "png", width = 10, height = 5, scale = 0.65)

################################################################################
#### Validation
################################################################################
# Load required data
validation <- read_rds("03_Data/03_Results/99_ModelValidation.rds")
dat_pref <- read_rds("03_Data/03_Results/99_ModelValidation(Data).rds")[[1]]
dat_rand <- read_rds("03_Data/03_Results/99_ModelValidation(Data).rds")[[2]]

# Put the k-fold cross validation data from observed and random preferences
# together
dat_pref$Group <- "Realized"
dat_rand$Group <- "Random"
dat <- rbind(dat_pref, dat_rand)

# We want to plot this data and add the information from the validation table
# too. Let's prepare a column that indicates the text that we want to plot on
# top of the data. Let's first round the values from the validation table
validation[, c(2:4)] <- round(validation[, c(2:4)], 2)

# Match the Grouping names of the validation table to the groups above
validation$Group <- as.factor(c("Realized", "Random"))

# Create a dataframe which we use to annotate the two facets
text <- data.frame(
    Group = c("Realized", "Random")
  , Text = c(
        paste0(
            "$\\bar{r}_s = "
          , validation$Mean[1]
          , "$, $95%-CI = "
          , validation$LCL[1]
          , "$, $"
          , validation$UCL[1], "$"
        )
      , paste0(
            "$\\bar{r}_s = "
          , validation$Mean[2]
          , "$, $95%-CI = "
          , validation$LCL[2]
          , "$, $"
          , validation$UCL[2], "$"
        )
    #   "$\\bar{r}_s = -0.45$, $95%-CI = -0.48$, $-0.42$"
    # , "$\\bar{r}_s = 0.04$, $95%-CI = 0.01$, $0.08$"
  )
)

# Reorder the factors in the Group variable
dat$Group <- factor(dat$Group, levels = c("Realized", "Random"))

library(msir)
loess <- dat %>%
  group_by(Group) %>%
  nest() %>%
  mutate(data = lapply(data, function(x){
    l <- loess.sd(x$Frequency ~ x$Rank)
    df <- data.frame(
        Loess = l$y
      , Rank  = l$x
      , Upper = l$upper
      , Lower = l$lower
    )
    return(df)
  })) %>%
  unnest(data)

# Plot the data
ggplot(loess, aes(x = Rank, y = Loess)) +
  geom_jitter(aes(x = Rank, y = Frequency)
    , data  = dat,
    , alpha = 0.15
    , size  = 1
  ) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper)
    , alpha    = 0.4
    , fill     = "darkorchid4"
    , col      = "darkorchid4"
    , linetype = "dashed"
  ) +
  geom_line(size = 1, col = "darkorchid4") +
  facet_wrap("Group") +
  theme_classic() +
  coord_capped_cart(left = "both", bottom = "both") +
  geom_text(data = text
    , mapping = aes(
        x = -Inf
      , y = -Inf
      , label = TeX(Text, output = "character")
    )
    , hjust   = -0.05
    , vjust   = -0.5
    , parse   = TRUE
    , size    = 3
  ) +
  ylab("Frequency")

# Plot the data
ggplot(data = dat, aes(x = Rank, y = Frequency)) +
  geom_jitter(alpha = 0.15, size = 1) +
  geom_smooth(method = "loess", se = T, level = 0.95, fill = "blue") +
  theme_classic() +
  facet_wrap("Group") +
  coord_capped_cart(left = "both", bottom = "both") +
  geom_text(data = text
    , mapping = aes(
        x = -Inf
      , y = -Inf
      , label = TeX(Text, output = "character")
    )
    , hjust   = -0.05
    , vjust   = -0.5
    , parse   = TRUE
    , size    = 3
  )

################################################################################
#### Interactions
################################################################################
# Create color ramp
pal <- colorRampPalette(magma(20))

# Check possible interactions
summary(best)

# Visualize them
visInt(best, xVar = "log_sl_", yVar = "Water", colorPalette = pal)
visInt(best, xVar = "log_sl_", yVar = "Woodland", colorPalette = pal)
visInt(best, xVar = "cos_ta_", yVar = "DistanceToWater", colorPalette = pal)
visInt(best, xVar = "cos_ta_", yVar = "HumansBuff5000", colorPalette = pal)

# Prepare visualizations using raster
int1 <- visInt2(best, yVar = "log_sl_", xVar = "Water")
int2 <- visInt2(best, yVar = "log_sl_", xVar = "Woodland")
int3 <- visInt2(best, yVar = "cos_ta_", xVar = "DistanceToWater")
int4 <- visInt2(best, yVar = "cos_ta_", xVar = "HumansBuff5000")

# Visualize using raster
par(mfrow = c(2, 2))
plot(int1, col = pal(20), ylab = "log_sl_", xlab = "Water")
plot(int2, col = pal(20), ylab = "log_sl_", xlab = "Woodland")
plot(int3, col = pal(20), ylab = "cos_ta_", xlab = "DistanceToWater")
plot(int4, col = pal(20), ylab = "cos_ta_", xlab = "HumansBuff5000")

################################################################################
#### Random Effects
################################################################################
# Rename model
mod <- best

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
  gather(key = Covariate, value = Mean, 2:8)
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
  , "Woodland", "HumanInfluence"
))

# Visualize. Note that I am transforming the variance using mean - 2 *
# sqrt(Variance). This was taken from here: https://stackoverflow.com/questions
# /13847936/plot-random-effects-from-lmer-lme4-package-using-qqmath-or-dotplot-
# how-to-mak
p <- ggplot(rfs, aes(x = Mean, y = id)) +
  geom_point() +
  facet_wrap("Covariate", nrow = 2) +
  geom_errorbarh(aes(
      xmin = Mean - 2 * sqrt(Variance)
    , xmax = Mean + 2 * sqrt(Variance)
  ), colour = "black", height = 0)
