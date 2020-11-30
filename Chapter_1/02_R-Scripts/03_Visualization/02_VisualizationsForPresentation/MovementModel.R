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
library(cowplot)    # To add image to ggplot
library(lemon)      # For nice axis labels
library(imager)     # To import images
library(viridis)    # For nice colors
library(raster)     # To create rasters (visInt)
library(ggdark)     # Dark ggplot theme

################################################################################
#### Prepare Data
################################################################################
# Old directory
wd <- "/media/david/My Passport/Backups/WildDogs/15. PhD/00_WildDogs"
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
      LCI = Coefficient - 1.96 * SE
    , UCI = Coefficient + 1.96 * SE
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
coeffs$Covariate[coeffs$Covariate == "Shrubs"] <- "Shrubs/Grassland"
coeffs$Covariate[coeffs$Covariate == "Trees"] <- "Woodland"
coeffs$Covariate[coeffs$Covariate == "Water"] <- "Water"
coeffs$Covariate[coeffs$Covariate == "DistanceToWater"] <- "DistanceToWater"
coeffs$Covariate[coeffs$Covariate == "cos_ta_"] <- "cos(ta)"
coeffs$Covariate[coeffs$Covariate == "log_sl_"] <- "log(sl)"
coeffs$Covariate[coeffs$Covariate == "HumansBuff5000"] <- "HumanInfluence"
coeffs$Covariate[coeffs$Covariate == "log_sl_:ActivityMainActivity"] <- "log(sl):MainActivity"
coeffs$Covariate[coeffs$Covariate == "log_sl_:Water"] <- "log(sl):Water"
coeffs$Covariate[coeffs$Covariate == "log_sl_:Trees"] <- "log(sl):Woodland"
coeffs$Covariate[coeffs$Covariate == "cos_ta_:DistanceToWater"] <- "cos(ta):DistanceToWater"
coeffs$Covariate[coeffs$Covariate == "cos_ta_:HumansBuff5000"] <- "cos(ta):HumanInfluence"
coeffs$Preference <- ifelse(coeffs$Coefficient > 0, "Preferred", "Avoided")
coeffs$Preference <- factor(coeffs$Preference, levels = c("Preferred", "Avoided"))

# Specify the order in which the coefficients should be plotted
order <- c(
      "Water"
    , "DistanceToWater"
    , "Woodland"
    , "Shrubs/Grassland"
    , "HumanInfluence"
    , "cos(ta):HumanInfluence"
    , "cos(ta):DistanceToWater"
    # , "log(sl):MainActivity"
    , "log(sl):Water"
    , "log(sl):Woodland"
    , "log(sl)"
    , "cos(ta)"
)

# Specify colors of axis labels
labcols <- c("black", "orange")[c(2, 2, 2, 1, 2, 1, 2, 1, 1, 1, 2, 2)]

# Let's ignore "main activity" for the moment
coeffs <- subset(coeffs, Covariate != "log(sl):MainActivity")

# Prepare plot with Covariates on the y-axis and the corresponding
# coefficients on the x-axis
p1 <- ggplot(data = coeffs, aes(y = Covariate, x = Coefficient, col = factor(Preference))) +
  geom_point(shape = 20, size = 2) +
  geom_errorbarh(aes(
      xmin = Coefficient - 1.96 * SE
    , xmax = Coefficient + 1.96 * SE
    , height = 0)
  ) +
  geom_text(aes(label = Significance, hjust = 0.5, vjust = -0.2), show.legend = F) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
  scale_y_discrete(limits = rev(order)) +
  # theme_cowplot() +
  xlim(c(-1.2, 1.2)) +
  coord_capped_cart(left = "both", bottom = "both") +
  labs(x = expression(beta*"-Coefficient")) +
  # scale_color_manual(values = c("#548CBE", "#FFA500")) +
  scale_color_manual(values = c("#5B9BD5", "orange")) +
  dark_theme_classic() +
  theme(
    #   panel.grid.minor = element_line(size = 0.0)
    # , panel.grid.major = element_line(size = 0.0)
    # , panel.border     = element_blank()
    # , axis.line        = element_line()
    # , axis.ticks       = element_line(colour = "black")
    , legend.title     = element_blank()
    # , axis.text.y      = element_text(color = rev(labcols))
  )

# Load a picture of a wild dog to the plot
# dog <- load.image("/home/david/ownCloud/University/15. PhD/General/Images/WildDog_Running.svg")

# Add it to the previous plot
# p <- ggdraw() +
#   draw_image(dog, x = 0.9, y = 0.7, hjust = 0.5, vjust = 0.5, scale = 0.2) +
#   draw_plot(p1)

# Store the plot
# CairoPDF("test.pdf", width = 8, height = 5, bg = "transparent")
ggsave("test.png", device = "png", width = 10, height = 5, scale = 0.65)

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
