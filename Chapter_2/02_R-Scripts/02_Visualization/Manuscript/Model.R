################################################################################
#### Plot Movement Model Results
################################################################################
# Description: Plot of the movement model from Hofmann et al. 2023

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Load required packages
library(tidyverse)    # For data wrangling
library(lubridate)    # To handle dates nicely
library(ggpubr)       # For nice plots
library(glmmTMB)      # To handle model results
library(lemon)        # For nice axis labels
library(latex2exp)    # For latex expressions
library(colorspace)   # To darken or lighten colors

# Suppress the scientific notion
options(scipen = 999)

################################################################################
#### Prepare Data
################################################################################
# Load the movement model
best <- readRDS("03_Data/02_CleanData/MovementModel.rds")

# Extract model results (note that newer versions of glmmTMB and TMB are not
# compatible with older models and we need to manually extract the values of
# interest)
best <- tibble(
    Covariate   = names(fixef(best)$cond)
  , Coefficient = fixef(best)$cond
  , SE          = sqrt(diag(vcov(best, full = T)))[1:length(Coefficient)]
)

# Calculate confidence intervals
coeffs <- best[-1, ] %>%
  mutate(
      LCI_90 = Coefficient - qnorm(1 - (1 - 0.90) / 2) * SE
    , UCI_90 = Coefficient + qnorm(1 - (1 - 0.90) / 2) * SE
    , LCI_95 = Coefficient - qnorm(1 - (1 - 0.95) / 2) * SE
    , UCI_95 = Coefficient + qnorm(1 - (1 - 0.95) / 2) * SE
    , LCI_99 = Coefficient - qnorm(1 - (1 - 0.99) / 2) * SE
    , UCI_99 = Coefficient + qnorm(1 - (1 - 0.99) / 2) * SE
  )

# Compute p-values
coeffs <- mutate(coeffs, pvalue = 2 * pnorm(-abs(Coefficient) / SE))

# Add stars indicating the significance
coeffs$Significance <- sapply(1:nrow(coeffs), function(x){
  if (coeffs$pvalue[x] <= 0.01){
    return("***")
  } else if (coeffs$pvalue[x] <= 0.05) {
    return("**")
  } else if (coeffs$pvalue[x] <= 0.1) {
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
  , replacement = "ln(sl)"
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
coeffs$Preference <- factor(coeffs$Preference, levels = c("Avoided", "Preferred"))

# Create a tidy dataframe for the confidence intervals
confs <- coeffs %>%
  dplyr::select(Covariate, Coefficient, LCI_90, UCI_90, LCI_95, UCI_95, LCI_99, UCI_99) %>%
  pivot_longer(!Covariate & !Coefficient) %>%
  separate(name, into = c("Limit", "Level"), sep = "_", convert = T) %>%
  spread(key = Limit, value = value) %>%
  mutate(alpha = case_when(
      Level == 90 ~ 0.50
    , Level == 95 ~ 0.75
    , Level == 99 ~ 1.00
  )) %>%
  mutate(size = case_when(
      Level == 90 ~ 2.0
    , Level == 95 ~ 1.0
    , Level == 99 ~ 0.3
  )) %>%
  mutate(Level = as.factor(Level))

# Join the information on preference
confs <- coeffs %>%
  dplyr::select(Covariate, Preference) %>%
  left_join(confs, ., by = "Covariate")

# Specify the order in which the coefficients should be plotted
order <- c(
      "Water"
    , "DistanceToWater '*m^0.5'"
    , "Woodland"
    , "Shrubs/Grassland"
    , "HumanInfluence"
    , "sl"
    , "cos(ta)"
    , "ln(sl)"
    , "cos(ta):sl"
    , "cos(ta):ln(sl)"
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
    , "ln(sl)"
    , "cos(ta):sl"
    , "cos(ta):ln(sl)"
    , "sl:LowActivity"
    , "sl:Water"
    , "sl:Woodland"
    , "sl:Shrubs/Grassland"
    , expression(sl:DistanceToWater^0.5)
    , "cos(ta):HumanInfluence"
    , expression(cos(ta):DistanceToWater^0.5)
)

# Prepare dataset for plotting confidence intervals
coeffs2 <- coeffs %>%
  dplyr::select(Covariate, Coefficient, Preference, LCI_90:UCI_99) %>%
  gather(key = confidence_level, value = value, LCI_90:UCI_99) %>%
  separate(col = confidence_level, into = c("Type", "Level"), sep = "_") %>%
  spread(key = Type, value = value) %>%
  mutate(Level = paste0(Level, "%"))

# Prepare plot with Covariates on the y-axis and the corresponding
# coefficients on the x-axis
p1 <- ggplot(data = coeffs, aes(y = Covariate, x = Coefficient, col = factor(Preference))) +
  geom_point(shape = 3, size = 2.5) +
  geom_errorbarh(
      aes(
        xmin = LCI
      , xmax = UCI
      , size = factor(Level)
    )
    , data = coeffs2
    , height = 0
    , alpha  = 0.5
  ) +
  geom_text(
      aes(label = Significance, hjust = 0.5, vjust = 0)
    , show.legend = F
  ) +
  geom_vline(
      xintercept = 0
    , color      = "darkgrey"
    , lty        = 2
    , lwd        = 0.3
  ) +
  scale_y_discrete(
      labels = rev(groups)
    , limits = rev(order)
  ) +
  theme_classic() +
  xlim(c(-1.3, 0.6)) +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
  ) +
  labs(x = expression(beta*"-Coefficient")) +
  scale_color_manual(
      name   = "Preference"
    , values = c("#5B9BD5", "orange")
  ) +
  scale_size_manual(
      name   = "Confidence Level"
    , values = c(2, 1, 0.3)
  ) +
  theme(
    , legend.position   = "bottom"
    , legend.margin     = margin(0, 50, 0, -20)
    , legend.box.margin = margin(-5, -10, -5, -10)
    , legend.text       = element_text(face = 3)
    , legend.title      = element_text(face = 3)
  ) +
  guides(
      colour = guide_legend(title.position = "top", title.hjust = 0.5)
    , size   = guide_legend(title.position = "top", title.hjust = 0.5, override.aes = list(colour = "#3CBB75FF"))
  )

# # Legacy Plot
# p1 <- ggplot(data = coeffs, aes(y = Covariate, x = Coefficient, col = factor(Preference))) +
#   geom_point(shape = 3, size = 2.5) +
#   geom_errorbarh(aes(
#       xmin = LCI_90
#     , xmax = UCI_90)
#     , height = 0, size = 2, alpha = 0.5
#   ) +
#   geom_errorbarh(aes(
#       xmin = LCI_95
#     , xmax = UCI_95)
#     , height = 0, size = 1, alpha = 0.75
#   ) +
#   geom_errorbarh(aes(
#       xmin = LCI_99
#     , xmax = UCI_99)
#     , height = 0, size = 0.3, alpha = 1
#   ) +
#   geom_text(
#       aes(label = Significance, hjust = 0.5, vjust = 0)
#     , show.legend = F
#   ) +
#   geom_vline(
#       xintercept = 0
# Add annotations
cols <- hcl.colors(n = 9, palette = "BuPu")
p2 <- p1 + annotate("rect"
    , xmin  = -1.30
    , xmax  = -1.20
    , ymin  = 12.75
    , ymax  = 17
    , alpha = 0.8
    , fill  = cols[6]
  ) + annotate("rect"
    , xmin  = -1.30
    , xmax  = -1.20
    , ymin  = 6.75
    , ymax  = 12.25
    , alpha = 0.8
    , fill  = cols[5]
  ) + annotate("rect"
    , xmin  = -1.30
    , xmax  = -1.20
    , ymin  = 1
    , ymax  = 6.25
    , alpha = 0.8
    , fill  = cols[4]
  ) + annotate(geom = "text"
    , x        = -1.25
    , y        = 15
    , label    = "HABITAT KERNEL"
    , color    = darken(cols[6])
    , angle    = 90
    , fontface = 3
    , size     = 2.5
  ) + annotate(geom = "text"
    , x        = -1.25
    , y        = 9.5
    , label    = "MOVEMENT KERNEL"
    , color    = darken(cols[5])
    , angle    = 90
    , fontface = 3
    , size     = 2.5
  ) + annotate(geom = "text"
    , x        = -1.25
    , y        = 3.5
    , label    = "INTERACTIONS"
    , color    = darken(cols[4])
    , angle    = 90
    , fontface = 3
    , size     = 2.5
  )

# Add lines to separate kernels
p3 <- p2 + annotate(geom = "segment"
    , x      = -1.30
    , xend   = 0.6
    , y      = 6.5
    , yend   = 6.5
    , colour = "gray80"
    , lty    = 1
    , lwd    = 0.3
  ) + annotate(geom = "segment"
    , x      = -1.30
    , xend   = 0.6
    , y      = 12.5
    , yend   = 12.5
    , colour = "gray80"
    , lty    = 1
    , lwd    = 0.3
  )

# Show the plot
p3

# Store the plot
ggsave("04_Manuscript/Figures/MovementModel.png"
  , width  = 7
  , height = 7
  , scale  = 0.75
  , device = png
)
