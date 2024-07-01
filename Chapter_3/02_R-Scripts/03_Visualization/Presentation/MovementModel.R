################################################################################
#### Movement Model Visualization
################################################################################
# Description: Results from the seasonal step selection models

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Load required packages
library(tidyverse)    # To wrangle data
library(kableExtra)   # To generate a nice table
library(lemon)        # For capped coordinate axes
library(colorspace)   # To darken colors
library(latex2exp)    # To use latex expressions
library(ggh4x)        # For nested plots
library(ggdark)       # For dark theme

# Load custom functions
source("02_R-Scripts/00_Functions.R")

################################################################################
#### Table
################################################################################
# Load results and do some cleaning
dat <- "/home/david/ownCloud/University/15. PhD/Chapter_3/03_Data/03_Results/MovementModels.rds" %>%
  read_rds() %>%
  # unnest(Models) %>%
  subset(Covariate != "(Intercept)") %>%
  rename(
      `z-value`  = zvalue
    , `p-value`  = pvalue
    , `Variance` = RandomVariance
    , `SD`       = RandomSD
  ) %>%
  subset(!is.na(SE)) %>%
  mutate(FittingCovariates = factor(FittingCovariates, levels = c("Static", "Dynamic"))) %>%
  mutate(Formula = factor(Formula, levels = c("Simple", "Full"), labels = c("Simple Formula", "Complex Formula")))

# Add stars indicating the significance
dat$Significance <- sapply(1:nrow(dat), function(x){
  if (dat$`p-value`[x] <= 0.01){
    return("***")
  } else if (dat$`p-value`[x] <= 0.05) {
    return("**")
  } else if (dat$`p-value`[x] <= 0.1) {
    return("*")
  } else {
    return("")
  }
})

# Reorder the columns
dat <- select(dat, FittingCovariates, Formula, Season, NumberRandomSteps, Covariate, Coefficient, SE, `z-value`, `p-value`, Significance, Variance, SD, everything())

# Calculate confidence intervals
coeffs <- dat %>%
  mutate(
      LCI_90 = Coefficient - qnorm(1 - (1 - 0.90) / 2) * SE
    , UCI_90 = Coefficient + qnorm(1 - (1 - 0.90) / 2) * SE
    , LCI_95 = Coefficient - qnorm(1 - (1 - 0.95) / 2) * SE
    , UCI_95 = Coefficient + qnorm(1 - (1 - 0.95) / 2) * SE
    , LCI_99 = Coefficient - qnorm(1 - (1 - 0.99) / 2) * SE
    , UCI_99 = Coefficient + qnorm(1 - (1 - 0.99) / 2) * SE
  )

# Create a tidy dataframe for the confidence intervals
coeffs <- coeffs %>%
  dplyr::select(FittingCovariates, Formula, Season, NumberRandomSteps, Covariate, Coefficient, SE, `z-value`, `p-value`, Significance, Variance, SD, LCI_90, UCI_90, LCI_95, UCI_95, LCI_99, UCI_99) %>%
  pivot_longer(!FittingCovariates & !Formula & !Season & !Covariate & !NumberRandomSteps & !Coefficient & !SE & !`z-value` & !`p-value` & !Significance & !Variance & !SD) %>%
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
  mutate(Level = as.factor(Level)) %>%
  mutate(Kernel = factor(case_when(
      Covariate %in% c(
          "Water"
        , "SqrtDistanceToWater"
        , "Shrubs"
        , "Trees"
        , "Humans"
      ) ~ "Habitat"
    , Covariate %in% c("sl"
      , "log_sl"
      , "cos_ta"
      , "cos_ta:sl"
      , "cos_ta:log_sl"
      , "sl:LightTypeDark"
      , "log_sl:LightTypeDark"
      , "sl:Temperature"
      , "log_sl:Temperature"
    ) ~ "Movement"
    , .default = "Interaction"
  ), levels = c("Movement", "Habitat", "Interaction"), labels = c("Movement Kernel", "Habitat-Selection Function", "Interaction"))) %>%
  mutate(Covariate = factor(Covariate
    , levels = c(
        "sl"
      , "log_sl"
      , "cos_ta"
      , "cos_ta:sl"
      , "cos_ta:log_sl"
      , "sl:LightTypeDark"
      , "log_sl:LightTypeDark"
      , "sl:Temperature"
      , "log_sl:Temperature"
      , "Humans"
      , "Trees"
      , "Shrubs"
      , "Water"
      , "SqrtDistanceToWater"
      , "SqrtDistanceToPans"
      , "sl:Trees"
      , "sl:Shrubs"
      , "sl:Water"
      , "cos_ta:SqrtDistanceToWater"
      , "cos_ta:Humans"
    )
    , labels = c(
        "sl"
      , "log(sl)"
      , "cos(ta)"
      , "cos(ta):sl"
      , "cos(ta):log(sl)"
      , "sl:Dark"
      , "log(sl):Dark"
      , "sl:\\underline{Temperature}"
      , "log(sl):\\underline{Temperature}"
      , "Humans"
      , "\\underline{Trees}"
      , "\\underline{Shrubs}"
      , "\\underline{Water}"
      , "\\underline{DistanceToWater}$^{0.5}$"
      , "\\underline{DistanceToPans}$^{0.5}$"
      , "sl:\\underline{Trees}"
      , "sl:\\underline{Shrubs}"
      , "sl:\\underline{Water}"
      , "cos(ta):\\underline{DistanceToWater}$^{0.5}$"
      , "cos(ta):Humans"
    )
  ))

################################################################################
#### Figure of Model Results
################################################################################
# Prepare the plot
dodge <- 0.6
p1 <- ggplot(data = subset(coeffs, NumberRandomSteps == 25), aes(x = Covariate, y = Coefficient, col = factor(Season), group = factor(Season))) +
  geom_hline(
      yintercept = 0
    , color      = "darkgrey"
    , lty        = 1
    , lwd        = 0.3
  ) +
  geom_errorbar(
      mapping  = aes(ymin = LCI, ymax = UCI, linewidth = factor(Level))
    , width    = 0
    , alpha    = 0.5
    , position = position_dodge(width = dodge)
  ) +
  geom_point(
      shape    = 3
    , size     = 1
    , position = position_dodge(width = dodge)
  ) +
  # geom_text(
  #     mapping     = aes(label = Significance)
  #   , vjust       = -0.25
  #   , show.legend = F
  #   , position    = position_dodge(width = dodge)
  #   , angle       = 90
  # ) +
  facet_nested(Formula + FittingCovariates ~ Kernel, scales = "free_x", space = "free_x") +
  dark_theme_awesome() +
  labs(y = expression(beta*"-Coefficient")) +
  scale_x_discrete(labels = TeX) +
  scale_color_viridis_d(begin = 0.3, name = "Season", direction = -1) +
  scale_linewidth_manual(
      name   = "Confidence Level"
    , values = c(1.5, 0.75, 0.3)
  ) +
  theme(
    , legend.position   = "bottom"
    , legend.margin     = margin(0, 50, 0, -20)
    , legend.box.margin = margin(-5, -10, -5, -10)
    , legend.text       = element_text(face  = 3)
    , legend.title      = element_text(face  = 3)
    , panel.grid.minor  = element_blank()
    , strip.background  = element_rect(fill = "gray30", color = NA)
    , plot.caption      = element_text(face = 3, color = "gray30")
    , panel.grid.major  = element_line(color = adjustcolor("white", alpha.f = 0.25))
    , axis.text.x       = element_text(angle = 45, hjust = 1)
    , axis.text.y       = element_text(size = 8)
    , axis.title.y      = element_text(angle = 90)
    , plot.background   = element_blank()
    , panel.background  = element_blank()
  ) +
  guides(
      colour    = guide_legend(title.position = "top", title.hjust = 0.5)
    , linewidth = guide_legend(title.position = "top", title.hjust = 0.5, override.aes = list(colour = "#3CBB75FF"))
  )

# Store the plot to file
ggsave("05_Presentation/MovementModel.png"
  , plot   = p1
  , device = png
  , bg     = "transparent"
  , width  = 7
  , height = 4
  , scale  = 1.3
)

################################################################################
#### Stability of Estimates with Regards to Number of Random Steps
################################################################################
# Plot showing how the estimates change for different number of random steps
dodge <- 0.8
p2 <- ggplot(data = coeffs, aes(x = Covariate, y = Coefficient, col = factor(NumberRandomSteps), group = factor(NumberRandomSteps))) +
  geom_hline(
      yintercept = 0
    , color      = "darkgrey"
    , lty        = 1
    , lwd        = 0.3
  ) +
  geom_errorbar(
      mapping  = aes(ymin = LCI, ymax = UCI, linewidth = factor(Level))
    , width    = 0
    , alpha    = 0.5
    , position = position_dodge(width = dodge)
  ) +
  geom_point(
      shape    = 3
    , size     = 1
    , position = position_dodge(width = dodge)
  ) +
  # geom_text(
  #     mapping     = aes(label = Significance)
  #   , vjust       = -0.25
  #   , show.legend = F
  #   , position    = position_dodge(width = dodge)
  #   , angle       = 90
  # ) +
  facet_nested(Formula + FittingCovariates + factor(Season) ~ Kernel, scales = "free_x", space = "free_x") +
  dark_theme_awesome() +
  labs(y = expression(beta*"-Coefficient")) +
  scale_x_discrete(labels = TeX) +
  scale_color_viridis_d(begin = 0.3, name = "Number of Random Steps", direction = -1) +
  scale_linewidth_manual(
      name   = "Confidence Level"
    , values = c(1.5, 0.75, 0.3)
  ) +
  theme(
    , legend.position   = "bottom"
    , legend.margin     = margin(0, 50, 0, -20)
    , legend.box.margin = margin(-5, -10, -5, -10)
    , legend.text       = element_text(face  = 3)
    , legend.title      = element_text(face  = 3)
    , panel.grid.minor  = element_blank()
    , strip.background  = element_rect(fill = "gray30", color = NA)
    , plot.caption      = element_text(face = 3, color = "gray30")
    , panel.grid.major  = element_line(color = adjustcolor("white", alpha.f = 0.25))
    , axis.text.x       = element_text(angle = 45, hjust = 1)
    , axis.text.y       = element_text(size = 8)
    , axis.title.y      = element_text(angle = 90)
    , plot.background   = element_blank()
    , panel.background  = element_blank()
  ) +
  guides(
      colour    = guide_legend(title.position = "top", title.hjust = 0.5)
    , linewidth = guide_legend(title.position = "top", title.hjust = 0.5, override.aes = list(colour = "#3CBB75FF"))
  )

# Store the plot to file
ggsave("05_Presentation/Stability.png"
  , plot   = p2
  , device = png
  , bg     = "transparent"
  , width  = 6
  , height = 5
  , scale  = 1.6
)
