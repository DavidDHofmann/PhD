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

dat %>%
  group_by(Covariate, FittingCovariates) %>%
  summarize(Coefficient = mean(Coefficient), .groups = "drop") %>%
  mutate(Coefficient = round(Coefficient, 1))

dat %>%
  subset(Covariate == "SqrtDistanceToWater") %>%
  group_by(FittingCovariates) %>%
  summarize(Coefficient = mean(Coefficient), .groups = "drop") %>%
  mutate(Coefficient = round(Coefficient, 1))

dat %>%
  subset(Covariate == "Water") %>%
  group_by(FittingCovariates, Season) %>%
  summarize(Coefficient = mean(Coefficient), .groups = "drop") %>%
  mutate(Coefficient = round(Coefficient, 1))

dat %>%
  subset(Covariate == "sl:LightTypeDark") %>%
  group_by(FittingCovariates, Season) %>%
  summarize(Coefficient = mean(Coefficient), .groups = "drop") %>%
  mutate(Coefficient = round(Coefficient, 1))

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

# Make four nice tables (simple + static, simple + dynamic, complex + static, complex + dynamic)
tab <- coeffs %>%
  subset(NumberRandomSteps == 25) %>%
  arrange(FittingCovariates, Formula, Season, Kernel, Covariate) %>%
  select(FittingCovariates, Formula, Season, Covariate, Coefficient, SE, `z-value`, `p-value`, Significance, Variance, SD) %>%
  distinct() %>%
  rename(Fitting = FittingCovariates) %>%
  mutate(
      Variance = ifelse(is.na(Variance), " - ", formatC(Variance, format = "f", flag = "0", digits = 3))
    , SD       = ifelse(is.na(SD), " - ", formatC(SD, format = "f", flag = "0", digits = 3))
  )
tab %>%
  subset(Fitting == "Static" & Formula == "Simple Formula") %>%
  select(-c(Fitting, Formula, `z-value`, SE, Variance)) %>%
  kbl(booktabs = T, format = "latex", escape = F, digits = 3, align = "l") %>%
    add_header_above(c("", "", "Fixed Effects" = 3, "Random Effects" = 1), align = "l") %>%
    collapse_rows(1, latex_hline = c("custom"), custom_latex_hline = 1) %>%
    writeLines("04_Manuscript/Figures/MovementModelStaticSimple.tex")
tab %>%
  subset(Fitting == "Static" & Formula == "Complex Formula") %>%
  select(-c(Fitting, Formula, `z-value`, SE, Variance)) %>%
  kbl(booktabs = T, format = "latex", escape = F, digits = 3, align = "l") %>%
    add_header_above(c("", "", "Fixed Effects" = 3, "Random Effects" = 1), align = "l") %>%
    collapse_rows(1, latex_hline = c("custom"), custom_latex_hline = 1) %>%
    writeLines("04_Manuscript/Figures/MovementModelStaticFull.tex")
tab %>%
  subset(Fitting == "Dynamic" & Formula == "Simple Formula") %>%
  select(-c(Fitting, Formula, `z-value`, SE, Variance)) %>%
  kbl(booktabs = T, format = "latex", escape = F, digits = 3, align = "l") %>%
    add_header_above(c("", "", "Fixed Effects" = 3, "Random Effects" = 1), align = "l") %>%
    collapse_rows(1, latex_hline = c("custom"), custom_latex_hline = 1) %>%
    writeLines("04_Manuscript/Figures/MovementModelDynamicSimple.tex")
tab %>%
  subset(Fitting == "Dynamic" & Formula == "Complex Formula") %>%
  select(-c(Fitting, Formula, `z-value`, SE, Variance)) %>%
  kbl(booktabs = T, format = "latex", escape = F, digits = 3, align = "l") %>%
    add_header_above(c("", "", "Fixed Effects" = 3, "Random Effects" = 1), align = "l") %>%
    collapse_rows(1, latex_hline = c("custom"), custom_latex_hline = 1) %>%
    writeLines("04_Manuscript/Figures/MovementModelDynamicFull.tex")

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
  theme_awesome() +
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
    , axis.text.x       = element_text(angle = 45, hjust = 1)
    , axis.text.y       = element_text(size = 8)
    , axis.title.y      = element_text(angle = 90)
  ) +
  guides(
      colour    = guide_legend(title.position = "top", title.hjust = 0.5)
    , linewidth = guide_legend(title.position = "top", title.hjust = 0.5, override.aes = list(colour = "#3CBB75FF"))
  )

# Store the plot to file
ggsave("04_Manuscript/Figures/MovementModel.png"
  , plot   = p1
  , device = png
  , bg     = "white"
  , width  = 5
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
  theme_awesome() +
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
    , panel.spacing.x   = unit(0, "lines")
    , axis.text.x       = element_text(angle = 45, hjust = 1)
    , axis.text.y       = element_text(size = 6)
    , axis.title.y      = element_text(angle = 90)
  ) +
  guides(
      colour    = guide_legend(title.position = "top", title.hjust = 0.5)
    , linewidth = guide_legend(title.position = "top", title.hjust = 0.5, override.aes = list(colour = "#3CBB75FF"))
  )

# Store the plot to file
ggsave("04_Manuscript/Figures/Stability.png"
  , plot   = p2
  , device = png
  , bg     = "white"
  , width  = 6
  , height = 5
  , scale  = 1.6
)
