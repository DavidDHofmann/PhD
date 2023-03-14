################################################################################
#### Plot of Model of Average Activity
################################################################################
# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_7"
setwd(wd)

# Load required packages
library(tidyverse)   # For data wrangling
library(broom)       # To clean model outputs
library(lemon)       # For capped axes
library(latex2exp)   # For latex expressions
library(grid)        # To add png to plot
library(egg)         # To add png to plot
library(png)         # To add png to plot

# Write this to file
coefs <- read_rds("03_Data/03_Results/99_StartActivityModelCoefficients.rds")

# Compute weights for each beta and run weighted regression to compute
# std.errors and confidence intervals
coefs_summary <- coefs %>%
  nest(Coefs = -c(ToD, term)) %>%
  mutate(Coefs = map(Coefs, function(x) {
    x$weight <- (1 / (x$std.error ** 2)) / sum(1 / (x$std.error ** 2))
    mod <- lm(estimate ~ 1, weights = weight, data = x)
    mod_summary <- tidy(mod)
    mod_summary <- mod_summary[, -1]
    mod_summary[, c("LCI_90", "UCI_90")] <- confint(mod, level = 0.90)
    mod_summary[, c("LCI_95", "UCI_95")] <- confint(mod, level = 0.95)
    mod_summary[, c("LCI_99", "UCI_99")] <- confint(mod, level = 0.99)
    return(mod_summary)
  })) %>% unnest(Coefs) %>% mutate(Significance = case_when(
      p.value <= 0.01 ~ "***"
    , p.value <= 0.05 ~ "**"
    , p.value <= 0.10 ~ "*"
    , TRUE ~ ""
  ))

# Prepare dataframe for confidence intervals
confs <- coefs_summary %>%
  select(ToD, term, Significance, LCI_90:UCI_99) %>%
  pivot_longer(LCI_90:UCI_99) %>%
  separate(name, into = c("Limit", "Level"), sep = "_", convert = T) %>%
  spread(key = Limit, value = value) %>%
  mutate(Level = as.factor(Level))

# Specify the order in which the coefficients should be plotted
order <- c(
    "(Intercept)"
  , "maxMoonlightIntensity"
  , "I(maxMoonlightIntensity^2)"
  # , "maxMoonlightIntensity:meanCloudCoverNight"
  , "maxMoonlightIntensity:maxMoonDelay"
  , "maxMoonlightIntensity:maxTemperature"
  , "maxTemperature"
  # , "meanCloudCoverNight"
  , "meanActX12"
  , "sin(2 * pi/365 * Day)"
  , "cos(2 * pi/365 * Day)"
)

# Prepare labels with "expressions"
labels <- c(
    "Intercept"
  , TeX("$_{max}MoonlightIntensity$")
  , TeX("$_{max}MoonlightIntensity^{0.5}$")
  # , TeX("$_{max}MoonlightIntensity$ x $_{mean}CloudCoverNight$")
  , TeX("$_{max}MoonlightIntensity$ x $_{max}MoonDelay$")
  , TeX("$_{max}MoonlightIntensity$ x $_{max}Temperature$")
  , TeX("$_{max}Temperature$")
  # , TeX("$_{mean}CloudCoverNight$")
  , TeX("$_{mean}Activity_{12 Hours}$")
  , TeX("$sin(\\frac{2\\pi}{365} * Day)$")
  , TeX("$cos(\\frac{2\\pi}{365} * Day)$")
)

# Plot
p1 <- ggplot() +
  geom_vline(xintercept = 0, lty = 2, size = 0.25, col = "gray70") +
  geom_errorbarh(
      aes(
      , y     = term
      , xmin  = LCI
      , xmax  = UCI
      , size  = Level
      , col   = ToD
      , group = ToD
    )
    , data     = confs
    , height   = 0
    , alpha    = 0.5
    # , position = position_dodge(width = 0.3)
  ) +
  geom_point(
      aes(
          x     = estimate
        , y     = term
        , col   = ToD
        , group = ToD
      )
    , data = coefs_summary
    , pch  = 18
  ) +
  geom_text(
    aes(
        x     = estimate
      , y     = term
      , label = Significance
      , hjust = 0.5
      , vjust = 0
      , angle = 0
      , col   = ToD
      , group = ToD
    )
    , data        = coefs_summary
    , show.legend = F
    # , position = position_dodge(width = 0.3)
  ) +
  scale_size_manual(
      name   = "Confidence Level"
    , values = c(2, 1, 0.3)
  ) +
  scale_color_viridis_d(begin = 0.2, end = 0.8) +
  scale_y_discrete(
      labels = rev(labels)
    , limits = rev(order)
  ) +
  guides(
      colour = "none"
    , size   = guide_legend(title.position = "top", title.hjust = 0.5, override.aes = list(colour = "#3CBB75FF"))
  ) +
  theme(
    , legend.position   = "bottom"
    , legend.margin     = margin(0, 50, 0, -20)
    , legend.box.margin = margin(-5, -10, -5, -10)
    , legend.text       = element_text(face  = 3)
    , legend.title      = element_text(face  = 3)
    , legend.key        = element_rect(fill  = "white")
    , axis.text.x       = element_text(angle = 45, vjust = 0.5)
    , panel.border      = element_blank()
    , axis.line         = element_line()
    , panel.background  = element_blank()
    , panel.grid.major  = element_line(color = "gray96")
    , panel.grid.minor  = element_line(color = "gray96")
    , plot.margin       = margin(1, 0, 0, 0, "cm")
  ) +
  xlab(expression(beta*"-Coefficient")) +
  ylab("Covariates")
  # coord_capped_cart(bottom = "both", left = "both")

# Store plot to file
ggsave("04_Manuscript/99_StartActivityModelCoefficients.png"
  , plot   = p1
  , width  = 4
  , height = 3
  , scale  = 1.5
  , bg     = "white"
)
