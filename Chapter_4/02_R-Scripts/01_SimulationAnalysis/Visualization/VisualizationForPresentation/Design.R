################################################################################
#### Design
################################################################################
# Description: Visualization(s) of the study design

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)        # To handle spatial data
library(tidyverse)     # For data wrangling
library(sf)            # For plotting
library(ggpubr)        # To extract legends
library(amt)           # To compute step metrics
library(pbmcapply)     # For multicore abilities
library(ggdark)        # For dark themes

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

# Load custom functions
source("02_R-Scripts/Functions.R")

################################################################################
#### Autocorrelation Scenarios
################################################################################
# Load observed movement data and covariates
dat <- "03_Data/Simulation.rds" %>%
  read_rds() %>%
  subset(Replicate == 2) %>%
  dplyr::select(AutocorrRange, Covariates, Movement)

# Prepare the data
dat_cov <- dat %>%
  dplyr::select(AutocorrRange, Covariates) %>%
  mutate(Covariates = map(Covariates, function(x) {
    as.data.frame(x, xy = T)
  })) %>%
  unnest(Covariates) %>%
  gather(key = covariate, value = value, 4:6) %>%
  nest(Data = -c(AutocorrRange, covariate))

# Create maps
dir.create("05_Presentation/99_Design", showWarnings = F)
lapply(1:nrow(dat_cov), function(x) {
  p <- ggplot(dat_cov$Data[[x]], aes(x = x, y = y, fill = value)) +
    geom_raster() +
    scale_fill_viridis_c(option = "viridis") +
    coord_sf() +
    theme_minimal() +
    scale_x_continuous(breaks = seq(0, 300, by = 100), limits = c(0, 300), expand = c(0, 0)) +
    scale_y_continuous(breaks = seq(0, 300, by = 100), limits = c(0, 300), expand = c(0, 0)) +
    theme_void() +
    theme(legend.position = "none")
  filename <- paste0(dat_cov$covariate[x], "_", dat_cov$AutocorrRange[x], ".png")
  ggsave(file.path("05_Presentation/99_Design", filename), bg = "white", scale = 2)
})

# Let's also store one legend
l <- ggplot(dat_cov$Data[[1]], aes(x = x, y = y, fill = value)) +
  geom_raster() +
  scale_fill_viridis_c(option = "viridis", "covariate value") +
  coord_sf() +
  theme_minimal() +
  scale_x_continuous(breaks = seq(0, 300, by = 100), limits = c(0, 300), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0, 300, by = 100), limits = c(0, 300), expand = c(0, 0)) +
  dark_theme_void() +
  theme(
      legend.position = "bottom"
    , legend.key.width = unit(1.5, "cm")
    , legend.key.height = unit(1, "cm")
  ) +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = 0.5))
l <- get_legend(l)

# Store
ggsave("05_Presentation/99_Design/CovariatesLegend.png"
  , plot   = l
  , scale  = 1.4
  , device = png
  , bg     = "transparent"
  , width  = 2.5
  , height = 0.75
)

################################################################################
#### Missingness Scenarios
################################################################################
# Simulate data
dat <- expand_grid(x = 1:20, y = 1:10, Missingness = c(0, 0.1, 0.2, 0.3, 0.4, 0.5))
dat <- nest(dat, Data = -Missingness)
dat <- mutate(dat, Data = map2(Missingness, Data, function(x, y) {
  index <- sample(nrow(y), size = nrow(y) * x)
  y$Obtained <- T
  y$Obtained[index] <- F
  return(y)
})) %>% unnest()

# Prepare plot
p  <- ggplot(dat, aes(x = x, y = y, color = Obtained)) +
  geom_point(size = 2.5, pch = 20) +
  facet_wrap(~ Missingness, ncol = 1) +
  scale_color_manual(values = c("gray30", "white"), name = "location obtained") +
  dark_theme_void() +
  theme(
      legend.position  = "bottom"
    , panel.spacing    = unit(1, "cm")
    , strip.background = element_blank()
    , plot.background  = element_blank()
    , panel.background = element_blank()
    , strip.text.x     = element_blank()
  ) +
  guides(color = guide_legend(title.position = "top", title.hjust = 0.5)) +
  coord_cartesian(clip = "off")

# Grab legend
l <- get_legend(p)
p <- p + theme(legend.position = "none")

# Store main plot
ggsave("05_Presentation/99_Design/Missingness.png"
  , plot   = p
  , scale  = 1.4
  , device = png
  , bg     = "transparent"
  , width  = 1.5
  , height = 4.5
)

# Store legend
ggsave("05_Presentation/99_Design/MissingnessLegend.png"
  , plot   = l
  , scale  = 1.4
  , device = png
  , bg     = "transparent"
  , width  = 1.2
  , height = 0.5
)

################################################################################
#### Forgiveness Scenarios
################################################################################
# Simulate a trajectory with a tendency to move to the right
set.seed(1234)
sim <- lapply(1:20, function(i) {
  xyt <- data.frame(
      x  = rnorm(1, mean = i, sd = 2)
    , y  = rnorm(1, mean = 0, sd = 1)
    , dt = rpois(1, lambda = 1.8) + 1
  )
}) %>% do.call(rbind, .)

# Compute from-to coordinates
sim$x_to <- lead(sim$x)
sim$y_to <- lead(sim$y)
sim <- sim[-nrow(sim), ]

# Assign different forgiveness values
sim <- tibble(
    Forgiveness = 1:5
  , Simulation  = lapply(Forgiveness, function(x) {sim})
) %>% unnest(cols = Simulation)

# Check if step-duration aligns with forgiveness
sim$Retained <- sim$dt <= sim$Forgiveness

# Visualize
p <- ggplot(data = sim) +
  geom_segment(mapping = aes(x = x, y = y, xend = x_to, yend = y_to, color = as.factor(dt), linetype = Retained)) +
  dark_theme_void() +
  scale_color_viridis_d(name = "step duration", direction = -1) +
  scale_linetype_manual(values = c(3, 1), name = "step retained") +
  facet_wrap(~ Forgiveness, ncol = 1) +
  theme(
      legend.position  = "bottom"
    , legend.box       = "vertical"
    , panel.spacing    = unit(1, "cm")
    , strip.background = element_blank()
    , strip.text.x     = element_blank()
    , panel.background = element_blank()
    , plot.background  = element_blank()
  ) +
  guides(
      color    = guide_legend(title.position = "top", title.hjust = 0.5)
    , linetype = guide_legend(title.position = "top", title.hjust = 0.5)
  )

# Grab legend
l <- get_legend(p)
p <- p + theme(legend.position = "none")

# Store main plot
ggsave("05_Presentation/99_Design/Forgiveness.png"
  , plot   = p
  , scale  = 0.8
  , device = png
  , bg     = "transparent"
  , width  = 2.25
  , height = 7
)

# Store legend
ggsave("05_Presentation/99_Design/ForgivenessLegend.png"
  , plot   = l
  , scale  = 1.4
  , device = png
  , bg     = "transparent"
  , width  = 1.8
  , height = 0.8
)

################################################################################
#### Model Approaches: Uncorrected
################################################################################
# Load simulations and rarify data
sims <- "03_Data/Simulation.rds" %>%
  read_rds() %>%
  pull(MovementFilename) %>%
  lapply(read_rds) %>%
  do.call(rbind, .) %>%
  rarifyData(missingness = 0.5) %>%
  computeBursts(max_duration = 5) %>%
  computeMetrics()

# Fit distributions to the different step durations
params <- lapply(1:5, function(x) {
  sub <- sims[sims$duration == x, ]
  fit_sl <- fit_distr(sub$sl, dist_name = "gamma")
  fit_ta <- fit_distr(sub$relta, dist_name = "vonmises")
  fit <- tibble(
      duration = x
    , shape    = fit_sl$params$shape
    , scale    = fit_sl$params$scale
    , kappa    = fit_ta$params$kappa
    , mu       = as.numeric(fit_ta$params$mu)
  )
  return(fit)
}) %>% do.call(rbind, .)

# Prepare density lines for each step duration
params$dens <- lapply(1:nrow(params), function(x) {
  dens <- tibble(
      sl    = seq(0, 20, length.out = 1000)
    , ta    = seq(-pi, pi, length.out = 1000)
    , pr_sl = dgamma(sl, shape = params$shape[x], scale = params$scale[x])
    , pr_ta = dvonmises(ta, kappa = params$kappa[x], mu = params$mu[x])
  )
  return(dens)
})

# Prepare a dataframe for plotting
toplot <- params %>%
  select(duration, dens) %>%
  unnest(dens) %>%
  mutate(used = duration == 1)

# Define limits
ylims_sl <- c(min(toplot$pr_sl), max(toplot$pr_sl))
ylims_ta <- c(min(toplot$pr_ta), max(toplot$pr_ta))

# Create a legend
legend <- ggplot(toplot, aes(x = sl, y = pr_sl, color = as.factor(duration))) +
  geom_line(linewidth = 1.3) +
  dark_theme_minimal() +
  scale_color_viridis_d(name = "Step-Duration", direction = -1) +
  theme(
    legend.position = "bottom", panel.grid.minor = element_blank()
  )
legend <- get_legend(legend)
legend <- ggarrange(legend)

# Store the legend
ggsave(
    filename = "05_Presentation/99_Design/ApproachLegend.png"
  , plot     = legend
  , width    = 3.5
  , height   = 0.3
  , bg       = "transparent"
  , device   = "png"
)

# Visualize densities
p1 <- ggplot(toplot, aes(x = sl, y = pr_sl, color = as.factor(duration), alpha = used)) +
  geom_line(linewidth = 1.3) +
  dark_theme_minimal() +
  scale_color_viridis_d(direction = -1) +
  scale_alpha_discrete(range = c(0.2, 1)) +
  scale_y_continuous(breaks = seq(0, 0.2, length.out = 3), limits = c(0, 0.3)) +
  theme(
      legend.position  = "none"
    , panel.grid       = element_blank()
    , panel.background = element_blank()
    , plot.background  = element_blank()
  ) +
  xlab("sl") +
  ylab("density")
p2 <- ggplot(toplot, aes(x = ta, y = pr_ta, color = as.factor(duration), alpha = used)) +
  geom_line(linewidth = 1) +
  dark_theme_minimal() +
  scale_color_viridis_d(direction = -1) +
  scale_alpha_discrete(range = c(0.2, 1)) +
  scale_x_continuous(
      breaks = c(-pi, -pi/2, 0, pi/2, pi)
    , labels = c(expression(-pi, -pi/2, 0, pi/2, pi))
  ) +
  scale_y_continuous(breaks = seq(0, 0.2, length.out = 3), limits = c(0, 0.3)) +
  theme(
      legend.position  = "none"
    , panel.grid       = element_blank()
    , panel.background = element_blank()
    , plot.background  = element_blank()
  ) +
  xlab("ta") +
  ylab("")

# Arrange them on a single plot and store them
p_uncorrected <- ggarrange(p1, p2, ncol = 2, align = "hv") + theme_minimal()
ggsave("05_Presentation/99_Design/ApproachUncorrected.png"
  , plot   = p_uncorrected
  , bg     = "transparent"
  , width  = 15
  , height = 4.5
  , device = png
  , scale  = 0.35
)

################################################################################
#### Model Approaches: Naive
################################################################################
# Exactly the same plot (we will later add a multiplier in the figure)
ggsave("05_Presentation/99_Design/ApproachNaive.png"
  , plot   = p_uncorrected
  , bg     = "transparent"
  , width  = 15
  , height = 4.5
  , device = png
  , scale  = 0.35
)

################################################################################
#### Model Approaches: Dynamic + Model
################################################################################
# Visualize densities
p1 <- ggplot(toplot, aes(x = sl, y = pr_sl, color = as.factor(duration))) +
  geom_line(linewidth = 1.3) +
  dark_theme_minimal() +
  scale_color_viridis_d(direction = -1) +
  scale_y_continuous(breaks = seq(0, 0.2, length.out = 3), limits = c(0, 0.3)) +
  theme(
      legend.position  = "none"
    , panel.grid       = element_blank()
    , panel.background = element_blank()
    , plot.background  = element_blank()
  ) +
  xlab("sl") +
  ylab("density")
p2 <- ggplot(toplot, aes(x = ta, y = pr_ta, color = as.factor(duration))) +
  geom_line(linewidth = 1) +
  dark_theme_minimal() +
  scale_color_viridis_d(direction = -1) +
  scale_x_continuous(
      breaks = c(-pi, -pi/2, 0, pi/2, pi)
    , labels = c(expression(-pi, -pi/2, 0, pi/2, pi))
  ) +
  scale_y_continuous(breaks = seq(0, 0.2, length.out = 3), limits = c(0, 0.3)) +
  theme(
      legend.position  = "none"
    , panel.grid       = element_blank()
    , panel.background = element_blank()
    , plot.background  = element_blank()
  ) +
  xlab("ta") +
  ylab("")

# Arrange them on a single plot and store them
p_dynamic <- ggarrange(p1, p2, ncol = 2, align = "hv") + theme_minimal()
ggsave("05_Presentation/99_Design/ApproachDynamic.png"
  , plot   = p_dynamic
  , bg     = "transparent"
  , width  = 15
  , height = 4.5
  , device = png
  , scale  = 0.35
)

################################################################################
#### Remaining Plots Done in LibreImpress
################################################################################
