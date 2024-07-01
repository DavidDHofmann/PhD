################################################################################
#### Dynamic Variables Plots
################################################################################
# Description: Visualizations to highlight how dynamic / seasonal the system is

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Load required packages
library(terra)        # To handle raster files
library(tidyverse)    # For data wrangling
library(lubridate)    # To handle timestamps
library(ggpubr)       # To arrange plots
library(mgcv)         # For smooth lines

# Load custom functions
source("02_R-Scripts/00_Functions.R")

################################################################################
#### Precipitation
################################################################################
# Load precipitation data
files <- dir(path = "03_Data/02_CleanData/00_Rainmaps", pattern = ".tif$", full.names = T)
dat <- lapply(files, rast) %>% rast()

# Compute mean precipitation in each layer
dat <- global(dat, mean)
dat$Timestamp <- ymd_h(rownames(dat))
rownames(dat) <- NULL

# Compute monthly means
means_prec <- dat %>%
  mutate(Year = year(Timestamp), Month = month(Timestamp)) %>%
  group_by(Year, Month) %>%
  summarize(Precipitation = sum(mean), .groups = "drop") %>%
  group_by(Month) %>%
  summarize(Precipitation = mean(Precipitation))

# Store for later
dat_prec <- dat

# Predict precipitation
mod        <- gamm(Precipitation ~ s(Month, bs = "cc", k = 12), data = means_prec)
newdata    <- data.frame(Month = seq(1, 12, length.out = 1000))
fit        <- predict(mod$gam, newdata, se.fit = T)
fit_prec   <- cbind(newdata, Precipitation = fit$fit, SE = fit$se)

################################################################################
#### Temperature
################################################################################
# Load precipitation data
files <- dir(path = "03_Data/02_CleanData/00_Tempmaps", pattern = ".tif$", full.names = T)
dat <- lapply(files, rast) %>% rast()

# Compute mean precipitation in each layer
dat <- global(dat, mean)
dat$Timestamp <- ymd_h(rownames(dat))
rownames(dat) <- NULL

# Compute monthly means
means_temp <- dat %>%
  mutate(Year = year(Timestamp), Month = month(Timestamp)) %>%
  group_by(Year, Month) %>%
  summarize(
      minTemperature = min(mean)
    , maxTemperature = max(mean)
    , .groups = "drop"
  ) %>%
  group_by(Month) %>%
  summarize(
      minTemperature = mean(minTemperature)
    , maxTemperature = mean(maxTemperature)
    , .groups = "drop"
  ) %>%
  pivot_longer(2:3, names_to = "TempClass", values_to = "Temperature") %>%
  mutate(TempClass = factor(TempClass, levels = c("minTemperature", "maxTemperature")))

# Store for later
dat_temp <- dat

# Predict minimum and maximum temperature
mod        <- gamm(Temperature ~ s(Month, bs = "cc", k = 5), data = subset(means_temp, TempClass == "minTemperature"))
newdata    <- expand.grid(Month = seq(1, 12, length.out = 1000))
fit        <- predict(mod$gam, newdata, se.fit = T)
fit1       <- cbind(newdata, minTemperature = fit$fit, SE = fit$se)
mod        <- gamm(Temperature ~ s(Month, bs = "cc", k = 5), data = subset(means_temp, TempClass == "maxTemperature"))
newdata    <- expand.grid(Month = seq(1, 12, length.out = 1000))
fit        <- predict(mod$gam, newdata, se.fit = T)
fit2       <- cbind(newdata, maxTemperature = fit$fit, SE = fit$se)
fit        <- cbind(fit1, maxTemperature = fit2$maxTemperature)
fit_temp   <- pivot_longer(fit, c(minTemperature, maxTemperature), names_to = "TempClass", values_to = "Temperature")
fit_temp   <- mutate(fit_temp, TempClass = factor(TempClass, levels = c("minTemperature", "maxTemperature")))

################################################################################
#### NDVI
################################################################################
# Load ndvi data
files <- dir(path = "03_Data/02_CleanData/00_NDVI", pattern = ".tif$", full.names = T)
dat <- lapply(files, rast) %>% rast()

# Compute mean precipitation in each layer
dat_ndvi <- global(dat, mean)
dat_ndvi$Timestamp <- ymd(rownames(dat_ndvi))
rownames(dat_ndvi) <- NULL

# Compute monthly means
means_ndvi <- dat_ndvi %>%
  mutate(Year = year(Timestamp), Month = month(Timestamp)) %>%
  group_by(Year, Month) %>%
  summarize(NDVI = mean(mean), .groups = "drop") %>%
  group_by(Month) %>%
  summarize(NDVI = mean(NDVI))

# Store for later
dat_ndvi <- dat

# Predict NDVI
mod        <- gamm(NDVI ~ s(Month, bs = "cc", k = 12), data = means_ndvi)
newdata    <- data.frame(Month = seq(1, 12, length.out = 1000))
fit        <- predict(mod$gam, newdata, se.fit = T)
fit_ndvi   <- cbind(newdata, NDVI = fit$fit, SE = fit$se)

################################################################################
#### Flood
################################################################################
# Load ndvi data
files <- dir(path = "03_Data/02_CleanData/00_Floodmaps/02_Resampled", pattern = ".tif$", full.names = T)
dat <- lapply(files, rast) %>% rast()

# We only want to plot relatively cloud free images. Identify the cloud cover in
# each image for this
clouds <- c()
for (i in 1:nlyr(dat)){
  clouds[i] <- sum(values(dat[[i]]) == 2)
  clouds[i] <- clouds[i] / ncell(dat[[i]])
}

# Let's identify the number of maps with a cloud coverage above 0.1
sum(clouds > 0.1)
wat <- vect("03_Data/02_CleanData/MajorWaters.gpkg")
wat <- wat[1, ]
wat <- buffer(wat, width = 15000)

# Keep only those images with a cloud coverage below 10%
indices <- which(clouds < 0.1)
dat <- dat[[indices]]
dat <- subst(dat, 2, NA)

# Compute extent of delta
dat_flood <- dat %>%
  mask(wat) %>%
  expanse(byValue = T, unit = "km") %>%
  as.data.frame() %>%
  subset(value == 1) %>%
  mutate(Timestamp = ymd(names(dat))) %>%
  select(Timestamp, Area = area)

# Compute monthly means
means_flood <- dat_flood %>%
  mutate(Year = year(Timestamp), Month = month(Timestamp)) %>%
  group_by(Year, Month) %>%
  summarize(Flood = mean(Area), .groups = "drop") %>%
  group_by(Month) %>%
  summarize(Flood = mean(Flood))

# Store for later
dat_flood <- dat

# Predict Flood
mod        <- gamm(Flood ~ s(Month, bs = "cc", k = 12), data = means_flood)
newdata    <- data.frame(Month = seq(1, 12, length.out = 1000))
fit        <- predict(mod$gam, newdata, se.fit = T)
fit_flood  <- cbind(newdata, Flood = fit$fit, SE = fit$se)

################################################################################
#### Store Data
################################################################################
# Let's store all of the data for later
save(
    means_prec
  , means_temp
  , means_ndvi
  , means_flood
  , dat_prec
  , dat_temp
  , dat_ndvi
  , dat_flood
  , fit_prec
  , fit_temp
  , fit_ndvi
  , fit_flood
  , file = "03_Data/03_Results/SeasonalCovariates.rds"
)

# Reload
load("03_Data/03_Results/SeasonalCovariates.rds")

################################################################################
#### Visualizations
################################################################################
# Data so we can overlay the seasons
seasons <- tibble(
    Season = c("Wet", "Dry", "Wet")
  , Start  = c(1, 4, 10)
  , End    = c(4, 10, 12)
)

# Visualize precipitation
p1 <- ggplot(fit_prec, aes(
      x    = Month
    , y    = Precipitation
    , ymin = Precipitation - 1.96 * SE
    , ymax = Precipitation + 1.96 * SE
  )) +
  # geom_rect(data = subset(seasons, Season == "Dry"), aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, group = Season), inherit.aes = F, alpha = 0.1) +
  geom_rect(data = subset(seasons, Season == "Wet"), aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, group = Season), inherit.aes = F, alpha = 0.1) +
  geom_ribbon(alpha = 0.2, fill = "purple") +
  geom_point(data = means_prec, aes(x = Month, y = Precipitation), inherit.aes = F, col = "purple") +
  geom_line(color = "purple") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(3, 6, 9, 12), labels = month.abb[c(3, 6, 9, 12)]) +
  xlab("") +
  ylab("Precipitation (mm)") +
  facet_wrap(~ "Precipitation") +
  theme(strip.background = element_rect(fill = "gray95", color = "transparent"))

# Visualize temperature
p2 <- ggplot(fit_temp, aes(
      x     = Month
    , y     = Temperature
    , ymin  = Temperature - 1.96 * SE
    , ymax  = Temperature + 1.96 * SE
    , color = TempClass
    , fill  = TempClass
  )) +
  # geom_rect(data = subset(seasons, Season == "Dry"), aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, group = Season), inherit.aes = F, alpha = 0.1) +
  geom_rect(data = subset(seasons, Season == "Wet"), aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, group = Season), inherit.aes = F, alpha = 0.1) +
  geom_ribbon(alpha = 0.2, color = NA) +
  geom_point(data = means_temp, aes(x = Month, y = Temperature, col = TempClass), inherit.aes = F) +
  geom_line() +
  theme_minimal() +
  scale_fill_manual(values = c("orange", "red"), name = "", labels = c("minimum", "maximum")) +
  scale_color_manual(values = c("orange", "red"), name = "", labels = c("minimum", "maximum")) +
  theme(
      legend.position = c(0.5, 0.1)
    , legend.direction = "horizontal"
    , legend.key.height = unit(0.1, "cm")
    , panel.grid.minor = element_blank()
  ) +
  theme() +
  scale_x_continuous(breaks = c(3, 6, 9, 12), labels = month.abb[c(3, 6, 9, 12)]) +
  ylim(c(0, 40)) +
  xlab("") +
  ylab("Temperature (°C)") +
  facet_wrap(~ "Temperature") +
  theme(strip.background = element_rect(fill = "gray95", color = "transparent"))

# Visualize ndvi
p3 <- ggplot(fit_ndvi, aes(
      x    = Month
    , y    = NDVI
    , ymin = NDVI - 1.96 * SE
    , ymax = NDVI + 1.96 * SE
  )) +
  # geom_rect(data = subset(seasons, Season == "Dry"), aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, group = Season), inherit.aes = F, alpha = 0.1) +
  geom_rect(data = subset(seasons, Season == "Wet"), aes(xmin = Start, xmax = End, ymin = 0.2, ymax = Inf, group = Season), inherit.aes = F, alpha = 0.1) +
  geom_ribbon(alpha = 0.2, fill = "darkgreen") +
  geom_point(data = means_ndvi, aes(x = Month, y = NDVI), inherit.aes = F, col = "darkgreen") +
  geom_line(color = "darkgreen") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(3, 6, 9, 12), labels = month.abb[c(3, 6, 9, 12)]) +
  ylim(c(0.2, 0.62)) +
  xlab("") +
  facet_wrap(~ "NDVI") +
  theme(strip.background = element_rect(fill = "gray95", color = "transparent"))

# Visualize the flood
p4 <- ggplot(fit_flood, aes(x = Month, y = Flood, ymin = Flood - 1.96 * SE, ymax = Flood + 1.96 * SE)) +
  # geom_rect(data = subset(seasons, Season == "Dry"), aes(xmin = Start, xmax = End, ymin = 0, ymax = Inf, group = Season), inherit.aes = F, alpha = 0.1) +
  geom_rect(data = subset(seasons, Season == "Wet"), aes(xmin = Start, xmax = End, ymin = 4500, ymax = Inf, group = Season), inherit.aes = F, alpha = 0.1) +
  geom_ribbon(alpha = 0.2, fill = "cornflowerblue") +
  geom_point(data = means_flood, aes(x = Month, y = Flood), inherit.aes = F, col = "cornflowerblue") +
  geom_line(color = "cornflowerblue") +
  theme_minimal() +
  theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(breaks = c(3, 6, 9, 12), labels = month.abb[c(3, 6, 9, 12)]) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ","), limits = c(4500, 7500)) +
  xlab("") +
  ylab(expression(paste("Flood Extent (", km^2, ")"))) +
  facet_wrap(~ "Flood Extent") +
  theme(strip.background = element_rect(fill = "gray95", color = "transparent"))

################################################################################
#### Put Plots Together
################################################################################
p <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, align = "hv", labels = "auto")
ggsave("04_Manuscript/Figures/SeasonalCovariates.png"
  , plot   = p
  , width  = 5
  , height = 3
  , bg     = "white"
  , scale  = 1.35
  , device = png
)
