################################################################################
#### Plot of the Covariates
################################################################################
# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(terra)
library(ggplot2)
library(ggpubr)

# Specify an area to plot
ext <- ext(21.74909, 24.30089, -20.64901, -18.14901)

# Load the covariates to plot
water <- rast("03_Data/02_CleanData/01_LandCover_WaterCoverAveraged_MERGED.tif")
dista <- rast("03_Data/02_CleanData/01_LandCover_DistanceToWaterAveraged_MERGED.tif")
human <- rast("03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluenceBuff_FACEBOOK.grd")
grass <- rast("03_Data/02_CleanData/01_LandCover_NonTreeVegetationAveraged_MODIS.tif")
trees <- rast("03_Data/02_CleanData/01_LandCover_TreeCoverAveraged_MODIS.tif")

# Crop to area of interest
water <- crop(water, ext)
dista <- crop(dista, ext)
human <- crop(human[[6]], ext)
grass <- crop(grass, ext)
trees <- crop(trees, ext)

# Remove outliers for nicer plotting
upper <- quantile(values(trees), 0.99, na.rm = TRUE)
lower <- quantile(values(trees), 0.01, na.rm = TRUE)
values(trees)[values(trees) > upper] <- upper
values(trees)[values(trees) < lower] <- lower

# Create list
covars <- list(water, dista, human, grass, trees)
names(covars) <- c("Water", "Distance", "Humans", "Grass", "Trees")
plots <- list()
for (i in 1:length(covars)) {
  cov <- as.data.frame(covars[[i]], xy = T)
  names(cov) <- c("x", "y", "value")
  plots[[i]] <- ggplot(cov, aes(x = x, y = y, fill = value)) +
    geom_raster() +
    scale_fill_viridis_c() +
    theme_void() +
    theme(legend.position = "none") +
    coord_sf()
}

# Store them
for (i in 1:length(plots)) {
  name <- paste0("06_Presentation/Covariate_", names(covars)[[i]], ".png")
  ggsave(plot = plots[[i]], filename = name)
}
