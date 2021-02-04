################################################################################
#### Plot Showing Behavior at Boundary
################################################################################
# Description: Plot of the source points / source areas

# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)
library(NLMR)
library(davidoff)

# Create random raster
r <- nlm_gaussianfield(ncol = 200, nrow = 200)
r <- r > 0.65
plot(r, col = c("transparent", "orange"))

# Get its extent
ext <- as(extent(r), "SpatialPolygons")

# Extend raster
r_ext <- extendRaster(r, c(-50, 250, -50, 250))

# Store plot
png("test.png", width = 1980, height = 1980, bg = "transparent")
plot(r_ext, col = c("transparent", "orange"), axes = F, legend = F, box = F)
plot(ext, border = "orange", add = T, lwd = 5)
dev.off()

# Create raster
n <- 15
r <- raster(ncol = n, nrow = n, vals = rbinom(n = n**2, size = 1, prob = 0.1), xmn = 0, xmx = n, ymn = 0, ymx = n)
png("test.png", width = 1980, height = 1980, bg = "transparent")
plot(r, col = c("gray20", "orange"), box = F, axes = F, legend = F)
plot(rasterToPolygons(r), add = T, border = "black", lwd = 8)
dev.off()
