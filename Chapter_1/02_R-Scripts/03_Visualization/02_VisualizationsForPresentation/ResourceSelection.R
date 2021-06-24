################################################################################
#### Resource Selection Example
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(raster)
library(tidyverse)
library(viridis)
library(NLMR)
library(rgeos)
library(adehabitatHR)
library(smoothr)
library(ggdark)
library(lemon)
library(sf)
library(ggpubr)

# Set a seed
set.seed(1234)

# Simulate a landscape
r <- nlm_gaussianfield(ncol = 100, nrow = 100)

# Create a circle in its middle
circle <- gCentroid(as(extent(r), "SpatialPolygons"))
circle <- gBuffer(circle, width = 30)
circle <- as(circle, "SpatialLines")

# Distribute a couple of random points on the circle
points <- as.data.frame(coordinates(spsample(circle, n = 20, type = "random")))
points$x <- points$x + rnorm(sd = 2, n = nrow(points))
points$y <- points$y + rnorm(sd = 2, n = nrow(points))

# Make a convex polygon around it
poly <- mcp(xy = SpatialPoints(points), percent = 100)

# Visualize
plot(r)
plot(poly, add = T)
points(points)

# Distribute points at high values
points <- spsample(poly, n = 1000, type = "random")
points$value <- raster::extract(r, points)
points$keep <- rbinom(n = nrow(points), size = 1, prob = points$value ** 4)
points <- points[points$keep == 1, ]
points$value <- NULL
points$keep <- NULL

# Distribute some random points
random <- spsample(poly, n = 200, type = "random")

# Put the points together
points$Type <- "Observed"
random$Type <- rep("Pseudo-Absence", length(random))
points <- rbind(points, random)

# Place a mcp around the observed dots
home <- mcp(SpatialPoints(coordinates(points)), percent = 100)
home <- smooth(home, method = "ksmooth")
home <- as(home, "SpatialLines")

# Prepare all for plotting
home <- st_as_sf(home)
rs <- as.data.frame(r, xy = T)
points <- as.data.frame(points, xy = F)

# Visualize
p1 <- ggplot(data = points, aes(x = x, y = y, col = Type)) +
  geom_sf(data = home, col = "white", inherit.aes = F, lwd = 1, lty = 2) +
  geom_point(pch = 20, size = 2) +
  scale_color_manual(values = c("orange", "gray30")) +
  dark_theme_classic() +
  coord_sf(xlim = c(0, 100), ylim = c(0, 100))
p2 <- ggplot(data = points, aes(x = x, y = y, col = Type)) +
  geom_raster(data = rs, aes(x = x, y = y, fill = layer, col = NA)) +
  # geom_sf(data = home, col = "white", inherit.aes = F, lwd = 1, lty = 2) +
  # geom_point(pch = 20, size = 2) +
  scale_color_manual(values = c("orange", "gray30")) +
  scale_fill_viridis() +
  dark_theme_classic() +
  coord_sf(xlim = c(0, 100), ylim = c(0, 100)) +
  labs(fill = "Suitability") +
  theme(
      legend.position  = "bottom"
    , legend.box       = "vertical"
  )
p <- ggarrange(p1, p2, legend = F)

# Store plots
ggsave(plot = p, "test1.png", scale = 1.2)
ggsave(plot = p2, "test2.png", scale = 1.2)

png("test1.png", width = 1080, height = 1080, pointsize = 30)
plot(r, col = "black", axes = F, box = F, legend = F)
plot(random, add = T, pch = 20, col = "gray30", cex = 0.75)
plot(points, add = T, pch = 20, col = "orange", cex = 0.75)
plot(home, add = T, border = "white", lty = 2, lwd = 2)
axis(1)
axis(2)
dev.off()

# Visualize map
png("test2.png", width = 1080, height = 1080, pointsize = 30)
plot(r, col = viridis(100), horizontal = T, box = F, axes = F)
plot(points, add = T, pch = 20, col = "orange", cex = 0.75)
plot(home, add = T, border = "white", lty = 2, lwd = 2)
axis(1)
axis(2)
dev.off()
