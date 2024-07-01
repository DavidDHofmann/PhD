################################################################################
#### Connectivity Models
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(gdistance)
library(terra)
library(ggplot2)

# Define source-target locations
pts <- rbind(
    c(0.2, 0.5)
  , c(0.8, 0.5)
)

# Set seed for reproducability
set.seed(12356)

################################################################################
#### Permeability Surface
################################################################################
# Simulate a permeability surface
n   <- 100
r   <- rast(ncol = n, nrow = n, xmin = 0, xmax = 1, ymin = 0, ymax = 1)
r[] <- abs(rnorm(ncell(r)))
r   <- focal(r, focalMat(r, d = 0.1, type = "circle"), fun = mean)
r   <- raster(r)
r   <- trim(r)
r   <- (r - minValue(r)) / (maxValue(r) - minValue(r))

################################################################################
#### Least-Cost Analysis
################################################################################
# Run least-cost model
trans <- transition(r, transitionFunction = mean, directions = 8)
trans <- geoCorrection(trans)
cost1 <- accCost(trans, fromCoords = pts[1, ])
cost2 <- accCost(trans, fromCoords = pts[2, ])
corr  <- cost1 + cost2
corr[corr > quantile(corr, 0.05)] <- NA
corr <- (corr - minValue(corr)) / (maxValue(corr) - minValue(corr))
corr[is.na(corr)] <- 1
corr <- corr * -1 + 1

################################################################################
#### Circuit Theory
################################################################################
# Run randomized shortest path
rand <- passage(trans, origin = pts[1, ], goal = pts[2, ], theta = 0.001)
rand[cellFromXY(rand, pts)] <- NA
rand <- (rand - minValue(rand)) / (maxValue(rand) - minValue(rand))
rand <- sqrt(rand)

################################################################################
#### Plots
################################################################################
# Scale all layers to values between 0 and 1
perm <- as.data.frame(r, xy = T)
corr <- as.data.frame(corr, xy = T)
rand <- as.data.frame(rand, xy = T)
perm$Type = "Permeability"
corr$Type = "Least-Cost Corridor"
rand$Type = "Conductance (Circuit Theory)"
names(corr)[3] <- names(rand)[3] <- names(perm)[3] <- "Value"

# Make points a data.frame
pts <- as.data.frame(pts)
names(pts) <- c("x", "y")

# Put them together
dat <- rbind(perm, corr, rand)
dat$Type <- factor(dat$Type, levels = c("Permeability", "Least-Cost Corridor", "Conductance (Circuit Theory)"))

# Visualize
p <- ggplot() +
  geom_raster(data = dat, aes(x = x, y = y, fill = Value)) +
  geom_point(data = pts, aes(x = x, y = y), size = 3, col = "white") +
  facet_wrap(~Type) +
  coord_sf() +
  scale_fill_viridis_c(breaks = c(0, 0.5, 1), labels = c("Low", "Medium", "High"), name = "") +
  ylim(c(0.3, 0.7)) +
  theme_minimal() +
  theme(
      axis.title.y      = element_text(angle = 0, vjust = 0.5)
    , legend.position   = "bottom"
    , legend.key.width  = unit(0.5, "cm")
    , legend.key.height = unit(0.2, "cm")
    , strip.background  = element_rect(color = "white", fill = "gray95")
    , legend.text       = element_text(size = 5)
  )

# Store
ggsave("/home/david/ownCloud/University/15. PhD/FinalThesis/Chapter_00/Figures/ConnectivityModels.png"
  , plot   = p
  , width  = 5
  , height = 1.8
  , bg     = "white"
  , scale  = 1.35
  , device = png
)
