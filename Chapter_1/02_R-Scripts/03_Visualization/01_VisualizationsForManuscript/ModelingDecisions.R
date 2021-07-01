################################################################################
#### Figures for Modeling Decisions
################################################################################
# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse) # To wrangle data
library(raster)    # To handle spatial data
library(ggpubr)    # To arrange multiple plots
library(lemon)     # For capped coordinates
library(sf)        # For plotting spatial things in ggplot
library(ggridges)  # For ridgeline plot
library(NLMR)      # To simulate landscapes
library(davidoff)  # Custom functions
library(png)       # To load png image
library(smoothr)   # To smooth spatial polygons

################################################################################
#### Number of Simulated Individuals
################################################################################
# Simulate some traversal data
dat <- expand_grid(
    Nindivs   = seq(0, 500, 10)
  , Replicate = 1:100
)
dat$MeanTraversalFrequency <- sapply(dat$Nindivs, function(x){
  traversals <- rbinom(n = x, prob = 0.4, size = 1)
  traversals <- mean(traversals)
  return(traversals)
})
dat$MeanTraversalFrequency[dat$Nindivs == 0] <- 0

# Simulate some checkpoints
check_r <- raster(ncol = 10, nrow = 10, xmn = 0, xmx = 10, ymn = 0, ymx = 10)
check_p <- as(check_r, "SpatialPolygons")
values(check_r) <- rbinom(size = 1, n = 100, prob = 0.2)
plot(check_r, col = c("white", "orange"))
plot(check_p, add = T, lwd = 0.5)

# Convert to sf
check_rsf <- as.data.frame(check_r, xy = T)
check_rsf$layer <- as.factor(check_rsf$layer)
check_psf <- st_as_sf(as(check_p, "SpatialLines"))

# Plot the checkpoints
p1 <- ggplot() +
  geom_raster(data = check_rsf, aes(x = x, y = y, fill = layer)) +
  geom_sf(data = check_psf, size = 0.2, alpha = 0.2) +
  scale_fill_manual(values = c("white", "black")) +
  theme_void() +
  theme(legend.position = "none")

# Plot the relative traversal frequency
p2 <- dat %>%
  group_by(Nindivs) %>%
  summarize(
      UCI = quantile(MeanTraversalFrequency, 0.975)
    , LCI = quantile(MeanTraversalFrequency, 0.025)
    , MeanTraversalFrequency = mean(MeanTraversalFrequency)
  ) %>%
  ggplot(aes(x = Nindivs, y = MeanTraversalFrequency)) +
  geom_ribbon(aes(
      ymin = LCI
    , ymax = UCI
  ), alpha = 0.5, fill = "orange", color = "orange") +
  geom_line() +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
  ) +
  scale_x_continuous(
      labels = function(x){format(x, big.mark = "'")}
    , breaks = c(0, 250, 500, 750, 1000)
  ) +
  scale_y_continuous(
      labels = function(x){format(x, big.mark = "'")}
    , breaks = c(0, 0.25, 0.5, 0.75, 1)
  ) +
  ylim(c(0, 0.9)) +
  xlab("# Trajectories") +
  ylab("Rel. Traversal Frequency") #+
  # annotation_custom(ggplotGrob(p1), xmin = 500, xmax = 1000, ymin = 0, ymax = 0.3)

################################################################################
#### Individual Variability
################################################################################
# Generate some arbitrary data to plot
x <- seq(-0.2, 0.2, by = 0.001)
y <- dnorm(x, mean = 0, sd = 0.05)
ci <- qnorm(c(0.025, 0.975), mean = 0, sd = 0.05)

# # Also sample some points
samp1 <- rnorm(20, mean = 0, sd = 0.05)
samp2 <- rnorm(20, mean = 0, sd = 0.05)
samp3 <- rnorm(20, mean = 0, sd = 0.05)

# Shift data to three different means
dat1 <- data.frame(Covariate = "x", Mean = -0.2, x = x - 0.2, y = y)
dat2 <- data.frame(Covariate = "y", Mean = 0.4, x = x + 0.4, y = y)
dat3 <- data.frame(Covariate = "z", Mean = -0.3, x = x - 0.3, y = y)
samp1 <- data.frame(Covariate = "x", Estimates = samp1 - 0.2)
samp2 <- data.frame(Covariate = "y", Estimates = samp2 + 0.4)
samp3 <- data.frame(Covariate = "z", Estimates = samp3 - 0.3)

# Put all together
dat <- rbind(dat1, dat2, dat3)
samp <- rbind(samp1, samp2, samp3)

# Calculate summaries
averaged <- data.frame(
    Covariate = c("x", "y", "z")
  , Mean = c(-0.2, 0.4, -0.3)
  , LCI = ci[1] + c(-0.2, 0.4, -0.3)
  , UCI = ci[2] + c(-0.2, 0.4, -0.3)
)

# Generate plot
p3 <- ggplot(dat, aes(x = x, y = Covariate, height = y, group = Covariate)) +
  geom_ridgeline(scale = 0.1, col = "orange", fill = "orange", alpha = 0.5) +
  geom_jitter(
      data        = samp
    , inherit.aes = F
    , position    = position_nudge(y = -0.1)
    , size        = 1
    , col         = "black"
    , alpha       = 0.2
    , aes(x = Estimates, y = Covariate)
  ) +
  geom_point(
      data        = averaged
    , inherit.aes = F
    , position    = position_nudge(y = -0.1)
    , size        = 2
    , aes(x = Mean, y = Covariate)
  ) +
  geom_errorbarh(
      data        = averaged
    , inherit.aes = F
    , position    = position_nudge(y = -0.1)
    , height      = 0
    , aes(y = Covariate, xmin = LCI, xmax = UCI)
  ) +
  geom_vline(xintercept = 0, lty = 2, col = "gray50") +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
    , xlim   = c(-0.5, 0.5)
  ) +
  labs(x = expression(beta*"-Coefficient"))

################################################################################
#### Dispersal Duration
################################################################################
# Sample some dispersal durations
set.seed(1234)
durs <- data.frame(Duration = rgamma(100, scale = 25, shape = 5))

# Define parametric distribution
x <- 0:500
y <- dgamma(x, scale = 25, shape = 5)
dat <- data.frame(x, y)

# Visualize
p4 <- ggplot(dat, aes(x = x, y = y)) +
  geom_area(col = "orange", fill = "orange", alpha = 0.5) +
  geom_rug(data = durs, inherit.aes = F, aes(x = Duration), col = "gray30") +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
  ) +
  geom_vline(xintercept = 500, lty = 2, col = "gray50") +
  geom_segment(
      aes(x = 440, xend = 490, y = 0.006, yend = 0.006)
    , arrow = arrow(length = unit(0.02, "npc"), type = "closed")
  ) +
  annotate("text"
    , x        = 295
    , y        = 0.006
    , fontface = 3
    , size     = 3
    , label    = "Simulate maximal\ndispersal duration\nfirst, then subsample from\nsimulated trajectories"
  ) +
  ylab("Density") +
  xlab("Observed Dispersal Duration (in steps)")

################################################################################
#### Boundary Behavior
################################################################################
# Load a picture of a wild dog that we want to add to the plot
dog <- readPNG("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/01_RawData/DAVID/WildDogTop.png")

# Simulate landscape
r <- nlm_gaussianfield(ncol = 200, nrow = 100)
r <- r > 0.65
plot(r, col = c("transparent", "orange"))

# Get its extent
ext <- as(extent(r), "SpatialPolygons")

# Extend raster
r_ext <- extendRaster(r, c(-20, 220, -20, 120))
plot(r_ext)

# Get the buffer
buffer <- as(extent(r_ext), "SpatialPolygons") - ext

# Generate a set of "random steps"
from <- data.frame(x = 100, y = 60)
to <- data.frame(x = c(100, 80, 180, 50, 120, 70), y = c(115, 90, 80, 50, 40, 10))
steps <- lapply(1:nrow(to), function(x){
  spLines(SpatialPoints(rbind(from, to[x, ])))
}) %>% do.call(rbind, .)

# Convert to sf
r_sf <- as.data.frame(r_ext, xy = T)
ext_sf <- st_as_sf(ext)
steps_sf <- st_as_sf(steps)
buffer_sf <- st_as_sf(buffer)

# Plot everything
p5 <- ggplot() +
  geom_raster(data = r_sf, aes(x = x, y = y, fill = layer), alpha = 0.5) +
  geom_sf(data = buffer_sf, col = "black", fill = NA, lty = 2) +
  geom_sf(data = steps_sf[2:nrow(steps_sf), ], col = "black") +
  geom_sf(data = steps_sf[1, ], col = "red") +
  geom_point(data = from, aes(x = x, y = y)) +
  scale_fill_manual(values = c("white", "orange")) +
  theme_void() +
  coord_sf(
      xlim = c(-20, 220)
    , ylim = c(-20, 120)
  ) +
  theme(legend.position = "none") +
  annotation_raster(dog, xmin = 77, xmax = 120, ymin = 30, ymax = 90) +
  annotate("label", x = 0, y = 110, label = "BUFFER", fontface = 3, alpha = 0.8, size = 2.5)

################################################################################
#### Selection of Source Points
################################################################################
# Generate three habitat patches
size <- 3.8
pol1 <- data.frame(
    x = size * c(0, 02, 05, 10, 20, 26, 20, 19, 18, 16, 05)
  , y = size * c(0, 10, 13, 07, 12, 08, 07, 04, 01, 03, 03)
)
pol2 <- data.frame(
    x = size * c(0, 00, 04, 05, 10, 15, 15, 16, 07, 04) + 80
  , y = size * c(0, 02, 04, 10, 10, 09, 05, 03, 04, 01) + 60
)
pol3 <- data.frame(
    x = size * c(4, 01, 00, 03, 04, 10, 12, 14, 16, 14, 15) + 140
  , y = size * c(0, 05, 06, 12, 10, 15, 08, 10, 08, 04, 02) + 000
)
pol1 <- Polygon(pol1)
pol2 <- Polygon(pol2)
pol3 <- Polygon(pol3)
pol1 <- Polygons(list(pol1), 1)
pol2 <- Polygons(list(pol2), 2)
pol3 <- Polygons(list(pol3), 3)

# Put them into a dataframe
pols <- SpatialPolygons(list(pol1, pol2, pol3))

# Smooth
pols <- smooth(pols, method = "ks")

# Get extent
ext <- extent(c(0, 200, 0, 100))
ext <- as(ext, "SpatialPolygons")

# Create buffer
buffer <- extent(c(-20, 220, -20, 120))
buffer <- as(buffer, "SpatialPolygons") - ext

# Distribute some randomly placed source points
pts <- spsample(pols, n = 500, type = "random")
pts_buffer <- spsample(buffer, n = 500, type = "random")

# Convert stuff to sf
pols_sf <- st_as_sf(pols)
ext_sf <- st_as_sf(ext)
buffer_sf <- st_as_sf(buffer)
pts_sf <- st_as_sf(pts)
pts_buffer_sf <- st_as_sf(pts_buffer)

# Prepare plot
p6 <- ggplot() +
  geom_sf(data = pts_sf, col = "orange", size = 0.5) +
  geom_sf(data = pts_buffer_sf, col = "orange", alpha = 0.5, size = 0.5) +
  geom_sf(data = pols_sf, col = "black", fill = "orange", alpha = 0.3) +
  geom_sf(data = buffer_sf, lty = 2, col = "black", fill = NA) +
  theme_void() +
  coord_sf(
      xlim = c(-20, 220)
    , ylim = c(-20, 120)
  ) +
  annotate("label", x = 40, y = 20, label = "PATCH A", fontface = 3, alpha = 0.8, size = 2.5) +
  annotate("label", x = 118, y = 85, label = "PATCH B", fontface = 3, alpha = 0.8, size = 2.5) +
  annotate("label", x = 170, y = 25, label = "PATCH C", fontface = 3, alpha = 0.8, size = 2.5) +
  annotate("label", x = 0, y = 110, label = "BUFFER", fontface = 3, alpha = 0.8, size = 2.5)

################################################################################
#### Combine Plots
################################################################################
# Add plot titles
p2 <- p2 +
  labs(title = "Decision 1", subtitle = "Number of Simulated Trajectories") +
  theme(plot.title = element_text(face = 2), plot.subtitle = element_text(face = 3))
p3 <- p3 +
  labs(title = "Decision 2", subtitle = "Handling Individual Variability") +
  theme(plot.title = element_text(face = 2), plot.subtitle = element_text(face = 3))
p4 <- p4 +
  labs(title = "Decision 3", subtitle = "Dispersal Duration") +
  theme(plot.title = element_text(face = 2), plot.subtitle = element_text(face = 3))
p5 <- p5 +
  labs(title = "Decision 4", subtitle = "Boundary Behavior") +
  theme(plot.title = element_text(face = 2), plot.subtitle = element_text(face = 3))
p6 <- p6 +
  labs(title = "Decision 5", subtitle = "Source Point Location") +
  theme(plot.title = element_text(face = 2), plot.subtitle = element_text(face = 3))

# Show plots individually
p2
p3
p4
p5
p6

# Put plots together
p7 <- ggarrange(p2, p3, p4, ncol = 3)
p8 <- ggarrange(p5, p6, ncol = 2)
p9 <- ggarrange(p7, p8, nrow = 2)

# Also add the checkpoints
p10 <- p9 +
  annotation_custom(ggplotGrob(p1), xmin = 0.2, xmax = 0.3, ymin = 0.58, ymax = 0.68) +
  annotate("text"
    , x        = 0.25
    , y        = 0.69
    , fontface = 3
    , size     = 2.3
    , label    = "Checkpoints"
  )

# Store them
ggsave("04_Manuscript/99_ModelingDecisions.png", plot = p10, height = 6, width = 10)
