################################################################################
#### Figures for the Graphical Abstract
################################################################################
# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse) # To wrangle data
library(raster)    # To handle spatial data
library(terra)     # To handle spatial data
library(rgeos)     # To manipulate spatial data
library(igraph)    # To create networks
library(NLMR)      # To create random landscapes
library(gdistance) # To calculate least-cost paths
library(pbmcapply) # For multicore abilities with progress bar
library(spatstat)  # To quickly rasterize lines
library(maptools)  # To quickly rasterize lines
library(viridis)   # For nice colors
library(sf)        # To plot spatial stuff with ggplot
library(ggpubr)    # To arrange multiple plots
library(ggnetwork) # To plot networks with ggplot
library(ggridges)  # For ridgeline plot
library(lemon)     # For capped coordinate system

################################################################################
#### Simulate Permeability Surface (could be any set of covariates)
################################################################################
# Create a permeability surface
set.seed(2)
perm <- nlm_gaussianfield(ncol = 100, nrow = 50)
names(perm) <- "permeability"

# Visualize
plot(perm)

# Place two source areas
source_areas <- SpatialPoints(rbind(
    c(20, 25)
  , c(80, 25)
))
source_areas$ID <- 1:length(source_areas)

# Create buffers
source_areas <- gBuffer(source_areas, byid = T, width = 6, quadsegs = 50)

# Visualize them
plot(perm)
plot(source_areas, add = T)
text(source_areas, "ID")

# Create transition matrix
t <- transition(perm, directions = 8, transitionFunction = mean)
t <- geoCorrection(t, type = "c")

# Calculate least-cost path between the two points
path <- shortestPath(t
  , origin  = gCentroid(source_areas[1, ])
  , goal    = gCentroid(source_areas[2, ])
  , output  = "SpatialLines"
)

# Buffer the path to create a corridor
corr <- gBuffer(path, width = 1.5, byid = F)

# Visualize them
plot(perm)
plot(corr, add = T)
plot(source_areas, add = T)
text(source_areas, "ID")

# Create paths between the two areas
createPath <- function(start, end, reverse = F){
  points1 <- coordinates(spsample(start, 20, type = "random"))
  points2 <- arrange(as.data.frame(coordinates(spsample(corr, 20, type = "random"))), x)
  points3 <- coordinates(spsample(end, 20, type = "random"))
  track <- rbind(points1, points2, points3)
  if (reverse){
    track <- track[nrow(track):1, ]
  }
  return(track)
}

# Try it
plot(perm)
plot(corr, add = T, border = "red", lty = 2, lwd = 2)
points(createPath(source_areas[1, ], source_areas[2, ], reverse = F), type = "o", pch = 16)
plot(perm)
plot(corr, add = T, border = "red", lty = 2, lwd = 2)
points(createPath(source_areas[1, ], source_areas[2, ], reverse = T), type = "o", pch = 16)

# Simulate several trajectories at each source location
nsims <- 5
sims <- lapply(1:nsims, function(x){

  # From source 1
  path <- createPath(source_areas[1, ], source_areas[2, ], reverse = F)
  path <- as.data.frame(path)
  path$StepNumber <- 1:nrow(path)
  path$SourceArea <- 1
  path$ID <- x
  path1 <- path

  # From source 2
  path <- createPath(source_areas[1, ], source_areas[2, ], reverse = T)
  path <- as.data.frame(path)
  path$StepNumber <- 1:nrow(path)
  path$SourceArea <- 2
  path$ID <- x
  path2 <- path

  # Put together and make IDs unique
  path <- rbind(path1, path2)

  # Return both
  return(path)

}) %>% do.call(rbind, .)

# Assign unique IDs
sims$ID <- group_indices(sims, SourceArea, ID)

# Arrange
sims <- arrange(sims, ID)

# Prepare paths as spatial lines
ids <- unique(sims$ID)
tracks <- lapply(1:length(ids), function(x){
  sub <- sims[sims$ID == ids[x], ]
  coordinates(sub) <- c("x", "y")
  track <- spLines(sub)
  track$SourceArea <- sub$SourceArea[1]
  track$ID <- sub$ID[1]
  return(track)
}) %>% do.call(rbind, .)

# Plot them
plot(tracks, col = tracks$SourceArea)

################################################################################
#### Heatmap
################################################################################
# Function to rasterize spatial lines (l) quickly onto a raster (r)
rasterizeSpatstat <- function(l, r){

  # Create im layer
  values(r) <- 0
  im <- as.im.RasterLayer(r)
  summed <- im

  # Prepare progress bar
  pb <- txtProgressBar(
      min     = 0
    , max     = length(l)
    , initial = 0
    , style   = 3
    , width   = 55
  )

  # Go through all lines and rasterize them
  for (y in 1:length(l)){
    line    <- as.psp(l[y, ], window = im)
    line    <- as.mask.psp(line)
    line_r  <- as.im.owin(line, na.replace = 0)
    summed  <- Reduce("+", list(summed, line_r))
    setTxtProgressBar(pb, y)
  }

  # Return heatmap as a raster
  return(raster(summed))

}

# Rasterize simulated trajctories
heatmap <- rasterizeSpatstat(
    l = as(tracks, "SpatialLines")
  , r = disaggregate(perm, fact = 2)
)

# Visualize them
plot(heatmap, col = viridis(50))

################################################################################
#### Betweenness
################################################################################
# Overlay the study area with a regular raster grid
regular <- raster(extent(perm), ncol = 100, nrow = 50)
values(regular) <- 1:ncell(regular)

# We want to determine all connections from one raster cell on the regular grid
# to another cell as per our simulations. For this, we'll identify the
# "visitation history" of each trajectory, i.e. the sequence of raster cells
# across which it moves.
visits <- data.frame(
    ID         = sims$ID
  , StepNumber = sims$StepNumber
  , x          = coordinates(sims)[, 1]
  , y          = coordinates(sims)[, 2]
  , cell       = raster::extract(regular, sims[, c("x", "y")])
)

# Function to retrieve the visitation history from a sequence of values
visitHist <- function(x, singlecount = F){
  transitions <- data.frame(from = lag(x), to = x) %>%
    group_by(from, to) %>%
    na.omit() %>%
    summarize(TotalConnections = n(), .groups = "drop")
  if (singlecount){
    transitions$TotalConnections = 1
  }
  return(transitions)
}

# Apply it to each simulated trajectory separately
visits <- nest(visits, data = -ID)
history <- pbmclapply(
    X                  = visits$data
  , ignore.interactive = T
  , mc.cores           = 1
  , FUN                = function(y){
      visitHist(y$cell, singlecount = T)
  }) %>%
  do.call(rbind, .) %>%
  group_by(from, to) %>%
  summarize(TotalConnections = sum(TotalConnections), .groups = "drop") %>%
  ungroup() %>%
  mutate(weight = mean(TotalConnections) / TotalConnections)

# Function to calculate betweenness
betwe <- function(network = NULL, raster = NULL){
  betweenness <- raster
  values(betweenness) <- betweenness(network)
  names(betweenness) <- "betweenness"
  return(betweenness)
}

# Coerce visits into a network (each cell in the regular raster serves as
# potential node/vertex)
vertices <- 1:ncell(regular)
net <- graph_from_data_frame(history, vertices = vertices)

# Calculate Betweenness
betweenness <- betwe(net, regular)

# Let's apply some transformations to make the graph look nicer
betweenness <- betweenness %>%
  disaggregate(method = "bilinear", fact = 2) %>%
  aggregate(fact = 2, fun = max) %>%
  sqrt()

################################################################################
#### Interpatch Connectivity
################################################################################
# Determine for each trajectory into which areas it moves
ints <- gIntersects(source_areas, tracks, byid = T)
visits <- cbind(tracks$SourceArea, ints) %>%
  as.data.frame() %>%
  gather(key = To, value = Reached, 2:ncol(.)) %>%
  setNames(c("From", "To", "Reached")) %>%
  subset(From != To & Reached == 1) %>%
  group_by(From, To) %>%
  count() %>%
  mutate(n = n / nsims)

# Let's create our own data to make it more interesting
visits <- data.frame(
    From         = c(1, 2)
  , To           = c(2, 1)
  , RelFrequency = c(0.8, 0.6)
  , Duration     = c(10, 15)
)

# Create network
net <- graph_from_data_frame(visits, vertices = source_areas$ID)
lay <- coordinates(gCentroid(source_areas, byid = T))

################################################################################
#### Visualize All Maps
################################################################################
# Visualize
par(mfrow = c(2, 2))
plot(perm)
plot(corr, add = T)
plot(source_areas, add = T)
text(source_areas, "ID")
plot(heatmap, col = viridis(50))
plot(source_areas, add = T, border = "white")
text(source_areas, "ID", col = "white")
plot(betweenness, col = magma(50))
plot(source_areas, add = T, border = "white")
text(source_areas, "ID", col = "white")
plot(perm, col = "white")
plot(net
  , layout      = lay
  , rescale     = F
  , add         = T
  , vertex.size = 500
  , edge.label  = round(E(net)$RelFrequency, 2)
  , edge.width  = E(net)$n
  , edge.curved = 0.2
  , edge.arrow.size = 0.5
)

################################################################################
#### Visualization with GGplot
################################################################################
# Convert objects for plotting with ggplot
heatmap_gg <- as.data.frame(heatmap, xy = T)
betweenness_gg <- as.data.frame(betweenness, xy = T)
source_areas_gg <- st_as_sf(source_areas)
source_areas_gg$Name <- c("A", "B")
net_gg <- ggnetwork(net, layout = lay, arrow.gap = 8, scale = F)

# Plot of heatmap
p1 <- ggplot() +
  geom_raster(
      data    = heatmap_gg
    , mapping = aes(x = x, y = y, fill = layer)
  ) +
  scale_fill_gradientn(
      # colours = viridis(100)
      colours = c("white", rev(magma(100)))
    , guide   = guide_colorbar(
      , title          = expression("Traversal Frequency")
      , show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.6, "cm")
      , barwidth       = unit(10.0, "cm")
    )
  ) +
  geom_sf(
      data        = source_areas_gg
    , fill        = NA
    , col         = "black"
    , lty         = 1
    , lwd         = 0.5
    , show.legend = "line"
  ) +
  geom_sf_text(
      data     = source_areas_gg
    , mapping  = aes(label = Name)
    , col      = "white"
    , nudge_y  = 9
    , fontface = 2
    , size     = 5
  ) +
  coord_sf(xlim = c(10, 90), ylim = c(10, 40)) +
  labs(
      x        = NULL
    , y        = NULL
    , fill     = NULL
  ) +
  theme_void() +
  theme(legend.position  = "none")

# Plot of betweenness
p2 <- ggplot() +
  geom_raster(
      data    = betweenness_gg
    , mapping = aes(x = x, y = y, fill = layer)
  ) +
  scale_fill_gradientn(
      # colours = magma(100)
      colours = c("white", rev(magma(100)))
    , guide   = guide_colorbar(
      , show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.6, "cm")
      , barwidth       = unit(10.0, "cm")
    )
  ) +
  geom_sf(
      data        = source_areas_gg
    , fill        = NA
    , col         = "black"
    , lty         = 1
    , lwd         = 0.5
    , show.legend = "line"
  ) +
  coord_sf(xlim = c(10, 90), ylim = c(10, 40)) +
  labs(
      x        = NULL
    , y        = NULL
    , fill     = NULL
  ) +
  theme_void() +
  theme(legend.position  = "none")

# Prepare a plot for interpatch connectivity (I'll do the rest in powerpoint)
p3 <- ggplot() +
  geom_raster(
      data    = betweenness_gg
    , mapping = aes(x = x, y = y, fill = layer)
  ) +
  scale_fill_gradientn(
      colours = "white"
    , guide   = guide_colorbar(
      , show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.6, "cm")
      , barwidth       = unit(10.0, "cm")
    )
  ) +
  geom_sf(
      data        = source_areas_gg
    , fill        = NA
    , col         = "black"
    , lty         = 1
    , lwd         = 0.5
    , show.legend = "line"
  ) +
  coord_sf(xlim = c(10, 90), ylim = c(10, 40)) +
  labs(
      x        = NULL
    , y        = NULL
    , fill     = NULL
  ) +
  theme_void() +
  theme(legend.position  = "none")

################################################################################
#### Plot of Arbitrary Movement Model
################################################################################
# Generate some arbitrary data to plot
x <- seq(-0.2, 0.2, by = 0.001)
y <- dnorm(x, mean = 0, sd = 0.05)
ci_90 <- qnorm(c(0.050, 0.950), mean = 0, sd = 0.05)
ci_95 <- qnorm(c(0.025, 0.975), mean = 0, sd = 0.05)
ci_99 <- qnorm(c(0.005, 0.995), mean = 0, sd = 0.05)

# Shift data to three different means
dat1 <- data.frame(Covariate = "v", Mean = -0.2, x = x - 0.2, y = y)
dat2 <- data.frame(Covariate = "w", Mean = 0.4, x = x + 0.4, y = y)
dat3 <- data.frame(Covariate = "x", Mean = -0.3, x = x - 0.3, y = y)
dat4 <- data.frame(Covariate = "y", Mean = -0.2, x = x - 0.2, y = y)
dat5 <- data.frame(Covariate = "z", Mean = 0.2, x = x + 0.2, y = y)

# Put all together
dat <- rbind(dat1, dat2, dat3, dat4, dat5)

# Calculate summaries
averaged <- data.frame(
    Covariate = c("v", "w", "x", "y", "z")
  , Mean = c(-0.2, 0.4, -0.3, -0.2, 0.2)
  , LCI_90 = ci_90[1] + c(-0.2, 0.4, -0.3, -0.2, 0.2)
  , UCI_90 = ci_90[2] + c(-0.2, 0.4, -0.3, -0.2, 0.2)
  , LCI_95 = ci_95[1] + c(-0.2, 0.4, -0.3, -0.2, 0.2)
  , UCI_95 = ci_95[2] + c(-0.2, 0.4, -0.3, -0.2, 0.2)
  , LCI_99 = ci_99[1] + c(-0.2, 0.4, -0.3, -0.2, 0.2)
  , UCI_99 = ci_99[2] + c(-0.2, 0.4, -0.3, -0.2, 0.2)
)

# Generate plot
p4 <- ggplot(dat, aes(x = x, y = Covariate, height = y, group = Covariate)) +
  geom_errorbarh(
      data        = averaged
    , inherit.aes = F
    , height      = 0
    , size        = 0.5
    , alpha       = 1
    , col         = "orange"
    , aes(y = Covariate, xmin = LCI_99, xmax = UCI_99)
  ) +
  geom_errorbarh(
      data        = averaged
    , inherit.aes = F
    , height      = 0
    , size        = 1
    , alpha       = 1
    , col         = "orange"
    , aes(y = Covariate, xmin = LCI_95, xmax = UCI_95)
  ) +
  geom_errorbarh(
      data        = averaged
    , inherit.aes = F
    , height      = 0
    , size        = 2
    , alpha       = 0.5
    , col         = "orange"
    , aes(y = Covariate, xmin = LCI_90, xmax = UCI_90)
  ) +
  geom_point(
      data        = averaged
    , inherit.aes = F
    , size        = 3
    , col         = "orange"
    , aes(x = Mean, y = Covariate)
  ) +
  geom_vline(xintercept = 0, lty = 2, col = "gray50") +
  theme_classic() +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
    , xlim   = c(-0.5, 0.5)
  ) +
  labs(x = expression(beta*"-Coefficient"))

# Show plots
p1 <- p1 +
  labs(title = "Heatmap", subtitle = "(Traversal Frequency)") +
  theme(plot.title = element_text(face = 2, hjust = 0.5), plot.subtitle = element_text(face = 3, hjust = 0.5))
p2 <- p2 +
  labs(title = "Betweenness", subtitle = "(Bottlenecks & Dispersal Corridors)") +
  theme(plot.title = element_text(face = 2, hjust = 0.5), plot.subtitle = element_text(face = 3, hjust = 0.5))
p3 <- p3 +
  labs(title = "Inter-Patch Connectivity", subtitle = "(Frequency & Duration of Movements between Patches)") +
  theme(plot.title = element_text(face = 2, hjust = 0.5), plot.subtitle = element_text(face = 3, hjust = 0.5))
p4 <- p4 +
  labs(title = "Integrated Step\nSelection Model") +
  theme(plot.title = element_text(face = 2, hjust = 0.5), plot.subtitle = element_text(face = 3, hjust = 0.5))

# Store the plots
ggsave("test.png", plot = p1)
ggsave("test2.png", plot = p2)
ggsave("test3.png", plot = p3)
ggsave("test4.png", plot = p4, height = 3, width = 4)
