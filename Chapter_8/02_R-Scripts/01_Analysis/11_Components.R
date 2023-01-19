################################################################################
#### Components
################################################################################
# Description: Use network analysis to identify "Connectivity Components"

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_8")

# Load required packages
library(terra)     # To handle spatial data
library(raster)    # To handle spatial data
library(tidyverse) # To wrangle data
library(pbmcapply) # To run stuff in parallel
library(Rcpp)      # To import C++ scripts
library(igraph)    # For network analysis

# Load some spatial features to plot
africa <- vect("03_Data/02_CleanData/Africa.shp")
water  <- vect("03_Data/02_CleanData/MajorWaters.shp")
refer  <- vect("03_Data/02_CleanData/ReferenceShape.shp")

# Load custom R functions
source("02_R-Scripts/00_Functions.R")

# Load custom C++ functions
sourceCpp("02_R-Scripts/00_Functions.cpp")

# Load the reference raster and derive an extent from it
ras <- rast("03_Data/02_CleanData/ReferenceRaster.tif")
ext <- as.polygons(ext(ras))
ext <- as(ext, "Spatial")

# Create hexagons that span the study area
# size    <- 0.05
size    <- 0.1
ver     <- spsample(ext, type = "hexagonal", cellsize = size)
ver     <- as.data.frame(ver)
ver     <- arrange(ver, desc(y), x)
ver     <- SpatialPoints(coords = ver, proj4string = CRS("+init=epsg:4326"))
hex     <- HexPoints2SpatialPolygons(ver, dx = size)
hex     <- vect(hex)
ver$ID  <- 1:length(ver)
hex$ID  <- 1:length(hex)
ver     <- vect(ver)

# Let's also rasterize the hexagons
hexr <- rasterize(hex, ras, "ID")

# Visualize
plot(hexr)
plot(hex, add = T)
plot(ext, add = T, border = "red")
# text(hex, label = "ID", cex = 0.5)
# text(ver, label = "ID", cex = 0.5, col = "blue")

################################################################################
#### Inter-Patch Connectivity between Hexagons
################################################################################
# Load dispersal simulations
sims <- read_rds("03_Data/03_Results/DispersalSimulation.rds")

# Keep only desired columns
sims <- sims[, c("x", "y", "TrackID", "StepNumber", "SourceArea", "FloodLevel")]

# For now, keep only 50 trajectories
# sims <- sims %>%
#   nest(Data = -c(FloodLevel, SourceArea, TrackID)) %>%
#   group_by(FloodLevel, SourceArea) %>%
#   sample_n(50) %>%
#   unnest(Data)

################################################################################
#### Function to Compute Betweenness
################################################################################
#' Calculate Betweenness for Simulations
#'
#' Function to compute betweenness on simulated data
#' @export
#' @param simulations \code{data.frame} containing the simulated data
#' @param raster \code{RasterLayer} based on which betweenness should be computed
#' @param steps numeric, how many steps should be considered?
#' @param area numeric, ID of the source areas that should be considered
#' @param flood character, one of "Min", "Mean", "Max"
#' @param messages boolean, should messages be printed during the rasterization
#' @param eps numeric, the interpolation distance (if desired) for point interpolation
#' @param filename character, the filename to which the raster should be stored (tempfile() by default)
#' @param mc.cores numeric, How many cores should be used?
#' @return \code{RasterLayer}
getHistory <- function(
      simulations = NULL      # Simulated trajectories
    , raster      = NULL      # Raster onto which we rasterize
    , steps       = 500       # How many steps should be considered
    , area        = NULL      # Simulations from which areas?
    , flood       = "Min"     # Which flood level?
    , messages    = T         # Print update messages?
    , eps         = NULL      # Interpolation distance
    , filename    = tempfile()
    , mc.cores    = detectCores() - 1
  ) {

  # # TESTING
  # simulations <- sims
  # raster      <- raster(hexr)
  # steps       <- 2000
  # area        <- NULL
  # flood       <- "Min"
  # messages    <- T
  # eps         <- 100
  # mc.cores    <- detectCores() / 2

  # Subset to corresponding data
  sub <- simulations[
    simulations$StepNumber <= steps &
    simulations$FloodLevel == flood, ]
  if (!is.null(area)) {
    sub <- sub[sub$SourceArea %in% area, ]
  }

  # Nest tracks by their IDs
  sub <- nest(sub, data = -TrackID)

  # Determine the visitation history of each path. Note that we will
  # "interpolate" each of the simulated steps. This will allow us to determine
  # cell-transitions at a much finer scale than if we would simply use the start
  # and endpoint of each step.
  if (messages) {
    cat("Computing visitation history...\n")
  }
  if (mc.cores > 1) {
      history <- pbmclapply(
          X                  = sub$data
        , ignore.interactive = T
        , mc.cores           = 1
        , FUN                = function(path) {
          if (!is.null(eps)) {
            path <- interpolatePath(path$x, path$y, eps = eps)
          }
          visits <- raster::extract(raster, path)
          visits <- visitHist(visits, singlecount = T)
          return(visits)
        })
    } else {
      if (messages) {
        pb <- txtProgressBar(min = 0, max = length(sub$data), style = 3)
      }
      history <- lapply(1:length(sub$data), function(x) {
        path <- sub$data[[x]]
        if (!is.null(eps)) {
          path <- interpolatePath(path$x, path$y, eps = eps)
        }
        visits <- raster::extract(raster, path)
        visits <- visitHist(visits, singlecount = T)
        if (messages) {
          setTxtProgressBar(pb, x)
        }
        return(visits)
      })
  }

  # Aggregate visitation histories across all paths
  history <- history %>%
    do.call(rbind, .) %>%
    group_by(from, to) %>%
    summarize(TotalConnections = sum(TotalConnections), .groups = "drop") %>%
    ungroup()

  # Return the history
  return(history)
}

# Get histories for minimum and maximum extent
design <- tibble(FloodLevel = c("Min", "Max"))
design$History <- lapply(design$FloodLevel, function(x) {
  hist <- getHistory(sims
    , raster   = raster(hexr)
    , steps    = 2000
    , area     = NULL
    , flood    = x
    , messages = T
    , eps      = 100
    , mc.cores = detectCores() - 1
  )
  return(hist)
})

# Create networks
design$Network <- lapply(design$History, function(x) {
  x$weight <- x$TotalConnections
  net <- graph_from_data_frame(x, vertices = ver$ID)
  return(net)
})

# Let's derive some centralization measures
design$Metrics <- lapply(design$Network, function(x) {
  degr <- centr_degree(x)$centralization
  betw <- centr_betw(x, directed = F)$centralization
  clos <- centr_clo(x)$centralization
  eige <- centr_eigen(x, directed = F)$centralization
  metrics <- data.frame(degr, betw, clos, eige)
  names(metrics) <- c("Degree", "Betweenness", "Closeness", "Eigen")
  return(metrics)
})

# Let's also prepare the five figure summary
design$FiveFigures <- lapply(design$Network, function(x) {
  nverts <- vcount(x)
  nedges <- ecount(x)
  dens <- edge_density(x)
  comp <- components(x)$no
  diam <- diameter(x, unconnected = T, directed = F)
  tran <- transitivity(x, type = "global")
  five <- data.frame(nverts, nedges, dens, comp, diam, tran)
  names(five) <- c("NVertices", "Nedges", "Density", "Components", "Diameter", "Transitivity")
  return(five)
})

# Finally, we can use different clustering algorithms to identify components
design$Components <- lapply(design$Network, function(x) {
  cl_greed <- cluster_fast_greedy(as.undirected(x, mode = c("collapse")))
  cl_louva <- cluster_louvain(as.undirected(x, mode = c("collapse")))
  cl_leide <- cluster_leiden(as.undirected(x, mode = c("collapse")))
  cl_walkt <- cluster_walktrap(as.undirected(x, mode = c("collapse")))
  cl_infom <- cluster_infomap(as.undirected(x, mode = c("collapse")))
  cl_label <- cluster_label_prop(as.undirected(x, mode = c("collapse")))
  # cl_eigen <- cluster_leading_eigen(as.undirected(net, mode = c("collapse")))
  # cl_sping <- cluster_spinglass(as.undirected(net, mode = c("collapse")))
  # cl_optim <- cluster_optimal(as.undirected(net, mode = c("collapse")))
  comps <- tibble(
      Algorithm = c("Greedy", "Louvain", "Leiden", "Walktrap", "Infomap", "LabelPropagation")
    , NGroups = c(
        length(unique(cl_greed$membership))
      , length(unique(cl_louva$membership))
      , length(unique(cl_leide$membership))
      , length(unique(cl_walkt$membership))
      , length(unique(cl_infom$membership))
      , length(unique(cl_label$membership))
    )
    , Membership = list(
        cl_greed$membership
      , cl_louva$membership
      , cl_leide$membership
      , cl_walkt$membership
      , cl_infom$membership
      , cl_label$membership
    )
    , Modularity = c(
        mean(cl_greed$modularity)
      , mean(cl_louva$modularity)
      , NA
      , mean(cl_walkt$modularity)
      , mean(cl_infom$modularity)
      , mean(cl_label$modularity)
    )
  )
  return(comps)
})

# Extract the groups
groups <- design %>%
  select(FloodLevel, Components) %>%
  unnest(Components) %>%
  unnest(Membership)

# Function to plot them
plotGroups <- function(algorithm = NULL, level = NULL) {
  sub <- subset(groups, Algorithm == algorithm & FloodLevel == level)
  cols <- sample(rainbow(length(unique(sub$Membership))))
  cols <- adjustcolor(cols, alpha.f = 0.5)
  plot(refer)
  plot(water, add = T, col = "cornflowerblue", border = NA)
  plot(hex, col = cols[as.factor(sub$Membership)], add = T, border = NA)
}

# Compare modularity of the different algorithms
design %>%
  select(FloodLevel, Components) %>%
  unnest(Components) %>%
  select(-c(NGroups, Membership)) %>%
  subset(Algorithm != "Leiden") %>%
  ggplot(aes(x = Algorithm, y = Modularity, fill = FloodLevel)) +
    geom_col(position = position_dodge())

# Try it
par(mfrow = c(2, 1))
plotGroups(algorithm = "Louvain", level = "Min")
plotGroups(algorithm = "Louvain", level = "Max")

par(mfrow = c(2, 1))
plotGroups(algorithm = "Infomap", level = "Min")
plotGroups(algorithm = "Infomap", level = "Max")

par(mfrow = c(2, 1))
plotGroups(algorithm = "LabelPropagation", level = "Min")
plotGroups(algorithm = "LabelPropagation", level = "Max")

par(mfrow = c(2, 1))
plotGroups(algorithm = "Greedy", level = "Min")
plotGroups(algorithm = "Greedy", level = "Max")

# Compare the metrics
design %>%
  select(FloodLevel, Metrics) %>%
  unnest(Metrics) %>%
  pivot_longer(Degree:Eigen, names_to = "Metric", values_to = "Value") %>%
  subset(Metric != "Closeness") %>%
  ggplot(aes(x = FloodLevel, y = Value)) +
    geom_col() +
    facet_wrap(~ Metric, scales = "free")

# Compare the five figure summaries
design %>%
  select(FloodLevel, FiveFigures) %>%
  unnest(FiveFigures) %>%
  pivot_longer(NVertices:Transitivity, names_to = "Metric", values_to = "Value") %>%
    ggplot(aes(x = FloodLevel, y = Value)) +
      geom_col() +
      facet_wrap(~ Metric, scales = "free")
