################################################################################
#### Network Metrics Example
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(igraph)
library(raster)
library(tidyverse)
library(viridis)
library(davidoff)

# Create a grid
r <- raster(ncol = 10, nrow = 10, vals = 1:100, xmn = 0, xmx = 10, ymn = 0, ymx = 10)

# Visualize it
plot(r, box = F, legend = F, axes = F)
text(r)

# Create example visitation histories
hist1 <- c(1, 2, 12, 11, 12, 1, 2, 12, 11, 12, 13, 22, 32, 42, 43, 44, 55
  , 65, 66, 76, 86, 87, 88, 98, 99, 100, 90, 89, 90, 100, 99, 89, 90, 89, 99)
hist2 <- c(10, 9, 8, 18, 17, 7, 6, 5, 15, 25, 34, 44, 55, 65, 66, 76, 86, 85
  , 84, 83, 93, 92, 82, 83, 93, 94, 83, 82, 83, 92)

# Function to retrieve visitation history
visitHist <- function(x){
  transitions <- data.frame(from = lag(x), to = x) %>%
    group_by(from, to) %>%
    na.omit() %>%
    summarize(TotalConnections = n())
  return(transitions)
}

# Apply function to our histories
hist1 <- visitHist(hist1)
hist2 <- visitHist(hist2)

# We will only keep unique transitions
hist1$TotalConnections <- 1
hist2$TotalConnections <- 1

# Put histories together
hist <- rbind(hist1, hist2) %>%
  group_by(from, to) %>%
  summarize(weight = sum(TotalConnections))

# Generate graph
graph <- graph_from_data_frame(hist, vertices = 1:100)
lay <- as.matrix(as.data.frame(r, xy = T)[, c("x", "y")])
is.weighted(graph)

# Calculate network metrics
netMet <- function(
    network   = NULL    # The network based on which the metrics are calculated
  , raster    = NULL    # The raster onto which the metrics are calculated
  , tempfile  = F       # Should the resulting raster go to a temporary file?
  , metrics = c("betweenness", "closeness", "degree")
  ){

    # Calculate the desired network metrics and put them as raster values
    result <- vector(mode = "list", length = 3)
    if ("betweenness" %in% metrics){
      betweenness <- raster
      values(betweenness) <- betweenness(network)
      names(betweenness) <- "betweenness"
      result[[1]] <- betweenness
    }
    if ("closeness" %in% metrics){
      closeness <- raster
      values(closeness) <- closeness(network)
      names(closeness) <- "closeness"
      result[[2]] <- closeness
    }
    if ("degree" %in% metrics){
      degree <- raster
      values(degree) <- degree(network)
      names(degree) <- "degree"
      result[[3]] <- degree
    }

    # Remove NULLs from the lsit
    result <- plyr::compact(result)

    # Put all into a stack
    result <- stack(result)

    # In case the user desires to store a temporary file, do so
    if (tempfile){
      result <- writeRaster(result, tempfile())
    }

    # Return the final stack
    return(result)
}

# Apply the function to our data
metrics <- netMet(network = graph, raster = r, metrics = c("betweenness", "degree"))

################################################################################
#### Visualizations
################################################################################
# Visualize graph
png("test1.png", width = 1080, heigh = 1080, pointsize = 30)
par(mar = c(1, 1, 1, 0))
plot(r
  , box    = F
  , legend = F
  , axes   = F
  , col = viridis(100, end = 0.6)
)
plot(
    graph
  , layout             = lay
  , add                = T
  , rescale            = F
  , vertex.label       = NA
  , vertex.color       = "transparent"
  , vertex.frame.color = NA
  , edge.color         = "white"
  , edge.arrow.size    = 0.7
  , edge.width         = E(graph)$weight * 1.5
)
text(r, col = colTrans("white", percent = 20))
dev.off()

# Visualize Betweenness
png("test2.png", width = 1080, heigh = 1080, pointsize = 30)
par(mar = c(1, 1, 1, 0))
plot(metrics[[1]]
  , box    = F
  , legend = F
  , axes   = F
  , col = viridis(100)
)
plot(
    graph
  , layout             = lay
  , add                = T
  , rescale            = F
  , vertex.label       = NA
  , vertex.color       = "transparent"
  , vertex.frame.color = NA
  , edge.color         = "white"
  , edge.arrow.size    = 0.7
)
text(metrics[[1]], col = colTrans("white", percent = 50))
dev.off()

# Visualize Degree
png("test3.png", width = 1080, heigh = 1080, pointsize = 30)
par(mar = c(1, 1, 1, 0))
plot(metrics[[2]]
  , box    = F
  , legend = F
  , axes   = F
  , col = viridis(100)
)
plot(
    graph
  , layout             = lay
  , add                = T
  , rescale            = F
  , vertex.label       = NA
  , vertex.color       = "transparent"
  , vertex.frame.color = NA
  , edge.color         = "white"
  , edge.arrow.size    = 0.7
)
text(metrics[[2]], col = colTrans("white", percent = 50))
dev.off()
