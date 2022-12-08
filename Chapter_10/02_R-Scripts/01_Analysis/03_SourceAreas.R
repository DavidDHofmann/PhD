################################################################################
#### Source Areas
################################################################################
# Description: Generation of source area polygons

# Clear R's brain
rm(list = ls())

# Change the working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_10")

# Load required packages
library(raster)         # To handle spatial data
library(terra)          # To handle spatial data

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load some landmarks according to which we will delineate the source areas
wat <- vect("03_Data/02_CleanData/MajorWaters.shp")
wat <- wat[c(1, 2, 4), ]
roa <- vect("03_Data/02_CleanData/Roads.shp")
riv <- vect("03_Data/02_CleanData/MajorRivers.shp")
riv <- riv[-2, ]

# Buffer slightly
wat_buff <- buffer(wat, width = 500)
wat_buff <- makeValid(wat_buff)

# Specify centroids for the different source areas
src <- data.frame(
    x  = c(22.7, 23.75, 23.4, 22.7, 21.8)
  , y  = c(-20.1, -19.4, -18.8, -18.3, -19.2)
  , ID = 1:5
)
src <- vect(as.matrix(src[, c("x", "y")]), crs = crs(wat), atts = src)
src <- buffer(src, width = 20000)

# Let's also generate target zones. One target zone is the center of the delta,
# the remaining targets are placed on an ellipse around the delta. We'll make
# use of an ellipse to generate said areas
cen <- cbind(22.8, -19.26)
ell <- ellipse(cen[1], cen[2], 1.4, 1.25, phi = 45, spatial = T)
crs(ell) <- crs(wat)

# Generate a set of lines to dissect the ellipse
# Center of ellipse
x_end <- 22.8 + sin(seq(0, 2 * pi, by = pi / 4)) * 2
y_end <- -19.26 + cos(seq(0, 2 * pi, by = pi / 4)) * 2
xy_end <- cbind(x_end, y_end)
lns <- lapply(1:nrow(xy_end), function(x) {
  vect(rbind(cen, xy_end[x, ]), type = "line")
})
lns <- do.call(rbind, lns)
crs(lns) <- crs(wat)

# Buffer the ellipse very slightly
out <- buffer(ell, width = 5000)
ell <- out - ell
crs(ell) <- crs(wat)

# Now dissect the ellipse using the lines
bor <- ell - buffer(lns, width = 0.1)
bor <- disagg(bor)
bor$ID <- c(14, 13, 12, 7, 8, 9, 10, 11)

# Finally, we want to get a polygon for chiefs island, as a target in the center
# of the delta
chi      <- out - wat_buff
chi      <- disagg(chi)
chi$Area <- expanse(chi)
chi      <- chi[chi$Area %in% sort(chi$Area, decreasing = T)[2], ]
chi$ID   <- 6

# Put everything together
areas <- rbind(src, chi, bor)
areas$x    <- NULL
areas$y    <- NULL
areas$Area <- NULL

# Assign IDs to the different areas
areas$Type <- c(rep("Main", 6), rep("Buffer", nrow(areas) - 6))

# Sort by ID
areas <- areas[order(areas$ID), ]

# Visualize all
plot(wat, col = "cornflowerblue", border = NA)
plot(riv, col = "cornflowerblue", border = NA, add = T)
plot(roa, col = "gray", add = T)
plot(lns, add = T, lty = 2)
plot(areas, col = adjustcolor(sample(rainbow(length(areas))), alpha.f = 0.2), add = T)
text(areas, "ID", pos = 1)
text(areas, "Type", pos = 3)

# Store them
writeVector(areas, "03_Data/02_CleanData/SourceAreas.shp", overwrite = T)

# Let's also store the cutlines
writeVector(lns, "03_Data/02_CleanData/Cutlines.shp", overwrite = T)

################################################################################
#### Session Information
################################################################################
# Store session information
session <- devtools::session_info()
readr::write_rds(session, file = "02_R-Scripts/99_SessionInformation/03_SourceAreas.rds")

# Print to terminal
cat("Done :) \n")
