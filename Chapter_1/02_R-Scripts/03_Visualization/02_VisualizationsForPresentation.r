############################################################
#### Preparation of all Plots for the Thesis
############################################################
# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
wd <- "C:/Users/david/switchdrive/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)
library(raster)
library(rgeos)
library(smoothr)
library(tmap)
library(Cairo)
library(igraph)
library(rgdal)
library(davidoff)
library(glmmTMB)
library(cowplot)
library(parallel)
library(animation)
library(viridis)
library(velox)
library(imager)       # To import images

############################################################
#### Random Network
############################################################
# Create random network
net <- erdos.renyi.game(n = 12, 20, type = "gnm")
plot(net, vertex.label = NA)

############################################################
#### Plot of Historic Range
############################################################
# Load map of africa
africa    <- shapefile("03_Data/02_CleanData/00_General_Africa")
africa2   <- shapefile("03_Data/02_CleanData/00_General_Africa")
historic  <- shapefile("03_Data/01_RawData/DAVID/HistoricRange")
dogs      <- shapefile("03_Data/02_CleanData/00_General_WildDogs_IUCN")
kaza      <- shapefile("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")

# Buffer slightly
africa <- gBuffer(africa, width = 0.01)

# Disaggregate
africa <- disaggregate(africa)

# Identify sizes of areas
africa$Size <- gArea(africa, byid = T)

# Keep only the largest two
africa <- subset(africa, Size %in% sort(africa$Size, decreasing = T)[1:2])

# Simplify and smoothen africa shape
africa <- gSimplify(africa, tol = 0.5)
africa <- smooth(africa, method = "ksmooth")

# Smoothen and crop historic range
historic <- smooth(historic, method = "ksmooth")
historic <- gIntersection(historic, africa, byid = F)

# Simplify and smoothen dog distribution
dogs <- gSimplify(dogs, tol = 0.2)
dogs <- smooth(dogs, method = "ksmooth")

# Create a buffered polygon of africa
africa2 <- gBuffer(africa, width = 100/111000)

# Simplify and smoothen kaza
kaza <- gSimplify(kaza, tol = 0.1)
kaza <- smooth(kaza, method = "ksmooth")
plot(dogs)

# Identify wild dog strongholds
strong <- rbind(disaggregate(dogs)[c(7, 19, 29), ])

# Prepare a map of Africa
p1 <- tm_shape(africa2) +
    tm_polygons(col = "gray70", border.col = "gray70", lwd = 2) +
  tm_shape(africa) +
    tm_polygons(
        col = "black"
      , lwd = 0.7
      , border.col = "gray70"
    ) +
  tm_shape(historic) +
    tm_polygons(
        col = "orange"
      , lwd = 0.7
      , border.col = "black"
      , alpha = 0.2
    ) +
  tm_shape(dogs) +
    tm_polygons(
        col           = "orange"
      , alpha         = 0.8
      , border.alpha  = 0
    ) +
  tm_shape(strong) +
    tm_polygons(
        col           = lighten("orange", 1.3)
      , border.alpha  = 0
    ) +
  tm_layout(
      asp         = 0.8
    , frame       = "white"
    , frame.lwd   = 3
    , legend.show = FALSE
    , bg.color    = "white"
)

# Prepare a map of Africa where the kaza is added
p2 <- tm_shape(africa) +
    tm_polygons(col = "gray70", border.col = "gray70", lwd = 2) +
  tm_shape(africa) +
    tm_polygons(
        col = "black"
      , lwd = 0.7
      , border.col = "gray70"
    ) +
  tm_shape(historic) +
    tm_polygons(
        col = "orange"
      , lwd = 0.7
      , border.col = "black"
      , alpha = 0.2
    ) +
  tm_shape(dogs) +
    tm_polygons(
        col           = "orange"
      , alpha         = 0.8
      , border.alpha  = 0
    ) +
  tm_shape(kaza) +
    tm_borders(
        col = "white"
      , lwd = 3
    ) +
  tm_layout(
      asp         = 0.8
    , frame       = "white"
    , frame.lwd   = 3
    , legend.show = FALSE
    , bg.color    = "white"
)

# Store the plot
CairoPDF("Test", width = 5.25, height = 6)
p2
dev.off()

############################################################
#### Plot of KAZA-TFCA (Large)
############################################################
# Load required data
africa2 <- "03_Data/02_CleanData/00_General_Africa.shp" %>% readOGR()
kaza <- "03_Data/02_CleanData/00_General_KAZA_KAZA.shp" %>% readOGR()
prot <- "03_Data/02_CleanData/02_LandUseTypes_Protected_PeaceParks.shp" %>% readOGR()
water <- "03_Data/02_CleanData/03_LandscapeFeatures_MajorWaters_GEOFABRIK.shp" %>% readOGR(.)

# Clip africa layer
africa <- gIntersection(africa2, africa, byid = T)

# Get the extent of the KAZA
kaza_ext <- as(extent(kaza), "SpatialPolygons")
crs(kaza_ext) <- CRS("+init=epsg:4326")

# Create labels for countries
labels_countries <- data.frame(
    x = c(16, 24, 11, 38, 34)
  , y = c(-9, -30, -26, -13, -23)
  , Label = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")
)
coordinates(labels_countries) <- c("x", "y")
crs(labels_countries) <- CRS("+init=epsg:4326")

# Create lines pointing towards these countries
l1 <- rbind(c(17, -10), c(19, -12)) %>% spLines()
l2 <- rbind(c(24, -29), c(24, -23)) %>% spLines()
l3 <- rbind(c(13, -25), c(17, -22)) %>% spLines()
l4 <- rbind(c(35, -13), c(30, -14)) %>% spLines()
l5 <- rbind(c(33, -22), c(30, -19)) %>% spLines()
lines_countries <- rbind(l1, l2, l3, l4, l5)
crs(lines_countries) <- crs(labels_countries)

# Prepare a plot of the KAZA. Create labels for countries first.
labels_countries2 <- data.frame(
    x = c(20.39, 23.94, 20.07, 25.69, 28.22)
  , y = c(-15.28, -19.94, -19.39, -15.22, -18.9)
  , Label = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")
)
coordinates(labels_countries2) <- c("x", "y")
crs(labels_countries2) <- CRS("+init=epsg:4326")

# Create labels for some geographical landmarks
labels_waters <- data.frame(
    x     = c(22.6, 23.7, 27.8, 25.6)
  , y     = c(-19, -18.2, -17, -20.7)
  , Label = c(
    "Okavango\nDelta", "Linyanti\nSwamp", "Lake\nKariba", "Makgadikgadi\nPans"
  )
)
coordinates(labels_waters) <- c("x", "y")
crs(labels_waters) <- CRS("+init=epsg:4326")

# Create labels for some national parks
labels_nationalparks <- data.frame(
    x = c(26.56, 28.61, 21.15, 25.87, 20.38, 23.58, 23.71, 24.51, 20.78)
  , y = c(-19.08, -17.05, -17.26, -14.66, -16.08, -21.4, -19.29, -18.65, -18.81)
  , Label = paste0(c(
      "Hwange", "Matusadona", "Luengue-Luiana", "Kafue", "Mavinga"
    , "Central Kalahari", "Moremi", "Chobe", "Khaudum"
  ), "\nNP")
)
coordinates(labels_nationalparks) <- c("x", "y")
crs(labels_nationalparks) <- CRS("+init=epsg:4326")

# Prepare a map of Africa
p1 <- tm_shape(africa) +
    tm_polygons(
        col = "black"
      , lwd = 0.7
      , border.col = "gray30"
    ) +
  tm_shape(kaza) +
    tm_polygons(
        col = "orange"
      , alpha = 0.5
      , border.col = "orange"
    ) +
  tm_shape(kaza_ext) +
    tm_borders(
        col = "white"
      , lty = 3
      , lwd = 1.5
  )

# Store plot
CairoPDF("Test")
p1
dev.off()

# Plot of kaza only
p2 <- tm_shape(prot) +
    tm_polygons(
      , col           = "black"
      , border.col    = NA
      , lwd           = 0
      , legend.show   = F
    ) +
  tm_shape(kaza, is.master = T) +
    tm_polygons(
        col = "orange"
      , alpha = 0.5
      , border.col = "orange"
    ) +
  tm_shape(africa) +
    tm_borders(
        col = "gray70"
    ) +
  tm_layout(bg.color = "black")

# Plot of kaza and protected areas
p3 <- tm_shape(prot) +
    tm_polygons(
      , col           = "gray40"
      , border.col    = NA
      , lwd           = 0
      , legend.show   = F
    ) +
  tm_shape(kaza, is.master = T) +
    tm_polygons(
        col = "orange"
      , alpha = 0.5
      , border.col = "orange"
    ) +
  tm_shape(africa) +
    tm_borders(
        col = "gray70"
    ) +
  tm_layout(bg.color = "black")

CairoPDF("Test1")
p2
dev.off()

CairoPDF("Test2")
p3
dev.off()

############################################################
#### Plot Dispersal Durations
############################################################
library(lubridate)

# Load Cutoff Dates
cut <- read_csv("03_Data/02_CleanData/00_General_Dispersers_Popecol_CutoffDates.csv")

# Add Two hours
cut$StartDate <- cut$StartDate + hours(2)
cut$EndDate <- cut$EndDate + hours(2)

# Calculate dispersal duration
durations <- cut %>%
  mutate(DispersalDuration = EndDate - StartDate) %>%
  group_by(DogName) %>%
  summarize(DispersalDuration = sum(DispersalDuration))

# Visualize
plot(density(as.numeric(durations$DispersalDuration)), xlab = "Days")
library(ggpubr)
ggdensity(wdata, x = "weight",
          fill = "#0073C2FF", color = "#0073C2FF",
          add = "mean", rug = TRUE)

################################################################################
#### Movement Model Results
################################################################################
# Old directory
wd <- "/media/david/My Passport/Backups/WildDogs/15. PhD/00_WildDogs"
setwd(wd)

# Load the movement model
move.mod <- "03_Data/03_Results/99_MovementModel.rds" %>% readRDS()

# Select best model
best <- move.mod$Model[[1]]

# Summary of the best model
summary(best)

# Calculate p-values
coeffs <- getCoeffs(best, pvalue = TRUE)[-1, ] %>%

  # Calculate confidence intervals
  mutate(
      LCI = Coefficient - 1.96 * SE
    , UCI = Coefficient + 1.96 * SE
  )

# Add stars indicating the significance
coeffs$Significance <- sapply(1:nrow(coeffs), function(x){
  if (coeffs$pvalue[x] <= 0.01){
    return("***")
  } else if (coeffs$pvalue[x] <= 0.05){
    return("**")
  } else if (coeffs$pvalue[x] <= 0.1){
    return("*")
  }
})

# Rename covariates
coeffs$Covariate[coeffs$Covariate == "Shrubs"] <- "Shrubs/Grassland"
coeffs$Covariate[coeffs$Covariate == "Trees"] <- "Trees"
coeffs$Covariate[coeffs$Covariate == "Water"] <- "Water"
coeffs$Covariate[coeffs$Covariate == "DistanceToWater"] <- "DistanceToWater"
coeffs$Covariate[coeffs$Covariate == "cos_ta_"] <- "cos(ta)"
coeffs$Covariate[coeffs$Covariate == "log_sl_"] <- "log(sl)"
coeffs$Covariate[coeffs$Covariate == "HumansBuff5000"] <- "HumanInfluence"
coeffs$Covariate[coeffs$Covariate == "log_sl_:ActivityMainActivity"] <- "log(sl):MainActivity"
coeffs$Covariate[coeffs$Covariate == "log_sl_:Water"] <- "log(sl):Water"
coeffs$Covariate[coeffs$Covariate == "log_sl_:Trees"] <- "log(sl):Trees"
coeffs$Covariate[coeffs$Covariate == "cos_ta_:DistanceToWater"] <- "cos(ta):DistanceToWater"
coeffs$Covariate[coeffs$Covariate == "cos_ta_:HumansBuff5000"] <- "cos(ta):HumanInfluence"
coeffs$Preference <- ifelse(coeffs$Coefficient > 0, "Preferred", "Avoided")
coeffs$Preference <- factor(coeffs$Preference, levels = c("Preferred", "Avoided"))

# Specify the order in which the coefficients should be plotted
order <- c(
      "Water"
    , "DistanceToWater"
    , "Trees"
    , "Shrubs/Grassland"
    , "HumanInfluence"
    , "cos(ta)"
    , "cos(ta):HumanInfluence"
    , "cos(ta):DistanceToWater"
    , "log(sl)"
    , "log(sl):MainActivity"
    , "log(sl):Water"
    , "log(sl):Trees"
)

# Specify colors of axis labels
labcols <- c("black", "orange")[c(2, 2, 2, 1, 2, 1, 2, 1, 1, 1, 2, 2)]

# Prepare plot with Covariates on the y-axis and the corresponding
# coefficients on the x-axis
p1 <- ggplot(data = coeffs, aes(y = Covariate, x = Coefficient, col = factor(Preference))) +
  geom_point(shape = 1, size = 2) +
  geom_errorbarh(aes(
      xmin = Coefficient - 1.96 * SE
    , xmax = Coefficient + 1.96 * SE
    , height = 0.2)
  ) +
  geom_text(aes(label = Significance, hjust = 0.5, vjust = -0.2), show.legend = F) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
  scale_y_discrete(limits = rev(order)) +
  theme_cowplot() +
  xlim(c(-2, 2)) +
  coord_capped_cart(left = "both", bottom = "both") +
  labs(x = expression(beta*"-Coefficient")) +
  scale_color_manual(values = c("orange", "black")) +
  theme(
      panel.grid.minor = element_line(size = 0.0)
    , panel.grid.major = element_line(size = 0.0)
    , panel.border     = element_blank()
    , axis.line        = element_line()
    , axis.ticks       = element_line(colour = "black")
    , legend.title     = element_blank()
    # , axis.text.y      = element_text(color = rev(labcols))
  )

# Load a picture of a wild dog to the plot
dog <- load.image("/home/david/ownCloud/University/15. PhD/General/Images/WildDog_Running.svg")

# Add it to the previous plot
p <- ggdraw() +
  draw_image(dog, x = 0.9, y = 0.7, hjust = 0.5, vjust = 0.5, scale = 0.2) +
  draw_plot(p1)

# Store the plot
CairoPDF("test.pdf", width = 8, height = 5, bg = "transparent")
p
dev.off()

############################################################
#### Movement Model Random Effects
############################################################
# Rename model
mod <- best

# Check out coefficients per individual
coeffs <- coef(mod)
ranefs <- ranef(mod, condVar = T)
coeffs$cond$id
ranefs$cond$id

# Note: ranef yields the difference between the individual specific effect and
# the mean level effect. coef, on the other hand, yields the individual specific
# effect. Thus, the followin two lines yield (approximately) the same
mean(coeffs$cond$id$cos_ta_) + ranefs$cond$id$cos_ta_[3]

# We now want to visualize the individual variation. There are two possibilities
# for this: lme4::dotplot() or a ggplot. The dotplot is easier, yet not
# customizable. Let's first do the dotplot, then recreate it in ggplot.
lme4:::dotplot.ranef.mer(ranef(mod)$cond)

# Maybe scalefree?
lme4:::dotplot.ranef.mer(ranef(mod)$cond, scales = list(x = list(relation = "free")))

# Prepare dataframe that we need to plot the same in ggplot
rfs <- ranefs$cond$id %>%
  rownames_to_column() %>%
  gather(key = Covariate, value = Mean, 2:8)
names(rfs)[1] <- "id"

# We need to add the conditional variance
condVar <- attributes(ranefs$cond$id)$condVar
names(condVar) <- attributes(ranefs$cond$id)$names
condVar <- as.data.frame(do.call(rbind, condVar))
names(condVar) <- attributes(ranefs$cond$id)$row.names
condVar <- rownames_to_column(condVar)
names(condVar)[1] <- "Covariate"
condVar <- gather(condVar, key = id, value = Variance, 2:17)

# Join data to rfs dataframe
rfs <- left_join(rfs, condVar)

# Rename stuff nicely
rfs$Covariate <- gsub(rfs$Covariate, pattern = "cos_ta_", replacement = "cos(ta)")
rfs$Covariate <- gsub(rfs$Covariate, pattern = "log_sl_", replacement = "log(sl)")
rfs$Covariate <- gsub(rfs$Covariate, pattern = "sl_", replacement = "sl")
rfs$Covariate <- gsub(rfs$Covariate, pattern = "HumansBuff5000", replacement = "HumanInfluence")

# Make covariates a factor
rfs$Covariate <- factor(rfs$Covariate, levels = c(
  "cos(ta)", "sl", "log(sl)", "Water", "DistanceToWater", "Shrubs"
  , "Trees", "HumanInfluence"
))

# Visualize. Note that I am transforming the variance using mean - 2 *
# sqrt(Variance). This was taken from here: https://stackoverflow.com/questions
# /13847936/plot-random-effects-from-lmer-lme4-package-using-qqmath-or-dotplot-
# how-to-mak
p <- ggplot(rfs, aes(x = Mean, y = id)) +
  geom_point() +
  facet_wrap("Covariate", nrow = 2) +
  geom_errorbarh(aes(
      xmin = Mean - 2 * sqrt(Variance)
    , xmax = Mean + 2 * sqrt(Variance)
  ), colour = "black", height = 0) +

################################################################################
#### Simulated Trajectories
################################################################################
# Load required data
source_areas  <- shapefile("03_Data/03_Results/99_SourceAreas2.shp")
source_points <- shapefile("03_Data/03_Results/99_SourcePoints2.shp")

# Load some trajectories
sims <- read_rds("03_Data/03_Results/99_DispersalSimulationSub.rds")

# Only keep those from the static source point
sims <- subset(sims, PointSampling == "Static")

# Create animation
ani.options(interval = .1, ani.width = 1980, ani.height = 1980)
saveVideo({
  for (i in 1:2000){
    trajs <- sims2tracks(sims, steps = i)
    plot(source_areas, col = "gray50", border = "gray20", bg = "black")
    plot(trajs, add = T, col = "orange", lwd = 2.5)
  }
}, video.name = "99_Simulations99.mp4")

################################################################################
#### Two Network Views
################################################################################
# Load required data
source_areas  <- shapefile("03_Data/03_Results/99_SourceAreas1.shp")
source_points <- shapefile("03_Data/03_Results/99_SourcePoints1.shp")
kaza <- shapefile("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")

# Option 1: Source point view
############################################################
#### WORK HERE
############################################################

# Option 2: Raster view
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")
r <- aggregate(r, fact = 50000 / 250, fun = "max")

# Coerce raster to network
layout  <- as.matrix(as.data.frame(r, xy = T)[, c(1, 2)])
graph <- make_empty_graph(n = length(vertices))

# Coerce raster to shapefile
r <- as(r, "SpatialPolygons")

# Visualize
plot(r, border = "gray50")
plot(kaza, add = T)
plot(graph, layout = layout, vertex.size = 12, vertex.label = NA, add = T, rescale = F)

############################################################
#### How to rasterize Stuff
############################################################
# Function to rasterize using velox
rasterizeVelox <- function(l, r){
  summed <- raster(r)
  summed <- setValues(r, 0)
  for (i in 1:length(l)){
    line <- l[i, ]
    line$value <- 1
    vx <- velox(r)
    vx$rasterize(line, field = "value", background = 0)
    vx <- vx$as.RasterLayer()
    summed <- calc(stack(summed, vx), sum)
    cat(i, "out of", length(l), "done...\n")
  }
  return(summed)
}

# Function to create random lines
randomLine <- function(nodes = 2){
  x <- runif(nodes)
  y <- runif(nodes)
  line <- spLines(cbind(x, y))
}

# Create a couple of lines
lines <- lapply(rep(10, 5), randomLine) %>% do.call(rbind, .)

# Add some info
lines$ID <- 1:length(lines)

# Plot the lines
plot(lines, col = viridis(length(lines)))

# Create an empty raster onto which we can rasterize
r <- raster(ncol = 20, nrow = 20, xmn = 0, xmx = 1, ymn = 0, ymx = 1, vals = 0)

# Rasterize lines
lines_r <- rasterizeVelox(lines, r)


# Visualize
plot(lines_r, col = "black", box = F, axes = F, legend = F)
plot(lines, col = viridis(length(lines)), add = T)

plot(lines_r, col = magma(20), box = F, axes = F, legend = F)
plot(lines, add = T, lwd = 0.4, col = "white")
text(lines_r)
