############################################################
#### Some Descriptive Statistics
############################################################
# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_0"
setwd(wd)

# Load required packages
library(raster)
library(rgdal)
library(RColorBrewer)
library(viridis)
library(adehabitatLT)
library(animation)
library(rosm)
library(lubridate)
library(tidyverse)
library(rgeos)
library(spatstat)
library(maptools)
library(GISTools)
library(adehabitatHR)
library(NLMR)
library(gdistance)

############################################################
#### Loading GPS Data
############################################################
# Load the original gps fixes
data <- "03_Data/02_CleanData/00_General_Dispersers_Popecol(Regular).csv" %>%
  read_csv()

# Maybe we want to plot the gps fixes. So make them spatial
gps_dots <- data
coordinates(gps_dots) <- c("x", "y")
crs(gps_dots) <- CRS("+init=epsg:4326")

# Let's also load and reproject the trajectories
gps_traj <- "03_Data/02_CleanData/00_General_Dispersers_Popecol(Regular).shp" %>%
  readOGR() %>%
  spTransform(CRS("+init=epsg:4326"))

############################################################
#### Preparing Bounding Boxes
############################################################
# Define a bounding box for the reduced extent
box1 <- extent(c(22, 27, -20.7, -17.7))

# Define a bounding box for which we calculated dynamic floodmaps
box2 <- extent(c(21.74909, 24.30089, -20.64901, -18.14901))

# Prepare a third bounding box for the extent of Maun
box3 <- extent(c(23.2, 24.1, -20.3, -19.8))

# Prepare bounding box for which we will do the social analysis
box4 <- extent(c(23.25, 24, -19.8, -19.05)) %>% as(., "SpatialPolygons")

############################################################
#### General Statistics
############################################################
# Because some individuals were not included in the analysis, I will remove them
data <- subset(data, !(DogName %in% c("Appalachia", "Pula", "Rattan", "Ripley")))

# Identify the first and last date
range(data$Timestamp)

# Identify the total number of dogs and the number of dispersing individuals
length(unique(data$DogName))
length(unique(data$DogName[data$State == "Disperser"]))

# Identify the total number of fixes and the number of fixes during dispersal
nrow(data)
nrow(data[data$State == "Disperser", ])

# Identify the number of fixes for each dog
noFixes <- data %>%
  group_by(DogName, State) %>%
  summarize(NoFixes = n()) %>%
  spread(State, NoFixes) %>%
  arrange(Disperser)

# Look at the table
noFixes

# Get the some summary statistics about the number of fixes and distance travelled
summary(noFixes$Disperser)
summary(noFixes$Resident)

# Calculate standard deviation
sd(noFixes$Disperser, na.rm = TRUE)
sd(noFixes$Resident, na.rm = TRUE)

# Prepare table with Coalition, #Days dispersing, #Fixes, DIstanceTravelled
info <- data %>%
  group_by(DogName) %>%
  nest() %>%
  mutate(NoFixesTotal = map(data, function(x){
    nrow(x)
  }) %>% do.call(rbind, .)) %>%
  mutate(NoFixesDispersal = map(data, function(x){
    subset(x, State == "Disperser") %>% nrow()
  }) %>% do.call(rbind, .)) %>%
  mutate(DaysDispersing = map(data, function(x){
    subset(x, State == "Disperser") %>%
    mutate(Date = as.Date(Timestamp)) %>%
    select(Date) %>%
    unique() %>%
    nrow()
  }) %>% do.call(rbind, .)) %>%
  subset(NoFixesDispersal > 0)

# Remove data
info$data <- NULL

# Look at the results
info

############################################################
#### Identify Dispersal Distances
############################################################
# Subset to dispersers, convert projection
disp <- gps_dots %>%
  subset(State == "Disperser") %>%
  spTransform(CRS("+init=epsg:32734")) %>%
  as.data.frame()

# Coerce to ltraj
disp <- as.ltraj(
    xy    = disp[, c("x", "y")]
  , date  = disp$Timestamp
  , id    = disp$DogName
)

# Calculate total distances per individual
distances <- lapply(disp, function(x){
  sum(x$dist, na.rm = T) / 1000
}) %>% do.call(rbind, .)

# Put back the dispersers' names
names <- lapply(disp, function(x){
  attr(x, "id")
}) %>% do.call(rbind, .)

# Put all together
distances <- data.frame(DogName = names, DistanceTravelled = distances)

# We are also intrested in straight line distances. For this, we only keep the
# first and last relocation of each individual
disp <- gps_dots %>%
  subset(State == "Disperser") %>%
  spTransform(CRS("+init=epsg:32734")) %>%
  as.data.frame() %>%
  group_by(DogName) %>%
  nest() %>%
  mutate(data = map(data, function(x){
    rbind(head(x, 1), tail(x, 1))
  })) %>%
  unnest()

# Coerce to ltraj
disp <- as.ltraj(
    xy    = disp[, c("x", "y")]
  , date  = disp$Timestamp
  , id    = disp$DogName
)

# Calculate total distances per individual
distances2 <- lapply(disp, function(x){
  sum(x$dist, na.rm = T) / 1000
}) %>% do.call(rbind, .)

# Put back the dispersers' names
names <- lapply(disp, function(x){
  attr(x, "id")
}) %>% do.call(rbind, .)

# Put all together
distances2 <- data.frame(DogName = names, StraightLineDistanceTravelled = distances2)

# Put all distances into one table
distances <- left_join(distances, distances2, by = "DogName")

# Join this data to our info table
info <- left_join(info, distances, by = "DogName")

# Look at the result
info

# Get some summaries
summary(info$DaysDispersing)
summary(info$DistanceTravelled)
summary(info$StraightLineDistanceTravelled)

mean(info$NoFixesTotal)
sd(info$NoFixesTotal)
mean(info$NoFixesDispersal)
sd(info$NoFixesDispersal)
mean(info$DaysDispersing)
sd(info$DaysDispersing)
mean(info$StraightLineDistanceTravelled)
sd(info$StraightLineDistanceTravelled)
mean(info$DistanceTravelled)
sd(info$DistanceTravelled)

############################################################
#### Identify Pack IDs
############################################################
# We also want to assign pack ids prior to dispersal of dispersers to our info
# table
packID <- subset(data, DogName %in% unique(data$DogName[data$State == "Disperser"])) %>%
  group_by(DogName) %>%
  nest() %>%
  mutate(data = map(data, function(x){
    subset(x, Timestamp < min(Timestamp[State == "Disperser"])) %>%
    tail(1)
  })) %>%
  unnest() %>%
  select(DogName, CurrentPack)

# Put that into the table
info <- left_join(info, packID, by = "DogName")

# Manually fill in the blanks
info$CurrentPack[info$DogName == "Stetson"] <- "MT"

# Look at the result
info

############################################################
#### Identify Sex
############################################################
# Load required data
sex <- "03_Data/01_RawData/POPECOL/CollarSettings.csv" %>%
  read_csv() %>%
  select(c("Dog Name", "Sex")) %>%
  setNames(c("DogName", "Sex")) %>%
  unique() %>%
  na.omit()

# Join this data to the info table
info <- left_join(info, sex, by = "DogName")

# Manually assign sex to Abrahms' dispersers
info$Sex[info$DogName %in% c("Lupe", "Scorpion", "Stetson")] <- "M"

# Round values
write_csv(info, "03_Data/03_Results/99_GPSInfo.csv")

############################################################
#### Plot of Study Area
############################################################
# Load the data to plot
africa  <- "03_Data/02_CleanData/00_General_Africa.shp" %>%
  shapefile()
kaza    <- "03_Data/02_CleanData/00_General_KAZA_KAZA.shp" %>%
  shapefile()
delta   <- "03_Data/02_CleanData/03_LandscapeFeatures_MajorWaters_GEOFABRIK" %>%
  shapefile()
prot    <- "03_Data/02_CleanData/02_LandUseTypes_Protected_PeaceParks(3Classes)" %>%
  shapefile()

# Plot the map of Africa with the study extent.
png("99_Africa_ESRI.png"
  , width   = 1980
  , height  = 1080
  , bg      = "transparent"
)
plot(africa
  , col     = viridis(1, begin = 0.2)
  , border  = "white"
)
plot(kaza
  , add     = TRUE
  , col     = viridis(1, begin = 0.6)
  , border  = viridis(1, begin = 0.6)
  , lwd     = 1
)
plot(africa
  , border  = "white"
  , add     = TRUE
)
plot(extent(kaza)
  , add = TRUE
  , col = viridis(1, begin = 0.6)
  , lwd = 7
)
scalebar(
    type  = "line"
  , below = "Km"
  , d     = 500
  , lwd   = 10
  , cex   = 5
  , col   = "white"
)
dev.off()

# Plot the KAZA
png("99_KAZA_KAZA.png"
  , width   = 1980
  , height  = 1080
  , bg      = "transparent"
)
plot(kaza
  , col = viridis(1, begin = 0.6)
)
plot(delta
  , add     = TRUE
  , col     = "cadetblue1"
  , border  = "cadetblue1"
)
plot(africa
  , border  = "white"
  , add     = TRUE
  , lwd     = 4
)
scalebar(
    type  = "line"
  , below = "Km"
  , d     = 150
  , lwd   = 10
  , cex   = 5
  , col   = "white"
)
dev.off()

# Plot the KAZA with National Parks
png("99_KAZA_KAZA(NationalParks).png"
  , width   = 1980
  , height  = 1080
  , bg      = "transparent"
)
plot(kaza
  , col = viridis(1, begin = 0.6)
)
plot(delta
  , add     = TRUE
  , col     = "cadetblue1"
  , border  = "cadetblue1"
)
plot(africa
  , border  = "white"
  , add     = TRUE
  , lwd     = 4
)
plot(crop(subset(prot, Desig == "National Park"), kaza)
  , col     = adjustcolor("green", alpha.f = 0.2)
  , border  = "white"
  , add     = TRUE
  , lwd     = 3
)
scalebar(
    type  = "line"
  , below = "Km"
  , d     = 150
  , lwd   = 10
  , cex   = 5
  , col   = "white"
)
dev.off()

############################################################
#### Plotting Trajectories
############################################################
# Plot all trajectories onto a bing map.
png("99_AllTrajectories(Indivs).png", width = 1920, height = 1080)

# Get the bing map
bmaps.plot(box1
  , type  = "Aerial"
  , key   = "Ai0A3oHx1c3rivWQfR2nAmSle5AVyo6RxjtzUFLkxTT2_qbxyNDSge0A5jxxnDgW"
  , stoponlargerequest = FALSE
)

# Add the trajectories onto the satellite map
gps_traj %>%

  # Aggregate them by individual
  aggregate(., by = "id") %>%

  # Transform to the crs of the bmap
  spTransform(., CRS("+init=epsg:3857")) %>%

  # Plot the tracks
  plot(.
    , col = rev(viridis(nrow(.)))
    , add = TRUE
    , lwd = 3
  )

# Store the plot
dev.off()

# We also want to plot all trajectories but with different colours for
# dispersers and residents
png("99_AllTrajectories(State).png", width = 1920, height = 1080)

# Get the bing map
bmaps.plot(bbox(gps_traj)
  , type  = "Aerial"
  , key   = "Ai0A3oHx1c3rivWQfR2nAmSle5AVyo6RxjtzUFLkxTT2_qbxyNDSge0A5jxxnDgW"
  , stoponlargerequest = FALSE
)

# Add the trajectories onto the satellite map
gps_traj %>%

  # Aggregate them by individual
  aggregate(., by = "state") %>%

  # Transform to the crs of the bmap
  spTransform(., CRS("+init=epsg:3857")) %>%

  # Plot the tracks
  plot(.
    , col = rev(viridis(nrow(.), begin = 0.3, end = 1))
    , add = TRUE
    , lwd = 3
  )

# Store the plot
dev.off()

############################################################
#### Water
############################################################
# Load data
water <- "03_Data/02_CleanData/01_LandCover_Water(Merged).tif" %>%
  stack()

# Calculate average water coverage
water <- calc(water, sum)

# Prepare Plot
png("99_Covariates_Water.png"
  , height = 1080
  , width  = 1080
)
plot(water
  , col    = viridis(20)
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
  , ext    = box2
)
dev.off()

############################################################
#### Trees
############################################################
# Load Data
trees <- "03_Data/02_CleanData/01_LandCover_TreeCover_MODIS.tif" %>%
  raster()

# Prepare Plot
png("99_Covariates_Trees.png"
  , height = 1080
  , width  = 1080
)
plot(log(trees)
  , col    = viridis(20)
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
  , ext    = box2
)
dev.off()

############################################################
#### Shrubs
############################################################
# Load Data
shrubs <- "03_Data/02_CleanData/01_LandCover_NonTreeVegetation_MODIS.tif" %>%
  raster()

# Prepare Plot
png("99_Covariates_Shrubs.png"
  , height = 1080
  , width  = 1080
)
plot(log(shrubs)
  , col    = viridis(20)
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
  , ext    = box2
)
dev.off()

############################################################
#### Protected
############################################################
# Load Data
protected <- "03_Data/02_CleanData/02_LandUseTypes_Protected_PeaceParks(1Class).tif" %>%
  shapefile()

# Prepare Plot
png("99_Covariates_Protected.png"
  , height = 1080
  , width  = 1080
)
plot(crop(protected, box2)
  , col = viridis(20)[15]
  , bg = viridis(20)[1]
)
dev.off()

############################################################
#### Humans
############################################################
# Load Data
humans <- "03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluence(Buffer5000).tif" %>%
  raster()

# Prepare Plot
png("99_Covariates_Humans.png"
  , height = 1080
  , width = 1080
)
plot(humans
  , col    = viridis(20)
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
  , ext    = box2
)
dev.off()

############################################################
#### Humans Detailed
############################################################
# Load the desired layers
farms1 <- "03_Data/02_CleanData/04_AnthropogenicFeatures_Agriculture_Globelands.tif" %>%
  raster() %>%
  crop(box1)
farms2 <- "03_Data/02_CleanData/04_AnthropogenicFeatures_Agriculture_Croplands.tif" %>%
  raster() %>%
  crop(box1)
farms3 <- "03_Data/02_CleanData/04_AnthropogenicFeatures_Farms_Gabriele.shp" %>%
  shapefile() %>%
  crop(box1)
roads <- "03_Data/02_CleanData/04_AnthropogenicFeatures_Roads_GEOFABRIK.shp" %>%
  shapefile()

# Replace 0s with NAs to allow stacking of multiple rasters
values(farms1)[values(farms1) == 0] <- NA
values(farms2)[values(farms2) == 0] <- NA

# Prepare a plot for the farms
png("99_Covariates_Humans(Farms).png"
  , height = 1080
  , width  = 1080
)
plot(farms1
  , col     = viridis(20)[20]
  , colNA   = viridis(20)[1]
  , legend  = FALSE
  , axes    = FALSE
  , ext     = box3
)
plot(farms2
  , col     = viridis(20)[20]
  , colNA   = viridis(20)[1]
  , legend  = FALSE
  , axes    = FALSE
  , add     = TRUE
  , ext     = box3
)
plot(farms3
  , col     = viridis(20)[20]
  , border  = viridis(20)[20]
  , add     = TRUE
)
dev.off()

# Prepare an empty background
bg <- farms1
values(bg) <- 0

# Prepare a plot for the roads
png("99_Covariates_Humans(Roads).png"
  , height = 1080
  , width  = 1080
)
plot(bg
  , col     = viridis(20)[1]
  , legend  = FALSE
  , axes    = FALSE
  , ext     = box3
)
plot(roads
  , col = "white"
  , add = TRUE
)
dev.off()

# Prepare plot for the merged product
png("99_Covariates_Humans(HumanDensity).png"
  , height = 1080
  , width  = 1080
)
plot(humans
  , col    = viridis(50)
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
  , ext    = box3
)
dev.off()

############################################################
#### Roads
############################################################
# Load Data
roads1 <- "03_Data/02_CleanData/04_AnthropogenicFeatures_DistanceToRoads.tif" %>%
  raster()
roads2 <- "03_Data/02_CleanData/04_AnthropogenicFeatures_Roads_GEOFABRIK.shp" %>%
  shapefile()

# Prepare a plot for the distance to roads map
png("99_Covariates_Roads.png"
  , height = 1080
  , width  = 1080
)
plot(roads1
  , col    = viridis(50)
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
  , ext    = box2
)
plot(crop(roads2, box2)
  , col = "white"
  , lwd = 1
  , add = TRUE
)
dev.off()

############################################################
#### Social Landscape
############################################################
# Load Data of social landscape
social <- "03_Data/02_CleanData/05_SocialFeatures_SocialLandscape(UtilisationDistributions).grd" %>%
  stack()

# Plot a couple of the landscapes to file
for (i in 30:39){
  png(paste0("99_Covariates_Social_No", i, ".png")
    , height = 1080
    , width  = 1080
  )
  plot(social[[i]]
    , col    = viridis(50)
    , legend = FALSE
    , axes   = FALSE
    , box    = FALSE
    , ext    = box4
  )
  dev.off()
}

# Put all together for better visibility
social <- calc(social, sum)

# Prepare a plot for one ud
png("99_Covariates_Social.png"
  , height = 1080
  , width  = 1080
)
plot(social
  , col    = viridis(50)
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
  , ext    = box2
)
plot(box4
  , border = "white"
  , add    = TRUE
  , lwd    = 4
)
dev.off()

# Prepare the plot for the reduced extent for better visibility
png("99_Covariates_Social(ReducedExtent).png"
  , height = 1080
  , width  = 1080
)
plot(social_sum
  , col    = viridis(50)
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
  , ext    = box4
)
plot(gps
  , col = "white"
  , pch = "."
  , add = TRUE
)
dev.off()

# Example of how we calculated social landscapes
set.seed(1234)

# Create example data
x1 <- rnorm(1:1000, mean = 0, sd = 0.3)
y1 <- rnorm(1:1000, mean = 0, sd = 0.1)
id1 <- rep("Pack1", 1000)

x2 <- rnorm(1:1000, mean = 1, sd = 0.2)
y2 <- rnorm(1:1000, mean = 1, sd = 0.4)
id2 <- rep("Pack2", 1000)

# Put everything into a dataframe
dat <- data.frame(
    x = c(x1, x2)
  , y = c(y1, y2)
  , id = c(id1, id2)
)

# Make the data spatial
coordinates(dat) <- c("x", "y")

# Plot the data
plot(dat)

# Calculate kernelUD
ud <- kernelUD(dat, h = "href", grid = 100, same4all = TRUE)

# Convert uds to rasters
ud1 <- raster(ud[[1]])
ud2 <- raster(ud[[2]])

# Plot the utilisation distributions
png("99_Covariates_Social(UD1).png"
  , height = 1080
  , width  = 1080
)
plot(ud1
  , col    = viridis(50)
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
)
plot(dat[1:1000, ]
  , col = "white"
  , pch = "."
  , add = TRUE
  , cex = 0.5
)
dev.off()
png("99_Covariates_Social(UD2).png"
  , height = 1080
  , width  = 1080
)
plot(ud2
  , col    = viridis(50)
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
)
plot(dat[1001:2000, ]
  , col = "white"
  , pch = "."
  , add = TRUE
  , cex = 0.5
)
dev.off()

# Calculate probability of not encountering pack i
values(ud1) <- 1 - values(ud1) / sum(values(ud1))
values(ud2) <- 1 - values(ud2) / sum(values(ud2))

# Calculate probability of encounter at all
ud <- 1 - ud1 * ud2

# Plot the result
png("99_Covariates_Social(MergedUDs).png"
  , height = 1080
  , width  = 1080
)
plot(ud
  , col    = viridis(50)
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
)
dev.off()

############################################################
#### The ORI Algorithm
############################################################
# Some illustrations that will help to explain the ORI algorithm. Let's load the
# masks
wet_mask <- "03_Data/01_RawData/MODIS/MCD43A4/02_Masks/MaskWater.shp" %>%
  shapefile()
dry_mask <- "03_Data/01_RawData/MODIS/MCD43A4/02_Masks/MaskDryland.shp" %>%
  shapefile()

png("99_ORI_Algorithm1.png"
  , height = 1080
  , width  = 1080
)
bmaps.plot(box2
  , type  = "Aerial"
  , key   = "Ai0A3oHx1c3rivWQfR2nAmSle5AVyo6RxjtzUFLkxTT2_qbxyNDSge0A5jxxnDgW"
  , stoponlargerequest = FALSE
)
wet_mask %>%
  spTransform(., CRS("+init=epsg:3857")) %>%
  plot(., add = TRUE, col = viridis(1, begin = 0.2))
dry_mask %>%
  spTransform(., CRS("+init=epsg:3857")) %>%
  plot(., add = TRUE, col = viridis(1, begin = 0.8))
dev.off()

# Create a histogram that illustrates the idea of dryland vs. wetland
x1 <- rpois(1000, lambda = 5)
x2 <- rpois(1000, lambda = 15)
x <- data.frame(water = x1, dryland = x2)
x <- gather(x, "Category", "Value", 1:2)
png("99_ORI_Algorithm3.png"
  , height = 1080
  , width  = 1080
)
ggplot(x, aes(Value, fill = Category)) +
  geom_density(alpha = 0.9) +
  scale_fill_viridis(discrete = T, begin = 0.2, end = 0.8, direction = -1) +
  theme_bw()
dev.off()

############################################################
#### Random Maps to Illustrate Resistance Mapping
############################################################
set.seed(12345)

# Specify raster size (ncol = nrow)
n <- 50

# Simulate Water
water <- nlm_gaussianfield(ncol = n, nrow = n, autocorr_range = 20)
water <- round(water)
plot(water, col = viridis(20, end = 0.8))

# Simulate Human density
humans <- nlm_gaussianfield(ncol = n, nrow = n)
plot(humans, col = viridis(20))

# Simulate Protected areas
protected <- nlm_mosaictess(ncol = n, nrow = n, germs = 5)
plot(protected, col = viridis(2, end = 0.8))

# Simulate a resistance map
res1 <- nlm_distancegradient(ncol = n, nrow = n, origin = c(1, 1, 1, 1))
res2 <- nlm_distancegradient(ncol = n, nrow = n, origin = c(n, n, n, n))
res3 <- nlm_gaussianfield(ncol = n, nrow = n)
res <- res1 + res2 + 0.3 * res3
plot(res, col = rev(viridis(20)))

# Plot and store all the random maps
png("99_RandMap_Water.png"
  , height = 1080
  , width  = 1080
)
plot(water
  , col    = viridis(2, end = 0.8)
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
)
dev.off()

png("99_RandMap_Humans.png"
  , height = 1080
  , width  = 1080
)
plot(humans
  , col    = viridis(5)
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
)
dev.off()

png("99_RandMap_Protected.png"
  , height = 1080
  , width  = 1080
)
plot(protected
  , col    = viridis(2, end = 0.8)
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
)
dev.off()

png("99_RandMap_Resistance.png"
  , height = 1080
  , width  = 1080
)
plot(res
  , col    = rev(viridis(5, end = 0.8))
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
)
dev.off()

############################################################
#### Show How Least Cost Paths and Corridors Work
############################################################
# Specify size of raster
n <- 100

# Simulate Gaussian Field
map1 <- nlm_gaussianfield(ncol = n, nrow = n)

# Simulate an area that is avoided in the center
map2 <- nlm_distancegradient(ncol = n, nrow = n, origin = c(n/2, n/2, n/2, n/2))

# Put the two maps together
map <- map1 + map2

# Rescale them to 0 and 1
map <- map / max(values(map))

# Plot the map
plot(map, col = viridis(50))

# Let's use this as our "premeability" map, where larger values indicate easier
# permeability. We now convert the raster to a transition layer
trans <- transition(map, transitionFunction = mean, directions = 16)

# Apply correction
trans <- geoCorrection(trans, type = "c", multpl = FALSE)

# Create points to connect
A <- data.frame(x = 5, y = 50) %>% SpatialPoints()
B <- data.frame(x = 95, y = 50) %>% SpatialPoints()
points <- rbind(A, B)

# Plot map and points together
png("99_LeastCostPath(ResistanceMap).png"
  , height = 1080
  , width  = 1080
)
plot(map
  , col    = viridis(50)
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
)
plot(points, add = TRUE, col = "white", pch = 19, cex = 5)
dev.off()

# Find least cost path visually
path <- shortestPath(trans
  , origin  = points[1, ]
  , goal    = points[2, ]
  , output  = "SpatialLines"
)

# Plot the least cost path
png("99_LeastCostPath(ShortestPath).png"
  , height = 1080
  , width  = 1080
)
plot(map
  , col    = viridis(50)
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
)
plot(points, add = TRUE, col = "white", pch = 19, cex = 5)
plot(path, add = TRUE, col = "white", lwd = 5)
dev.off()

# We can also calculate the cost map from the two locations
A_cost <- accCost(trans, points[1, ])
B_cost <- accCost(trans, points[2, ])

# Plot the costs from point A
png("99_LeastCostPath(CostA).png"
  , height = 1080
  , width  = 1080
)
plot(A_cost
  , col    = rev(viridis(50))
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
)
plot(points[1, ], add = TRUE, col = "white", pch = 19, cex = 5)
dev.off()

# Plot the costs from point B
png("99_LeastCostPath(CostB).png"
  , height = 1080
  , width  = 1080
)
plot(B_cost
  , col    = rev(viridis(50))
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
)
plot(points[2, ], add = TRUE, col = "white", pch = 19, cex = 5)
dev.off()

# Put the maps together
cost <- A_cost + B_cost

# Plot the cost map
png("99_LeastCostPath(CostAB).png"
  , height = 1080
  , width  = 1080
)
plot(cost
  , col    = rev(viridis(50))
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
)
plot(points, add = TRUE, col = "white", pch = 19, cex = 5)
dev.off()

# Identify a desired percentile
quant <- quantile(cost, probs = 0.1, na.rm = TRUE)

# Reclassify values above the quantile to get
values(cost) <- ifelse(values(cost) > quant, NA, values(cost))

# Plot the resulting corridor
png("99_LeastCostPath(CostABPercentile).png"
  , height = 1080
  , width  = 1080
  , bg = "transparent"
)
plot(cost
  , col    = rev(viridis(50))
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
)
plot(points, add = TRUE, col = "white", pch = 19, cex = 5)
dev.off()

# Compare Least Cost Path and Least Cost Corridor
png("99_LeastCostPath(CostPathCorridor).png"
  , height = 1080
  , width  = 1080
  , bg = "transparent"
)
plot(cost
  , col    = rev(viridis(50))
  , legend = FALSE
  , axes   = FALSE
  , box    = FALSE
)
plot(path, add = TRUE, col = "white", lwd = 5)
plot(points, add = TRUE, col = "white", pch = 19, cex = 5)
dev.off()
