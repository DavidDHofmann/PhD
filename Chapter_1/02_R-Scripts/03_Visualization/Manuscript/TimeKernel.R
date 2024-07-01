################################################################################
#### Temporal Dispersal Kernel
################################################################################
# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)    # For data wrangling
library(raster)       # To handle spatial data
library(terra)        # To handle spatial data
library(rgdal)        # To read shapefiles
library(davidoff)     # Custom functions
library(pbmcapply)    # For running stuff in parallel
library(viridis)      # For nice colors
library(adehabitatHR) # For calculating HRs
library(sf)           # To plot spatial stuff with ggplot
library(tmaptools)    # To download satellite imagery
library(RStoolbox)    # To plot RGB rasters
library(ggpubr)       # To arrange plots
library(ggspatial)    # To plot north arrows & scale bars

################################################################################
#### Prepare Data
################################################################################
# Load simulated dispersal data
sims <- read_rds("03_Data/03_Results/99_DispersalSimulation.rds")
# sims <- read_rds("03_Data/03_Results/99_DispersalSimulationSub.rds")

# Load protected areas
prot <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")

# Create SpatialPoints for first location of each trajectory
first <- sims %>%
  dplyr::select("x", "y", "TrackID") %>%
  group_by(TrackID) %>%
  slice(1) %>%
  SpatialPointsDataFrame(
      coords      = cbind(.[["x"]], .[["y"]])
    , proj4string = CRS("+init=epsg:4326")
  )

# Assess from which protected area each trajectory left
first$Origin <- as.character(over(first, prot)$Name)
first <- first@data[, c("TrackID", "Origin")]

# Join information to simulated tracks
sims <- left_join(sims, first, by = "TrackID")
sims$Origin[sims$Area == "Buffer"] <- "Buffer"

# Remove tracks that leave from an unknown origin (should only be few)
sims <- sims[!is.na(sims$Origin), ]

# We now rasterize trajectories leaving from desired source locations after
# different numbers of steps. In contrast to creating heatmaps, we only want
# binary maps indicating through which areas the dispersers have moved. This can
# efficiently be achieved using the terra::rasterize function. Let's thus setup
# the study design through which we will loop.
design <- as_tibble(
  expand.grid(
      StepNumber       = seq(2, 2000, length.out = 20)
    , Origin           = c("Moremi", "Hwange")
    , stringsAsFactors = F
  )
)

# Nest by origin (so that we don't need to duplicate data all the time)
design <- design %>% nest(Steps = -Origin)

# Loop through the design, rasterize trajectories and asses "95% home range" to
# indicate areas visited after x number of steps
design$Kernel <- lapply(1:nrow(design), function(x){

  # For simpler indexing, extract the different number of steps
  steps <- as.vector(as.matrix(design$Steps[[x]]))

  # Subset simulations to desired source location
  sims_sub <- sims[sims$Origin == design$Origin[x], ]

  # Define extent for which we will rasterize. We don't need to consider the
  # full extent here
  ext <- ext(c(
      min(sims_sub$x)
    , max(sims_sub$x)
    , min(sims_sub$y)
    , max(sims_sub$y)
  )) + c(-1, +1, -1, +1) * metersToDegrees(1000)

  # Prepare reference raster with a desired resolution
  r <- rast(ext, res = metersToDegrees(1000))

  # Loop through different number of steps and rasterize them, then compute home
  # ranges around the visited pixels
  maps <- pbmclapply(steps, mc.cores = 1, ignore.interactive = T, function(y){

    # Subset simulations to desired number of steps
    sub <- sims_sub[sims_sub$StepNumber <= y, ]

    # Create trajectories and convert to terra::vector
    sub_track <- sims2tracks(sub, messages = F)
    sub_track <- vect(sub_track)

    # Rasterize trajectories onto reference raster
    rasterized <- terra::rasterize(sub_track, r, background = 0, touches = T)

    # Use visited pixels to create a utilization distribution and compute home
    # ranges
    pts <- rasterToPoints(raster(rasterized), spatial = T, fun = function(x){x > 0})
    pts <- as(pts, "SpatialPoints")
    suppressWarnings(ud <- kernelUD(pts))
    suppressWarnings(hr <- getverticeshr(ud, percent = 95))

    # Return the home range
    return(hr)

  }) %>% do.call(rbind, .)

  # Indicate the number of steps corresponding to each homerange
  maps$Steps <- steps

  # Return the homeranges
  return(maps)

})

################################################################################
#### Visualizations
################################################################################
# Load additional shapefiles
kaza    <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")
africa  <- readOGR("03_Data/02_CleanData/00_General_Africa_ESRI.shp")

# Create labels for countries
labels_countries <- data.frame(
    x = c(20, 23, 20, 26, 26.3)
  , y = c(-17, -21, -19, -17, -18.2)
  , Label = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")
)
coordinates(labels_countries) <- c("x", "y")
crs(labels_countries) <- CRS("+init=epsg:4326")

# Convert objects to sf for plotting with ggplot
kernel_moremi <- st_as_sf(design$Kernel[[1]])
kernel_hwange <- st_as_sf(design$Kernel[[2]])
prot          <- st_as_sf(prot)
kaza          <- st_as_sf(kaza)
africa        <- st_as_sf(africa)
labels_countries <- st_as_sf(labels_countries)

# Invert order
kernel_moremi <- kernel_moremi[rev(1:nrow(kernel_moremi)), ]
kernel_hwange <- kernel_hwange[rev(1:nrow(kernel_hwange)), ]
kernel_moremi$NationalPark <- "Moremi"
kernel_hwange$NationalPark <- "Hwange"
kernels <- rbind(kernel_moremi, kernel_hwange)
kernels$NationalPark <- as.factor(as.character(kernels$NationalPark))

# Download a satellite images for a background map
sat1 <- read_osm(kernel_moremi, type = "bing")
sat1 <- as(sat1, "Raster")
sat1 <- projectRaster(sat1, crs = CRS("+init=epsg:4326"), method = "ngb")
sat2 <- read_osm(kernel_hwange, type = "bing")
sat2 <- as(sat2, "Raster")
sat2 <- projectRaster(sat2, crs = CRS("+init=epsg:4326"), method = "ngb")

# Crop the image
sat <- crop(sat, extent(sat) + c(0, 0, 0, -4))

# Plot background for first kernel
p1 <- ggplot() +
  ggRGB(sat1, r = 1, g = 2, b = 3, ggLayer = T) +
  geom_sf(
      data = subset(prot, Name %in% c("Moremi", "Hwange"))
    , lwd  = 0.2
    , col  = "white"
    , fill = colTrans("white", percent = 80)
  ) +
  geom_sf(
      data = subset(prot, Name %in% c("Moremi"))
    , lwd  = 0.2
    , col  = "red"
    , fill = colTrans("red", percent = 80)
  ) +
  geom_sf(
      data = africa
    , lwd  = 0.5
    , col  = "white"
    , fill = NA
    , lty  = 2
  ) +
  geom_sf_text(
      data    = labels_countries
    , mapping = aes(label = Label)
    , col     = "white"
    , size    = 4
  ) +
  coord_sf(
      crs  = 4326
    , xlim = extent(sat1)[1:2]
    , ylim = extent(sat1)[3:4]
    , expand = F
  ) +
  theme_minimal() +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 0.5
    , height     = unit(0.15, "cm")
    , text_col   = "white"
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(1.2, "cm"),
    , width    = unit(1, "cm"),
    , style    = north_arrow_fancy_orienteering(
          fill      = c("white", "white")
        , line_col  = NA
        , text_col  = "white"
        , text_size = 12
      )
  )

# Plot first kernel
p2 <- ggplot() +
  geom_sf(
      data    = kernel_moremi
    , mapping = aes(fill = Steps)
    , col     = "black"
    , lwd     = 0.1
  ) +
  scale_fill_viridis_c(
      direction = -1
    , guide   = guide_colorbar(
      , title          = "Number of Steps"
      , show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.6, "cm")
      , barwidth       = unit(2.0, "cm")
    )
  ) +
  geom_sf(
      data = subset(prot, Name %in% c("Moremi", "Hwange"))
    , lwd  = 0.2
    , col  = "white"
    , fill = colTrans("white", percent = 80)
  ) +
  geom_sf(
      data = subset(prot, Name == "Moremi")
    , col  = "red"
    , fill = colTrans("red", 60)
    , lwd  = 0.4
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = extent(sat1)[1:2]
    , ylim   = extent(sat1)[3:4]
    , expand = F
  ) +
  theme_minimal() +
  theme(
      legend.position  = "none"
    , legend.box       = "vertical"
    , panel.grid.major = element_line(color = colTrans("gray70", 80), size = 0.2)
    , panel.grid.minor = element_blank()
  ) +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 0.5
    , height     = unit(0.15, "cm")
    , text_col   = "black"
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(1.2, "cm"),
    , width    = unit(1, "cm"),
    , style    = north_arrow_fancy_orienteering(
          fill      = c("black", "black")
        , line_col  = NA
        , text_col  = "black"
        , text_size = 12
      )
  )

# Plot background for second kernel
p3 <- ggplot() +
  ggRGB(sat2, r = 1, g = 2, b = 3, ggLayer = T) +
  geom_sf(
      data = subset(prot, Name %in% c("Moremi", "Hwange"))
    , lwd  = 0.2
    , col  = "white"
    , fill = colTrans("white", percent = 80)
  ) +
  geom_sf(
      data = subset(prot, Name %in% c("Hwange"))
    , lwd  = 0.2
    , col  = "red"
    , fill = colTrans("red", percent = 80)
  ) +
  geom_sf(
      data = africa
    , lwd  = 0.5
    , col  = "white"
    , fill = NA
    , lty  = 2
  ) +
  geom_sf_text(
      data    = labels_countries
    , mapping = aes(label = Label)
    , col     = "white"
    , size    = 4
  ) +

  coord_sf(
      crs  = 4326
    , xlim = extent(sat2)[1:2]
    , ylim = extent(sat2)[3:4]
    , expand = F
  ) +
  theme_minimal() +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 0.5
    , height     = unit(0.15, "cm")
    , text_col   = "white"
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(1.2, "cm"),
    , width    = unit(1, "cm"),
    , style    = north_arrow_fancy_orienteering(
          fill      = c("white", "white")
        , line_col  = NA
        , text_col  = "white"
        , text_size = 12
      )
  )

# Plot second kernel
p4 <- ggplot() +
  geom_sf(
      data    = kernel_hwange
    , mapping = aes(fill = Steps)
    , col     = "black"
    , lwd     = 0.1
  ) +
  scale_fill_viridis_c(
      direction = -1
    , guide   = guide_colorbar(
      , title          = "Number of Steps"
      , show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.6, "cm")
      , barwidth       = unit(2.0, "cm")
    )
  ) +
  geom_sf(
      data = subset(prot, Name %in% c("Moremi", "Hwange"))
    , lwd  = 0.2
    , col  = "white"
    , fill = colTrans("white", percent = 80)
  ) +
  geom_sf(
      data = subset(prot, Name == "Hwange")
    , col  = "red"
    , fill = colTrans("red", 60)
    , lwd  = 0.4
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = extent(sat2)[1: 2]
    , ylim   = extent(sat2)[3: 4]
    , expand = F
  ) +
  theme_minimal() +
  theme(
      legend.position  = "none"
    , legend.box       = "vertical"
    , panel.grid.major = element_line(color = colTrans("gray70", 80), size = 0.2)
    , panel.grid.minor = element_blank()
  ) +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 0.5
    , height     = unit(0.15, "cm")
    , text_col   = "black"
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(1.2, "cm"),
    , width    = unit(1, "cm"),
    , style    = north_arrow_fancy_orienteering(
          fill      = c("black", "black")
        , line_col  = NA
        , text_col  = "black"
        , text_size = 12
      )
  )

# Put plots together
p <- ggarrange(p1, p2, p3, p4
  , labels = c("a1", "a2", "b1", "b2")
)

# Prepare a plot from which we can grab the legend
legend <- ggplot() +
  geom_sf(
      data    = kernel_hwange
    , mapping = aes(fill = Steps)
    , col     = "black"
    , lwd     = 0.1
  ) +
  scale_fill_viridis_c(
      direction = -1
    , guide   = guide_colorbar(
      , title          = "Number of Steps"
      , show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.6, "cm")
      , barwidth       = unit(10.0, "cm")
    )
  ) +
  theme(
      legend.position  = "bottom"
    , legend.box       = "vertical"
    , panel.grid.major = element_line(color = colTrans("gray70", 80), size = 0.2)
    , panel.grid.minor = element_blank()
  )
legend <- get_legend(legend)

# Add legend
p <- ggarrange(p, legend, nrow = 2, heights = c(10, 1))

# Store
ggsave("04_Manuscript/99_TimeKernel.png", plot = p, scale = 1.3)
