################################################################################
#### Visualize Connections Between National Parks
################################################################################
# Description: In this script, we visualize the interpatch connectivity

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)      # For data wrangling
library(raster)         # To handle spatial data
library(igraph)         # For network analysis
library(rgdal)          # To load spatial data
library(sf)             # To plot spatial stuff with ggplot
library(ggspatial)      # To add scale bars etc to plots
library(rgeos)          # To calculate areas
library(ggnetwork)      # To plot network using ggplot
library(viridis)        # For nice colors
library(davidoff)       # Custom functions
library(cowplot)        # Grab ggplot legend
library(gtable)         # To combine ggplots
library(grid)           # To combine ggplots
library(gridExtra)      # To combine ggplots
library(ggpubr)         # To arrange multiple plots
library(ggdark)         # For dark ggplot themes

################################################################################
#### Prepare Data
################################################################################
# Load data on interpatch connectivity (essentially the from-to connections that
# each trajectory makes)
visits <- read_rds("03_Data/03_Results/99_InterpatchConnectivity.rds")
nsims_count <- read_rds("03_Data/03_Results/99_NumberSimulatedCountry.rds")
nsims_parks <- read_rds("03_Data/03_Results/99_NumberSimulatedPark.rds")

# Summarize connections by park. Compute average number of steps required to
# make connection, it's standard deviation, as well as how often the connection
# has been made.
visits_parks <- visits %>%
  group_by(From, To) %>%
  summarize(
      MeanStepNumber = mean(StepNumber)
    , SDStepNumber   = sd(StepNumber)
    , Frequency      = n()
    , .groups        = "drop"
  )

# Summarize connections by country. Compute average number of steps required to
# make connection, it's standard deviation, as well as how often the connection
# has been made.
visits_count <- visits %>%
  group_by(FromCountry, ToCountry) %>%
  summarize(
      MeanStepNumber = mean(StepNumber)
    , SDStepNumber   = sd(StepNumber)
    , Frequency      = length(unique(TrackID))
    , .groups        = "drop"
  )

# How many tracks leaving from each country reached another country?
reached_countries <- subset(visits, FromCountry != ToCountry) %>%
  dplyr::select(TrackID, FromCountry) %>%
  unique() %>%
  group_by(FromCountry) %>%
  summarize(SuccessfullSimulations = n()) %>%
  left_join(nsims_count, by = c("FromCountry")) %>%
  mutate(RelativeSuccess = SuccessfullSimulations / NumberSimulations)

# Use the number of initiated individuals to calculate the relative frequency
visits_parks$Simulations <- nsims_parks$NumberSimulations[match(visits_parks$From, nsims_parks$SourceArea)]
visits_count$Simulations <- nsims_count$NumberSimulations[match(visits_count$FromCountry, nsims_count$FromCountry)]
visits_parks$RelFrequency <- visits_parks$Frequency / visits_parks$Simulations
visits_count$RelFrequency <- visits_count$Frequency / visits_count$Simulations

# Show the results
visits_parks
visits_count

# Fill empty combinations (make implicit 0s explicit)
combs_parks <- expand_grid(
    From = unique(c(visits_parks$From, visits_parks$To))
  , To   = unique(c(visits_parks$From, visits_parks$To))
) %>% arrange(From, To)
combs_count <- expand_grid(
    FromCountry = unique(c(visits_count$FromCountry, visits_count$FromCountry))
  , ToCountry   = unique(c(visits_count$ToCountry, visits_count$ToCountry))
)

# Left_join
visits_parks <- left_join(combs_parks, visits_parks, by = c("From", "To"))
visits_count <- left_join(combs_count, visits_count, by = c("FromCountry", "ToCountry"))

# Replace NAs with 0s
visits_parks[is.na(visits_parks)] <- 0
visits_count[is.na(visits_count)] <- 0

# Replace the diagnoal (i.e. where from = to) with NAs
visits_parks[, 3:ncol(visits_parks)][visits_parks$From == visits_parks$To, ] <- NA
visits_count[, 3:ncol(visits_count)][visits_count$FromCountry == visits_count$ToCountry, ] <- NA

# Standard deviation of 0 does not make sense here. Replace it
visits_parks$SDStepNumber[visits_parks$SDStepNumber == 0] <- NA

# Show difference in successfull dispersers for some parks
visits_parks$FromPark <- visits$FromPark[match(visits_parks$From, visits$From)]
visits_parks$ToPark <- visits$ToPark[match(visits_parks$To, visits$To)]
subset(visits_parks, FromPark %in% c("Moremi", "Chobe") & ToPark %in% c("Moremi", "Chobe"))
visits_parks$FromPark <- NULL
visits_parks$ToPark <- NULL

################################################################################
#### Matrix Plots
################################################################################
# Visualize data on a matrix plot
p1 <- ggplot(visits_count, aes(x = ToCountry, y = FromCountry)) +
  geom_tile(aes(fill = Frequency), col = "gray20", lwd = 0.2) +
  geom_text(aes(label = round(Frequency, 2))) +
  coord_fixed() +
  scale_fill_distiller(palette = 5, na.value = "black") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev) +
  theme_minimal()
p2 <- ggplot(visits_count, aes(x = ToCountry, y = FromCountry)) +
  geom_tile(aes(fill = RelFrequency), col = "gray20", lwd = 0.2) +
  geom_text(aes(label = round(RelFrequency, 2))) +
  coord_fixed() +
  scale_fill_distiller(palette = 5, na.value = "black") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev) +
  theme_minimal()
p3 <- ggplot(visits_count, aes(x = ToCountry, y = FromCountry)) +
  geom_tile(aes(fill = MeanStepNumber), col = "gray20", lwd = 0.2) +
  geom_text(aes(label = round(MeanStepNumber, 2))) +
  coord_fixed() +
  scale_fill_distiller(palette = 5, na.value = "black") +
  scale_x_discrete(position = "top") +
  scale_y_discrete(limits = rev) +
  theme_minimal()
ggarrange(p1, p2, p3)

################################################################################
#### Prepare Additional Data for Network Plots
################################################################################
# Load additional shapefiles and raster for the background
kaza    <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")
africa  <- readOGR("03_Data/02_CleanData/00_General_Africa_ESRI.shp")
prot    <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")

# Simplify Protection zones
prot$Desig <- as.character(prot$Desig)
prot$Desig[prot$Desig == "Forest Reserve"] <- "Protected"

# Make nicer names
prot$Desig[prot$Desig == "National Park"] <- "National Parks"
prot$Desig[prot$Desig == "Protected"] <- "Protected Areas"
prot$Desig <- as.factor(as.character(prot$Desig))

# Get the extent of the KAZA
kaza_ext <- as(extent(kaza), "SpatialPolygons")
crs(kaza_ext) <- CRS("+init=epsg:4326")

# Rename
kaza$Name <- "KAZA-TFCA Borders"

# Prepare country labels
labels_countries <- data.frame(
    x = c(20.39, 23.94, 20.07, 25.99, 28.22)
  , y = c(-15.28, -21.80, -19.39, -14.52, -18.9)
  , Label = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")
)
coordinates(labels_countries) <- c("x", "y")
crs(labels_countries) <- CRS("+init=epsg:4326")

# Create labels for some national parks
labels_nationalparks <- data.frame(
    x = c(26.56, 28.61, 21.15, 25.87, 20.38, 23.58, 23.21, 24.51, 20.78, 22.63, 27.92, 28.54)
  , y = c(-19.08, -17.05, -17.26, -15.25, -16.08, -21.4, -19.29, -18.65, -18.81, -14.54, -17.76, -20.53)
  , Label = paste0(c(
      "Hwange", "Matusadona", "Luengue-Luiana", "Kafue", "Mavinga"
    , "Central Kalahari", "Moremi", "Chobe", "Khaudum", "Liuwa Plains"
    , "Chizarira", "Matobo"
  ), " NP")
)
coordinates(labels_nationalparks) <- c("x", "y")
crs(labels_nationalparks) <- CRS("+init=epsg:4326")

# Convert objects to sf
kaza                 <- st_as_sf(kaza)
africa               <- st_as_sf(africa)
prot                 <- st_as_sf(prot)
labels_countries     <- st_as_sf(labels_countries)
labels_nationalparks <- st_as_sf(labels_nationalparks)

# Convert reference raster to dataframe
r <- raster("03_Data/02_CleanData/00_General_Raster.tif")
r <- as.data.frame(r, xy = T)

################################################################################
#### Network Plot for National Parks
################################################################################
# Load protected areas
np <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")
np <- subset(np, Desig == "National Park")
np$ID <- 1:nrow(np)

# Remove self-loops and zero-links
visits <- subset(visits_parks, From != To & Frequency != 0)

# Create a network
net <- graph_from_data_frame(
    d        = visits
  , vertices = unique(np$ID)
  , directed = T
)

# Prepare layouts
lay <- coordinates(gCentroidWithin(np))

# Prepare networks for ggplotting with ggplot
net_p <- ggnetwork(net, layout = lay, arrow.gap = 0.1, scale = F)

# Make nice labels
net_p$Label <- paste0(net_p$FromName, " NP")

# Prepare color palette
pal <- colorRampPalette(plasma(100, begin = 1, end = 0.4))

# Plot of kaza
p1a <- ggplot() +
  geom_sf(
      data        = kaza
    , col         = "gray80"
    , fill        = NA
    , lty         = 1
    , lwd         = 1
    , show.legend = F
  ) +
  geom_sf(
      data        = africa
    , col         = "gray80"
    , fill        = NA
    , lty         = 2
    , lwd         = 0.5
    , show.legend = F
  ) +
  geom_sf_text(
      data     = labels_countries
    , mapping  = aes(label = Label)
    , col      = "white"
    , fontface = 2
    , size     = 5
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = c(min(r$x), max(r$x))
    , ylim   = c(min(r$y), max(r$y))
    , expand = F
  ) +
  labs(
      x        = NULL
    , y        = NULL
    , fill     = NULL
    # , title    = "Interpatch Connectivity"
    # , subtitle = "In Relation to Dispersal Duration"
  ) +
  guides(
      size  = guide_legend(title.position = "top", order = 2)
  ) +
  dark_theme_classic() +
  scale_fill_manual(values = c("gray30", "gray10")) +
  theme(
      legend.position       = "bottom"
    , legend.box            = "horizontal"
    , legend.title.align    = 0.5
    , panel.border          = element_rect(colour = "white", fill        = "transparent", size = 0.7)
    , axis.line             = element_blank()
    , legend.title          = element_text(size   = 10)
    , legend.text           = element_text(size   = 8)
    , legend.margin         = margin(c(0, 0, 0, 0))
    , legend.key            = element_blank()
    , plot.background       = element_rect(fill   = "transparent", color = NA)
    , panel.grid.major      = element_blank()
    , panel.grid.minor      = element_blank()
    , legend.background     = element_rect(fill   = "transparent")
    , legend.box.background = element_rect(fill   = "transparent", color = "transparent")
  ) +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 1
    , height     = unit(0.15, "cm")
    , bar_cols   = c("white", "white")
    , text_col   = "white"
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(1.5, "cm"),
    , width    = unit(1.2, "cm"),
    , style    = north_arrow_fancy_orienteering(
          fill      = c("white", "white")
        , line_col  = NA
        , text_col  = "white"
        , text_size = 12
      )
  )

# Plot of kaza and national parks
p1b <- ggplot() +
  geom_sf(
      data        = prot
    , mapping     = aes(fill = Desig)
    , col         = "black"
    , lwd         = 0.1
    , show.legend = F
  ) +
  geom_sf(
      data        = kaza
    , col         = "gray80"
    , fill        = NA
    , lty         = 1
    , lwd         = 1
    , show.legend = F
  ) +
  geom_sf(
      data        = africa
    , col         = "gray80"
    , fill        = NA
    , lty         = 2
    , lwd         = 0.5
    , show.legend = F
  ) +
  geom_sf_text(
      data     = labels_countries
    , mapping  = aes(label = Label)
    , col      = "white"
    , fontface = 2
    , size     = 5
  ) +
  geom_sf_text(
      data     = labels_nationalparks
    , mapping  = aes(label = Label)
    , nudge_y  = 0.3
    , fontface = 3
    , size     = 3
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = c(min(r$x), max(r$x))
    , ylim   = c(min(r$y), max(r$y))
    , expand = F
  ) +
  labs(
      x        = NULL
    , y        = NULL
    , fill     = NULL
    # , title    = "Interpatch Connectivity"
    # , subtitle = "In Relation to Dispersal Duration"
  ) +
  guides(
      size  = guide_legend(title.position = "top", order = 2)
  ) +
  dark_theme_classic() +
  scale_fill_manual(values = c("gray30", "gray10")) +
  theme(
      legend.position       = "bottom"
    , legend.box            = "horizontal"
    , legend.title.align    = 0.5
    , panel.border          = element_rect(colour = "white", fill        = "transparent", size = 0.7)
    , axis.line             = element_blank()
    , legend.title          = element_text(size   = 10)
    , legend.text           = element_text(size   = 8)
    , legend.margin         = margin(c(0, 0, 0, 0))
    , legend.key            = element_blank()
    , plot.background       = element_rect(fill   = "transparent", color = NA)
    , panel.grid.major      = element_blank()
    , panel.grid.minor      = element_blank()
    , legend.background     = element_rect(fill   = "transparent")
    , legend.box.background = element_rect(fill   = "transparent", color = "transparent")
  ) +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 1
    , height     = unit(0.15, "cm")
    , bar_cols   = c("white", "white")
    , text_col   = "white"
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(1.5, "cm"),
    , width    = unit(1.2, "cm"),
    , style    = north_arrow_fancy_orienteering(
          fill      = c("white", "white")
        , line_col  = NA
        , text_col  = "white"
        , text_size = 12
      )
  )

# Main Plot
p1 <- ggplot() +
  geom_sf(
      data        = prot
    , mapping     = aes(fill = Desig)
    , col         = "black"
    , lwd         = 0.1
    , show.legend = F
  ) +
  geom_sf(
      data        = kaza
    , col         = "gray80"
    , fill        = NA
    , lty         = 1
    , lwd         = 1
    , show.legend = F
  ) +
  geom_sf(
      data        = africa
    , col         = "gray80"
    , fill        = NA
    , lty         = 2
    , lwd         = 0.5
    , show.legend = F
  ) +
  geom_edges(
      data      = net_p
    , mapping   = aes(
        x    = x
      , y    = y
      , xend = xend
      , yend = yend
      , size = RelFrequency
      , col  = MeanStepNumber
    )
    , curvature = 0.2
    , arrow     = arrow(length = unit(6, "pt"), type = "closed", angle = 10)
  ) +
  geom_sf_text(
      data     = labels_countries
    , mapping  = aes(label = Label)
    , col      = "white"
    , fontface = 2
    , size     = 5
  ) +
  scale_size_area(
      name     = "Relative Frequency"
    , max_size = 1
  ) +
  scale_color_gradientn(
      colors  = pal(100)
    , guide   = guide_colorbar(
        title          = "Duration (Steps)"
      , show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.6, "cm")
      , barwidth       = unit(3.0, "cm")
      , order = 1
    )
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = c(min(r$x), max(r$x))
    , ylim   = c(min(r$y), max(r$y))
    , expand = F
  ) +
  labs(
      x        = NULL
    , y        = NULL
    , fill     = NULL
    # , title    = "Interpatch Connectivity"
    # , subtitle = "In Relation to Dispersal Duration"
  ) +
  guides(
      size  = guide_legend(title.position = "top", order = 2)
  ) +
  dark_theme_classic() +
  scale_fill_manual(values = c("gray30", "gray10")) +
  theme(
      legend.position       = "bottom"
    , legend.box            = "horizontal"
    , legend.title.align    = 0.5
    , panel.border          = element_rect(colour = "white", fill        = "transparent", size = 0.7)
    , axis.line             = element_blank()
    , legend.title          = element_text(size   = 10)
    , legend.text           = element_text(size   = 8)
    , legend.margin         = margin(c(0, 0, 0, 0))
    , legend.key            = element_blank()
    , plot.background       = element_rect(fill   = "transparent", color = NA)
    , panel.grid.major      = element_blank()
    , panel.grid.minor      = element_blank()
    , legend.background     = element_rect(fill   = "transparent")
    , legend.box.background = element_rect(fill   = "transparent", color = "transparent")
  ) +
  annotation_scale(
      location   = "bl"
    , width_hint = 0.2
    , line_width = 1
    , height     = unit(0.15, "cm")
    , bar_cols   = c("white", "white")
    , text_col   = "white"
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(1.5, "cm"),
    , width    = unit(1.2, "cm"),
    , style    = north_arrow_fancy_orienteering(
          fill      = c("white", "white")
        , line_col  = NA
        , text_col  = "white"
        , text_size = 12
      )
  )

# Plot for the separate legend of the dots
p2 <- ggplot() +
  geom_point(
      data    = net_p %>% dplyr::select(x, y, name, Simulations) %>% distinct()
    , mapping = aes(x = x, y = y)
    , col     = "black"
    , size    = 0.1
  ) +
  geom_point(
      data    = net_p %>% dplyr::select(x, y, name, Simulations) %>% distinct() %>% na.omit()
    , mapping = aes(x = x, y = y, size = Simulations, color = Simulations)
    , col     = "orange"
  ) +
  geom_sf_text(
      data     = labels_nationalparks
    , mapping  = aes(label = Label)
    , nudge_y  = 0.3
    , fontface = 3
    , size     = 3
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = c(min(r$x), max(r$x))
    , ylim   = c(min(r$y), max(r$y))
    , expand = F
  ) +
  labs(
      x        = NULL
    , y        = NULL
    , fill     = NULL
    # , title    = "Interpatch Connectivity"
    # , subtitle = "In Relation to Dispersal Duration"
  ) +
  guides(
    size = guide_legend(title = "Number of Simulations", title.position = "top")
  ) +
  dark_theme_classic() +
  theme(
      legend.position       = "bottom"
    , legend.box            = "horizontal"
    , legend.title.align    = 0.5
    , legend.margin         = margin(c(0, 0, 13, 0))
    , legend.title          = element_text(size = 10),
    , legend.text           = element_text(size = 8)
    , legend.key            = element_blank()
    , panel.background      = element_rect(fill = "transparent")
    , plot.background       = element_rect(fill = "transparent", color = NA)
    , panel.grid.major      = element_blank()
    , panel.grid.minor      = element_blank()
    , legend.background     = element_rect(fill = "transparent")
    , legend.box.background = element_rect(fill = "transparent", color = "transparent")
    , panel.border          = element_rect(colour = "white", fill        = "transparent", size = 0.7)
    , axis.line             = element_blank()
  )

# Extract legends from all plots
legend1 <- get_legend(p1)
legend2 <- get_legend(p2)

# Remove the original legend from main plot(s)
p1a <- p1a + theme(legend.position = "none")
p1b <- p1b + theme(legend.position = "none")
p3 <- p1 + theme(legend.position = "none")

# Add circles to main plot
g1 <- ggplotGrob(p3)
g2 <- ggplotGrob(p2)
g2 <- gtable_filter(g2, "panel")
pos <- c(subset(g1$layout, grepl("panel", g1$layout$name), select = t:r))
p4 <- gtable_add_grob(g1, g2, t = pos$t, l = pos$l)
p4 <- ggplotify::as.ggplot(p4)

# Put legends together
legends <- grid.arrange(legend1, legend2, nrow = 1, widths = c(2, 1))
legends <- gtable_add_padding(legends, unit(c(0, 1, 0, 0.3), "cm"))

# Put them below the first plot
p1a <- arrangeGrob(p1a, legends, heights = c(10, 1))
p1b <- arrangeGrob(p1b, legends, heights = c(10, 1))
p5 <- arrangeGrob(p4, legends, heights = c(10, 1))
p1a <- ggplotify::as.ggplot(p1a)
p1b <- ggplotify::as.ggplot(p1b)
p5 <- ggplotify::as.ggplot(p5)

# Visualize
p1a
p1b
p5

# Store
ggsave("InterpatchConnectivity1.png", plot = p1a)
ggsave("InterpatchConnectivity2.png", plot = p1b)
ggsave("InterpatchConnectivity3.png", plot = p5)
