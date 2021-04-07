################################################################################
#### Plot of Study Area
################################################################################
# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

################################################################################
#### Preparing Data
################################################################################
# Load required packages
library(rgdal)        # To load and store spatial data
library(raster)       # To manipulate raster data
library(rgeos)        # To manipulate spatial data
library(terra)        # To handle raster data
library(tidyverse)    # For data wrangling
library(davidoff)     # Custom functions
library(elevatr)      # To create nice hillshade map
library(sf)           # To plot spatial objects with ggplot
library(ggspatial)    # For north arrow and scale bar
library(cowplot)      # To grab legends and arrange multiple plots

# Load required data
africa  <- readOGR("03_Data/02_CleanData/00_General_Africa_ESRI.shp")
kaza    <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")
dogs    <- readOGR("03_Data/02_CleanData/00_General_WildDogs_IUCN.shp")
prot    <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")
water   <- rast("03_Data/02_CleanData/01_LandCover_WaterCoverAveraged_MERGED.tif")

# Aggregate for coarser resolution
water <- water %>%
  aggregate(fun = max, fact = 2) %>%
  raster()

# Identify isoltated pixels
clump <- clump(water)
clump_freq <- clump %>%
  freq() %>%
  data.frame() %>%
  subset(count <= 40)

# Remove them
rcl <- data.frame(old = clump_freq$value, new = NA)
water <- reclassify(clump, rcl = rcl)
water <- !is.na(water)

# Simplify Protection zones
prot$Desig[prot$Desig == "Forest Reserve"] <- "Protected"
prot$Desig <- as.factor(as.character(prot$Desig))

# Remove small islands for plotting
africa <- subset(africa, !(ID %in% c(27:41, 689, 690)))

# Get the extent of the KAZA
kaza_ext <- as(extent(kaza), "SpatialPolygons")
crs(kaza_ext) <- CRS("+init=epsg:4326")

# Create labels for countries
labels_countries <- data.frame(
    x = c(16, 24, 11, 44, 35)
  , y = c(-8, -30, -27, -14, -24)
  , Label = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")
)
coordinates(labels_countries) <- c("x", "y")
crs(labels_countries) <- CRS("+init=epsg:4326")

# Create lines pointing towards these countries
l1 <- rbind(c(17, -10), c(19, -12)) %>% spLines()
l2 <- rbind(c(24, -29), c(24, -23)) %>% spLines()
l3 <- rbind(c(13, -25), c(17, -22)) %>% spLines()
l4 <- rbind(c(37, -14), c(29, -15)) %>% spLines()
l5 <- rbind(c(33, -22), c(30, -19)) %>% spLines()
lines_countries <- rbind(l1, l2, l3, l4, l5)
crs(lines_countries) <- crs(labels_countries)

# Visualize
plot(labels_countries)
plot(lines_countries, add = T)
plot(lines_countries[4], add = T, col = "red")
plot(africa, add = T)

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
    x     = c(22.6, 23.7, 27.8, 25.6, 22.9, 27.8, 27.25)
  , y     = c(-19, -18.2, -17, -20.7, -15.0, -14.4, -15.58)
  , Label = c(
    "Okavango\nDelta", "Linyanti\nSwamp", "Lake\nKariba", "Makgadikgadi\nPans"
    , "Barotse\nFloodplain", "Lukanga\nSwamp", "Kafue\nFlats"
  )
)
coordinates(labels_waters) <- c("x", "y")
crs(labels_waters) <- CRS("+init=epsg:4326")

# Visualize
plot(water)
plot(labels_waters, add = T, col = "red")

# Create labels for some national parks
labels_nationalparks <- data.frame(
    x = c(26.56, 28.61, 21.15, 25.87, 20.38, 23.58, 23.71, 24.51, 20.78, 22.63, 27.92, 28.54)
  , y = c(-19.08, -17.05, -17.26, -14.66, -16.08, -21.4, -19.29, -18.65, -18.81, -14.54, -17.76, -20.53)
  , Label = paste0(c(
      "Hwange", "Matusadona", "Luengue-Luiana", "Kafue", "Mavinga"
    , "Central Kalahari", "Moremi", "Chobe", "Khaudum", "Liuwa Plains"
    , "Chizarira", "Matobo"
  ), "\nNP")
)
coordinates(labels_nationalparks) <- c("x", "y")
crs(labels_nationalparks) <- CRS("+init=epsg:4326")

# Visualize
plot(subset(prot, Desig == "National Park"))
plot(labels_nationalparks, add = T, col = "red")

# Prepare hillshade for Africa
terrain_africa <- africa %>%
  get_elev_raster(z = 4) %>%
  crop(africa) %>%
  mask(mask = africa, updatevalue = NA) %>%
  calc(., function(x){x ** 2}) %>%
  terrain(opt = c("slope", "aspect"))
hill_africa <- hillShade(
    aspect = terrain_africa$aspect
  , slope  = terrain_africa$slope
)

# Visualize it
plot(hill_africa, col = gray(80:100/100))
plot(africa, add = T, border = "gray80")

# Prepare hillshade for KAZA
ext <- as(extent(kaza) + c(-1, 1, -1, 1), "SpatialPolygons")
crs(ext) <- CRS("+init=epsg:4326")
terrain_kaza <- kaza %>%
  get_elev_raster(z = 6) %>%
  crop(ext) %>%
  mask(mask = ext, updatevalue = NA) %>%
  calc(., function(x){x ** 2}) %>%
  terrain(opt = c("slope", "aspect"))
hill_kaza <- hillShade(
    aspect = terrain_kaza$aspect
  , slope  = terrain_kaza$slope
)

# Visualize it
plot(hill_kaza, col = gray(80:100/100))
plot(kaza, add = T, border = "gray80")

# Remove undesired variables
rm(clump, clump_freq, ext, rcl, terrain_africa, terrain_kaza)
gc()

# Reduce resolution for now
water <- aggregate(water, fact = 5, fun = max)
hill_africa <- aggregate(hill_africa, fact = 5, fun = mean)
hill_kaza <- aggregate(hill_kaza, fact = 5, fun = mean)

# Convert everything to be compatible with ggplot
hill_africa           <- as.data.frame(hill_africa, xy = T)
hill_kaza             <- as.data.frame(hill_kaza, xy = T)
water                 <- as.data.frame(water, xy = T)
africa                <- st_as_sf(africa)
dogs                  <- st_as_sf(dogs)
kaza                  <- st_as_sf(kaza)
prot                  <- st_as_sf(prot)
lines_countries       <- st_as_sf(lines_countries)
labels_countries      <- st_as_sf(labels_countries)
labels_countries2     <- st_as_sf(labels_countries2)
labels_waters         <- st_as_sf(labels_waters)
labels_nationalparks  <- st_as_sf(labels_nationalparks)

################################################################################
#### Plot of Africa
################################################################################
# Plot of Africa
p1 <- ggplot() +
  geom_raster(
      data        = hill_africa
    , mapping     = aes(x = x, y = y, fill = layer)
    , show.legend = F
    , alpha       = 0.3
  ) +
  scale_fill_gradientn(
      colours  = gray(60:100/100)
    , na.value = "transparent"
  ) +
  geom_sf(
      data  = africa
    , fill  = NA
    , color = "gray80"
    , lwd   = 0.4
  ) +
  geom_sf(
      data  = dogs
    , fill  = "cornflowerblue"
    , color = NA
    , alpha = 0.8
  ) +
  geom_sf(
      data  = kaza
    , fill  = NA
    , color = "black"
    , lwd   = 1
  ) +
  geom_sf_text(
      data    = labels_countries
    , mapping = aes(label = Label)
  ) +
  geom_sf(
      data = lines_countries
    , col  = "black"
  ) +
  coord_sf(
    crs = 4326
  ) +
  labs(
      x = NULL
    , y = NULL
  ) +
  theme_void() +
  theme(
      panel.background = element_blank()
    , panel.border     = element_rect(colour = "black", fill = NA, size = 1)
    , plot.margin = unit(c(0, 0, 0, 0), "null")
  )

################################################################################
#### Plot of KAZA-TFCA
################################################################################
# Plot of KAZA
p2 <- ggplot() +
  geom_raster(
      data        = hill_kaza
    , mapping     = aes(x = x, y = y, alpha = layer)
    , show.legend = F
  ) +
  scale_alpha(
    range = c(0.1, 0)
  ) +
  geom_sf(
      data        = prot
    , aes(fill    = factor(Desig))
    , col         = "#6ba36b"
    , lwd         = 0.1
    , show.legend = F
  ) +
  scale_fill_manual(
    values = c("#70ab70", "#d9f0d3")
  ) +
  annotate(
      geom = "raster"
    , x    = water$x
    , y    = water$y
    , fill = scales::colour_ramp(c("transparent", "#96d0ff"))(water$layer)
  ) +
  geom_sf(
      data = africa
    , col  = "gray80"
    , fill = NA
    , lwd  = 0.2
    , lty  = 1
  ) +
  geom_sf(
      data = kaza
    , col  = "black"
    , fill = NA
    , lwd  = 1
    , lty  = 1
  ) +
  geom_sf_text(
      data     = labels_countries2
    , mapping  = aes(label = Label)
    , size     = 5
    , fontface = 2
  ) +
  geom_sf_text(
      data     = labels_waters
    , mapping  = aes(label = Label)
    , col      = "#064886"
    , size     = 2
    , fontface = 3
  ) +
  geom_sf_text(
      data     = labels_nationalparks
    , mapping  = aes(label = Label)
    , col      = "#004800"
    , size     = 2
    , fontface = 3
  ) +
  coord_sf(
      crs    = 4326
    , xlim   = c(xmin(extent(kaza)) - 0.2, xmax(extent(kaza)) + 0.2)
    , ylim   = c(ymin(extent(kaza)) - 0.2, ymax(extent(kaza)) + 0.2)
    , expand = F
  ) +
  labs(
      x = NULL
    , y = NULL
  ) +
  theme(
      panel.background = element_blank()
    , panel.border     = element_rect(colour = "black", fill = NA, size = 1)
    , plot.margin = unit(c(0, 0, 0, 0), "null")
  ) +
  annotation_scale(
      location   = "bl"
    , line_width = 0.5
    , height     = unit(0.15, "cm")
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(0.75, "cm"),
    , width    = unit(0.75, "cm"),
    , style    = north_arrow_orienteering(
          fill      = c("black", "black")
        , line_col  = "black"
        , text_col  = "black"
        , text_size = 9
      )
  )

################################################################################
#### Common Legend
################################################################################
# Create an arbitrary object with the number of levels for which we want to
# display the colors
prot$Random <- factor(rep_len(c(
    "National Park"
  , "Protected"
  , "Water"
), nrow(prot)), levels = c("Water", "National Park", "Protected"))

# Prepare legend
p3 <- ggplot() +
  geom_sf(
      data     = prot
    , aes(fill = factor(Random))
  ) +
  geom_sf(
      data        = africa
    , mapping     = aes(col = "Country Borders")
    , fill        = NA
    , show.legend = "line"
  ) +
  geom_sf(
      data        = kaza
    , mapping     = aes(col = "KAZA-TFCA")
    , fill        = NA
    , show.legend = "line"
  ) +
  scale_fill_manual(
    values = c("Water" = "cornflowerblue", "National Park" = "#70ab70", "Protected" = "#d9f0d3")
    , guide = guide_legend(
      override.aes = list(
        linetype = c("blank", "blank", "blank")
      )
    )
  ) +
  scale_color_manual(
      values = c("gray80", "black")
    , guide = guide_legend(
        override.aes = list(
          linetype = c(1, 1)
        )
      )
  ) +
  theme(
      legend.title     = element_blank()
    , legend.spacing.y = unit(0, "cm")
  )

# Extract legend
legend <- get_legend(p3)

################################################################################
#### Putting Plots together
################################################################################
p4 <- ggdraw(xlim = c(0, 22.5), ylim = c(0, 13.5)) +
  draw_plot(p1, x = 1, y = 1, width = 5.5, height = 7) +
  draw_plot(p2, x = 8, y = 1, width = 14, height = 12) +
  draw_plot(legend, x = 1, y = 8.5, width = 5.5, height = 4)
ggsave("test.pdf", plot = p4, width = 22.5 / 3, height = 13.5 / 3)

ggplot() +
  coord_equal(xlim = c(0, 20), ylim = c(0, 10), expand = F) +
  annotation_custom(
    ggplotGrob(p1), xmin = 0, xmax = 5, ymin = 0, ymax = 0.5
  ) +
  annotation_custom(
    ggplotGrob(p2), xmin = 0.3, xmax = 1, ymin = 0, ymax = 1
  )


x1 <- rnorm(100)
x2 <- 2 * x1
x3 <- seq(-pi, pi, 0.01)
y <- sin(x3)
p1 <- ggplot(data.frame(x1), aes(x = x1)) + geom_boxplot()
p2 <- ggplot(data.frame(x2), aes(x = x2)) + geom_histogram()
p3 <- ggplot(data.frame(x3, y), aes(x = x3, y = y)) + geom_line()

gt <- arrangeGrob(
    grobs         = list(p1, p2, p3)
  , layout_matrix = rbind(
      c(3, 3, 2, 2, 2, 2)
    , c(1, 1, 2, 2, 2, 2)
    , c(1, 1, 2, 2, 2, 2)
  )
)
grid.draw(gt)

ggdraw() +
  draw_plot(p1, x = 0, y = 0) +
  draw_plot(p2, x = 0.4, y = 0) +
  draw_plot(p3, x = 0, y = 0.5)


library(grid)
dummy_grob <- function(id)  {
  grobTree(rectGrob(gp=gpar(fill=id, alpha=0.5)), textGrob(id))
}
gs <- lapply(1:3, dummy_grob)
grid.arrange(ncol=4, grobs=gs,
               top="top\nlabel", bottom="bottom\nlabel",
               left="left\nlabel", right="right\nlabel")
grid.rect(gp=gpar(fill=NA))
plots <- gs
plots <- list(p1, p2, legend)
plots <- list(p1, p2, p3)
gt <- arrangeGrob(
    grobs         = plots
  , layout_matrix = rbind(
      c(3, 3, 2, 2, 2, 2)
    , c(1, 1, 2, 2, 2, 2)
    , c(1, 1, 2, 2, 2, 2)
  )
)
grid.draw(gt)
p1 <- ggplot() +
  geom_sf(data = africa) +
  theme_minimal()
p2 <- ggplot() +
  geom_sf(data = kaza) +
  theme_minimal()
p3 <- ggplot() +
  geom_sf(data = labels_countries) +
  theme_minimal()





library(egg)
plot_grid(p1, p2, rel_widths = c(0.2, 0.8), align = "v", axis = "b")
test <- grid.arrange(
    grobs = list(p1, p2, legend)
  , layout_matrix = rbind(
      c(3, 2, 2, 2)
    , c(1, 2, 2, 2)
  )
)
plot(test)
ggsave("hello.png")
plot_grid(legend, p2, p1)
ggarrange(p1, p2, ncol = 2, heights = c(1, 4))
test <- arrangeGrob(
    grobs = list(p1, p2)
  , layout_matrix = rbind(
      c(NA, 2)
    , c(1, 2)
  )
)
grid.newpage()
grid.draw(test)
library(grid)

plot(test)

library(ggpubr)
library(egg)
install.packages("egg")
egg::ggarrange(p1, p2)
library(gridExtra)
grid.arrange(
    grobs = list(p1, p2)
  , widths = c(1, 3)
  , layout_matrix = rbind(
      c(NA, 2)
    , c(1, 2)
  )
)

################################################################################
#### LEGACY
################################################################################
# Load required packages
library(rgdal)        # To load and store spatial data
library(raster)       # To manipulate raster data
library(rgeos)        # To manipulate spatial data
library(terra)        # To handle raster data
library(tmap)         # For beautiful maps
library(tidyverse)    # For data wrangling
library(davidoff)     # Custom functions
library(grid)         # To arrange multiple plots
library(Cairo)        # To store plots
library(elevatr)      # To create nice hillshade map

################################################################################
#### Plot of the Study Area
################################################################################
# Load required data
africa  <- readOGR("03_Data/02_CleanData/00_General_Africa_ESRI.shp")
kaza    <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")
dogs    <- readOGR("03_Data/02_CleanData/00_General_WildDogs_IUCN.shp")
prot    <- readOGR("03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")
water   <- rast("03_Data/02_CleanData/01_LandCover_WaterCoverAveraged_MERGED.tif")

# Clean up the water file for nicer illustration
water <- aggregate(water, fun = max, fact = 2)
water <- raster(water)
clump <- clump(water)
clump_freq <- data.frame(freq(clump))
clump_freq <- subset(clump_freq, count <= 40)
rcl <- data.frame(old = clump_freq$value, new = NA)
water <- reclassify(clump, rcl = rcl)
water <- !is.na(water)

# Remove 0 from water layer
water[water == 0] <- NA
plot(water)

# Simplify Protection zones
prot$Desig[prot$Desig == "Forest Reserve"] <- "Protected"
prot$Desig <- as.factor(as.character(prot$Desig))

# Remove small islands for plotting
africa <- subset(africa, !(ID %in% c(27:41, 689, 690)))

# Get the extent of the KAZA
kaza_ext <- as(extent(kaza), "SpatialPolygons")
crs(kaza_ext) <- CRS("+init=epsg:4326")

# Create labels for countries
labels_countries <- data.frame(
    x = c(16, 24, 11, 44, 35)
  , y = c(-8, -30, -27, -14, -24)
  , Label = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")
)
coordinates(labels_countries) <- c("x", "y")
crs(labels_countries) <- CRS("+init=epsg:4326")

# Create lines pointing towards these countries
l1 <- rbind(c(17, -10), c(19, -12)) %>% spLines()
l2 <- rbind(c(24, -29), c(24, -23)) %>% spLines()
l3 <- rbind(c(13, -25), c(17, -22)) %>% spLines()
l4 <- rbind(c(37, -14), c(29, -15)) %>% spLines()
l5 <- rbind(c(33, -22), c(30, -19)) %>% spLines()
lines_countries <- rbind(l1, l2, l3, l4, l5)
crs(lines_countries) <- crs(labels_countries)

# Visualize
plot(labels_countries)
plot(lines_countries, add = T)
plot(lines_countries[4], add = T, col = "red")
plot(africa, add = T)

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
    x     = c(22.6, 23.7, 27.8, 25.6, 22.9, 27.8, 27.25)
  , y     = c(-19, -18.2, -17, -20.7, -15.0, -14.4, -15.58)
  , Label = c(
    "Okavango\nDelta", "Linyanti\nSwamp", "Lake\nKariba", "Makgadikgadi\nPans"
    , "Barotse\nFloodplain", "Lukanga\nSwamp", "Kafue\nFlats"
  )
)
coordinates(labels_waters) <- c("x", "y")
crs(labels_waters) <- CRS("+init=epsg:4326")

# Visualize
plot(water)
plot(labels_waters, add = T, col = "red")

# Create labels for some national parks
labels_nationalparks <- data.frame(
    x = c(26.56, 28.61, 21.15, 25.87, 20.38, 23.58, 23.71, 24.51, 20.78, 22.63, 27.92, 28.54)
  , y = c(-19.08, -17.05, -17.26, -14.66, -16.08, -21.4, -19.29, -18.65, -18.81, -14.54, -17.76, -20.53)
  , Label = paste0(c(
      "Hwange", "Matusadona", "Luengue-Luiana", "Kafue", "Mavinga"
    , "Central Kalahari", "Moremi", "Chobe", "Khaudum", "Liuwa Plains"
    , "Chizarira", "Matobo"
  ), "\nNP")
)
coordinates(labels_nationalparks) <- c("x", "y")
crs(labels_nationalparks) <- CRS("+init=epsg:4326")

# Visualize
plot(subset(prot, Desig == "National Park"))
plot(labels_nationalparks, add = T, col = "red")

# Create elevation map for africa
elev_africa <- get_elev_raster(africa, z = 4)
elev_africa <- crop(elev_africa, africa)
elev_africa <- mask(elev_africa, africa, updatevalue = NA)
terrain_africa <- terrain(elev_africa ** 1.6, opt = c("slope", "aspect"))
hill_africa <- hillShade(aspect = terrain_africa$aspect, slope = terrain_africa$slope)

# Visualize it
plot(hill_africa, col = gray(80:100/100))
plot(africa, add = T, border = "gray80")

# Create elevation map for kaza
ext <- as(extent(kaza) + c(-1, 1, -1, 1), "SpatialPolygons")
crs(ext) <- CRS("+init=epsg:4326")
elev_kaza <- get_elev_raster(kaza, z = 6)
elev_kaza <- crop(elev_kaza, ext)
elev_kaza <- mask(elev_kaza, ext, updatevalue = NA)
terrain_kaza <- terrain(elev_kaza ** 1.6, opt = c("slope", "aspect"))
hill_kaza <- hillShade(aspect = terrain_kaza$aspect, slope = terrain_kaza$slope)

# Visualize it
plot(hill_kaza, col = gray(80:100/100))
plot(africa, add = T, border = "gray80")
plot(kaza, add = T, border = "gray80")

# Prepare a map of Africa
p1 <- tm_shape(hill_africa) +
  tm_raster(
      style       = "cont"
    , palette     = gray(60:100/100)
    , alpha       = 0.2
    , legend.show = F
  ) +
    tm_shape(dogs) +
  tm_polygons(
      col        = "cornflowerblue"
    , border.col = "cornflowerblue"
  ) +
  tm_shape(africa, is.master = T) +
    tm_borders(
        col        = "gray80"
      , lwd        = 0.5
    ) +
  tm_shape(kaza) +
    tm_borders(
        col = "black"
      , lty = 1
      , lwd = 1.5
    ) +
  tm_shape(lines_countries) +
    tm_lines(
        col = "gray30"
      , lty = 1
      , lwd = 1
    ) +
  tm_shape(labels_countries) +
    tm_text(
        text = "Label"
      , size = 0.6
      , col  = "black"
    ) +
  tm_layout(
      asp         = 0.8
    , frame       = "black"
    , frame.lwd   = 3
    , legend.show = F
  ) +
  tm_credits("a"
    , position  = c("right", "top")
    , size      = 1.5
    , col       = "black"
    , fontface  = "bold"
)

# Plot of the KAZA
p2 <- tm_shape(hill_kaza) +
  tm_raster(
      style       = "cont"
    , palette     = gray(60:100/100)
    , alpha       = 0.2
    , legend.show = F
  ) +
  tm_shape(prot) +
    tm_polygons(
        col         = "Desig"
      , palette     = c("#70ab70", "#d9f0d3")
      , lwd         = 0.1
      , border.col  = "#6ba36b"
      , legend.show = F
    ) +
  tm_shape(water) +
    tm_raster(
        palette     = "#96d0ff"
      , legend.show = F
    ) +
  tm_shape(kaza, is.master = T) +
    tm_borders(
        col = "black"
      , lty = 1
      , lwd = 2
    ) +
  tm_shape(africa) +
    tm_borders(
        col = "gray50"
      , lwd = 0.5
    ) +
  tm_shape(labels_countries2) +
    tm_text("Label"
      , col       = "black"
      , fontface  = 2
      , size      = 1.5
    ) +
  tm_shape(labels_waters) +
    tm_text("Label"
      , fontface  = 3
      , size      = 0.5
      , col       = "#064886"
    ) +
  tm_shape(labels_nationalparks) +
    tm_text("Label"
      , col       = "#004800"
      , fontface  = 3
      , size      = 0.5
    ) +
  tm_graticules(
      n.y                 = 5
    , n.x                 = 5
    , labels.inside.frame = FALSE
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_layout(
    , frame                   = "gray20"
    , frame.lwd               = 3
    , asp                     = 1.2
    , legend.outside          = TRUE
    , legend.outside.position = "left"
    , legend.stack            = "vertical"
    , legend.text.size        = 0.8
  ) +
  tm_scale_bar(
        position  = c("right", "bottom")
      , text.size = 0.5
      , text.col  = "black"
      , width     = 0.125
  ) +
  tm_credits("b"
    , position = c("left", "top")
    , size     = 1.5
    , col      = "black"
    , fontface = "bold"
  ) +
  tm_compass(
      color.dark  = "black"
    , color.light = "black"
    , text.color  = "black"
    , position    = c("left", "bottom")
)

################################################################################
#### Prepare a common Legend
################################################################################
# We now prepare a legend that the two plots share. However, we have to assign
# to one of the two plots. Let's use p2 for this
p2 <- p2 + tm_add_legend(
    type    = "fill"
  , labels  = c(
      "Wild Dog Populations"
    , "Water Bodies"
    , "Protected Areas"
    , "National Parks (NP)"
  )
  , col = c(
      "cornflowerblue"
    , "#96d0ff"
    , "#d9f0d3"
    , "#70ab70"
  )) + tm_add_legend(
    type    = "line"
  , labels  = c("KAZA-TFCA Borders", "Country Borders")
  , col     = c("black", "gray50")
  , lty     = c(1, 1)
  , lwd     = c(2, 1)
) + tm_layout(legend.frame.lwd = 2, legend.text.size = 1.05)

tmap_save(
    tm        = p2
  , filename  = "04_Manuscript/99_StudyArea.png"
  , insets_tm = p1
  , insets_vp = viewport(0.164, 0.341, width = 0.56, height = 0.56)
  , width     = 9
  , height    = 5.25
  , scale     = 1.1
)

# # Combine plots and store them
# CairoPDF("04_Manuscript/99_StudyArea.pdf", width = 9, height = 5.25, bg = "transparent")
# p2
# print(p1, vp = viewport(0.158, 0.341, width = 0.56, height = 0.56))
# dev.off()
#
# # Store as png as well
# png("04_Manuscript/99_StudyArea.png", width = 1080, height = 720, pointsize = 20)
# p2
# print(p1, vp = viewport(0.158, 0.341, width = 0.56, height = 0.56))
# dev.off()
