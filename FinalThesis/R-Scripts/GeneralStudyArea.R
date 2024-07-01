################################################################################
#### General Study Area
################################################################################
# Description: Plot of the general study area

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/FinalThesis"
setwd(wd)

# Load required packages
library(terra)        # To handle raster files
library(raster)       # To handle raster files
library(sf)           # To handle spatial data
library(tidyverse)    # For data wrangling
library(ggspatial)    # To add scale bars and north arrows
library(colorspace)   # To lighten or darken colors
library(smoothr)      # To smoothen geometries
library(tidyterra)    # To handle spatial data in a tidy way

################################################################################
#### Load Required data
################################################################################
# Load map of africa and create a simplified version of it
africa      <- vect("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/00_General_Africa_ESRI.shp")
africa_simple      <- buffer(africa, width = 1000)
africa_simple      <- aggregate(africa_simple)
africa_simple      <- disagg(africa_simple)
africa_simple$Area <- expanse(africa_simple, unit = "km")
africa_simple      <- filter(africa_simple, Area %in% sort(Area, decreasing = T)[1:2])

# Also load and merge the present and historic distribution of wild dogs
dogs1 <- vect("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/00_General_WildDogs_IUCN.shp")[1, ]
dogs2 <- vect("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/01_RawData/DAVID/HistoricRange.shp")
dogs2 <- smooth(dogs2)
dogs2 <- crop(dogs2, aggregate(africa))
dogs  <- rbind(dogs1, dogs2)
values(dogs) <- data.frame(Type = c("Present", "Historic"))
rm(dogs1, dogs2)

# Load an average water map and simplify it
# water <- rast("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/01_LandCover_WaterCoverAveraged_MERGED.tif")
# water <- aggregate(water, fun = max, fact = 2)
# water <- raster(water)
# clump <- clump(water)
# clump_freq <- data.frame(freq(clump))
# clump_freq <- subset(clump_freq, count <= 40)
# rcl <- data.frame(old = clump_freq$value, new = NA)
# water <- reclassify(clump, rcl = rcl)
# water <- !is.na(water)
# water <- as.data.frame(water, xy = T)
# water <- subset(water, layer)
# rm(clump, clump_freq, rcl)
water  <- vect("/home/david/ownCloud/University/15. PhD/Chapter_2/03_Data/02_CleanData/MajorWaters.gpkg")

# Load additional data we would like to plot
kaza   <- vect("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/00_General_KAZA_KAZA.shp")
prot   <- vect("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/02_LandUse_Protected_PEACEPARKS.shp")
river  <- vect("/home/david/ownCloud/University/15. PhD/Chapter_2/03_Data/02_CleanData/MajorRivers.gpkg")
roads  <- vect("/home/david/ownCloud/University/15. PhD/Chapter_2/03_Data/02_CleanData/Roads.gpkg")
main   <- vect("/home/david/ownCloud/University/15. PhD/General/Supervision & Collaboration/CovariatesForRaph/Covariates/StudyArea.gpkg")
vills  <- read_sf("/home/david/ownCloud/University/15. PhD/Chapter_2/03_Data/02_CleanData/Villages.gpkg")
vills  <- cbind(st_drop_geometry(vills), st_coordinates(vills)) %>%
  rename(x = X, y = Y) %>%
  mutate(place = gsub(place, pattern = "Village", replacement = "Villages")) %>%
  mutate(place = gsub(place, pattern = "City", replacement = "Cities"))
disp <- read_csv("/home/david/ownCloud/University/15. PhD/Chapter_2/03_Data/02_CleanData/Dispersers.csv")
disp <- subset(disp, State == "Disperser")

# Reorder the levels of national parks
prot$Desig <- factor(prot$Desig, levels = c("National Park", "Forest Reserve", "Protected"))

# Three different extents
ext1 <- as.polygons(ext(prot), crs = "epsg:4326")
ext2 <- as.polygons(ext(c(22, 27, -20.7, -17.9)), crs = "epsg:4326")

# Convert all to sf
# africa    <- st_as_sf(africa)
# kaza      <- st_as_sf(kaza)
# dogs      <- st_as_sf(dogs) %>% arrange(Type)
# main      <- st_as_sf(main)
# ext1_sf   <- st_as_sf(ext1)
# ext2_sf   <- st_as_sf(ext2)

################################################################################
#### Plot 1: Africa
################################################################################
# Prepare plot
p1 <- ggplot() +
  geom_sf(data = africa, fill = "gray90", col = "white", linewidth = 0.3) +
  geom_sf(data = dogs, aes(fill = Type, col = Type)) +
  geom_sf(data = kaza, fill = "red", col = "red", alpha = 0.1, linewidth = 0.4) +
  # geom_sf(data = ext1, fill = "transparent", col = "black", linewidth = 0.5, lty = 2) +
  geom_sf(data = ext2, fill = "orange", col = "orange", alpha = 0.3, linewidth = 0.4) +
  coord_sf(xlim = c(-18, 48), ylim = c(-32, 35)) +
  # scale_fill_manual(values = c(adjustcolor("orange", alpha.f = 0.2), "orange")) +
  # scale_color_manual(values = c(adjustcolor("orange", alpha.f = 0.2), "orange")) +
  scale_fill_manual(values = c(adjustcolor("gray50", alpha.f = 0.2), "gray50")) +
  scale_color_manual(values = c(adjustcolor("gray50", alpha.f = 0.2), "gray50")) +
  theme_minimal() +
  theme(
      legend.position = "none"
    , panel.grid      = element_blank()
    , axis.ticks      = element_blank()
    , axis.text       = element_blank()
    , axis.title      = element_blank()
  )

# ################################################################################
# #### Plot 2: KAZA
# ################################################################################
# # Country labels
# labels_countries <- data.frame(
#     x = c(20.39, 23.94, 20.07, 25.69, 28.22)
#   , y = c(-15.28, -19.94, -19.39, -15.22, -18.9)
#   , Label = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")
# )
#
# # Create labels for some geographical landmarks
# labels_waters <- data.frame(
#     x     = c(22.6, 23.7, 27.8, 25.6, 22.9, 27.8, 27.25)
#   , y     = c(-19, -18.2, -17, -20.7, -15.0, -14.4, -15.58)
#   , Label = c(
#     "Okavango\nDelta", "Linyanti\nSwamp", "Lake\nKariba", "Makgadikgadi\nPans"
#     , "Barotse\nFloodplain", "Lukanga\nSwamp", "Kafue\nFlats"
#   )
# )
#
# # Create labels for some national parks
# labels_nationalparks <- data.frame(
#     x = c(26.56, 28.61, 21.15, 25.87, 20.38, 23.58, 23.71, 24.51, 20.78, 22.63, 27.92, 28.54)
#   , y = c(-19.08, -17.05, -17.26, -14.66, -16.08, -21.4, -19.29, -18.65, -18.81, -14.54, -17.76, -20.53)
#   , Label = paste0(c(
#       "Hwange", "Matusadona", "Luengue-Luiana", "Kafue", "Mavinga"
#     , "Central Kalahari", "Moremi", "Chobe", "Khaudum", "Liuwa Plains"
#     , "Chizarira", "Matobo"
#   ), "\nNP")
# )
#
# # Prepare plot
# p2 <- ggplot() +
#   geom_sf(data = prot, aes(fill = Desig), col = NA, alpha = 0.7) +
#   geom_sf(data = kaza, fill = "transparent", col = "black", linewidth = 0.5) +
#   geom_sf(data = ext1, fill = "transparent", col = "black", linewidth = 0.5, lty = 2) +
#   # geom_sf(data = ext2, fill = "transparent", col = "black", linewidth = 0.5, lty = 2) +
#   geom_raster(data = water, aes(x = x, y = y), fill = "cornflowerblue") +
#   geom_sf(data = river, col = "cornflowerblue", linewidth = 0.1) +
#   geom_text(
#       data     = labels_waters
#     , mapping  = aes(x = x, y = y, label = Label)
#     , col      = darken("cornflowerblue", 0.4)
#     , fontface = 3
#     , size     = 3
#   ) +
#   geom_text(
#       data     = labels_nationalparks
#     , mapping  = aes(x = x, y = y, label = Label)
#     , col      = darken("green", 0.8)
#     , fontface = 3
#     , size     = 3
#   ) +
#   geom_text(
#       data     = labels_countries
#     , mapping  = aes(x = x, y = y, label = Label)
#     , col      = "black"
#     , fontface = 2
#   ) +
#   geom_sf(data = ext2, fill = "orange", col = "orange", alpha = 0.3, linewidth = 0.4) +
#   scale_fill_brewer(palette = "Greens", direction = -1) +
#   coord_sf(
#       xlim   = c(xmin(ext1), xmax(ext1))
#     , ylim   = c(ymin(ext1), ymax(ext1))
#     , expand = F
#   ) +
#   annotation_scale(
#       location   = "br"
#     , width_hint = 0.2
#     , line_width = 0.5
#     , height     = unit(0.15, "cm")
#   ) +
#   annotation_north_arrow(
#       location = "br"
#     , height   = unit(1.5, "cm"),
#     , width    = unit(1.2, "cm"),
#     , pad_y    = unit(0.8, "cm")
#     , style    = north_arrow_fancy_orienteering(
#           fill      = c("black", "black")
#         , line_col  = NA
#         , text_col  = "black"
#         , text_size = 12
#       )
#   ) +
#   theme_minimal() +
#   theme(
#       legend.position = "none"
#     , axis.title      = element_blank()
#   )

################################################################################
#### Plot 3: BPC Area
################################################################################
# Create country labels
labels_countries <- data.frame(
    x     = c(25.2, 26.4, 25.7, 21.5, 22.5)
  , y     = c(-19.3, -18.2, -17.6, -17.6, -18.0)
  , Label = c("BOTSWANA", "ZIMBABWE", "ZAMBIA", "ANGOLA", "NAMIBIA")
)

# Create labels for some geographical landmarks
labels_waters <- data.frame(
    x     = c(22.6, 23.7, 22.7, 25.3)
  , y     = c(-19, -18.1, -20.4, -20.5)
  , Label = c("Okavango\nDelta", "Linyanti\nSwamp", "Lake\nNgami", "Makgadikgadi\nPans")
)

# Create labels for some national parks and source areas
labels_nationalparks <- data.frame(
    x     = c(26.3, 22.35, 23.67, 24.51, 24.7, 23.35)
  , y     = c(-19.08, -17.70, -19.35, -18.65, -20.4, -21.3)
  , Label = paste0(c("Hwange", "Luengue-Luiana", "Moremi", "Chobe", "Nxai Pan", "Central Kalahari"), "\nNP")
)

# Create labels for rivers
labels_rivers <- data.frame(
    x     = c(20.65, 23.18, 23.15, 23.9)
  , y     = c(-17.85, -18.65, -20.35, -20.25)
  , Label = c("Okavango", "Selinda", "Nhabe", "Boteti")
)

# Prepare plot
p2 <- ggplot() +
  geom_sf(data = prot, aes(fill = Desig), col = NA, alpha = 0.7) +
  # geom_sf(data = ext1, fill = "transparent", col = "black", linewidth = 0.5, lty = 2) +
  geom_sf(data = ext2, fill = "transparent", col = "black", linewidth = 0.5) +
  # geom_raster(data = water, aes(x = x, y = y), fill = "cornflowerblue") +
  geom_sf(data = water, fill = "cornflowerblue", col = "cornflowerblue") +
  geom_sf(data = river, col = "cornflowerblue", linewidth = 0.1) +
  geom_sf(data = roads, col = "gray70", lwd = 0.2) +
  geom_sf(data = africa, col = "black", fill = NA, lwd = 0.3) +
  geom_point(data = disp, mapping = aes(x = x, y = y), size = 0.1, col = "gray30") +
  geom_path(data = disp, mapping = aes(x = x, y = y, group = ID), linewidth = 0.1, col = "gray30") +
  geom_sf(data = main, fill = "orange", col = "orange", alpha = 0.3, linewidth = 0.4) +
  geom_text(
      data     = labels_waters
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = darken("cornflowerblue", 0.4)
    , fontface = 3
    , size     = 3
  ) +
  geom_text(
      data     = labels_rivers
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = darken("cornflowerblue", 0.4)
    , fontface = 3
    , size     = 1.5
    , angle    = c(-30, 40, 20, -32)
  ) +
  geom_text(
      data     = subset(labels_nationalparks, Label != "Moremi\nNP")
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = darken("green", 0.8)
    , fontface = 3
    , size     = 3
  ) +
  geom_text(
      data     = labels_countries
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = "black"
    , fontface = 2
  ) +
  coord_sf(
      xlim   = c(xmin(ext2), xmax(ext2))
    , ylim   = c(ymin(ext2), ymax(ext2))
    , expand = F
  ) +
  scale_fill_brewer(palette = "Greens", direction = -1, name = "") +
  annotation_scale(
      location   = "br"
    , width_hint = 0.2
    , line_width = 0.5
    , height     = unit(0.15, "cm")
  ) +
  annotation_north_arrow(
      location = "br"
    , height   = unit(1.5, "cm"),
    , width    = unit(1.2, "cm"),
    , pad_y    = unit(0.8, "cm")
    , style    = north_arrow_fancy_orienteering(
          fill      = c("black", "black")
        , line_col  = NA
        , text_col  = "black"
        , text_size = 12
      )
  ) +
  theme_minimal() +
  theme(
      legend.position = "none"
    , axis.title      = element_blank()
  )

################################################################################
#### Store all the Plots to File
################################################################################
ggsave(
    filename = "/home/david/ownCloud/University/15. PhD/FinalThesis/Chapter_00/Figures/StudyArea1.png"
  , plot     = p1
  , device   = png
  , bg       = "white"
)
# ggsave(
#     filename = "/home/david/ownCloud/University/15. PhD/FinalThesis/Chapter_00/Figures/StudyArea2.png"
#   , plot     = p2
#   , device   = png
#   , bg       = "white"
#   , scale    = 1.3
#   , width    = 7
#   , height   = 6
# )
ggsave(
    filename = "/home/david/ownCloud/University/15. PhD/FinalThesis/Chapter_00/Figures/StudyArea2.png"
  , plot     = p2
  , device   = png
  , bg       = "white"
  , scale    = 1.3
  , width    = 7
  , height   = 4
)
