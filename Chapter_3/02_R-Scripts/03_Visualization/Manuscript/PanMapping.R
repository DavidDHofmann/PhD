################################################################################
#### Plots to Exemplify our Pan Mapping Approach
################################################################################
# Clear R's brain
rm(list = ls())

# Load required packages
library(rgdal)         # To handle spatial data
library(sf)            # To handle spatial data
library(raster)        # To handle spatial data
library(terra)         # To handle spatial data
library(lubridate)     # To handle dates
library(tidyverse)     # For data wrangling
library(ggridges)      # For ridgelines in ggplot
library(rpart)         # To train a classifier
library(randomForest)  # To train a classifier
library(rpart.plot)    # To visualize classifier
library(caret)         # To assess variable importance
library(ggpubr)        # To arrange multiple ggplots
library(ggstance)      # To plot pointranges
library(tidyterra)     # To use dplyr with terra objects
library(ggspatial)     # For annotation scale and arrow
library(grid)          # To add pseudofacette
library(gtable)        # To add pseudofacette
library(ggh4x)         # To get nested wraps
library(viridis)       # Access to the true viridis color scale

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Load custom functions
source("02_R-Scripts/00_Functions.R")

################################################################################
#### Sentinel Tiles
################################################################################
# Reload the moving windows
wind <- read_rds("03_Data/02_CleanData/Windows.rds")

# Extract the different objects that we want to plot
windows <- wind %>%
  pull(Window) %>%
  do.call(rbind, .) %>%
  st_as_sf()
gps <- wind %>%
  pull(GPS) %>%
  do.call(rbind, .)
tiles <- wind %>%
  pull(Tiles) %>%
  do.call(rbind, .) %>%
  st_as_sf() %>%
  select(-Overlap) %>%
  distinct()
africa <- read_sf("03_Data/02_CleanData/Africa.gpkg")

# Create country labels
labels_countries <- data.frame(
    x     = c(25.2, 27.3, 25.8, 21.5, 21.5)
  , y     = c(-21.0, -19.5, -17.2, -17.3, -18.1)
  , Label = c("BOTSWANA", "ZIMBABWE", "ZAMBIA", "ANGOLA", "NAMIBIA")
)

# Define a highlight color
col <- "cornflowerblue"

# Visualize all
p1 <- ggplot() +
  geom_sf(
      data      = africa
    , fill      = NA
    , color     = "black"
    , linewidth = 0.1
  ) +
  geom_sf(
      data  = tiles
    , fill  = "gray50"
    , color = "gray90"
    , alpha = 0.15
  ) +
  geom_sf_text(
      data    = tiles
    , mapping = aes(label = TileID)
    , color   = "gray50"
    , size    = 2
  ) +
  geom_sf(
      data      = windows
    , fill      = col
    , color     = col
    , alpha     = 0.2
    , linewidth = 0.1
  ) +
  geom_point(
      data    = gps
    , mapping = aes(x = x, y = y)
    , col     = "black"
    , size    = 0.1
  ) +
  geom_text(
      data     = labels_countries
    , mapping  = aes(x = x, y = y, label = Label)
    , col      = "black"
    , fontface = 2
  ) +
  theme_awesome() +
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
      legend.position = "bottom"
    , panel.border    = element_rect(colour = "gray30", fill = NA, size = 0.5)
  ) +
  xlab("") +
  ylab("") +
  coord_sf(
      crs    = 4326
    , xlim   = c(20.5, 28)
    , ylim   = c(-22, -17)
    , expand = F
  )

# Generate a custom legend
df <- africa %>%
  slice(1:2) %>%
  mutate(Name = factor(
      c("Sentinel 2 Tiles", "Moving Windows")
    , levels = c("Sentinel 2 Tiles", "Moving Windows"))
  )

# Legend item for moving windows
gps$Name <- "Moving Windows"
l1 <- ggplot() +
  geom_point(
      data    = gps[1, ]
    , mapping = aes(x = x, y = y, fill = Name, color = Name)
    , shape   = 21
    , size    = 5
    , alpha   = 0.2
  ) +
  scale_fill_manual(values = col, name = "") +
  scale_color_manual(values = col, name = "") +
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_blank())

# Legend item for tiles
gps$Name <- "Sentinel 2 Tiles"
l2 <- ggplot() +
  geom_point(
      data    = gps[1, ]
    , mapping = aes(x = x, y = y, fill = Name, color = Name)
    , shape   = 22
    , size    = 5
    , alpha   = 0.15
  ) +
  scale_fill_manual(values = c("gray50"), name = "") +
  scale_color_manual(values = c("gray90"), name = "") +
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_blank())

# Legend item for GPS dispersal points
gps$Name <- "GPS Data during Dispersal"
l3 <- ggplot() +
  geom_point(
      data    = gps[1, ]
    , mapping = aes(x = x, y = y, color = Name)
  ) +
  scale_color_manual(values = c("black"), name = "") +
  theme_void() +
  theme(legend.position = "bottom", legend.title = element_blank())

# Extract the legends and put them together
l1 <- get_legend(l1)
l2 <- get_legend(l2)
l3 <- get_legend(l3)
p2 <- ggarrange(l2, l1, l3, nrow = 1)

# Adjust margins
p1 <- p1 + theme(plot.margin = unit(c(5.5, 5.5, 1.0, 5.5), "pt"))
p2 <- p2 + theme(plot.margin = unit(c(1.0, 5.5, 5.5, 5.5), "pt"))
p <- ggarrange(p1, p2, heights = c(0.95, 0.05), nrow = 2)

# Store them separately
ggsave(
    plot     = p
  , filename = "04_Manuscript/Figures/MovingWindows.png"
  , device   = png
  , bg       = "white"
  , width    = 6
  , height   = 4.5
  , scale    = 1.2
)

################################################################################
#### Reflectances
################################################################################
# Load results from pan mapping
validation <- read_rds("03_Data/03_Results/PanMapping.rds")

# Load some things to add to the plots
africa <- st_read("03_Data/02_CleanData/Africa.gpkg")
classe <- st_read("03_Data/02_CleanData/TrainingClasses.gpkg")

# Function to add a pseudofacette to our plots
addPseudofacette <- function(p, label) {

  # Convert to grob and add another row
  z <- ggplotGrob(p)
  z <- gtable_add_rows(z, z$height[7], pos = 6)

  # Check the layout
  # gtable_show_layout(z)

  # Span strip
  z <- gtable_add_grob(z, list(
      rectGrob(gp = gpar(col = NA, fill = "gray95", size = .5))
    , textGrob(label = label, gp = gpar(cex = .75, col = "black"))
    ), t = 7, l = 5, b = 7, r = 17, name = c("a", "b"))

  # Add small gap between strips - below row 6
  z <- gtable_add_rows(z, unit(2/10, "line"), 7)

  # Return it
  return(as_ggplot(z))
}

# Slightly adjust naming of the columns
validation <- mutate(validation, Data = map(Data, function(x) {
  oldnames <- colnames(x)
  newnames <- str_replace_all(oldnames, "\\d+", function(m) sprintf("%02d", as.integer(m)))
  newnames[2:length(newnames)] <- toupper(newnames[2:length(newnames)])
  names(x) <- newnames
  return(x)
}))

# Define colors
cols <- viridis(3, begin = 0.5, direction = -1)

# Histograms of Landsat & Sentinel Reflectances
p1 <- validation$Data[validation$Satellite == "Landsat"][[1]] %>%
  pivot_longer(B01:BEST, names_to = "Band", values_to = "Reflectance") %>%
  ggplot(aes(x = Reflectance, y = Class, col = Class, fill = Class)) +
    geom_density_ridges(alpha = 0.5) +
    facet_wrap("Band", scales = "free") +
    theme_awesome() +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    theme(
        axis.text.x       = element_blank()
      , axis.text.y       = element_blank()
      , axis.title.y      = element_blank()
      , legend.position   = "bottom"
      , legend.title      = element_blank()
      , legend.key.height = unit(0.1, "cm")
      , panel.grid.minor  = element_blank()
    )
p2 <- validation$Data[validation$Satellite == "Sentinel"][[1]] %>%
  pivot_longer(B01:BEST, names_to = "Band", values_to = "Reflectance") %>%
  ggplot(aes(x = Reflectance, y = Class, col = Class, fill = Class)) +
    geom_density_ridges(alpha = 0.5) +
    facet_wrap("Band", scales = "free") +
    theme_awesome() +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    theme(
        axis.text.x       = element_blank()
      , axis.text.y       = element_blank()
      , axis.title.y      = element_blank()
      , legend.position   = "bottom"
      , legend.title      = element_blank()
      , legend.key.height = unit(0.1, "cm")
      , panel.grid.minor  = element_blank()
    )

# Add a pseudofacettes with the title
p1 <- addPseudofacette(p1, "Landsat: Reflectance Profiles")
p2 <- addPseudofacette(p2, "Sentinel: Reflectance Profiles")

################################################################################
#### Decision Tree for CART
################################################################################
# Plot of Decision Trees
dend <- subset(validation, Model == "Cart" & Satellite == "Landsat") %>%
    pull(ModelObject) %>%
    first() %>%
    dendro_data() %>%
    stretchDendro(0.18)
p3 <- ggplot() +
  geom_segment(data = dend$segments, aes(x = x, xend = xend, y = y, yend = yend)) +
  geom_label(data = dend$leaf_labels, aes(x = x, y = y, label = label, fill = label), alpha = 0.5, vjust = 1, label.size = NA, size = 2.5) +
  geom_label(data = dend$labels, aes(x = x, y = y, label = label), vjust = 0.5, size = 3, color = "black", label.size = NA) +
  geom_text(data = dend$labels[1, ], aes(x = x, y = y, label = "yes"), nudge_x = -1, nudge_y = -0.04, fontface = 3, size = 2.5) +
  geom_text(data = dend$labels[1, ], aes(x = x, y = y, label = "no"), nudge_x = 1, nudge_y = -0.04, fontface = 3, size = 2.5) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(trans = "log") +
  facet_wrap(~ "Landsat: CART Decision Tree") +
  theme_awesome() +
  theme(
      panel.grid      = element_blank()
    , axis.text       = element_blank()
    , axis.title.x    = element_blank()
    , axis.title.y    = element_blank()
    , legend.position = "none"
  )

# Use the dendro package to create a gg-plottable object
dend <- subset(validation, Model == "Cart" & Satellite == "Sentinel") %>%
    pull(ModelObject) %>%
    first() %>%
    dendro_data() %>%
    stretchDendro(0.05)
p4 <- ggplot() +
  geom_segment(data = dend$segments, aes(x = x, xend = xend, y = y, yend = yend)) +
  geom_label(data = dend$leaf_labels, aes(x = x, y = y, label = label, fill = label), alpha = 0.5, vjust = 1, label.size = NA, size = 2.5) +
  geom_label(data = dend$labels, aes(x = x, y = y, label = label), vjust = 0.5, size = 3, color = "black", label.size = NA) +
  geom_text(data = dend$labels[1, ], aes(x = x, y = y, label = "yes"), nudge_x = -1, nudge_y = -0.04, fontface = 3, size = 2.5) +
  geom_text(data = dend$labels[1, ], aes(x = x, y = y, label = "no"), nudge_x = 1, nudge_y = -0.04, fontface = 3, size = 2.5) +
  scale_fill_manual(values = cols) +
  scale_y_continuous(trans = "log") +
  facet_wrap(~ "Sentinel: CART Decision Tree") +
  theme_awesome() +
  theme(
      panel.grid      = element_blank()
    , axis.text       = element_blank()
    , axis.title.x    = element_blank()
    , axis.title.y    = element_blank()
    , legend.position = "none"
  )

################################################################################
#### Variable Importance for RandomForest
################################################################################
p5 <- subset(validation, Model == "RandomForest" & Satellite == "Landsat") %>%
  pull(ModelObject) %>%
  first() %>%
  varImp() %>%
  as.data.frame() %>%
  mutate(Band = rownames(.)) %>%
  ggplot(aes(x = reorder(Band, Overall), y = Overall)) +
    geom_segment(aes(x = reorder(Band, Overall), xend = reorder(Band, Overall), y = 0, yend = Overall), col = "gray50") +
    geom_point() +
    theme_awesome() +
    coord_flip() +
    xlab("Band / Index") +
    ylab("Importance") +
    facet_wrap(~ "Landsat: Random Forest Variable Importance") +
    theme(axis.title.y = element_text(angle = 90))
p6 <- subset(validation, Model == "RandomForest" & Satellite == "Sentinel") %>%
  pull(ModelObject) %>%
  first() %>%
  varImp() %>%
  as.data.frame() %>%
  mutate(Band = rownames(.)) %>%
  ggplot(aes(x = reorder(Band, Overall), y = Overall)) +
    geom_segment(aes(x = reorder(Band, Overall), xend = reorder(Band, Overall), y = 0, yend = Overall), col = "gray50") +
    geom_point() +
    theme_awesome() +
    coord_flip() +
    xlab("Band / Index") +
    ylab("Importance") +
    facet_wrap(~ "Sentinel: Random Forest Variable Importance") +
    theme(axis.title.y = element_text(angle = 90))

# Put the plots together
p <- ggarrange(p1, p2, p3, p4, p5, p6, labels = c("a1", "a2", "b1", "b2", "c1", "c2"), nrow = 3, ncol = 2)
ggsave(
    plot     = p
  , filename = "04_Manuscript/Figures/Reflectances.png"
  , bg       = "white"
  , height   = 6
  , width    = 5
  , scale    = 2
  , device   = png
)

################################################################################
#### Validation
################################################################################
# Define highlight a color
col <- hcl.colors(10, palette = "viridis")[8]

# Visualize confusion matrices
p6 <- validation %>% select(Satellite, Model, Confusion) %>%
  unnest(Confusion) %>%
  ggplot(aes(x = observed, y = predicted, fill = Freq)) +
    geom_tile(colour = "black") +
    geom_text(aes(label = Freq), col = "black", size = 2) +
    facet_nested_wrap(~ Satellite + Model, ncol = 4) +
    scale_fill_gradient(low = "white", high = adjustcolor(col, alpha.f = 0.3)) +
    theme_minimal() +
    theme(
        legend.position  = "none"
      , panel.grid       = element_blank()
      , strip.background = element_rect(fill = "grey95", color = "white")
    ) +
    coord_equal() +
    xlab("Observed") +
    ylab("Predicted")
p7 <- validation %>%
  select(Satellite, Model, Specificity, Sensitivity, Accuracy) %>%
  pivot_longer(3:5, names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(y = Metric, x = Value, xmin = 0, xmax = Value)) +
    geom_vline(xintercept = 1, lty = 2, col = "gray50", linewidth = 0.5) +
    geom_pointrangeh(position = position_dodgev(height = 0.5), fatten = 2) +
    theme_minimal() +
    xlab("Score") +
    ylab("Metric") +
    facet_nested_wrap(~ Satellite + Model, ncol = 4) +
    theme(
        legend.position  = "none"
      , panel.grid.minor = element_blank()
      , strip.background = element_rect(fill = "grey95", color = "white")
    )

# Put the plots together
p <- ggarrange(p6, p7, nrow = 2, labels = "auto", align = "hv")

# Store them
ggsave(
    plot     = p
  , filename = "04_Manuscript/Figures/ClassificationValidation.png"
  , height   = 5
  , width    = 8
  , bg       = "white"
  , scale    = 1.2
  , device   = png
)
