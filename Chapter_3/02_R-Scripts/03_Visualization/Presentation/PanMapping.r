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
library(ggdark)        # Access to dark themes

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Load results from pan mapping
validation <- read_rds("03_Data/03_Results/99_PanMapping.rds")
print(validation)

# Load some things to add to the plots
africa <- st_read("03_Data/02_CleanData/00_General_Africa.shp")
classe <- st_read("03_Data/02_CleanData/00_General_TrainingClasses.shp")

################################################################################
#### Reflectances
################################################################################
cols <- c("orange", "#ECCB27", "cornflowerblue")

# Histograms of Landsat & Sentinel Reflectances
p1 <- validation$Data[validation$Satellite == "Landsat"][[1]] %>%
  pivot_longer(B1:best, names_to = "Band", values_to = "Reflectance") %>%
  ggplot(aes(x = Reflectance, y = Class, col = Class, fill = Class)) +
    geom_density_ridges(alpha = 0.3) +
    facet_wrap("Band", scales = "free") +
    dark_theme_minimal() +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    ggtitle("Landsat Reflectance Profiles") +
    xlab("") +
    ylab("") +
    theme(
        panel.background = element_blank()
      , plot.background  = element_blank()
      , panel.grid.major = element_blank()
      , panel.grid.minor = element_blank()
      , axis.text.x      = element_blank()
      , axis.text.y      = element_blank()
      , legend.title     = element_blank()
      , legend.position  = "bottom"
    )

p2 <- validation$Data[validation$Satellite == "Sentinel"][[1]] %>%
  pivot_longer(B1:best, names_to = "Band", values_to = "Reflectance") %>%
  ggplot(aes(x = Reflectance, y = Class, col = Class, fill = Class)) +
    geom_density_ridges(alpha = 0.3) +
    facet_wrap("Band", scales = "free") +
    dark_theme_minimal() +
    scale_fill_manual(values = cols) +
    scale_color_manual(values = cols) +
    ggtitle("Sentinel Reflectance Profiles") +
    xlab("") +
    ylab("") +
    theme(
        panel.background = element_blank()
      , plot.background  = element_blank()
      , panel.grid.major = element_blank()
      , panel.grid.minor = element_blank()
      , axis.text.x      = element_blank()
      , axis.text.y      = element_blank()
      , legend.title     = element_blank()
      , legend.position  = "bottom"
    )

# Store the plots
p3 <- ggarrange(p1, p2, nrow = 1, common.legend = T, legend = "bottom")
ggsave("05_Presentation/Reflectances.png", plot = p3, width = 10, height = 5, scale = 1.5)

################################################################################
#### Validation
################################################################################
# Visualize confusion matrices
p4 <- validation %>% select(Satellite, Model, Confusion) %>%
  unnest(Confusion) %>%
  ggplot(aes(x = observed, y = predicted)) +
    geom_tile(fill = "orange", aes(alpha = Freq), col = "gray70", size = 0.4) +
    geom_text(aes(label = Freq), col = "white", size = 2) +
    facet_wrap(~ Satellite + Model, ncol = 4) +
    dark_theme_minimal() +
    coord_equal() +
    xlab("Observed") +
    ylab("Predicted") +
    theme(
        legend.position  = "none"
      , panel.background = element_blank()
      , plot.background  = element_blank()
      , panel.grid       = element_blank()
    )
ggsave("05_Presentation/ReflectanceValidation.png", plot = p4, width = 10, height = 5, scale = 0.8)
