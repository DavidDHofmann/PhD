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

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
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
# Histograms of Landsat & Sentinel Reflectances
p1 <- validation$Data[validation$Satellite == "Landsat"][[1]] %>%
  pivot_longer(B1:best, names_to = "Band", values_to = "Reflectance") %>%
  ggplot(aes(x = Reflectance, y = Class, col = Class, fill = Class)) +
    geom_density_ridges() +
    facet_wrap("Band", scales = "free") +
    theme_minimal() +
    theme(
        axis.text.x = element_blank()
      , axis.text.y = element_blank()
    ) +
    scale_fill_viridis_d(alpha = 0.5) +
    scale_color_viridis_d() +
    ggtitle("Landsat Reflectance Profiles") +
    theme(legend.position = "bottom")
p2 <- validation$Data[validation$Satellite == "Sentinel"][[1]] %>%
  pivot_longer(B1:best, names_to = "Band", values_to = "Reflectance") %>%
  ggplot(aes(x = Reflectance, y = Class, col = Class, fill = Class)) +
    geom_density_ridges() +
    facet_wrap("Band", scales = "free") +
    theme_minimal() +
    theme(
        axis.text.x = element_blank()
      , axis.text.y = element_blank()
    ) +
    scale_fill_viridis_d(alpha = 0.5) +
    scale_color_viridis_d() +
    ggtitle("Sentinel Reflectance Profiles") +
    theme(legend.position = "bottom")

################################################################################
#### Variable Importance
################################################################################
# Assess and plot variable importance
p3 <- subset(validation, Model == "RandomForest" & Satellite == "Landsat") %>%
  pull(ModelObject) %>%
  first() %>%
  varImp() %>%
  as.data.frame() %>%
  mutate(Band = rownames(.)) %>%
  ggplot(aes(x = reorder(Band, Overall), y = Overall)) +
    geom_segment(aes(x = reorder(Band, Overall), xend = reorder(Band, Overall), y = 0, yend = Overall), col = "gray50") +
    geom_point(col = "orange") +
    theme_minimal() +
    coord_flip() +
    xlab("Band / Index") +
    ylab("Importance") +
    ggtitle("Landsat 8 Variable Importance")
p4 <- subset(validation, Model == "RandomForest" & Satellite == "Sentinel") %>%
  pull(ModelObject) %>%
  first() %>%
  varImp() %>%
  as.data.frame() %>%
  mutate(Band = rownames(.)) %>%
  ggplot(aes(x = reorder(Band, Overall), y = Overall)) +
    geom_segment(aes(x = reorder(Band, Overall), xend = reorder(Band, Overall), y = 0, yend = Overall), col = "gray50") +
    geom_point(col = "orange") +
    theme_minimal() +
    coord_flip() +
    xlab("Band / Index") +
    ylab("Importance") +
    ggtitle("Sentinel 2 Variable Importance")

# Put the plots together
p <- ggarrange(p1, p2, p3, p4, labels = "auto")
ggsave(
    plot     = p
  , filename = "/home/david/ownCloud/University/15. PhD/Chapter_2/04_Manuscript/99_Reflectances.png"
  , bg       = "white"
  , scale    = 1.2
)

# Alternative to visualize variable importance
p5 <- validation %>%
  select(Satellite, Model, Varimp) %>%
  unnest(Varimp) %>%
  ggplot(aes(y = reorder(Band, Overall), x = Overall, xmin = 0, xmax = Overall, colour = Model)) +
    geom_pointrangeh(position = position_dodgev(height = 0.5), fatten = 1.5) +
    theme_minimal() +
    xlab("Band / Index") +
    ylab("Importance") +
    facet_wrap("Satellite", scales = "free") +
    theme(legend.position = "bottom") +
    scale_color_manual(values = c("cornflowerblue", "orange"))

################################################################################
#### Validation
################################################################################
# Visualize confusion matrices
p6 <- validation %>% select(Satellite, Model, Confusion) %>%
  unnest(Confusion) %>%
  ggplot(aes(x = observed, y = predicted, fill = Freq)) +
    geom_tile(colour = "black") +
    geom_text(aes(label = Freq), col = "black", size = 2) +
    facet_wrap(~ Satellite + Model, ncol = 4) +
    scale_fill_gradient(low = "white", high = adjustcolor("orange", alpha.f = 0.3)) +
    theme_minimal() +
    theme(legend.position = "none", panel.grid = element_blank()) +
    coord_equal() +
    xlab("Observed") +
    ylab("Predicted")
p7 <- validation %>%
  select(Satellite, Model, Specificity, Sensitivity, Accuracy) %>%
  pivot_longer(3:5, names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(y = Metric, x = Value, xmin = 0, xmax = Value)) +
    geom_vline(xintercept = 1, lty = 2, col = "gray50", lwd = 0.5) +
    geom_pointrangeh(position = position_dodgev(height = 0.5), fatten = 1.5, col = "orange") +
    theme_minimal() +
    xlab("") +
    ylab("Metric") +
    facet_wrap(~ Satellite + Model, ncol = 4) +
    theme(legend.position = "bottom") +
    scale_color_manual(values = c("cornflowerblue", "orange")) +
    scale_x_continuous(breaks = c(0, 0.5, 1))

# Put the plots together
p <- ggarrange(p6, p7, nrow = 2)

# Store them
ggsave(
    plot     = p
  , filename = "/home/david/ownCloud/University/15. PhD/Chapter_2/04_Manuscript/99_ClassificationValidation.png"
  , bg       = "white"
)
