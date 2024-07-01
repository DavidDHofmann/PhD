################################################################################
#### Distance Plot
################################################################################
# Description: Visualization of the distance to AOI

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Load required packages
library(tidyverse)    # For wrangling data
library(colorspace)   # To change colors

# Reload the metrics
metrics <- read_rds("03_Data/03_Results/DistanceAOI.rds")
metrics <- mutate(metrics, LWR = Number - SE, UPR = Number + SE)
metrics <- mutate(metrics, FloodLevel = paste0(FloodLevel, "-Flood"))
metrics <- mutate(metrics, FloodLevel = factor(FloodLevel, levels = c("Min-Flood", "Max-Flood")))

# Visualize
cols <- c("orange", "cornflowerblue")
p <- ggplot(metrics, aes(x = FloodLevel, y = Number, ymin = LWR, ymax = UPR, fill = FloodLevel)) +
  geom_col(color = "white", alpha = 0.8) +
  geom_errorbar(width = 0.2, aes(color = FloodLevel)) +
  facet_grid(~ Name) +
  scale_color_manual(name = "Flood-Level", values = lighten(cols, 0.5)) +
  scale_fill_manual(
      values = cols
    , guide  = guide_legend(
      , title          = "Flood-Level"
      , show.limits    = T
      , title.position = "bottom"
      , title.hjust    = 0.5
      , label.position = "bottom"
      , ticks          = F
      , nrow           = 1
    )
  ) +
  theme(
      legend.position   = "none"
    , legend.key.height = unit(0.2, 'cm')
    , legend.key.width  = unit(1, 'cm')
    , panel.grid.major  = element_line(color = "gray95")
    , panel.grid.minor  = element_blank()
    , axis.text         = element_text(size  = 6)
    , axis.title        = element_text(size  = 9)
    , panel.background  = element_blank()
    , strip.background  = element_rect(fill = "gray95", color = "transparent")
  ) +
  xlab("") +
  ylab("Density")

# Store it
ggsave("04_Manuscript/Figures/HWCDifferenceAOI.png"
  , bg     = "white"
  , width  = 4
  , height = 3
  , scale  = 1.2
  , plot   = p
  , device = png
)
