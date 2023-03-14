################################################################################
#### Distance Plot
################################################################################
# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_8"
setwd(wd)

# Load required packages
library(tidyverse)
library(colorspace)

# Reload the metrics
metrics <- read_rds("03_Data/03_Results/DistanceAOI.rds")
metrics <- mutate(metrics, LWR = Number - SE, UPR = Number + SE)
metrics <- mutate(metrics, FloodLevel = paste0(FloodLevel, "-Flood"))

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
  ) +
  xlab("") +
  ylab("Density")

# Store it
ggsave("04_Manuscript/99_HWCDifferenceAOI.png"
  , bg     = "white"
  , width  = 4
  , height = 3
  , scale  = 1.2
  , plot   = p
)
