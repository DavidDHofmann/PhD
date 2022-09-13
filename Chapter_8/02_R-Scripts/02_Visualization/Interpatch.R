################################################################################
#### Table of Inter-Patch Connectivity
################################################################################
# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_8"
setwd(wd)

# Load required packages
library(tidyverse)    # To wrangle data
library(kableExtra)   # To generate a nice table
library(sf)           # To handle spatial data
library(igraph)       # For network analysis
library(ggnetwork)    # To plot networks

# Load results on interpatch connectivity
dat <- read_rds("03_Data/03_Results/BootstrappedInterpatchConnectivity.rds")

# Load source areas
source <- read_sf("03_Data/02_CleanData/SourceAreas.shp")
source <- st_make_valid(source)

# Prepare centroids
lay <- st_point_on_surface(source)
lay <- st_coordinates(lay)

# Let's write a function that takes a vector of source and target areas and
# returns a network plot of the desired metric
networkPlot <- function(sources, targets) {

  # Subset to desired data
  sub <- select(dat, SourceArea, CurrentArea, Freq, StepNumber, FloodLevel)
  sub <- subset(sub, SourceArea %in% sources & CurrentArea %in% targets)
  sub <- subset(sub, SourceArea != CurrentArea)

  # Generate a network
  net <- graph_from_data_frame(
      d        = sub
    , vertices = unique(source$ID)
    , directed = T
  )

  # Prepare networks for ggplotting with ggplot
  net_p <- suppressWarnings(ggnetwork(net, layout = lay, arrow.gap = 0.1, scale = F))
  net_p <- subset(net_p, !is.na(FloodLevel))

  # Plot
  ggplot() +
    geom_sf(data = source) +
    geom_edges(
        data      = net_p
      , mapping   = aes(
          x    = x
        , y    = y
        , xend = xend
        , yend = yend
        , size = Freq
        , col  = StepNumber
      )
      , curvature = 0.2
      , arrow     = arrow(length = unit(6, "pt"), type = "closed", angle = 10)
    ) +
    scale_size_area(
        name     = "Frequency"
      , max_size = 1
    ) +
    scale_color_viridis_c(
        name   = "StepNumber"
      , option = "magma"
    ) +
    facet_wrap(~ FloodLevel) +
    theme_minimal()

}

# Apply it
networkPlot(1:5, 6)
networkPlot(6, 1:5)
networkPlot(6, 7:14)
networkPlot(1:6, 1:6)
networkPlot(1:6, 7:14)
networkPlot(1:14, 1:14)
networkPlot(1:6, 14:14)

################################################################################
#### Matrix Plots
################################################################################
dat %>%
  select(FloodLevel, SourceArea, CurrentArea, Freq) %>%
  ggplot(aes(x = as.factor(SourceArea), y = as.factor(CurrentArea))) +
    geom_tile(aes(fill = Freq)) +
    geom_text(aes(label = round(Freq, 2)), size = 2) +
    coord_fixed() +
    scale_fill_viridis_c(option = "magma") +
    scale_x_discrete(position = "top") +
    theme_minimal() +
    facet_wrap(~ FloodLevel)
dat %>%
  select(FloodLevel, SourceArea, CurrentArea, StepNumber) %>%
  ggplot(aes(x = as.factor(SourceArea), y = as.factor(CurrentArea))) +
    geom_tile(aes(fill = StepNumber)) +
    geom_text(aes(label = round(StepNumber)), size = 2) +
    coord_fixed() +
    scale_fill_viridis_c(option = "magma") +
    scale_x_discrete(position = "top") +
    theme_minimal() +
    facet_wrap(~ FloodLevel)

################################################################################
#### Immigration & Emigration
################################################################################
# Count immigration and emigration in both scenarios
emig <- subset(dat, SourceArea %in% 1:5 & CurrentArea %in% 6:9) %>%
  select(SourceArea, CurrentArea, Freq, FloodLevel) %>%
  pivot_wider(names_from = FloodLevel, values_from = Freq) %>%
  select(Max, Min) %>%
  colSums()
immig <- subset(dat, SourceArea %in% 6:9 & CurrentArea %in% 1:5) %>%
  select(SourceArea, CurrentArea, Freq, FloodLevel) %>%
  pivot_wider(names_from = FloodLevel, values_from = Freq) %>%
  select(Max, Min) %>%
  colSums()
cbind(emig, immig)

################################################################################
#### Frequency
################################################################################
# Reshape the data
dat_arr_freq <- dat %>%
  select(-c(StepNumber, StepNumberSE)) %>%
  mutate(Combined = paste(round(dat$Freq), "$\\pm$", sprintf("%.2f", round(dat$FreqSE, 2)))) %>%
  select(-c(Freq, FreqSE)) %>%
  pivot_wider(names_from = SourceArea, values_from = Combined) %>%
  arrange(CurrentArea, desc(FloodLevel)) %>%
  select(To = CurrentArea, everything()) %>%
  select(-FloodLevel)

# Set the diagonals to "-"
for (i in 1:nrow(dat_arr_freq)) {
  for (j in 2:ncol(dat_arr_freq)) {
    if (dat_arr_freq$To[i] == names(dat_arr_freq)[j]) {
      dat_arr_freq[i, j] <- "-"
    }
  }
}

# Prepare the table
kbl(dat_arr_freq, booktabs = T, format = "latex", escape = F, caption = "Frequency") %>%
  add_header_above(c("", "From" = 9)) %>%
  collapse_rows(1, latex_hline = c("custom"), custom_latex_hline = 1) %>%
  kable_styling(latex_options = c("scale_down")) %>%
  writeLines("04_Manuscript/IPC_Frequency.tex")

################################################################################
#### Duration
################################################################################
# Reshape the data
dat_arr_dura <- dat %>%
  select(-c(Freq, FreqSE)) %>%
  mutate(Combined = paste(round(dat$StepNumber), "$\\pm$", sprintf("%.2f", round(dat$StepNumberSE, 2)))) %>%
  select(-c(StepNumber, StepNumberSE)) %>%
  pivot_wider(names_from = SourceArea, values_from = Combined) %>%
  arrange(CurrentArea, desc(FloodLevel)) %>%
  select(To = CurrentArea, everything()) %>%
  select(-FloodLevel)

# Set the diagonals to "-"
for (i in 1:nrow(dat_arr_dura)) {
  for (j in 2:ncol(dat_arr_dura)) {
    if (dat_arr_dura$To[i] == names(dat_arr_dura)[j]) {
      dat_arr_dura[i, j] <- "-"
    }
  }
}

# Prepare the table
kbl(dat_arr_dura, booktabs = T, format = "latex", escape = F, caption = "Duration") %>%
  add_header_above(c("", "From" = 9)) %>%
  collapse_rows(1, latex_hline = c("custom"), custom_latex_hline = 1) %>%
  kable_styling(latex_options = c("scale_down")) %>%
  writeLines("04_Manuscript/IPC_Duration.tex")
