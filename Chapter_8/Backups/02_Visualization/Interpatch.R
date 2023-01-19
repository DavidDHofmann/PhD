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
library(ggspatial)    # For scalebar and north arrow
library(igraph)       # For network analysis
library(ggnetwork)    # To plot networks
library(terra)        # To handle raster files
library(RColorBrewer) # For custom colors
library(ggpubr)       # To arrange multiple plots

# Load results on interpatch connectivity
dat <- read_rds("03_Data/03_Results/BootstrappedInterpatchConnectivity.rds")

# Load source areas and the two watermaps
areas <- read_sf("03_Data/02_CleanData/SourceAreas.shp")
areas <- st_make_valid(areas)
water  <- rast("03_Data/02_CleanData/WaterCover.tif")[[c("min", "max")]]
water  <- crop(water, ext(21.5, 24.3, -20.7, -18))
water  <- as.data.frame(water, xy = T) %>% pivot_longer(3:4
  , names_to  = "FloodLevel"
  , values_to = "Water"
)
water$FloodLevel <- ifelse(water$FloodLevel == "min", "Min", "Max")
water$FloodLevel <- factor(water$FloodLevel, levels = c("Min", "Max"))
dat$FloodLevel   <- factor(dat$FloodLevel, levels = c("Min", "Max"))

# Prepare centroids
lay <- st_point_on_surface(areas)
lay <- st_coordinates(lay)

# Generate area labels
labels_areas <- st_coordinates(st_point_on_surface(areas))
labels_areas <- cbind(labels_areas, st_drop_geometry(areas))

# Create facet labels
lab <- data.frame(
    FloodLevel = factor(c("Min", "Max"), levels = c("Min", "Max"))
  , Label      = c("a", "b")
  , x = c(21.5, 21.5)
  , y = c(-17.75, -17.75)
)

# Let's write a function that takes a vector of source and target areas and
# returns a network plot of the desired metric
networkPlot <- function(sources, targets, reverse = F, curvature = 0.2, labels = F) {

  # Subset to desired data
  sub <- select(dat, SourceArea, CurrentArea, Freq, StepNumber, FloodLevel)
  subi <- subset(sub, SourceArea %in% sources & CurrentArea %in% targets)
  if (reverse) {
    subi2 <- subset(sub, SourceArea %in% targets & CurrentArea %in% sources)
    subi <- rbind(subi, subi2)
  }
  sub <- subi
  sub <- subset(sub, SourceArea != CurrentArea)

  # Generate a network
  net <- graph_from_data_frame(
      d        = sub
    , vertices = unique(areas$ID)
    , directed = T
  )

  # Prepare networks for ggplotting with ggplot
  net_p <- suppressWarnings(ggnetwork(net, layout = lay, arrow.gap = 0.1, scale = F))
  net_p <- subset(net_p, !is.na(FloodLevel))
  net_p$FloodLevel <- factor(net_p$FloodLevel, levels = c("Min", "Max"))

  # Plot
  ggplot() +
    geom_raster(data = water, aes(x = x, y = y, fill = as.factor(Water))) +
    # geom_sf(
    #     data  = subset(areas, Type == "Main")
    #   , col   = "orange"
    #   , fill  = "orange"
    #   , alpha = 0.2
    # ) +
    # geom_sf(
    #     data  = subset(areas, Type == "Buffer")
    #   , col   = "purple"
    #   , fill  = "purple"
    #   , alpha = 0.2
    # ) +
    geom_sf(
        data  = areas
      , col   = "gray"
      , fill  = "gray"
      , alpha = 0.2
    ) +
    geom_point(
        data     = subset(labels_areas, Type == "Buffer")
      , mapping  = aes(x = X, y = Y)
      , col      = "gray"
      , size     = 4
    ) +
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
      , curvature = curvature
      , arrow     = arrow(length = unit(6, "pt"), type = "closed", angle = 10)
    ) +
    geom_text(
        data     = labels_areas
      , mapping  = aes(x = X, y = Y, label = ID)
      , col      = "black"
      , fontface = 3
      , size     = 2
    ) +
    {
      if (labels) {
        geom_text(
            data     = lab
          , aes(x    = x, y = y, label = Label)
          , fontface = 2
          , size     = 5
          , col      = "black"
        )
      }
    } +
    scale_fill_manual(
        values = c("transparent", "cornflowerblue")
      , name   = ""
    ) +
    scale_size_area(
        name     = "Frequency"
      , max_size = 1
    ) +
    scale_color_gradientn(
      colors = rev(brewer.pal(11, "RdYlGn"))
    ) +
    facet_wrap(~ FloodLevel) +
    guides(
        fill  = "none"
      , size  = guide_legend(order = 2)
    ) +
    theme(
        panel.background = element_blank()
      , panel.border     = element_rect(colour = "black", fill = NA, size = 1)
    ) +
    annotation_scale(
        location   = "bl"
      , width_hint = 0.2
      , line_width = 0.5
      , height     = unit(0.15, "cm")
    ) +
    annotation_north_arrow(
        location = "br"
      , height   = unit(1.5, "cm"),
      , width    = unit(1.2, "cm"),
      , style    = north_arrow_fancy_orienteering(
            fill      = c("black", "black")
          , line_col  = NA
          , text_col  = "black"
          , text_size = 12
        )
    ) +
    xlab("") +
    ylab("")
}

# Apply it
# networkPlot(1:5, 6)
# networkPlot(6, 1:5)
# networkPlot(6, 1:5, reverse = T)
# networkPlot(6, 7:14, curvature = 0)
# networkPlot(1:6, 1:6)
# networkPlot(1:6, 7:14)
# networkPlot(1:14, 1:14)
# networkPlot(2, c(1, 3), reverse = T)

# Let's store one of the plots
p1 <- networkPlot(6, 1:5, reverse = T, labels = T)

# Let's generate the same plot for any other source area
p2 <- list()
for (i in 1:6) {
  p2[[i]] <- networkPlot(i, base::setdiff(1:6, i), reverse = T)
}

# Put the plots together
p3 <- ggarrange(p2[[1]], p2[[2]], p2[[3]], ncol = 1, labels = "auto")
p4 <- ggarrange(p2[[4]], p2[[5]], p2[[6]], ncol = 1, labels = "auto")
p5 <- ggarrange(p3, p4, ncol = 2)

# And repeat the same for the borders
p6 <- list()
for (i in 1:6) {
  p6[[i]] <- networkPlot(i, 7:14, curvature = 0)
}

# Put the plots together
p7 <- ggarrange(p6[[1]], p6[[2]], p6[[3]], ncol = 1, labels = "auto")
p8 <- ggarrange(p6[[4]], p6[[5]], p6[[6]], ncol = 1, labels = "auto")
p9 <- ggarrange(p7, p8, ncol = 2)

################################################################################
#### Store Plots
################################################################################
ggsave("04_Manuscript/99_Interpatch.png"
  , plot   = p1
  , bg     = "white"
  , height = 4
  , width  = 8
  , scale  = 1
)
ggsave("04_Manuscript/99_IPCMain.png"
  , plot   = p5
  , bg     = "white"
  , height = 6
  , width  = 8
  , scale  = 1.8
)
ggsave("04_Manuscript/99_IPCBuffer.png"
  , plot   = p9
  , bg     = "white"
  , height = 6
  , width  = 8
  , scale  = 1.8
)

################################################################################
#### Matrix Plots
################################################################################
# dat %>%
#   select(FloodLevel, SourceArea, CurrentArea, Freq) %>%
#   ggplot(aes(x = as.factor(SourceArea), y = as.factor(CurrentArea))) +
#     geom_tile(aes(fill = Freq)) +
#     geom_text(aes(label = round(Freq, 2)), size = 2) +
#     coord_fixed() +
#     scale_fill_gradientn(
#       colors = brewer.pal(11, "RdYlGn")
#     ) +
#     scale_x_discrete(position = "top") +
#     theme_minimal() +
#     facet_wrap(~ FloodLevel)
# dat %>%
#   select(FloodLevel, SourceArea, CurrentArea, StepNumber) %>%
#   ggplot(aes(x = as.factor(SourceArea), y = as.factor(CurrentArea))) +
#     geom_tile(aes(fill = StepNumber)) +
#     geom_text(aes(label = round(StepNumber)), size = 2) +
#     coord_fixed() +
#     scale_fill_gradientn(
#       colors = rev(brewer.pal(11, "RdYlGn"))
#     ) +
#     scale_x_discrete(position = "top") +
#     theme_minimal() +
#     facet_wrap(~ FloodLevel)

# Number of successful dispersers
p1 <- dat %>%
  mutate(Label = paste(round(dat$Freq), "%+-%", sprintf("%.2f", round(dat$FreqSE, 2)))) %>%
  select(FloodLevel, SourceArea, CurrentArea, Label, Freq) %>%
  mutate(FloodLevel = factor(FloodLevel, levels = c("Max", "Min"))) %>%
  ggplot(aes(x = as.factor(SourceArea), y = as.factor(FloodLevel), fill = Freq)) +
    geom_tile(col = "black") +
    geom_text(aes(label = Label), size = 2, parse = T) +
    scale_x_discrete(position = "top") +
    scale_fill_gradientn(
        colors  = adjustcolor(brewer.pal(11, "RdYlGn"), alpha.f = 0.75)
      , guide   = guide_colorbar(
        , title          = "Frequency"
        , show.limits    = T
        , title.position = "bottom"
        , title.hjust    = 0.5
        , ticks          = F
        , barheight      = unit(0.2, "cm")
        , barwidth       = unit(5, "cm")
      )
    ) +
    facet_grid("CurrentArea", switch = "y") +
    theme_minimal() +
    xlab("From") +
    ylab("To") +
    theme(
        strip.placement   = "outside"
      , strip.text.y.left = element_text(angle = 0)
      , axis.title.x      = element_text(angle = 0, vjust = 0.5, size = 5)
      , axis.title.y      = element_text(angle = 0, vjust = 0.5, size = 5)
      , axis.text.y       = element_text(size = 5)
      , legend.position   = "none"
      , legend.box        = "vertical"
      , panel.grid.minor  = element_blank()
      , panel.grid.major  = element_blank()
      , panel.spacing.y   = unit(0.05, "lines")
    )

# Number of steps required
p2 <- dat %>%
  mutate(Label = paste(round(dat$StepNumber), "%+-%", sprintf("%.2f", round(dat$StepNumberSE, 2)))) %>%
  select(FloodLevel, SourceArea, CurrentArea, Label, StepNumber) %>%
  mutate(FloodLevel = factor(FloodLevel, levels = c("Max", "Min"))) %>%
  ggplot(aes(x = as.factor(SourceArea), y = as.factor(FloodLevel), fill = StepNumber)) +
    geom_tile(col = "black") +
    geom_text(aes(label = Label), size = 2, parse = T) +
    scale_x_discrete(position = "top") +
    scale_fill_gradientn(
        colors  = adjustcolor(rev(brewer.pal(11, "RdYlGn")), alpha.f = 0.75)
      , guide   = guide_colorbar(
        , title          = "StepNumber"
        , show.limits    = T
        , title.position = "bottom"
        , title.hjust    = 0.5
        , ticks          = F
        , barheight      = unit(0.2, "cm")
        , barwidth       = unit(5, "cm")
      )
    ) +
    facet_grid("CurrentArea", switch = "y") +
    theme_minimal() +
    xlab("From") +
    ylab("") +
    theme(
        strip.placement   = "outside"
      , strip.text.y.left = element_text(angle = 0)
      , axis.title.x      = element_text(angle = 0, vjust = 0.5, size = 5)
      , axis.title.y      = element_text(angle = 0, vjust = 0.5, size = 5)
      , axis.text.y       = element_text(size = 5)
      , legend.position   = "none"
      , legend.box        = "vertical"
      , panel.grid.minor  = element_blank()
      , panel.grid.major  = element_blank()
      , panel.spacing.y   = unit(0.05, "lines")
    )

# Arrange plots and store them
p3 <- ggarrange(p1, p2, labels = "auto")
# p3 <- p3 + coord_fixed(ratio = 0.6)
ggsave("04_Manuscript/99_IPCTable.png"
  , plot   = p3
  , bg     = "white"
  , scale  = 1.5
  , width  = 6
  , height = 4
)

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
  add_header_above(c("", "From" = 6)) %>%
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
  add_header_above(c("", "From" = 6)) %>%
  collapse_rows(1, latex_hline = c("custom"), custom_latex_hline = 1) %>%
  kable_styling(latex_options = c("scale_down")) %>%
  writeLines("04_Manuscript/IPC_Duration.tex")
