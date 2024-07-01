################################################################################
#### IPC Summary Table
################################################################################
# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_2"
setwd(wd)

# Load required packages
library(tidyverse)    # To wrangle data
library(kableExtra)   # To generate a nice table
library(ggpubr)       # To arrange multiple plots
library(colorspace)   # To manipulate colors

# Reload results on interpatch connectivity
dat <- read_rds("03_Data/03_Results/InterpatchConnectivityBootstrapped.rds")

# Function for plotting
plotTable <- function(data, xlab = "Source Area", ylab = "Current Area", scalename = "Metric", rev = F, days = T) {

  # Some people feel strong about having the dispersal duration in days, so here
  # we go
  if (days) {
    data$DispersalDuration   <- data$DispersalDuration / 5
    data$DispersalDurationSD <- data$DispersalDurationSD / 5
  }

  # Prepare labels
  toplot <- data %>%
    mutate(Label = paste(round(Metric), "%+-%", sprintf("%.2f", round(MetricSE, 2)))) %>%
    mutate(Label = ifelse(CurrentArea == SourceArea & SourceArea != "Overall" & CurrentArea != "Overall", NA, Label)) %>%
    select(FloodLevel, SourceArea, CurrentArea, Label, Metric) %>%
    mutate(FloodLevel = factor(FloodLevel, levels = c("Max", "Min"))) %>%
    mutate(CurrentArea = factor(CurrentArea, levels = c(as.character(1:14), "Overall"))) %>%
    mutate(SourceArea = factor(SourceArea, levels = c(as.character(1:6), "Overall")))

  # Identify minimum and maximum values
  minimum <- min(toplot$Metric[!is.infinite(toplot$Metric)], na.rm = T)
  maximum <- max(toplot$Metric[!is.infinite(toplot$Metric)], na.rm = T)
  mapping <- scales::rescale(seq(minimum, maximum, length.out = 20) ** 2)
  colors  <- hcl.colors(11, "RdYlGn", rev = rev)

  # Prepare the plot
  p <- ggplot(toplot, aes(x = 1, y = FloodLevel, fill = Metric)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = Label), size = 2, parse = T) +
    scale_fill_gradientn(
        colors   = adjustcolor(colors, alpha.f = 0.75)
      , na.value = "gray95"
      , values   = mapping
      , breaks   = seq(minimum, maximum, length.out = 3)
      , limits   = c(minimum, maximum)
      , labels   = c("Low", "Medium", "High")
      , guide    = guide_colorbar(
        , title          = scalename
        , show.limits    = T
        , title.position = "bottom"
        , title.hjust    = 0.5
        , ticks          = F
        , barheight      = unit(0.2, "cm")
        , barwidth       = unit(5, "cm")
      )
    ) +
    facet_grid(CurrentArea ~ SourceArea, switch = "y") +
    theme(
        legend.position  = "bottom"
      , strip.background = element_rect(fill = "gray95")
      , strip.placement  = "outside"
      , axis.text.y      = element_text(size = 5)
      , axis.text.x      = element_blank()
      , axis.ticks.x     = element_blank()
      , panel.grid       = element_blank()
      , panel.background = element_rect(fill = "transparent")
      , panel.spacing    = unit(0.05, "cm")
    ) +
    xlab(xlab) +
    ylab(ylab) +
    scale_x_continuous(position = "top")

  # Return the plot
  return(p)
}

################################################################################
#### Dispersal Success and Duration
################################################################################
# Get data to plot dispersal success
subdat <- dat %>%
  subset(Type == "Dispersal") %>%
  select(FloodLevel, SourceArea, CurrentArea, Metric = DispersalSuccess, MetricSE = DispersalSuccessSD)

# Plot it
p1 <- plotTable(subdat
  , scalename = "Number of Trajectories"
  , rev       = F
  , xlab      = "Source Area"
  , ylab      = "Target Area"
)

# Get data to plot dispersal durations
subdat <- dat %>%
  subset(Type == "Dispersal") %>%
  select(FloodLevel, SourceArea, CurrentArea, Metric = DispersalDuration, MetricSE = DispersalDurationSD)

# Plot it
p2 <- plotTable(subdat
  , scalename = "Dispersal Duration (days)"
  , rev       = T
  , xlab      = "Source Area"
  , ylab      = "Target Area"
)

################################################################################
#### Egression Success and Duration
################################################################################
# Get data to plot egression success
subdat <- dat %>%
  subset(Type == "Egression") %>%
  select(FloodLevel, SourceArea, CurrentArea, Metric = DispersalSuccess, MetricSE = DispersalSuccessSD)

# Plot it
p3 <- plotTable(subdat
  , scalename = "Number of Trajectories"
  , rev       = F
  , xlab      = "Source Area"
  , ylab      = "Egression Zone"
)

# Get data to plot egression durations
subdat <- dat %>%
  subset(Type == "Egression") %>%
  select(FloodLevel, SourceArea, CurrentArea, Metric = DispersalDuration, MetricSE = DispersalDurationSD)

# Plot it
p4 <- plotTable(subdat
  , scalename = "Egression Duration (days)"
  , rev       = T
  , xlab      = "Source Area"
  , ylab      = "Egression Zone"
)

# Combine plots
p  <- ggarrange(p1, p3, p2, p4, ncol = 2, nrow = 2, labels = c("a1", "a2", "b1", "b2"))
ggsave(plot = p, filename = "04_Manuscript/Figures/InterpatchConnectivityTable.png", scale = 2, height = 5, width = 6, device = png, bg = "transparent")

################################################################################
#### Emigration / Immigration Bar Chart
################################################################################
# Prepare data for immigration and emigration plots
dat_sub <- dat %>%
  select(FloodLevel, SourceArea, CurrentArea, DispersalSuccess, DispersalSuccessSD) %>%
  subset(SourceArea %in% 1:6 & CurrentArea %in% 1:6) %>%
  subset(SourceArea != CurrentArea) %>%
  mutate(FloodLevel = factor(FloodLevel, levels = c("Min", "Max")))
dat_emi <- select(dat_sub
    , FocalArea     = SourceArea
    , OtherArea     = CurrentArea
    , FloodLevel    = FloodLevel
    , Freq          = DispersalSuccess
    , FreqSE        = DispersalSuccessSD
  ) %>%
  mutate(Type = "Emigration") %>%
  mutate(LWR = Freq - 1.96 * FreqSE, UPR = Freq + 1.96 * FreqSE)
dat_imi <- select(dat_sub
    , FocalArea     = CurrentArea
    , OtherArea     = SourceArea
    , FloodLevel    = FloodLevel
    , Freq          = DispersalSuccess
    , FreqSE        = DispersalSuccessSD
  ) %>%
  mutate(Type = "Immigration") %>%
  mutate(LWR = Freq - 1.96 * FreqSE, UPR = Freq + 1.96 * FreqSE)
dat_sub <- rbind(dat_emi, dat_imi) %>%
  mutate(FocalArea = paste0("Focal Area ", FocalArea)) %>%
  mutate(FloodLevel = factor(paste0(FloodLevel, "-Flood"), levels = c("Min-Flood", "Max-Flood")))

# Prepare the plot
cols <- hcl.colors(n = 6, palette = "viridis")
p <- dat_sub %>%
  ggplot(aes(x = FloodLevel, y = Freq, ymin = LWR, ymax = UPR, fill = as.factor(OtherArea))) +
    geom_bar(stat = "identity"
      , position  = position_dodge(width = 0.9)
      , col       = "white"
      , linewidth = 0.2
    ) +
    geom_errorbar(aes(col = as.factor(OtherArea))
      , position    = position_dodge(width = 0.9)
      , width       = 0.3
      , show.legend = F
    ) +
    # facet_grid(Type ~ FocalArea, switch = "x") +
    facet_grid(Type ~ FocalArea) +
    # theme_minimal() +
    # scale_color_manual(values = darken(cols, 0.2)) +
    scale_color_manual(values = lighten(cols, 0.2)) +
    scale_fill_manual(
        values = cols
      , guide  = guide_legend(
        , title          = "Target-Area / Source-Area"
        , show.limits    = T
        , title.position = "bottom"
        , title.hjust    = 0.5
        , label.position = "bottom"
        , ticks          = F
        , nrow           = 1
      )
    ) +
    theme(
        legend.position   = "bottom"
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
    ylab("Frequency")

# Store it
ggsave("04_Manuscript/Figures/ImmigrationEmigration.png"
  , plot   = p
  , bg     = "white"
  , scale  = 1.5
  , width  = 6
  , height = 3
  , device = png
)
