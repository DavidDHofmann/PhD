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
library(igraph)       # For network analysis
library(ggnetwork)    # To plot networks
library(RColorBrewer) # For custom colors
library(ggpubr)       # To arrange multiple plots
library(colorspace)   # To lighten or darken colors

# Load results on interpatch connectivity
dat <- read_rds("03_Data/03_Results/BootstrappedInterpatchConnectivity.rds")

# Prepare data for immigration and emigration plots
dat_sub <- dat %>%
  select(FloodLevel, SourceArea, CurrentArea, Freq, FreqSE) %>%
  subset(CurrentArea %in% 1:6) %>%
  subset(SourceArea != CurrentArea) %>%
  mutate(FloodLevel = factor(FloodLevel, levels = c("Min", "Max")))
dat_emi <- select(dat_sub
    , FocalArea     = SourceArea
    , OtherArea     = CurrentArea
    , FloodLevel    = FloodLevel
    , Freq          = Freq
    , FreqSE        = FreqSE
  ) %>%
  mutate(Type = "Emigration\n(into)") %>%
  mutate(LWR = Freq - 1.96 * FreqSE, UPR = Freq + 1.96 * FreqSE)
dat_imi <- select(dat_sub
    , FocalArea     = CurrentArea
    , OtherArea     = SourceArea
    , FloodLevel    = FloodLevel
    , Freq          = Freq
    , FreqSE        = FreqSE
  ) %>%
  mutate(Type = "Immigration\n(from)") %>%
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
    ) +
    xlab("") +
    ylab("Frequency")

# Store it
ggsave("04_Manuscript/99_ImmigrationEmigration.png"
  , plot   = p
  , bg     = "white"
  , scale  = 1.5
  , width  = 6
  , height = 3
)

################################################################################
#### Matrix Plots
################################################################################
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
      , breaks  = seq(min(dat$Freq), max(dat$Freq), length.out = 3)
      , labels  = c("Low", "Medium", "High")
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
      , legend.position   = "bottom"
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
      , breaks  = seq(min(dat$StepNumber), max(dat$StepNumber), length.out = 3)
      , labels  = c("Low", "Medium", "High")
      , guide   = guide_colorbar(
        , title          = "Step Number"
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
      , legend.position   = "bottom"
      , legend.box        = "vertical"
      , panel.grid.minor  = element_blank()
      , panel.grid.major  = element_blank()
      , panel.spacing.y   = unit(0.05, "lines")
    )

# Arrange plots
p3 <- ggarrange(p1, p2, labels = "auto", ncol = 2)

# Store
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

# ################################################################################
# #### Backup
# ################################################################################
# # Add second legend
# df <- data.frame(x = c(0, 1), y = c(0, 1), col = factor(c("Source Area", "Emigration Zone"), levels = c("Source Area", "Emigration Zone")))
# l1 <- ggplot(df, aes(x = x, y = y, col = col)) +
#   geom_point(size = 5, alpha = 0.2) +
#   scale_color_manual(
#       values = c("orange", "purple")
#     , guide  = guide_legend(
#       , title          = ""
#       , label.position = "bottom"
#   )) +
#   theme(
#       legend.position = "bottom"
#     , legend.key      = element_blank()
#   )
# l1 <- get_legend(l1)
# p4 <- p3 + annotation_custom(
#     grob = l1
#     , xmin = 0.25
#     , xmax = 0.75
#     , ymin = 0.0
#     , ymax = 0.1
# )
#
# # Add some highlights
# highlight1 <- ggplot(data.frame(x = 1:6, y = rep(1, 6)), aes(x = x, y = y)) +
#   geom_point(size = 10, alpha = 0.2, col = rep("orange", 6)) +
#   theme_void() +
#   coord_cartesian(clip = "off")
# highlight2 <- ggplot(data.frame(x = rep(1, 14), y = 1:14), aes(x = x, y = y)) +
#   geom_point(size = 10, alpha = 0.2, col = c(rep("purple", 8), rep("orange", 6))) +
#   theme_void() +
#   coord_cartesian(clip = "off")
# p5 <- p4 + annotation_custom(
#     grob = ggplotGrob(highlight1)
#   , xmin = 0.123
#   , xmax = 0.465
#   , ymin = 0.93
#   , ymax = 1.0
#   ) + annotation_custom(
#     grob = ggplotGrob(highlight1)
#   , xmin = 0.123 + 0.488
#   , xmax = 0.465 + 0.5
#   , ymin = 0.93
#   , ymax = 1.0
#   ) + annotation_custom(
#     grob = ggplotGrob(highlight2)
#   , xmin = 0
#   , xmax = 0.095
#   , ymin = 0.115
#   , ymax = 0.958
#   ) + annotation_custom(
#     grob = ggplotGrob(highlight2)
#   , xmin = 0 + 0.488
#   , xmax = 0.095 + 0.488
#   , ymin = 0.115
#   , ymax = 0.958
# )
