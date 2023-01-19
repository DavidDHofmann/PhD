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

# Load results on interpatch connectivity
dat <- read_rds("03_Data/03_Results/BootstrappedInterpatchConnectivity.rds")

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

# Arrange plots
p3 <- ggarrange(p1, p2, labels = "auto")

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
