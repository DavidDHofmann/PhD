################################################################################
#### Valid Steps
################################################################################
# Description: Plot of number of valid steps

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)     # For data wrangling
library(ggpubr)        # To arrange plots
library(ggdark)        # For dark themes

# Set working directory
setwd("/home/david/ownCloud/University/15. PhD/Chapter_4")

# Load simulated data
dat <- read_csv("03_Data/ValidSteps.csv")
dat <- subset(dat, Forgiveness == 1)

# Compute summary stats from the simulations
summaries <- dat %>%
  group_by(Missingness, Forgiveness) %>%
  summarize(
      LCI         = quantile(NumberSteps, 0.025)
    , UCI         = quantile(NumberSteps, 0.975)
    , NumberSteps = mean(NumberSteps)
    , .groups     = "drop"
  )

# Example trajectory
pts <- data.frame(rbind(
    c(0, 0)
  , c(5, 3)
  , c(10, -1)
  , c(15, 4)
  , c(12, 7)
  , c(7, 8)
))
names(pts) <- c("x", "y")

# Rescale
pts$x <- pts$x / (max(pts$x) / max(summaries$Missingness))
pts$y <- pts$y / (max(pts$y) / max(summaries$NumberSteps))

# Plot one
p1 <- ggplot(summaries, aes(x = Missingness, y = NumberSteps / 1000, ymin = LCI / 1000, ymax = UCI / 1000)) +
  geom_ribbon(alpha = 0.5, fill = "orange", color = "orange", linewidth = 0.2) +
  geom_line(size = 0.4, col = "white") +
  dark_theme_minimal() +
  ylab("Fraction of usable steps") +
  theme(
      legend.position  = "bottom"
    , panel.grid.minor = element_blank()
    , panel.grid.major = element_line(linewidth = 0.2, color = "gray50")
    , panel.background = element_blank()
    , plot.background  = element_blank()
  )

# Plot two
p2 <- ggplot(pts, aes(x = x, y = y)) +
  geom_point(col = "orange", alpha = 0.5, size = 5) +
  dark_theme_minimal() +
  theme(
      legend.position  = "bottom"
    , panel.grid.minor = element_blank()
    , panel.grid.major = element_line(linewidth = 0.2, color = "gray50")
    , panel.background = element_blank()
    , plot.background  = element_blank()
  )

# Arrange plots
p <- ggarrange(p2, p1, nrow = 2)

# Store plot
ggsave("05_Presentation/PlotsForFieberg.png"
  , plot   = p
  , width  = 6
  , height = 5
  , bg     = "transparent"
  , scale  = 1
)
