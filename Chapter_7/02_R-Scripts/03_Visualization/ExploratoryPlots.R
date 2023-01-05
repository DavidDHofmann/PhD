################################################################################
#### Some Exploratory plots
################################################################################

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_7"
setwd(wd)

# Load required packages
library(tidyverse)   # For data wrangling
library(lemon)       # For capped axes
library(lubridate)   # To handle dates
library(hms)         # To handle times
library(latex2exp)   # For latex expressions
library(scales)      # To format times on axes
library(png)         # To import pngs
library(grid)        # To plot pngs
library(egg)         # To add png to plot
library(ggpubr)      # To arrange multiple plots

################################################################################
#### Raw Data
################################################################################
# Load activity data and add some time variables
dat <- read_csv("03_Data/02_CleanData/ActivityDataCovariates.csv")
dat$Month <- month(dat$Timestamp)
dat$Year  <- year(dat$Timestamp)
dat$Hour  <- hour(dat$TimestampRounded)
dat$Time  <- as_hms(dat$Timestamp)

# Aggregate activity by day and see if there are differences between the dogs
p1 <- dat %>%
  group_by(DogID, Date) %>%
  summarize(AverageDailyActivity = mean(Act), n = n(), .groups = "drop") %>%
  subset(n > 250) %>%
  ggplot(aes(x = DogID, y = AverageDailyActivity, fill = DogID, col = DogID)) +
    geom_boxplot(alpha = 0.3) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90), legend.position = "none") +
    scale_fill_viridis_d() +
    scale_color_viridis_d()

# Plot activity against month, regardless of the year
p2 <- ggplot(dat, aes(y = Time, x = yday(Timestamp), z = Act)) +
  stat_summary_2d(fun = "mean", bins = 100) +
  scale_fill_viridis_c(
      option   = "magma"
    , guide    = guide_colorbar(
        title          = "Activity"
      , show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.6, "cm")
      , barwidth       = unit(10.0, "cm")
    )
  ) +
  theme_minimal() +
  theme(
      legend.position = "bottom"
  ) +
  scale_y_time(breaks = date_breaks ("4 hours"), labels = label_time("%H:%M")) +
  xlim(c(0, 365)) +
  xlab("Day of Year")

# Plot activity against moonlight intensity
p3 <- ggplot(dat, aes(y = Time, x = maxMoonlightIntensity, z = Act)) +
  stat_summary_2d(fun = "mean", bins = 100) +
  scale_fill_viridis_c(
      option   = "magma"
    , guide    = guide_colorbar(
        title          = "Activity"
      , show.limits    = T
      , title.position = "top"
      , title.hjust    = 0.5
      , ticks          = T
      , barheight      = unit(0.6, "cm")
      , barwidth       = unit(10.0, "cm")
    )
  ) +
  theme_minimal() +
  theme(
      legend.position = "bottom"
    , axis.title.y    = element_blank()
    , axis.ticks.y    = element_blank()
    , axis.text.y     = element_blank()
  ) +
  scale_y_time(breaks = date_breaks ("4 hours"), labels = label_time("%H:%M")) +
  xlab(TeX("$_{max}MoonlightIntensity$"))

# Put plots together
p4 <- ggpubr::ggarrange(p2, p3
  , ncol          = 2
  , common.legend = T
  , legend        = "bottom"
  , align         = "hv"
)

# Store plots file
ggsave("04_Manuscript/99_AverageDailyActivity.png"
  , plot   = p1
  , width  = 6
  , height = 3
  , scale  = 1.5
  , bg     = "white"
)
ggsave("04_Manuscript/99_ActivityMoonlightIntensity.png"
  , plot   = p4
  , width  = 6
  , height = 6
  , scale  = 1.5
  , bg     = "white"
)

################################################################################
#### Aggregated Data
################################################################################
# Load aggregated activity data
dat_byphase <- read_csv("03_Data/02_CleanData/ActivityDataCovariatesAggregated.csv")
dat_byphase <- mutate(dat_byphase
  , ToD  = factor(ToD, levels = c("Morning", "Evening", "Night"))
)

# Define classes
levels <- c("Very Low", "Low", "Medium", "High", "Very High")
classes <- data.frame(
    MoonBinNumeric = 1:5
  , MoonClass      = factor(levels, levels = levels)
)

# Let's bin the nights by moon illumination
dat_byphase <- dat_byphase %>%
  mutate(MoonBin = cut(maxMoonlightIntensity, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), include.lowest = T)) %>%
  mutate(MoonBinNumeric = as.numeric(MoonBin)) %>%
  left_join(., classes, by = "MoonBinNumeric")

# Visualize the relationship between activity and moon-illumination
p1 <- dat_byphase %>%
  group_by(MoonClass, ToD) %>%
  summarize(
      n       = n()
    , mean    = round(mean(meanAct))
    , .groups = "drop"
  ) %>%
  pivot_longer(n:mean, names_to = "Metric", values_to = "Value") %>%
  ggplot(aes(x = MoonClass, y = Metric, label = Value)) +
    geom_text(size = 2) +
    facet_wrap(~ToD) +
    theme_void() +
    theme(axis.text.y = element_text())

p2 <- ggplot(dat_byphase, aes(x = MoonClass, y = meanAct, col = ToD, fill = ToD)) +
  geom_boxplot(alpha = 0.2, show.legend = F) +
  theme_minimal() +
  facet_wrap(~ ToD) +
  ylab("Average Activity") +
  xlab("Maximum Moon Illumination at Night") +
  scale_color_viridis_d(begin = 0.2, end = 0.8) +
  scale_fill_viridis_d(begin = 0.2, end = 0.8, alpha = 0.2) +
  theme(axis.text.x = element_text(angle = 45))

# Load pngs to overlay
imgs <- c(
    "04_Manuscript/TimeOfDay_Morning.png"
  , "04_Manuscript/TimeOfDay_Evening.png"
  , "04_Manuscript/TimeOfDay_Night.png"
  ) %>% lapply(function(x) {
     readPNG(x) %>% rasterGrob(., width = 0.3, hjust = 1, vjust = 0)
})

# Put all into its own dataframe
pngs <- tibble(
    ToD  = unique(dat_byphase$ToD)
  , grob = imgs
  , x    = 5.75
  , y    = 380
)

# Add them to the plot
p3 <- p2 + geom_custom(aes(data = grob, x = x, y = y)
  , data     = pngs
  , grob_fun = identity
)

# Put plots together
p4 <- ggpubr::ggarrange(p3, p1, nrow = 2, heights = c(0.75, 0.25), align = "hv")

# Save to file
# Store plot to file
ggsave("04_Manuscript/99_ActivityByMoonlight.png"
  , plot   = p4
  , width  = 6
  , height = 4
  , scale  = 1.5
  , bg     = "white"
)
