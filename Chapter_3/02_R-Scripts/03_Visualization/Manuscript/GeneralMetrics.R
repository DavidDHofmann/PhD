################################################################################
#### A Few General Metrics on the Dispersal Events
################################################################################
# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Load required packages
library(tidyverse)    # For plotting and data wrangling
library(lubridate)    # To handle dates
library(hms)          # To handle times
library(ggpubr)       # To arrange multiple plots

# Load custom functions
source("02_R-Scripts/00_Functions.R")

# Load the raw and cleaned gps data (for some metrics we'll use the raw data,
# for others we'll use the cleaned one)
dat_raw   <- read_csv("03_Data/02_CleanData/Dispersers.csv")
dat_clean <- read_rds("03_Data/02_CleanData/SSFExtracted.rds") %>% subset(case == 1)

# Create overview plot
p1 <- gpsOverview(dat_raw, id = "ID", timestamp = "Timestamp")
p2 <- gpsOverview(dat_clean, id = "ID", timestamp = "Timestamp")
ggarrange(p1, p2, ncol = 1)

################################################################################
#### Name and Sex of each Individual
################################################################################
# Sex
p1 <- dat_clean %>%
  select(ID, Sex) %>%
  distinct() %>%
    ggplot(aes(x = 1, y = ID, shape = Sex)) +
      geom_point(size = 5) +
      scale_shape_manual(values = c("\u2640", "\u2642")) +
      theme_awesome() +
      theme(
          legend.position  = "none"
        , panel.grid.minor = element_blank()
        , axis.text.x      = element_blank()
        , axis.title.x     = element_blank()
      ) +
      facet_wrap(~ "Sex")

################################################################################
#### Number of Fixes by Season
################################################################################
# Calculate the number of fixes per dog
number_fixes <- dat_clean %>%
  count(ID, SeasonClimate, name = "NumberOfFixes")
mean_number_fixes <- mean(number_fixes$NumberOfFixes)

# Prepare a plot
p2 <- ggplot(number_fixes, aes(x = NumberOfFixes, y = ID, fill = ID, color = ID, alpha = SeasonClimate)) +
  geom_col() +
  geom_vline(xintercept = mean_number_fixes, linetype = "22", col = "gray30") +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  scale_x_continuous(labels = function(x) {format(x, big.mark = ",")}) +
  scale_alpha_manual(values = c(0.5, 1)) +
  xlab("Number of Fixes") +
  facet_wrap(~ "Number of Fixes") +
  theme_awesome() +
  theme(
      legend.position  = "none"
    , panel.grid.minor = element_blank()
    , axis.text.y      = element_blank()
    , axis.title.y     = element_blank()
  )

################################################################################
#### Average Daily Distance Dispersed
################################################################################
# Compute the average daily distance dispersed
average_distance <- dat_clean %>%
  mutate(Date = as.Date(Timestamp)) %>%
  group_by(ID, Date) %>%
  summarize(Distance = mean(sl) * 5, .groups = "drop") %>%
  group_by(ID) %>%
  summarize(MeanDistance = mean(Distance), SDDistance = sd(Distance))
mean_average_distance <- mean(average_distance$MeanDistance)

p3 <- ggplot(average_distance, aes(x = MeanDistance, y = ID, fill = ID, color = ID)) +
  geom_col(alpha = 0.75) +
  geom_errorbarh(aes(xmin = MeanDistance - SDDistance, xmax = MeanDistance + SDDistance), height = 0) +
  geom_vline(xintercept = mean_average_distance, linetype = "22", col = "gray30") +
  scale_fill_viridis_d() +
  scale_color_viridis_d() +
  xlab("Average Daily Distance Dispersed (km)") +
  theme_awesome() +
  theme(
      legend.position  = "none"
    , panel.grid.minor = element_blank()
    , axis.text.y      = element_blank()
    , axis.title.y     = element_blank()
  ) +
  facet_wrap(~ "Daily Distance")

################################################################################
#### Arrange the Plots
################################################################################
ggarrange(p1, p2, p3, nrow = 1, align = "h")

mean_sls <- dat_clean %>%
  nest(Data = -c(ID)) %>%
  mutate(Data = map(Data, function(x) {
    x$NewBurst = lag(x$Night) != x$Night
    x$NewBurst[1] <- T
    x$Burst <- cumsum(x$NewBurst)
    x$NewBurst <- NULL
    return(x)
  })) %>%
  unnest(Data) %>%
  group_by(ID, Burst, Night) %>%
  summarize(
      Moonlit = mean(meanMoonlight) > 0.2
    , sl      = mean(sl)
    , .groups = "drop"
  ) %>%
  mutate(
    Time = case_when(
        Night & Moonlit ~ "Moonlit Night"
      , Night & !Moonlit ~ "Dark Night"
      , .default = "Day"
    )
  ) %>%
  group_by(ID, Time) %>%
  summarize(
    , MeanStepLength = mean(sl)
    , SDStepLength = sd(sl)
    , .groups = "drop"
  )
mean_sls %>%
  group_by(Time) %>%
  summarize(SD = sd(MeanStepLength), MeanStepLength = mean(MeanStepLength)) %>%
  ggplot(aes(x = Time, y = MeanStepLength, ymin = MeanStepLength - SD, ymax = MeanStepLength + SD)) +
    geom_col() +
    geom_errorbar()

# Function to create an overview of a gps dataset
gpsOverview <- function(data, id, timestamp) {

  data <- dat
  id <- "ID"
  timestamp <- "Timestamp"

  # Rename columns
  data$ID <- data[, id, drop = T]
  data$Timestamp <- data[, timestamp, drop = T]

  # Count the total number of individuals
  totalinds <- length(unique(data$ID))

  # Average daily and nightly distances. We'll need to identify bursts of "days"
  # and "nights" for this (we can't go by date as a night spans multiple dates).
  mean_sls <- dat_clean %>%
    nest(Data = -c(ID)) %>%
    mutate(Data = map(Data, function(x) {
      x$NewBurst = lag(x$Night) != x$Night
      x$NewBurst[1] <- T
      x$Burst <- cumsum(x$NewBurst)
      x$NewBurst <- NULL
      return(x)
    })) %>%
    unnest(Data) %>%
    # group_by(ID, Burst, Night) %>%
    group_by(ID, Night) %>%
    summarize(
        # Moonlit = mean(meanMoonlight) > 0.2
        sl      = mean(sl)
      , .groups = "drop"
    ) #%>%
    # mutate(
    #   Time = case_when(
    #       Night & Moonlit ~ "Moonlit Night"
    #     , Night & !Moonlit ~ "Dark Night"
    #     , .default = "Day"
    #   )
    # ) %>%
    # mutate(Time = factor(Time, levels = c("Dark Night", "Moonlit Night", "Day"))) %>%
    # group_by(ID, Time) %>%
    # summarize(
    #   , MeanStepLength = mean(sl)
    #   , SDStepLength = sd(sl)
    #   , .groups = "drop"
    # )

  # Means across individuals
  mean_sls2 <- mean_sls %>%
    group_by(Time) %>%
    summarize(MeanStepLength = mean(MeanStepLength))

  p5 <- ggplot(mean_sls, aes(x = sl, y = ID, fill = ID, col = ID)) +
    # geom_vline(data = mean_sls2, aes(xintercept = MeanStepLength), col = "black", lty = "22") +
    geom_col() +
    # geom_errorbarh(aes(xmin = MeanStepLength - SDStepLength, xmax = MeanStepLength + SDStepLength), height = 0) +
    facet_grid(~ Night) +
    scale_fill_viridis_d() +
    scale_color_viridis_d() +
    theme_awesome() +
    theme(
        axis.text.y      = element_blank()
      , axis.title.y     = element_blank()
      , legend.position  = "none"
      , panel.grid.minor = element_blank()
    )

  # Euclidean dispersal distance
  library(terra)
  test <- dat_raw %>%
    subset(State == "Disperser") %>%
    nest(Data = -ID) %>%
    mutate(Distance = map_dbl(Data, function(x) {
      fila <- x[c(1, nrow(x)), ]
      fila <- vect(fila, geom = c("x", "y"), crs = "epsg:4326")
      dist <- distance(fila, unit = "km")[1]
      return(dist)
    }))
  print(test, n = 30)
  subset(dat_raw, ID == "MadameChing") %>%
    ggplot(aes(x = x, y = y)) +
      geom_point() +
      geom_path() +
      coord_sf()

  # Put the plots together
  p <- ggarrange(p1, p2, p5, nrow = 1, align = "h", widths = c(1, 0.5, 2))
  p <- annotate_figure(p, top = text_grob(paste0("Total Number of Individuals: ", totalinds)))
  return(p)

}
