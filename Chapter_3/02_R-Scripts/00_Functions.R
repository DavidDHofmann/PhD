################################################################################
#### Function to Expand the Segments of a Dendrogram
################################################################################
# Function to stretch a dendrogram
stretchDendro <- function(dend, to = 0.2) {

  # Find segments that end in a leaf
  dend$segments$y <- sapply(1:nrow(dend$segments), function(i) {
    sub <- subset(dend$leaf_labels, x == dend$segments$x[i] & y == dend$segments$y[i])
    y <- ifelse(nrow(sub) > 0, to, dend$segments$y[i])
    return(y)
  })
  dend$segments$yend <- sapply(1:nrow(dend$segments), function(i) {
    sub <- subset(dend$leaf_labels, x == dend$segments$xend[i] & y == dend$segments$yend[i])
    y <- ifelse(nrow(sub) > 0, to, dend$segments$yend[i])
    return(y)
  })

  # Adjust leafs
  dend$leaf_labels$y <- to

  # Return stretched data
  return(dend)
}

################################################################################
#### Function to Compute Moonlight Statistics
################################################################################
# Function to calculate moonlight intensity (adapted from the moonlit package)
moonlightIntensity <- function(lon, lat, date, e) {

  # Unique combinations of input values for which to calculate moonlight stats
  if (length(lon) == length(lat) & length(lat) == length(date)) {
      d1 <- data.frame(date, lon, lat)
    } else {
      d1 <- expand.grid(date = unique(date), lon = unique(lon), lat = unique(lat))
  }

  # Calculate metrics using suncalc
  moonPosition <- suncalc::getMoonPosition(data = d1)
  moonIllum    <- suncalc::getMoonIllumination(date = date)
  sunPosition  <- suncalc::getSunlightPosition(data = d1)
  phaseAngle   <- acos((2 * moonIllum$fraction) - 1)

  # Prepare output table
  night                    <- d1
  night                    <- cbind(night, moonIllum[match(night$date, moonIllum$date), -1])
  night$moonAlt            <- moonPosition$altitude
  night$sunAlt             <- sunPosition$altitude
  night$moonAltDegrees     <- (180 / pi) * night$moonAlt
  night$sunAltDegrees      <- (180 / pi) * night$sunAlt
  night$moonIllum          <- moonIllum$fraction
  night$moonIllumCorrected <- (night$moonAlt > 0) * night$moonIllum
  night$phase              <- moonIllum$phase
  night$distance           <- moonPosition$distance
  rownames(night) <- NULL

  # Adding a "night" field - below 6 degrees equals to nautical twilight and
  # darker
  night$night <- night$sunAltDegrees < (0)

  ##############################################################################
  #### Modelling moonlight intensity on the ground
  ##############################################################################
  # Model provides an iluminance value for a given location and time. Outcome of
  # the model is a scalar quantity and represents a proportion of reference
  # illumination value. The reference value, in this case, is luminance of full
  # moon at 384400 km (average distance) in zenith. Model includes 4 variables,
  # each providing a correction factor of 0-1:
  # 1) visual extinction
  # 2) distance to the moon relative to mean distance
  # 3) phase effect
  # 4) Angle of incidence of moonlight
  # Final value is obtained by multiplying all 4 values.

  ##############################################################################
  #### 1) Visual Extinction
  ##############################################################################
  # 1) Visual extinction of the stars depending on the angle of the atmosphere
  # equations provided in: Green, Daniel W. E. 1992. "Magnitude Corrections for
  # Atmospheric Extinction." International Comet Quarterly 14: 55-59 Average
  # extinction coefficients (magnitude per air mass) are as follows:
  # At sea level: 0.28
  # At 500m asl: 0.24
  # At 1000m asl: 0.21
  # At 2000m asl: 0.16
  # Default value in the model is for 1000m asl and, if needed, should be
  # changed below:

  # Calculating thikness of the atmosphere (airmass value) at given elevation.
  # Note that elevation is above the horizon, so we need to use 90-elevation
  # from suncalc model. might show false results for moon altitude <0, but these
  # always have value 0 in end model so it doesn't matter. moon altitude must be
  # in radians!
  night$airMass <- 1 / (cos(pi / 2 - night$moonAlt) +
    0.025 * 2.71827 ** (-11 * cos(pi / 2 - night$moonAlt)))

  # Difference in magnitudes - in zenith the absortion reduces brightess by 0.24
  # magnitude, or ~20%, a minimal value possible. We will calculate the
  # difference between extinction at given angle and extinction at zenith and
  # transform it to proportion
  night$airMassMagnitude <- (night$airMass - 1) * e

  # Proportion
  night$atmExtinction <- 10 ** (-0.4 * night$airMassMagnitude)

  ##############################################################################
  #### 2) Distance to the Moon
  ##############################################################################
  # 2) Distance to the moon. The illumination is proportional to 1/(distance^2)
  # (inverse-square law, Meeus 1991), therefore relative illumination:
  night$distanceIllum = 1 / (night$distance / 384400) ** 2

  ##############################################################################
  #### 3) Phase Effect
  ##############################################################################
  # 3) phase effect. Relationship between phase angle (Sun-Moon-observer) and
  # normalised, disc-integrated brightness of the moon Data extracted from a
  # plot and interpolation function created for further use in the model Source:
  # Buratti, Bonnie J., John K. Hillier, and Michael Wang. 1996. "The Lunar
  # Opposition Surge: Observations by Clementine." Icarus 124 (December):
  # 490-99. https://doi.org/10.1006/icar.1996.0225.
  # Data points

  # phase angle in degrees
  phaseAngle <- c(
      0.367186581, 0.561328356, 0.755753776, 1.017911051, 1.343778816
    , 1.638020082, 1.9321973, 2.201542706, 2.446200869, 2.690859033, 2.960111111
    , 3.304228275, 3.798169858, 4.566540124, 5.659518624, 6.926950784
    , 8.294080819, 9.760921815, 11.2775818, 12.86918679, 14.51056309
    , 16.15193939, 17.79332234, 19.43470363, 21.07608991, 22.71746788
    , 24.55846678, 26.59940044, 28.6403448, 30.68133202, 32.72256159
    , 34.76382063, 36.805081, 38.84640163, 40.88817885, 42.93016095
    , 44.97243897, 47.01472501, 49.05701374, 51.09936005, 53.14188578
    , 55.18442088, 57.22695599, 59.2695018, 61.3123181, 63.35526696
    , 65.39821581, 67.44116333, 69.48411219, 71.52711193, 73.57035938
    , 75.61364299, 77.65692793, 79.7002102, 81.74354871, 83.78698764
    , 85.8304172, 87.87385613, 89.91732854, 91.961129, 94.00502721
    , 96.04892543, 98.09282096, 100.1367821, 102.1807767, 104.2247874
    , 106.2687927, 108.3127874, 110.3568141, 112.4010042, 114.4452425
    , 116.4894862, 118.5337607, 120.5781851, 122.6226792, 124.6673768
    , 126.7121119, 128.7568778, 130.8016611, 132.8463627, 134.8910482
    , 136.9357699, 138.9805559, 141.0253124, 143.0700395, 145.1147611
    , 147.1594494, 149.2040787, 151.2487227, 153.2933882, 155.3381112
    , 157.382853, 159.4275827, 161.4723968, 163.517239, 165.5620973
    , 167.6069529, 169.6517964, 171.6966668, 173.7414943, 175.786378
    , 177.8311961, 179.5020253
  )

  # Converting phase angle to % illuminated
  discIlluminated <- (1 + cos(phaseAngle * pi / 180)) / 2

  # Disc-integrated brighness
  phaseBrightness <- c(
      0.981495601, 0.957317682, 0.93441936, 0.911719223, 0.887147319
    , 0.864543232, 0.841650203, 0.819234668, 0.797948808, 0.776662948
    , 0.753826384, 0.731224944 ,0.709518906, 0.688329272, 0.669028531
    , 0.649231859 , 0.629197232, 0.608983672, 0.588516023, 0.568644051
    , 0.54830288, 0.52796171, 0.50765056, 0.487331904, 0.467035764
    , 0.446702099, 0.426908396, 0.409070349, 0.391280626, 0.373684202
    , 0.357181123, 0.340810937, 0.324446792, 0.308354473, 0.294321992
    , 0.281213721, 0.269440419, 0.25770336, 0.245978382, 0.234513149
    , 0.223857355, 0.213243844, 0.202630334, 0.192065148, 0.18272016
    , 0.17397319, 0.165226219, 0.156473208, 0.147726238, 0.139208809
    , 0.131808889, 0.124572064, 0.11734128, 0.110098415, 0.103109254
    , 0.096573137, 0.089994736, 0.083458618, 0.077073516, 0.072168356
    , 0.067704159, 0.063239962, 0.058763683, 0.054583393, 0.050554118
    , 0.04659733, 0.042616379, 0.038587104, 0.034702802, 0.031555452
    , 0.028625563, 0.025719836, 0.022953043, 0.020862794, 0.019086656
    , 0.018228687, 0.017539854, 0.016989954, 0.016518582, 0.015678735
    , 0.0147664, 0.014017161, 0.013557871, 0.012965687, 0.01224061
    , 0.011491371, 0.010591118, 0.009425079, 0.008325486, 0.007322543
    , 0.006579345, 0.005920715, 0.005207719, 0.004875281, 0.004669694
    , 0.004536595, 0.004391414, 0.004191868, 0.004113134, 0.003841101
    , 0.003822773, 0.003508456, 0.003536722
  )

  # Spline interpolation function
  brightness <- stats::splinefun(discIlluminated, phaseBrightness
    , method = "fmm"
    , ties   = mean
  )

  # Correction factor for model
  night$phaseEffect <- brightness(night$moonIllum)

  ##############################################################################
  #### 4) Angle of Incidence
  ##############################################################################
  # 4) angle of incidence. Correction factor is sine of the angle of moon's
  # altitude (Austin et al. 1976)
  night$incidence <- sin(night$moonAlt)

  ##############################################################################
  #### Total correction factor
  ##############################################################################
  # Assign values of 0 if not night and 0 if moon not visible then multiplied
  # variables 1:4 turns out there are some occasional NAs because atmExtinction
  # is infinite, so will have to correct for that by making a conditional
  # version of the model output
  night$illuminationModelComponents <- (night$night > 0) *
    (night$moonAlt > 0 ) * night$atmExtinction * night$distanceIllum *
    night$phaseEffect * night$incidence
  night$illuminationModel <- ifelse((night$moonAlt <= 0 | night$night <= 0)
    , yes = 0
    , no  = night$illuminationModelComponents
  )
  night$moonlightModel <- night$illuminationModel

  # Alternatively calculating values for moonlight also during the night
  night$moonlight24hComponents <- (night$moonAlt > 0 ) * night$atmExtinction *
    night$distanceIllum * night$phaseEffect * night$incidence
  night$moonlightModel24 <- ifelse((night$moonAlt <= 0)
    , yes = 0
    , no  = night$moonlight24hComponents
  )

  # Twilight illumination from empirical data
  # https://www.jstor.org/stable/44612241 Sun Position and Twilight Times for
  # Driver Visibility Assessment; Duane D. MacInnis, Peter B. Williamson and
  # Geoffrey P. Nielsen; SAE Transactions; Vol. 104, Section 6: JOURNAL OF
  # PASSENGER CARS: Part 1 (1995), pp. 759-783 Values calulated from the paper
  # for angles -18 to 5 degrees, by 0.25 degree, imported to R and spline
  # function applied
  solarAngle <- c(
      5, 4.75, 4.5, 4.25, 4, 3.75, 3.5, 3.25, 3, 2.75, 2.5, 2.25, 2, 1.75, 1.5
    , 1.25, 1, 0.75, 0.5, 0.25, 0, -0.25, -0.5, -0.75, -1, -1.25, -1.5, -1.75
    , -2, -2.25, -2.5, -2.75, -3, -3.25, -3.5, -3.75, -4, -4.25, -4.5, -4.75
    , -5, -5.25, -5.5, -5.75, -6, -6.25, -6.5, -6.75, -7, -7.25, -7.5, -7.75
    , -8, -8.25, -8.5, -8.75, -9, -9.25, -9.5, -9.75, -10, -10.25, -10.5
    , -10.75, -11, -11.25, -11.5, -11.75, -12, -12.25, -12.5, -12.75, -13
    , -13.25, -13.5, -13.75, -14, -14.25, -14.5, -14.75, -15, -15.25, -15.5
    , -15.75, -16, -16.25, -16.5, -16.75, -17, -17.25, -17.5, -17.75, -18
  )
  twilightIllumination <- c(
    4499.368645, 4251.166227, 4010.398231, 3776.215801, 3547.994929,
    3325.316934, 3107.949504, 2895.827638, 2689.034018, 2487.778469,
    2292.376367, 2103.225979, 1920.784922, 1745.546056, 1578.01326,
    1418.677674, 1267.99505, 1126.364901, 994.1121318, 871.4718063,
    758.577575, 655.4541802, 562.0142693, 478.0595544, 404.556034,
    338.5269901, 281.0455199, 231.5470057, 189.3612856, 153.7599131,
    123.995336, 99.33178815, 79.06826738, 62.55436811, 49.1999767,
    38.47993901, 29.93480689, 23.16869331, 17.84514167, 13.68176809,
    10.44428227, 7.911259933, 5.993606008, 4.527348005, 3.411591605,
    2.566090598, 1.927671576, 1.447053311, 1.086102387, 0.815522737,
    0.612949604, 0.461403805, 0.348056383, 0.263253638, 0.199755997,
    0.152149406, 0.116393841, 0.089479462, 0.069166322, 0.053788354,
    0.04210643, 0.033198579, 0.026378238, 0.02113351, 0.017082103,
    0.01393796, 0.01148655, 0.009566583, 0.008056477, 0.006931022,
    0.00596833, 0.00519104, 0.004556611, 0.004033257, 0.003596965,
    0.003229411, 0.002916468, 0.002647147, 0.002412827, 0.002206692,
    0.002023321, 0.001858385, 0.001708419, 0.001570657, 0.001442901,
    0.001323428, 0.001210913, 0.001104369, 0.001003101, 0.000906659,
    0.000814804, 0.000727464, 0.000644703
  )

  # Spline interpolation function
  twilightIlluminationFunction <- splinefun(
      x = solarAngle
    , y = twilightIllumination
  )

  # The model can only be interpolated to -18 degrees, so need to assign 0 to
  # everything below -18 degrees
  night$twilightModel <- twilightIlluminationFunction(night$sunAltDegrees)
  night$twilight      <- night$sunAltDegrees > -18
  night$twilightModel <- night$twilightModel * night$twilight

  # Combining twilight illumination with lunar illumination
  night$illumination <- night$illuminationModel * 0.32 + night$twilightModel

  # Generating output data frame
  d1$night          <- night$night
  d1$sunAltDegrees  <- night$sunAltDegrees
  d1$moonlightModel <- night$moonlightModel
  d1$twilightModel  <- night$twilightModel
  d1$illumination   <- night$illumination
  d1$moonPhase      <- night$moonIllum

  # This is an obsolete column that calculates moonlight intensity during the
  # day. It was needed for some tests in the past, it is now obsolete as it has
  # no use in ecological studies. It is left here in case it is needed in the
  # future
  d1$moonlightModel24 <- night$moonlightModel24
  d1$moonAltDegrees   <- night$moonAltDegrees

  # Return the dataframe
  return(d1)
}

################################################################################
#### Function to Compute Moonlight Summary
################################################################################
# Function to compute summary statistics about the moonlight intensity at night
moonlightSummary <- function(lat, lon, date, e, t) {

  # Nightly statistics are calculated for the FOLLOWING night!!!
  date <- as.Date(date)

  # Unique combinations through wich we need to loop
  if (length(lon) == length(lat) & length(lat) == length(date)) {
      data <- data.frame(date, lon, lat)
    } else {
      data <- expand.grid(date = unique(date), lon = unique(lon), lat = unique(lat))
  }

  # Determine sunrise and sunset for each of the dates
  data$sunset_last  <- as.POSIXct(NA)
  data$sunrise      <- as.POSIXct(NA)
  data$sunset       <- as.POSIXct(NA)
  data$sunrise_next <- as.POSIXct(NA)
  data <- split(data, data[, c("lon", "lat")])
  data <- lapply(data, function(x) {
    x$sunset_last  <- getSunlightTimes(x$date - days(1), lon = unique(x$lon), lat = unique(x$lat), keep = "sunset")$sunset
    x$sunrise      <- getSunlightTimes(x$date, lon = unique(x$lon), lat = unique(x$lat), keep = "sunrise")$sunrise
    x$sunset       <- getSunlightTimes(x$date, lon = unique(x$lon), lat = unique(x$lat), keep = "sunset")$sunset
    x$sunrise_next <- getSunlightTimes(x$date + days(1), lon = unique(x$lon), lat = unique(x$lat), keep = "sunrise")$sunrise
    return(x)
  })
  data <- do.call(rbind, data)
  rownames(data) <- NULL

  # Span the sequence of times to go through
  nights <- lapply(1:nrow(data), function(x) {
    night <- seq(from = data$sunset[x], to = data$sunrise_next[x], by = t)
    night <- data.frame(
        date         = data$date[x]
      , lon          = data$lon[x]
      , lat          = data$lat[x]
      , sunset_last  = data$sunset_last[x]
      , sunrise      = data$sunrise[x]
      , sunset       = data$sunset[x]
      , sunrise_next = data$sunrise_next[x]
      , timestamp    = night
    )
    return(night)
  })
  nights <- do.call(rbind, nights)

  # Compute moonlight statistics
  nights_stats <- moonlightIntensity(
      date = nights$timestamp
    , lon  = nights$lon
    , lat  = nights$lat
    , e    = e
  )
  names(nights_stats)[1] <- "timestamp"
  nights <- cbind(date = nights$date, nights_stats)

  # Aggregate by date and coordinates
  nights <- nights %>%
    group_by(date, lon, lat) %>%
    summarize(
        meanMoonlightIntensity = mean(moonlightModel)
      , maxMoonlightIntensity  = max(moonlightModel)
      , minMoonlightIntensity  = min(moonlightModel)
      , meanMoonPhase          = mean(moonPhase)
      , maxMoonPhase           = max(moonPhase)
      , minMoonPhase           = min(moonPhase)
      , maxMoonTime            = timestamp[moonlightModel == max(moonlightModel)][1]
      , .groups                = "drop"
    )

  # Bind with original data
  nights <- left_join(data, nights, by = c("date", "lon", "lat"))

  # return it
  return(nights)
}

################################################################################
#### Von Mises Distribution Functions
################################################################################
# Function to determine the pdf of a mixed von mises distribution
dvonmises <- function(x, kappa, mu){
  exp(kappa * cos(x - mu)) / (2 * pi * besselI(kappa, nu = 0))
}

# Function to randomly sample from a mixed von mises distribution
rvonmises <- function(n, kappa, mu, by = 0.01){
  x <- seq(-pi, +pi, by = by)
  probs <- dvonmises(x, kappa = kappa, mu = mu)
  random <- sample(x, size = n, prob = probs, replace = T)
  return(random)
}

################################################################################
#### GGplot Basic Theme
################################################################################
# Custom ggplot theme
theme_awesome <- function() {
    theme_minimal() +
    theme(
        axis.title.y      = element_text(angle = 0, vjust = 0.5)
      , legend.position   = "bottom"
      , legend.key.width  = unit(0.5, "cm")
      , legend.key.height = unit(0.2, "cm")
      , strip.background  = element_rect(color = "white", fill = "gray95")
      , legend.text       = element_text(size = 5)
    )
}

# Custom ggplot theme
dark_theme_awesome <- function() {
    dark_theme_minimal() +
    theme(
        axis.title.y      = element_text(angle = 0, vjust = 0.5)
      , legend.position   = "bottom"
      , legend.key.width  = unit(0.5, "cm")
      , legend.key.height = unit(0.2, "cm")
      , strip.background  = element_rect(fill = "gray10", color = NA)
      , legend.text       = element_text(size = 5)
      , panel.background  = element_blank()
      , plot.background   = element_blank()
    )
}

################################################################################
#### Function to Normalize a Terra Raster
################################################################################
# Function to normalize a raster
normalizeRaster <- function(r) {
  minmax <- minmax(r)
  r <- (r - minmax["min", ]) / (minmax["max", ] - minmax["min", ])
  return(r)
}

################################################################################
#### Function to Summarize GPS Sampling
################################################################################
# Function to summarize the sampling frequency
gpsSampling <- function(data, id, timestamp) {

  # Rename columns
  data$ID <- data[, id, drop = T]
  data$Timestamp <- data[, timestamp, drop = T]

  # Calculate summaries
  smry <- data %>%
    nest(GPS = -ID) %>%
    mutate(SamplingFrequency = map(GPS, function(x) {
      smry <- x %>%
        arrange(Timestamp) %>%
        mutate(dt = difftime(Timestamp, lag(Timestamp), units = "hours")) %>%
        summarize(
            min        = min(dt, na.rm = T)
          , Quantile25 = quantile(dt, 0.25, na.rm = T)
          , Quantile50 = quantile(dt, 0.5, na.rm = T)
          , Quantile75 = median(dt, na.rm = T)
          , mean       = mean(dt, na.rm = T)
          , max        = max(dt, na.rm = T)
        )
      return(smry)
    })) %>% dplyr::select(ID, SamplingFrequency) %>% unnest(SamplingFrequency)
  return(smry)
}

################################################################################
#### Function to Generate a GPS Overview
################################################################################
# Function to create an overview of a gps dataset
gpsOverview <- function(data, id, timestamp) {

  # Rename columns
  data$ID <- data[, id, drop = T]
  data$Timestamp <- data[, timestamp, drop = T]

  # Count the total number of individuals
  totalinds <- length(unique(data$ID))

  # Number of fixes per dog
  p1 <- data %>%
    count(ID, name = "NumberOfFixes") %>%
    ggplot(aes(x = NumberOfFixes, y = ID, fill = ID, color = ID)) +
      geom_col(alpha = 0.75, width = 0.75) +
      scale_fill_viridis_d() +
      scale_color_viridis_d() +
      scale_x_continuous(labels = function(x) {format(x, big.mark = ",")}) +
      theme_awesome() +
      theme(legend.position = "none", panel.grid.minor = element_blank()) +
      xlab("Number of Fixes") +
      facet_wrap(~ "Number of Fixes")

  # Hour of the day
  p2 <- data %>%
    mutate(TimestampRounded = round_date(Timestamp, "1 hour")) %>%
    mutate(Hour = hour(TimestampRounded)) %>%
    count(ID, Hour, name = "NumberOfFixes") %>%
    ggplot(aes(x = Hour, y = ID, fill = ID, color = ID)) +
      geom_point() +
      scale_fill_viridis_d() +
      scale_color_viridis_d() +
      theme_awesome() +
      theme(
          axis.text.y     = element_blank()
        , axis.title.y    = element_blank()
        , legend.position = "none"
      ) +
      xlab("Hour of the Day") +
      facet_wrap(~ "Time of the Day")

  # Sampling Frequency
  threshold <- 25
  p3 <- data %>%
    nest(GPS = -ID) %>%
    mutate(GPS = map(GPS, function(x) {
      x <- x %>%
        arrange(Timestamp) %>%
        mutate(dt = as.numeric(difftime(Timestamp, lag(Timestamp), units = "hours")))
      return(x)
    })) %>%
    unnest(GPS) %>%
    subset(!is.na(dt) & dt < threshold) %>%
    ggplot(aes(x = dt, y = ID, fill = ID, color = ID)) +
      geom_boxplot(alpha = 0.75) +
      stat_summary(geom = "point", fun = median, pch = 3, col = "black", alpha = 0.75) +
      scale_fill_viridis_d() +
      scale_color_viridis_d() +
      theme_awesome() +
      theme(
          axis.text.y     = element_blank()
        , axis.title.y    = element_blank()
        , legend.position = "none"
      ) +
      xlab("Duration between Subsequent\nGPS Fixes (Hours)") +
      facet_wrap(~ "Time Lag")

  # Put the plots together
  p <- ggarrange(p1, p2, p3, nrow = 1, align = "h")
  p <- annotate_figure(p, top = text_grob(paste0("Total Number of Individuals: ", totalinds)))
  return(p)

}

################################################################################
#### Function to summarize numeric data from dataframe
################################################################################
# Function to generate a summary over numeric variables in a dataframe
summarizeData <- function(x) {
  smry <- x %>%
    select(where(is.numeric)) %>%
    pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
    group_by(Variable) %>%
    summarize(
        Min    = min(Value, na.rm = T)
      , Median = median(Value, na.rm = T)
      , Max    = max(Value, na.rm = T)
      , Mean   = mean(Value, na.rm = T)
      , SD     = sd(Value, na.rm = T)
      , `NA`   = sum(is.na(Value))
    )
  return(smry)
}

################################################################################
#### Function to standardize (scale or normalize) values
################################################################################
# Function to standardize (scale or normalize) values
standardize <- function(x, operation = c("scale", "normalize")) {
  match.arg(operation)
  center <- ifelse(operation == "scale", mean(x, na.rm = T), min(x, na.rm = T))
  scale  <- ifelse(operation == "scale", sd(x, na.rm = T), diff(range(x, na.rm = T)))
  vals_standardized <- (x - center) / scale
  attr(vals_standardized, "standardization")    <- operation
  attr(vals_standardized, "standardize:center") <- center
  attr(vals_standardized, "standardize:scale")  <- scale
  return(vals_standardized)
}

################################################################################
#### Function to backtransform scaled or normalized values
################################################################################
# Function to backtransform scaled or normalized values
backtransform <- function(x) {
  if (is.null(attr(x, "standardization"))) {
    stop("This data has not been standardized")
  }
  center <- attr(x, "standardize:center")
  scale  <- attr(x, "standardize:scale")
  vals_backtransformed <- x * scale + center
  attr(vals_backtransformed, "standardize:center") <- NULL
  attr(vals_backtransformed, "standardize:scale")  <- NULL
  return(vals_backtransformed)
}

################################################################################
#### Function to determine Flood Season
################################################################################
# Function to determine the season within which a date falls
getHerbivoreSeason <- function(x) {
    DS <- as.Date("2020-01-01", format = "%Y-%m-%d") # Start Dispersed
    CS <- as.Date("2020-07-01", format = "%Y-%m-%d") # Start Concentrated

    # Convert dates from any year to 2020 dates
    d <- as.Date(strftime(x, format = "2020-%m-%d"))

    # Identify season
    season <- ifelse(d >= DS & d < CS, "Dispersed", "Concentrated")

    # Return the season
    return(season)
}

################################################################################
#### Function to determine Rainy Season
################################################################################
# Function to determine the season within which a date falls
getRainySeason <- function(x) {
    DS <- as.Date("2020-04-15", format = "%Y-%m-%d") # Start Dry
    WS <- as.Date("2020-10-15", format = "%Y-%m-%d") # Start Wet

    # Convert dates from any year to 2020 dates
    d <- as.Date(strftime(x, format = "2020-%m-%d"))

    # Identify season
    season <- ifelse(d >= DS & d < WS, "Dry", "Wet")

    # Return the season
    return(season)
}

################################################################################
#### Function to List Available MODIS Data (adapted from luna package)
################################################################################
# Function to list available modis data. This is a slightly improved version
# from the luna package that not only returns the filename, but the full url
# including some info
modSearch <- function(
    product
  , start_date
  , end_date
  , aoi
  , version   = "006"
  , download  = FALSE
  , server    = "LPDAAC_ECS"
  , limit     = 100000
  , ...) {

  # Testing
#   product    = "MOD44B"
#  product    = "MCD43A4"
#   start_date = "2011-01-01"
  # end_date   = "2011-03-01"
#   aoi        = extent(c(20, 28, -22, -16))
#   server     = "LPDAAC_ECS"
#   version    = "006"
#   limit <- 100000

  # Some checks
	if(missing(product)) stop("provide a product name")
	if(missing(start_date)) stop("provide a start_date")
	if(missing(end_date)) stop("provide an end_date")
	if(missing(aoi)) stop("provide an area of interest")

  # Subset to products of interest
	h <- luna:::.humanize()
	h <- h[h$short_name == product, ]
	pp <- h[h$version == version & h$provider == server, ]

  # Check if the desired products are available
  if (nrow(pp) < 1) {
    if (nrow(h) < 1) {
      stop("The requested product is not available for this through this function")
    } else {
      cat("Options for this product:\n")
      print(head(h, 10))
      cat("\n")
      stop("The requested product is not available for this product version or server")
    }
    } else if (nrow(pp) > 1) {
      warning("Multiple sources available, using first one")
      print(pp)
      pp <- pp[1, ]
  }

  # Find product URLs
	results <- luna:::.searchGranules(
      product    = product
    , start_date = start_date
    , end_date   = end_date
    , extent     = aoi
    , limit      = limit
    , version    = version
  )
	urls <- unique(results[, "Online Access URLs"])

  # Makre results readable
  info <- modInfo(urls)
  res <- cbind(info, url = urls)
  return(res)
}

################################################################################
#### Function to Download MODIS Data (adapted from luna package)
################################################################################
# Function to download files from MODIS urls
modDownload <- function(url, username, password, path = getwd(), overwrite = F) {
  path <- luna:::.getPath(path)
  if(missing(username)) stop("provide a username")
  if(missing(password)) stop("provide a password")
  ff <- luna:::.cmr_download(url, path, username, password, overwrite)
  ff <- file.path(path, basename(url))
  return(ff)
}

################################################################################
#### Function to extract information from modis filenames
################################################################################
# Function to extract information from modis filenames
modInfo <- function(urlname) {
  name <- basename(urlname)
  info <- strsplit(name, split = "\\.")
  info <- lapply(info, rbind)
  info <- do.call(rbind, info)
  info <- as.data.frame(info, stringsAsFactors = F)
  names(info) <- c("Product", "AcquisitionDate", "Tile", "Version", "ProductionDate", "Format")
  info$AcquisitionDate <- substr(info$AcquisitionDate, start = 2, stop = nchar(info$AcquisitionDate))
  info$AcquisitionDate <- as.Date(info$AcquisitionDate, "%Y%j")
  info$ProductionDate <- as.Date(as.POSIXct(info$ProductionDate, format = "%Y%j%H%M%S"))
  return(info)
}

################################################################################
#### Function to Determine the Corrected Sentinel 2 Names
################################################################################
# Function to determine the name of the corrected Sentinel 2 product
correctedName <- function(x) {
  corr <- gsub(basename(x), pattern = "MSIL1C", replacement = "MSIL2A")
  return(corr)
}

################################################################################
#### Function to Compute Normalized Difference Index
################################################################################
# Function to compute the normalized difference (nd) index of two bands
nd <- function(img, band_x, band_y) {
  x <- img[[band_x]]
  y <- img[[band_y]]
  nd <- (x - y) / (x + y)
  return(nd)
}

################################################################################
#### Function to calculate Distances on a Raster Efficiently
################################################################################
distanceTo <- function(x, value = 1, retainNAs = F) {

  # Get a mask of NAs if so desired
  if (retainNAs) {
    maski <- is.na(x)
  }

  # Project input layer and compute distances
  d <- x
  d <- project(d, "+init=epsg:32734", method = "near")
  d <- d == value
  d <- as.numeric(d)
  d <- subst(d, 0, NA)
  d <- distance(d)

  # Unproject and crop to original extent
  d <- project(d, x, method = "near", threads = T)
  d <- crop(d, x)

  # Add NAs back
  if (retainNAs) {
    d <- mask(d, maski, maskvalue = T)
  }
  return(d)
}

################################################################################
#### Function to Resample GPS Data to a Desired Resolution
################################################################################
# Function to resample fixes to a desired resolution
resampleFixes <- function(data, hours, start, tol = 0.5) {

  # Identify the first date at which a fix was taken and set the time to the
  # desired start time
  first <- range(data$Timestamp)[1] %>%
    update(., hour = start, min = 0, sec = 0)

  # Identify the last date at which a fix was taken and set the end time to 24
  last <- range(data$Timestamp)[2] %>%
    update(., hour = 24, min = 0, sec = 0)

  # Prepare a range of dates for which we would expect to find data according
  # to the specified sampling scheme
  dates <- seq(first, last, by = paste0(hours, " hours")) %>%
    as.data.frame() %>%
    set_names("Timestamp")

  # For each Timestamp we now identify the closest fix
  closest <- sapply(1:nrow(dates), function(x) {
    index <- which.min(abs(dates$Timestamp[x] - data$Timestamp))[1]
    close <- as.numeric(abs(dates$Timestamp[x] - data$Timestamp[index]), units = "hours") <= tol

    # In case the fix is close enough, return its index
    if (close) {
      return(index)
    } else {
      return(NA)
    }
  })

  # Return respective fixes
  closest <- na.omit(closest)
  if (length(closest) > 0) {
    return(data[closest, ])
  } else {
    return(NA)
  }
}

################################################################################
#### Function to Compute Bursts
################################################################################
# Function to compute bursts
computeBursts <- function(data) {
  bursted <- data %>%
    mutate(dt_ = Timestamp - lag(Timestamp)) %>%
    mutate(dt_ = as.numeric(dt_, units = "hours")) %>%
    ungroup() %>%
    mutate(NewBurst = ifelse(
      dt_ > 8.25 |
      (dt_ > 4.25 & hour(TimestampRounded) != 15) |
      is.na(dt_), yes = 1, no = 0)
    ) %>%
    mutate(burst_id = cumsum(NewBurst)) %>%
    dplyr::select(-c(NewBurst, dt_))
  return(bursted)
}

################################################################################
#### Function to Compute Step Metrics
################################################################################
# Function to compute the absolute turning angle
absTA <- function(dx, dy) {
  absta <- atan2(dy, dx)
  absta <- (absta - pi / 2) * (-1)
  absta <- ifelse(absta < 0, 2 * pi + absta, absta)
  return(absta)
}

# Function to compute the relative turning angle
relTA <- function(absta) {
  relta <- (absta[-1] - absta[-length(absta)])
  relta <- c(NA, relta)
  relta <- ifelse(relta > +pi, relta - 2 * pi, relta)
  relta <- ifelse(relta < -pi, 2 * pi + relta, relta)
  return(relta)
}

# Function to calculate the new absolute turning angle, given a random relative
# turning angle
newabsTA <- function(absta, ta) {
  absta <- absta + ta
  absta[absta > 2 * pi] <-
    absta[absta > 2 * pi] - 2 * pi
  absta[absta < 0] <-
    absta[absta < 0] + 2 * pi
  return(absta)
}

# # Function to compute step metrics (data needs to be projected)
# stepMetrics <- function(data) {
#
#     # Get relevant data
#     x <- data$x
#     y <- data$y
#     timestamp <- data$Timestamp
#
#     # Compute distances moved in x and y direction
#     dx <- c(x[-1], NA) - x
#     dy <- c(y[-1], NA) - y
#
#     # Calculate step length
#     sl <- sqrt(dx ** 2 + dy ** 2)
#
#     # Compute absolute turn angle
#     absta <- absTA(dx, dy)
#
#     # Compute relative turn angle
#     relta <- relTA(absta)
#
#     # Compute step duration
#     if (!is.null(timestamp)) {
#       dt <- difftime(lead(timestamp), timestamp, units = "hours")
#     }
#
#     # Put metrics into data.frame
#     metrics <- data.frame(x_to = lead(x), y_to = lead(y), sl = sl, absta = absta, relta = relta)
#     if (!is.null(timestamp)) {
#       metrics$dt <- dt
#     }
#
#     # Combine with original data
#     data <- cbind(data, metrics)
#     data <- data[-nrow(data), ]
#
#     # Return the metrics
#     return(data)
# }
# Function to compute step metrics (data needs to be projected)
stepMetrics <- function(data) {

    # Get relevant data
    x <- data$x
    y <- data$y
    timestamp <- data$Timestamp

    # Compute distances moved in x and y direction
    dx <- c(x[-1], NA) - x
    dy <- c(y[-1], NA) - y

    # Calculate step length
    sl <- sqrt(dx ** 2 + dy ** 2)

    # Compute absolute turn angle
    absta <- absTA(dx, dy)

    # Compute relative turn angle
    relta <- relTA(absta)

    # Put metrics into data.frame
    metrics <- data.frame(
        x_to  = lead(x)
      , y_to  = lead(y)
      , sl    = sl
      , absta = absta
      , relta = relta
    )
    if (!is.null(timestamp)) {
      metrics$dt <- difftime(lead(timestamp), timestamp, units = "hours")
    }

    # Return the metrics
    return(metrics)
}


################################################################################
#### Function to Generate Random Steps
################################################################################
# Function to generate random steps
randomSteps <- function(data, n_rsteps, scale, shape, slmax = NULL, slunit = c("m", "km")) {

  # step lengths are only allowed to be in meters or kilometers
  match.arg(slunit)

  # Generate a new column that indicates that the steps are "observed" steps
  data$case <- 1

  # Generate step ids
  data$step_id <- 1:nrow(data)

  # Cannot work with steps that have no turning angle, so remove them
  data <- subset(data, !is.na(relta))

  # Create a new dataframe into which we can put alternative/random steps
  rand <- data[rep(1:nrow(data), each = n_rsteps), ]

  # Indicate that these steps are random steps (case = 0)
  rand$case <- 0

  # Sample random turning angles
  rand$sl <- rgamma(n = nrow(rand)
    , scale = scale
    , shape = shape
  )

  # If a maximum step length is provided, truncate the sampled step length if
  # necessary
  if (!is.null(slmax)) {
    rand$sl <- ifelse(rand$sl > slmax, slmax, rand$sl)
  }

  # Sample random step lengths
  rand$relta_new <- runif(n = nrow(rand)
    , min = -pi
    , max = +pi
  )

  # Adjust step length to meters if needed
  if (slunit == "km") {
    rand$sl <- rand$sl * 1000
    data$sl <- data$sl * 1000
  }

  # Calculate new "absolute" turning angle
  rand$relta_diff <- rand$relta_new - rand$relta
  rand$absta <- rand$absta + rand$relta_diff
  rand$absta <- ifelse(rand$absta < 0, 2 * pi + rand$absta, rand$absta)
  rand$absta <- ifelse(rand$absta > 2 * pi, rand$absta - 2 * pi, rand$absta)
  rand$relta <- rand$relta_new

  # Remove undesired stuff
  rand$relta_new <- NULL
  rand$relta_diff <- NULL

  # Put steps together
  all <- rbind(data, rand)
  all <- arrange(all, ID, step_id, desc(case))

  # Calculate new endpoints
  all$x_to <- all$x + sin(all$absta) * all$sl
  all$y_to <- all$y + cos(all$absta) * all$sl

  # Adjust step length to km if needed
  if (slunit == "km") {
    all$sl <- all$sl / 1000
  }

  # Return the final dataframe
  return(all)

}

################################################################################
#### Function to Check which Points are inside an Extent
################################################################################
# Function to Check which Points are inside an Extent. Extent values have to be
# given as list.
pointsInside <- function(xy, extent) {

  # Check if x value is within boundaries
  xy[, 1] > extent$xmin & xy[, 1] < extent$xmax &

    # Check if y value is within boundaries
    xy[, 2] > extent$ymin & xy[, 2] < extent$ymax
}

################################################################################
#### Function to Reproject Coordinates
################################################################################
# Reproject coordinates to another projection
reprojectCoords <- function(xy, from = NULL, to = NULL) {
  xy <- as.matrix(xy)
  xy <- suppressWarnings(terra::project(xy, from = from, to = to))
  xy[is.nan(xy)] <- NA
  return(xy)
}

################################################################################
#### Function to Project Coordinates to Local CRS
################################################################################
# Function to Project Coordinates to Local CRS
projectCoords <- function(dat) {
  dat[, c("x", "y")] <- reprojectCoords(
      xy   = cbind(dat$x, dat$y)
    , from = "epsg:4326"
    , to   = "epsg:32734"
  )
  if (all(c("x_to", "y_to") %in% names(dat))) {
    dat[, c("x_to", "y_to")] <- reprojectCoords(
        xy   = cbind(dat$x_to, dat$y_to)
      , from = "epsg:4326"
      , to   = "epsg:32734"
    )
  }
  return(dat)
}

################################################################################
#### Function to Un-Project Coordinates from Local CRS
################################################################################
# Function to Project Coordinates to Local CRS
unprojectCoords <- function(dat) {
  if (all(c("x", "y") %in% names(dat))) {
    dat[, c("x", "y")] <- reprojectCoords(
        xy   = cbind(dat$x, dat$y)
      , from = "epsg:32734"
      , to   = "epsg:4326"
    )
  }
  if (all(c("x_to", "y_to") %in% names(dat))) {
    dat[, c("x_to", "y_to")] <- reprojectCoords(
        xy   = cbind(dat$x_to, dat$y_to)
      , from = "epsg:32734"
      , to   = "epsg:4326"
    )
  }
  return(dat)
}

################################################################################
#### Function to Project Extent to Local CRS
################################################################################
# Function to project an extent
projectExte <- function(e) {
  exte <- data.frame(x = xmin(e), y = ymin(e), x_to = xmax(e), y_to = ymax(e))
  exte <- projectCoords(exte)
  exte <- ext(c(exte$x, exte$x_to, exte$y, exte$y_to))
  return(exte)
}

################################################################################
#### Function to Un-Project Extent from Local CRS
################################################################################
# Function to unproject an extent
unprojectExte <- function(e) {
  exte <- data.frame(x = xmin(e), y = ymin(e), x_to = xmax(e), y_to = ymax(e))
  exte <- unprojectCoords(exte)
  exte <- ext(c(exte$x, exte$x_to, exte$y, exte$y_to))
  return(exte)
}

################################################################################
#### Function to Interpolate Between Points
################################################################################
# Function to interpolate coordinates between two points
interpolatePoints <- function(x1, x2, y1, y2, t1 = NULL, t2 = NULL, by = 1) {

  # Calculate length of line between points
  length <- sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)

  # Calculate how many segments we need
  nsegs <- max(ceiling(length / by), 1)

  # Interpolate between points
  x <- seq(x1, x2, length.out = nsegs + 1)
  y <- seq(y1, y2, length.out = nsegs + 1)

  # Return
  if (!is.null(t1) & !is.null(t2)) {
    timestamp <- seq(t1, t2, length.out = nsegs + 1)
    result <- data.frame(x, y, timestamp)
  } else {
    result <- data.frame(x, y)
  }
  return(result)
}

################################################################################
#### Function to Interpolate a Track
################################################################################
# Function to interpolate a track
interpolateTrack <- function(trk, by = 1) {

  # Switch to step-representation
  trk_stp <- data.frame(
      step_number  = 1:nrow(trk)
    , timestamp    = trk$timestamp
    , timestamp_to = c(trk$timestamp[2:nrow(trk)], NA)
    , x            = trk$x
    , x_to         = c(trk$x[2:nrow(trk)], NA)
    , y            = trk$y
    , y_to         = c(trk$y[2:nrow(trk)], NA)
  )
  trk_stp <- trk_stp[-nrow(trk_stp), ]

  # Interpolate
  trk_inter <- lapply(seq_len(nrow(trk_stp)), function(i) {
    inter <- interpolatePoints(
        x1 = trk_stp$x[i]
      , x2 = trk_stp$x_to[i]
      , y1 = trk_stp$y[i]
      , y2 = trk_stp$y_to[i]
      , t1 = trk_stp$timestamp[i]
      , t2 = trk_stp$timestamp_to[i]
      , by = by
    )
    inter$interpolated              <- T
    inter$interpolated[1]           <- F
    inter$interpolated[nrow(inter)] <- F
    inter$step_number        <- trk_stp$step_number[i]
    inter$step_number_within <- seq_len(nrow(inter))
    return(inter)
  })

  # Clean up and return
  trk_inter <- do.call(rbind, trk_inter)
  return(trk_inter)
}

################################################################################
#### Function to split a raster
################################################################################
# Function to split a raster
splitRas <- function(file, s, buffer = 0, overwrite = T, cores = 1, outdir = tempdir()) {

  # Testing
  r <- rast(file)

  # Identify width and height of raster (in pixels)
  x <- ncol(r)
  y <- nrow(r)

  # Generate split windows
  t <- s - 1
  srcwin <- expand.grid(x = 0:t * x / s, y = 0:t * y / s)
  srcwin$xto <- srcwin$x + x / s
  srcwin$yto <- srcwin$y + y / s

  # Add buffer
  srcwin[, 1] <- srcwin[, 1] - buffer
  srcwin[, 2] <- srcwin[, 2] - buffer
  srcwin[, 3] <- srcwin[, 3] + buffer
  srcwin[, 4] <- srcwin[, 4] + buffer

  # Check if values lie outside study area
  srcwin[, 1] <- ifelse(srcwin[, 1] < 0, 0,  srcwin[, 1])
  srcwin[, 2] <- ifelse(srcwin[, 2] < 0, 0,  srcwin[, 2])
  srcwin[, 3] <- ifelse(srcwin[, 3] > x, x,  srcwin[, 3])
  srcwin[, 4] <- ifelse(srcwin[, 4] > y, y,  srcwin[, 4])

  # Change the "to coordinates" to widths
  srcwin$x_width <- srcwin$xto - srcwin$x
  srcwin$y_width <- srcwin$yto - srcwin$y
  srcwin$xto <- NULL
  srcwin$yto <- NULL
  srcwin <- as.matrix(srcwin)

  # Prepare output names
  outnames <- substr(basename(file), start = 1, stop = nchar(basename(file)) - 4)
  outnames <- paste0(outnames, "_Tile_", sprintf("%02d", 1:(s**2)), ".tif")
  outnames <- file.path(outdir, outnames)

  # Check if some of the files already exist
  exists <- file.exists(outnames)

  # Crop to the different tiles
  invisible(mclapply(1:nrow(srcwin), mc.cores = cores, function(x) {
    if (overwrite | !exists[x]) {
      gdal_translate(
          src_dataset = file
        , dst_dataset = outnames[x]
        , srcwin      = srcwin[x, ]
      )
    }
  }))

  # Return the filenames
  return(outnames)
}

################################################################################
#### Function to Run GLMMTMB like Muff & Fieberg
################################################################################
# Fit glmmTMB model
glmm_clogit <- function(formula, data) {

  # Prepare Model call but do not fit
  model <- glmmTMB(formula
    , family  = poisson()
    , data    = data
    , doFit   = FALSE
  )

  # Set the variance of the intercept artificially high
  model$parameters$theta[1] <- log(1e6)

  # Tell glmmTMB not to change the first entry of the vector of variances and
  # give all other variances another indicator to make sure they can be freely
  # estimated
  nvarparm <- length(model$parameters$theta)
  model$mapArg <- list(theta = factor(c(NA, 1:(nvarparm - 1))))

  # Fit the model
  model <- glmmTMB:::fitTMB(model)

  # Return the model
  return(model)
}

################################################################################
#### Function to Prepare Model Calls for Desired Covariates
################################################################################
# Prepare model call
writeForm <- function(x, slope = T, stepmets = T) {

  # Identify the fixed effects structure
  fixed <- paste(x, collapse = " + ")

  if (slope) {

    # Create random slope for each fixed effect
    rando <- paste("(0 + ", x, "|id)", collapse = " + ")

    # Put the random and fixed effects together
    combi <- paste(fixed, rando, sep = " + ")

  } else {

    combi <- fixed

  }

  # Prepare basic iSSF model formula (or just SSF)
  if (stepmets) {
    form <- (case ~
      + cos_ta
      + sl
      + log_sl
      + (1|step_id)
      + (0 + cos_ta|id)
      + (0 + sl|id)
      + (0 + log_sl|id)
    )
  } else {
    form <- (case ~ + (1|step_id))
  }


  # Update formula with the desired covariates
  form <- update(form, paste("~ . + ", combi))

  # Return the final formula
  return(form)
}

################################################################################
#### Function to Extract Coefficients from a glmmTMB result
################################################################################
# Extract glmmTMB Coefficients
getCoeffs <- function(x, zvalue = F, pvalue = F, ranefs = F) {

  # Check out coefficients
  coeffs <- summary(x)$coefficients$cond %>%

    # Coerce the data to a dataframe
    as.data.frame() %>%

    # Keep the row names as own column
    mutate(Covariate = row.names(.)) %>%

    # Rename columns more nicely
    rename(
        Coefficient = Estimate
      , SE          = `Std. Error`
      , zvalue      = `z value`
      , pvalue      = `Pr(>|z|)`
    ) %>%

    # Arrange the variables nicely
    dplyr::select(Covariate, Coefficient, SE, zvalue, pvalue)

  # Add or remove info
  if (!zvalue) {
    coeffs <- dplyr::select(coeffs, -zvalue)
  }
  if (!pvalue) {
    coeffs <- dplyr::select(coeffs, -pvalue)
  }
  if (ranefs) {
    coeffs <- VarCorr(x)$cond %>%
      as.data.frame() %>%
      rename(`(Intercept)` = X.Intercept.) %>%
      pivot_longer(1:ncol(.), names_to = "Covariate", values_to = "RandomVariance") %>%
      mutate(RandomSD = sqrt(RandomVariance)) %>%
      left_join(coeffs, ., by = "Covariate")
  }

  # Return the final dataframe
  return(coeffs)
}

################################################################################
#### Function to Find Best Arrangement of a certain Number of Plots
################################################################################
# Function to find best arrangement of plots
bestMfrow <- function(n) {
  factors <- factors <- (1:n)[n %% (1:n) == 0]
  best_combination <- NULL
  min_difference <- Inf
  for (i in factors) {
    j <- n / i
    difference <- abs(i - j)
    if (difference < min_difference) {
      best_combination <- c(i, j)
      min_difference <- difference
    }
  }
  i <- best_combination[1]
  j <- best_combination[2]
  return(c(i, j))
}

################################################################################
#### Plot Method for SpatRasterCollection
################################################################################
# Method to plot one of our spat raster collections
plot.SpatRasterCollection <- function(x, ...) {
  if (!inherits(x, "SpatRasterCollection")) {
    stop("Input must be a SpatRasterCollection object")
  }
  par(mfrow = bestMfrow(length(x)))
  for (i in 1:length(x)) {
    plot(x[i], main = names(x)[i], ...)
  }
}

################################################################################
#### Function to Prepare an ISSF model (or multiple models) for Prediction
################################################################################
# Function to prepare a model so it can be used for prediction
prepareModel <- function(form, model) {

  # Remove intercept and create vector from coefficients
  coeffs        <- subset(model, Covariate != "(Intercept)")
  coffee        <- coeffs$Coefficient
  names(coffee) <- coeffs$Covariate

  # Get standard errors
  coeffs        <- subset(model, Covariate != "(Intercept)")
  errors        <- coeffs$SE
  names(errors) <- coeffs$Covariate

  # Prepare simplified model formula (no ranefs)
  form    <- terms(form)
  newform <- c("Intercept", "step_id", "0 +")
  newform <- paste(newform, collapse = "|")
  newform <- grep(newform, attr(form, "term.labels"))
  newform <- drop.terms(form, newform, keep.response = F)
  newform <- formula(newform)

  # Return the coefficients and the formula
  result <- list(
      betas   = coffee
    , ses     = errors
    , formula = newform
  )
  return(result)
}

# Function to prepare multiple models for prediction
prepareModels <- function(form, models) {

  # Go through the models and prepare them
  seasons <- unique(models$Season)
  coefs <- lapply(seasons, function(x) {

    # Remove intercept and create vector from coefficients
    coeffs        <- subset(models, Season == x & Covariate != "(Intercept)")
    coffee        <- coeffs$Coefficient
    names(coffee) <- coeffs$Covariate
    return(coffee)
  })
  names(coefs) <- seasons

  # Add standard errors if needed
  se <- lapply(seasons, function(x) {

    # Remove intercept and create vector from coefficients
    coeffs        <- subset(models, Season == x & Covariate != "(Intercept)")
    errors        <- coeffs$SE
    names(errors) <- coeffs$Covariate
    return(errors)
  })
  names(se) <- seasons

  # Prepare simplified model formula (no ranefs)
  form    <- terms(form)
  newform <- c("Intercept", "step_id", "0 +")
  newform <- paste(newform, collapse = "|")
  newform <- grep(newform, attr(form, "term.labels"))
  newform <- drop.terms(form, newform, keep.response = F)
  newform <- formula(newform)

  # Return the coefficients and the formula
  result <- list(
      betas   = coefs
    , ses     = se
    , formula = newform
  )
  return(result)
}


################################################################################
#### Function for Cross-Validation Procedure
################################################################################
# Function for Cross-Validation Procedure
crossVal <- function(formula, dat_train, dat_valid, random = FALSE) {

  # Train the model using the training data and prepare it for prediction
  mod <- glmm_clogit(formula, data = dat_train)
  mod <- getCoeffs(mod)
  mod <- prepareModel(formula, mod)

  # Prepare model matrix and predict selectio nscores
  modeldat         <- model.matrix(mod$formula, data = dat_valid)[, -1, drop = F]
  dat_valid$Scores <- exp(modeldat %*% mod$betas)

  # In case we dont want to randomize, we keep all steps
  validation_pref <- dat_valid %>%
    group_by(., step_id) %>%
    mutate(., Rank = order(order(Scores, decreasing = TRUE))) %>%
    subset(., case == 1) %>%
    dplyr::select(id, step_id, step_id_within, case, Rank) %>%
    mutate(Preferences = "True") %>%
    ungroup()

  # In case we want to randomize, sample one of the random steps (and drop the
  # observed step)
  validation_rand <- dat_valid %>%
    subset(., case == 0) %>%
    group_by(., step_id) %>%
    mutate(., Rank = order(order(Scores, decreasing = TRUE))) %>%
    sample_n(1) %>%
    dplyr::select(id, step_id, step_id_within, case, Rank) %>%
    mutate(Preferences = "Randomized") %>%
    ungroup()

  # Bind them
  valid <- rbind(validation_pref, validation_rand)
  return(valid)

}

# ################################################################################
# #### Function to Prepare Covariates for Prediction
# ################################################################################
# # Function that reduces both the covariate layers and the lookup table to only
# # the timestamps for which they are actually needed
# prepareCovariates <- function(covars, lookup, timestamps, type, readall = F) {
#
#   # Open a connection to all covariates
#   covars_sub <- subset(covars, Type == type) %>%
#     mutate(Raster = map(Filename, stack))
#
#   # Subset the lookup table to the entries for which we need to predict. Also
#   # replace the layerindices with actual layernames
#   lookup_sub <- subset(lookup
#     , Type      == type
#     & Timestamp %in% timestamps
#     ) %>%
#     group_by(Covariate) %>%
#     mutate(Layername = map2_chr(Covariate, Layerindex, function(x, y) {
#       names(covars_sub$Raster[[which(covars_sub$Covariate == x)]])[y]
#     })) %>%
#     ungroup()
#
#   # Go through the covariate rasters and drop all layers that are not needed, as
#   # they dropped from the lookup table. Load the remaining ones into memory.
#   covars_sub$Raster <- lapply(seq_len(nrow(covars_sub)), function(x) {
#     keep <- lookup_sub %>%
#       subset(Covariate == covars_sub$Covariate[x]) %>%
#       pull(Layername) %>%
#       unique()
#     r <- covars_sub$Raster[[x]][[keep]]
#     if (readall) {
#       r <- readAll(r)
#     }
#     return(r)
#   })
#
#   # Put covariates and lookup table into a list
#   result <- list(Covariates = covars_sub, Lookup = lookup_sub)
#   return(result)
#
# }

################################################################################
#### Several Functions to Deal with Covariate Layers Efficiently
################################################################################
# Class containing covariate data
setClass("covariates"
  , slots = c(
      covariates = "data.frame"
    , lookup     = "data.frame"
  )
)

# Constructor function to create a covariate object
covariates <- function(covariates, lookup) {
  object <- new("covariates"
    , covariates = covariates
    , lookup     = lookup
  )
  object@covariates$Raster <- vector("list", nrow(object@covariates))
  return(object)
}

# Method to get the extent
covariatesExtent <- function(covariates) {

  # Check if rasters have been loaded or not
  load <- !all(sapply(covariates@covariates$Raster, class) == "NULL")
  if (!load) {
    stop("Cannot get an extent from a covariates object where rasters haven't been loaded\n")
  }

  # Find the extent (given by the smalles expanse)
  extes <- lapply(covariates@covariates$Raster, function(x) {
    ext(x)
  })
  exte <- ext(c(
      xmin = max(sapply(extes, xmin))
    , xmax = min(sapply(extes, xmax))
    , ymin = max(sapply(extes, ymin))
    , ymax = min(sapply(extes, ymax))
  ))
  return(exte)
}

# Print method
setMethod("show", "covariates", function(object) {
  nams <- unique(object@covariates$Covariate)
  datr <- as.character(range(object@lookup$Timestamp))
  typs <- unique(object@covariates$Type)
  load <- !all(sapply(object@covariates$Raster, class) == "NULL")
  if (load) {
      inme <- all(sapply(object@covariates$Raster, inMemory))
      layr <- sum(sapply(object@covariates$Raster, nlyr))
    } else {
      inme <- F
      layr <- NA
  }
  file <- nrow(object@covariates)
  look <- nrow(object@lookup)
  cat("\n")
  cat("------------------------------------\n")
  cat("--------- covariate object ---------\n")
  cat("------------------------------------\n")
  cat("From:     :", datr[1], "\n")
  cat("To        :", datr[2], "\n")
  cat("Types     :", typs, "\n")
  cat("Loaded    :", load, "\n")
  cat("Memory    :", inme, "\n")
  cat("Files     :", file, "\n")
  cat("Layers    :", layr, "\n")
  cat("Lookup    :", look, "\n")
  cat("------------------------------------\n")
  cat("Covariates:\n -", paste(nams, collapse = "\n - "), "\n")
  cat("------------------------------------\n")
})

# Function to drop layers that are not needed anymore
dropLayers <- function(covariates) {
  covariates@covariates$Raster <- lapply(seq_len(nrow(covariates@covariates)), function(i) {
    keep <- subset(covariates@lookup
      , Type == covariates@covariates$Type[i]
      & Covariate == covariates@covariates$Covariate[i]
    ) %>% pull(Layername) %>% unique()
    keep <- covariates@covariates$Raster[[i]][[keep]]
    return(keep)
  })
  return(covariates)
}

# Function to load covariates as rasters
loadCovariates <- function(covariates, readall = F) {

  # Class check
  if (class(covariates) != "covariates") {
    stop("Please provide an object of the class 'covariates'\n")
  }
  load <- !all(sapply(covariates@covariates$Raster, class) == "NULL")
  if (load) {
    stop("Covariates are already loaded and will not be loaded again\n")
  }

  # Load layers and extract layernames, so we can associate them with layer
  # indices
  covariates@covariates$Raster <- lapply(covariates@covariates$Filename, function(x) {
    cov <- rast(x)
    if (readall) {
      cov <- setValues(cov, cov[])
    }
    return(cov)
  })
  covariates@lookup <- covariates@covariates %>%
    mutate(Names = map(Raster, function(x) {
      df <- data.frame(Layername = names(x), Layerindex = 1:nlyr(x))
    })) %>%
    dplyr::select(Type, Covariate, Names) %>%
    unnest(Names) %>%
    left_join(covariates@lookup, ., by = c("Type", "Covariate", "Layerindex"))

  # Drop layers that are not needed
  covariates <- dropLayers(covariates)

  # Return
  return(covariates)

}

# Method to subset
setMethod("subset", signature(x = "covariates"), function(x, ...) {

  # Check if raster layers are loaded
  loaded <- !all(sapply(x@covariates$Raster, class) == "NULL")

  # Full join
  all <- left_join(x@covariates, x@lookup, by = c("Covariate", "Type"))

  # Need to be careful with scopes
  condition <- substitute(...)
  indices   <- eval(condition, all, parent.frame())
  all       <- all[indices, ]

  # Split again
  x@covariates <- distinct(all[, c("Type", "Covariate", "Filename", "Raster")], Type, Covariate, Filename, .keep_all = T)
  if (loaded) {
      x@lookup <- distinct(all[, c("Timestamp", "Covariate", "Type", "Layerdate", "Layerindex", "Layername")])
    } else {
      x@lookup <- distinct(all[, c("Timestamp", "Covariate", "Type", "Layerdate", "Layerindex")])
  }

  # If layers are loaded, drop the ones that are not needed anymore
  if (loaded) {
    x <- dropLayers(x)
  }
  return(x)
})

# Function to extract covariates
extractCovariates <- function(covariates, xy) {
  extracted <- lapply(covariates@covariates$Raster, function(x) {
    extr <- terra::extract(x, xy)
  }) %>% do.call(cbind, .) %>% as.data.frame()
  names(extracted) <- covariates@covariates$Covariate
  return(extracted)
}

# Wrapper function to prepare covariates
prepareCovariates <- function(covars, timestamps, type, readall = F) {
  covars_sub <- subset(covars, Type %in% type & Timestamp %in% timestamps)
  covars_sub <- loadCovariates(covars_sub, readall = readall)
  return(covars_sub)
}

# Wrapper to simplify access to a specific covariate stack
getCovariates <- function(covars, covariate, timestamps, type) {
  covars_sub <- subset(covars, Covariate == covariate)
  covars_sub <- prepareCovariates(covars_sub, timestamps, type)
  covars_sub <- covars_sub@covariates$Raster[[1]]
  return(covars_sub)
}

# Function to scale covariates given a scaling table
scaleCovars <- function(covars, scaling) {
  scaled <- lapply(1:ncol(covars), function(x) {
    scale(covars[, x]
      , center  = scaling$center[scaling$Covariate == names(covars)[x]]
      , scale   = scaling$scale[scaling$Covariate == names(covars)[x]]
    )
  })
  scaled <- do.call(cbind, scaled)
  scaled <- as.data.frame(scaled)
  names(scaled) <- names(covars)
  return(scaled)
}
