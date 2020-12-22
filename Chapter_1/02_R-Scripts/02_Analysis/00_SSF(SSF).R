################################################################################
#### Step Selection Function - Generation of Random Steps
################################################################################
# Description: In this script we coerce our gps data to steps and generate
# random steps for the (time varying) integrated step selection analysis.

# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load packages
library(tidyverse)    # For data wrangling
library(davidoff)     # Custom functions
library(lubridate)    # To handle dates
library(amt)          # To coerce gps fixes to steps
library(raster)       # To handle spatial data
library(rgdal)        # To handle spatial data
library(pbmcapply)    # For multicore abilities

# Set seed for reproducability
set.seed(1234)

# Big loop over the two methods
methods <- c("iSSF", "TiSSF")
for (i in 1:length(methods)){

  # Select method
  method <- methods[i]

################################################################################
#### Data Cleaning
################################################################################

  # Load the gps data
  data <- read_csv("03_Data/02_CleanData/00_General_Dispersers_POPECOL.csv")

  # We only need dispersers' data
  data <- subset(data, State == "Disperser")

  # Let's create a timestamp that is rounded to the nearest hour
  data$TimestampRounded <- round_date(data$Timestamp, "1 hour")

  # Make sure there are no duplicates!
  table(duplicated(data[, c("DogName", "TimestampRounded")]))
  table(hour(data$TimestampRounded))

  # Count number of dogs
  length(unique(data$DogName))

  # Nest data by dog
  data <- data %>% group_by(DogName) %>% nest()

  # Resample to two hours. Can't use 1 hourly fixes because step speeds are
  # unproportionally large with 1 hour.
  data$data <- suppressMessages(
    pbmclapply(1:nrow(data)
      , ignore.interactive = T
      , mc.cores           = detectCores() - 1
      , FUN                = function(x){
        resFix2(data$data[[x]], hours = 2, start = 1, tol = 0.5)
      }
    )
  )
  data <- unnest(data)

  # In case we're using iSSF, we will only keep the fixes that were collected on
  # our pre-defined fixrate schedule. Otherwise, if using TiSSF, we're going to
  # keep 2-hourly fixes.
  if (method == "iSSF"){
    data <- subset(data, hour(TimestampRounded) %in% c(3, 7, 15, 19, 23))
    } else {
    data <- subset(data, hour(TimestampRounded) %in% seq(1, 23, 2))
  }

  # We're going to convert the fixes to steps, yet we want to indicate a new
  # burst if a step takes longer than 4.25 hours (8.25 hours between 07:00 and
  # 15:00). Let's thus calculate the timelags.
  data <- data %>%
    group_by(DogName) %>%
    mutate(dt_ = Timestamp - lag(Timestamp)) %>%
    mutate(dt_ = as.numeric(dt_, units = "hours")) %>%
    ungroup()

  # Indicate when a new burst starts (a new burst also starts whenever dt_ = NA)
  if (method == "iSSF"){
    data <- data %>%
      mutate(NewBurst = ifelse(
        dt_ > 8.25 |
        (dt_ > 4.25 & hour(TimestampRounded) != 15) |
        is.na(dt_), yes = 1, no = 0)
      ) %>%
      mutate(BurstID = cumsum(NewBurst)) %>%
      dplyr::select(-c(NewBurst))
    } else {
    data <- data %>%
      mutate(NewBurst = ifelse(
          dt_ > 8.25 |
          is.na(dt_), yes = 1, no = 0)
      ) %>%
      mutate(BurstID = cumsum(NewBurst)) %>%
      dplyr::select(-c(NewBurst))
  }

  # We can only work with bursts that contain at least three fixes
  data <- data %>%
    group_by(DogName, BurstID) %>%
    nest() %>%
    mutate(Nrow = map(data, nrow) %>% do.call(rbind, .)) %>%
    subset(Nrow >= 3) %>%
    unnest(data) %>%
    ungroup() %>%
    dplyr::select(-c(Nrow, dt_))

  # Create bursts that are unique
  data$BurstID <- data %>% group_indices(DogName, BurstID)

  # Now we can coerce the data to proper steps. We'll trick amt::make_track by
  # using the burstID as animalID
  tracks <- data %>%
    make_track(.
      , .x      = x
      , .y      = y
      , .t      = Timestamp
      , id      = BurstID
      , crs     = CRS("+init=epsg:4326")
      , State   = State
      , DogName = DogName
    ) %>%

    # Transform the tracks to utm
    transform_coords(CRS("+init=epsg:32734")) %>%

    # Nest the tracks
    nest(-"id") %>%

    # Turn to a step representation (look up the amt vignette for details)
    mutate(data = map(data, function(x){
      x %>%

        # The function creates steps from the fixes in each tibble row. The
        # option "keep_cols" allows to keep the time of day (tod) column that we
        # added
        steps(keep_cols = "start") %>%

        # Transform the difftime column to a numeric column to avoid that there
        # are heterogeneous units
        mutate(., dt_ = as.numeric(dt_, units = "hours"))
    })) %>%

    # Unnest the tibble
    unnest(data) %>%

    # Multiply turning angles with negative one (for some reason the package
    # calculates the turning angles conuterclockwise but we want them clockwise)
    mutate(ta_ = ta_ * (-1)) %>%

    # Add a column indicating the absolute turning angle (important for the
    # ellipses). We can use the function we created above.
    mutate(absta_ = absAngle(.)) %>%

    # We can only work with steps for which we have a turning angle. Let's get
    # rid of any steps where the turning angle is NA
    filter(!is.na(ta_))

  # Generate unique step id and indicate that steps are observed (case_) steps
  tracks <- tracks %>% mutate(
      step_id_ = 1:nrow(.)
    , case_    = 1
  )

  # Add a row indicating the method
  tracks$method <- method

  # Indicate if a fix was taken during time of activity or inactivity (we define
  # inactivity as anything between 07:00 and 15:00)
  tracks <- tracks %>% mutate(
    inactive =
      hour(round_date(t1_, "1 hour")) >= 7 &
      hour(round_date(t2_, "1 hour")) <= 15 &
      hour(round_date(t1_, "1 hour")) <
      hour(round_date(t2_, "1 hour"))
  )

  # Use the "inactive" indicator to calculate corrected step times (i.e.
  # assuming that 8 hours correspond to 4 hours during the time of inactivity)
  tracks <- tracks %>% mutate(dtcorr_ = ifelse(inactive, dt_ / 2, dt_))

  # Use the corrected step times to calculate the (corrected) step speed. Note
  # that the step speed is almost identical to the (normalized) step length
  # under iSSF. Yet under TiSSF the step lengths are not comparable because of
  # irregular sampling.
  tracks <- tracks %>% mutate(speed_ = sl_ / dtcorr_)

################################################################################
#### Generation of Random Steps
################################################################################
  # We can't work with 0 step lengths or step speeds
  tracks$sl_[tracks$sl_ == 0] <- 1
  tracks$speed_[tracks$speed_ == 0] <- 1

  # Fit a Gamma distribution to the step speed (again, fitting a gamma to step
  # speeds or fitting a gamma to (normalized) step lengths yields the same under
  # iSSF)
  sl <- fit_distr(tracks$speed_, "gamma")

  # Let's visualize the fit
  x <- seq(0, max(tracks$speed_), 1)
  y <- dgamma(x, scale = sl$params$scale, shape = sl$params$shape)
  hist(tracks$speed_, freq = F, breaks = 100)
  lines(y ~ x, col = "red")

  # Write the distribution to file
  if (method == "TiSSF"){
      write_rds(sl, "03_Data/03_Results/99_GammaDistribution(TiSSF).rds")
    } else {
      write_rds(sl, "03_Data/03_Results/99_GammaDistribution(iSSF).rds")
  }

  # Define the number of random steps
  nsteps <- 24

  # Now we use the fitted gamma distribution to sample step lengths. For the
  # turning angles we use a uniform distribution. Note that we need to scale the
  # step lengths depending on the (corrected) step time.
  randomSteps <- pbmclapply(
      X                  = 1:nrow(tracks)
    , mc.cores           = detectCores() - 1
    , ignore.interactive = T
    , FUN                = function(x){

    # Draw random turning angles
    ta_new <- runif(nsteps, min = -pi, max = pi)

    # Check out the step duration (corrected)
    duration <- tracks$dtcorr_[x]

    # Draw random step lengths with the fitted parameters. Note that I'm going
    # to adjust the scale parameter so that the distribution matches the
    # (corrected) step duration
    sl_new <- rgamma(nsteps
      , scale = sl$params$scale * duration
      , shape = sl$params$shape
    )

    # Put the step lengths and turning angles into a new dataframe and calculate
    # the new endpoints of each random step. We also indicate that the steps are
    # control steps (i.e. 0)
    randomSteps <- tracks[rep(x, nsteps), ] %>%
      mutate(.
        , absta_  = absta_ + (ta_new - ta_)
        , ta_     = ta_new
        , sl_     = sl_new
        , case_   = 0
        , x2_     = x1_ + sin(absta_) * sl_
        , y2_     = y1_ + cos(absta_) * sl_
        , speed_  = sl_ / dtcorr_
      )

    # Return the sampled steps
    return(randomSteps)

  }) %>% do.call(rbind, .)

  # We need to make sure that the absolute turning angle ranges from 0 to 2 * pi
  randomSteps$absta_[randomSteps$absta_ > 2 * pi] <-
    randomSteps$absta_[randomSteps$absta_ > 2 * pi] - 2 * pi
  randomSteps$absta_[randomSteps$absta_ < 0] <-
    randomSteps$absta_[randomSteps$absta_ < 0] + 2 * pi

  # Merge the dataframes of the observed and random steps
  ssf <- rbind(tracks, randomSteps) %>%
    arrange(step_id_, desc(case_)) %>%
    transform(case_ = as.logical(case_))

  # Prepare a filename
  if (method == "iSSF"){
    filename <- "00_General_Dispersers_POPECOL(iSSF)"
    } else {
    filename <- "00_General_Dispersers_POPECOL(TiSSF)"
  }

  # Store to file
  write_csv(ssf, paste0("03_Data/02_CleanData/", filename, ".csv"))

  # Coerce the steps to spatial lines. We can use the function we defined
  # earlier for this
  lines <- lineTrack(ssf, CRS("+init=epsg:32734"))

  # Transform the lines to WGS84
  lines <- spTransform(lines, CRS("+init=epsg:4326"))

  # Coerce the duration column to a format that can be stored (numeric rather
  # than difftime)
  lines$dt_ <- as.numeric(lines$dt_)

  # Store the lines
  writeOGR(lines
    , dsn       = "03_Data/02_CleanData"
    , layer     = filename
    , driver    = "ESRI Shapefile"
    , overwrite = TRUE
  )

  # End big loop
  cat("Data preparation for", method, "finished\n")

}
