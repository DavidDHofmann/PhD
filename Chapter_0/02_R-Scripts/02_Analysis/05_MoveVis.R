############################################################
#### Create Movement Animation
############################################################
# Clear R's brain
rm(list = ls())

# Set working directory
setwd("/home/david/ownCloud/University/14. FS 19/Masterarbeit/03_Data/01_SmallExtent")

# Load packages
library(tidyverse)
library(moveVis)
library(raster)
library(move)
library(lubridate)
library(viridis)

# Let's first check the maps that we can get using the package, as well as the formats that it suggests to use
get_maptypes()
suggest_formats()

# # Load a background map for plotting (I think we will actually use another map)
# map <- raster("01_LandCoverClasses30_Globeland.tif")

# Load GPS tracks
tracks <- read_csv("00_General_Dispersers_Popecol(Regular).csv") %>%

  # Coerce the stuff to a dataframe because the move package can't handle tibbles
  as.data.frame() %>%

  # We want to keep only those individuals that eventually dispersed
  subset(., DogName %in% unique(.$DogName[.$State == "Disperser"])) %>%

  # We will not keep Amacuro since her track is not nice to look at and because it will make the animation way longer than necessary
  subset(., DogName != "Amacuro")

# For each disperser we want to keep 30-Days of the pre-dispersal phase (if available) and the dispersal phase itself. We do not keep the gps locations after dispersal. To achieve this we could use our cutoff table or we simply identify the earliest and latest date at which each dog was a "Disperser". Let's go for the second approach.
first <- aggregate(Timestamp ~ DogName, data = subset(tracks, State == "Disperser"), FUN = min) %>% set_names(., c("DogName", "FirstDate"))
last <- aggregate(Timestamp ~ DogName, data = subset(tracks, State == "Disperser"), FUN = max) %>% set_names(., c("DogName", "LastDate"))

# Put the dates together into a dataframe
dates <- left_join(first, last)

# Make some changes
tracks <- dates %>%

  # Subtract 30 days from the first date
  mutate(., FirstDate = FirstDate - days(30)) %>%

  # Add 30 days to the last date
  mutate(., LastDate = LastDate + days(30)) %>%

  # Join the table with the tracks data
  left_join(tracks, ., by = "DogName") %>%

  # Keep only the fixes that fall into the period between the First and LastDate
  subset(., Timestamp >= FirstDate & Timestamp <= LastDate)

# Prepare a bounding box which we will use for all animations
bbox <- extent(c(22, 27, -20.7, -17.7)) %>% as(., "SpatialPolygons")
crs(bbox) <- CRS("+init=epsg:4326")

############################################################
#### Animation for All Dispersers
############################################################
# We want to plot all dispersers at the same time. To achieve this, we need to manipulate their timestamps. In principle we move all dates to the same starting point. However, due to the fact that we do not have 30 days of pre-dispersal data for all dispersing individuals we need to make sure that the day at which they start to disperse align. We can achieve this by identifying and adding the time difference between the start of the dispersal, to an arbitrary date (let's say 01.01.2019) for each dog
dates$Difference <- as.POSIXct("2019-01-01 00:00:00", tz = "UTC") - dates$FirstDate

# Let's join the time differences that we want to add to each dog
allDisp <- left_join(tracks, dates[, c("DogName", "Difference")], by = "DogName")

# Now add the time difference to the timestamp of each gps fix
allDisp$Timestamp <- allDisp$Timestamp + allDisp$Difference

# We can use the data to prepare an animation now. First, we coerce the data to a proper move object
move <- move(
    x       = allDisp$x
  , y       = allDisp$y
  , time    = allDisp$Timestamp
  , data    = allDisp
  , animal  = allDisp$DogName
  , proj    = CRS("+init=epsg:4326")
)

# As you can see the final object is a movestack, where each layer is a move object for one individual
show(move)

# Then we interpolate the gps fixes to a desired resolution.
m <- align_move(move, res = 2, digit = 0, unit = "hours")

# To have some nicer colour for our plots we use the colours provided by the viridis package. We need to get as many colours as we have individuals
cols <- viridis(length(unique(m@trackId)))

# We also want to add a counter that indicates the day. We cannot really use the date since we changed the dates artificially
day <- unique(timestamps(m)) %>% as.Date() %>% sort()
day <- day - day[1]
day <- as.numeric(day)
day <- as.character(paste0("Day: ", sprintf("%03d", day)))

# Now we prepare the frames for the animation
frames <- frames_spatial(
    m
  , trace_colour  = "black"
  , path_colours  = cols
  , path_legend   = FALSE
  , path_size     = 4
  , tail_size     = 1
  , tail_colour   = "white"
  , ext           = extent(bbox)
  , map_service   = "mapbox"
  , map_type      = "satellite"
  , map_token     = "pk.eyJ1IjoiZG9keDkiLCJhIjoiY2p3dnltejJjMGR4YjN5bXp0ZjA2ZXBzMCJ9.4hirgQ-1SfJ2KHI7SR54cQ") %>%
  add_progress() %>%
  add_northarrow(colour = "white", size = 3) %>%
  add_scalebar(colour = "white", distance = 50) %>%
  add_labels(x = "Longitude", y = "Latitude") %>%
  add_text(day, x = 22.5, y = -17.9, type = "label", size = 8)

# Look at a desired frame to make sure everything looks as desired. Note that the objects will all be misplaced somehow. This is because scaling will be changed when we specify the width and height manually.
frames[[50]]

# Store the animation
animate_frames(
    frames
  , width     = 1920
  , height    = 1080
  , out_file  = "DispersalEvents.mov"
  , overwrite = TRUE
)

############################################################
#### Animations for Everest, Abel and MadameChing
############################################################
# We also want to prepare three individual animations. We chose the individuals to be Everest, Abel and Madame Ching, since their tracks are rather nice to look at. In addition to their tracks we want to visualize the number of days since dispersal, as well as the cumulative distance travelled. This is rather challenging because we need to calculate these two metrics for every frame. Anyways, let's first create a dataframe in which we can store information that will be used in a loop.
names <- c("Everest", "Abel", "MadameChing")

# Loop through the individuals and prepare animations
for (i in 1:length(names)){

  # Subset to the desired individual
  sub <- subset(tracks, DogName == names[i])

  # Coerce the track to a move object
  move <- move(
      x       = sub$x
    , y       = sub$y
    , time    = sub$Timestamp
    , data    = sub
    , proj    = CRS("+init=epsg:4326")
  )

  # Align (interpolate) the track.
  m <- align_move(move, res = 1, digit = 0, unit = "hours")

  # Since we only work with one individual we don't need the movestack. Let's "unnest" the stack.
  m <- split(m)

  # This turned the stack into a list and we need to extract the move object from it
  m <- m[[1]]

  # Let's create a vector that indicates during which fixes the individual was dispersing
  dispersing <- m$time >= dates$FirstDate[dates$DogName == names[i]] & m$time <= dates$LastDate[dates$DogName == names[i]]

  # Using this we can create a colour vector that we will use in the animation. This will allow us to plot the dispersal phase in a different colour
  cols <- rep("green", length(dispersing))
  cols[dispersing] <- "red"

  # The colours need to be assigned to the move object
  m$colour <- cols

  # Now that we interpolated the track we can calculate the desired metrics, i.e. the days since start of dispersal, as well as the cumulative distance travelled so far. Let's first get the number of days since dispersal.
  days <- as.Date(m$time) - as.Date(dates$FirstDate[dates$DogName == names[i]])

  # Any gps fix that is taken before or after dispersal should not count as dispersal day. Lets make them zero.
  days[!dispersing] <- 0

  # Actually we want that the counter keeps the last dispersal day until the animation is finished. Let's replace all after dispersal days to the maximal value
  days[m$time > dates$LastDate[dates$DogName == names[i]]] <- max(days)

  # Now we also calculate the distance. We can calculate it easily by converting the move object to an ltraj object, which automatically calculates the distance
  dist <- m %>%

    # Transform the track to utm
    spTransform(., CRS("+init=epsg:32734")) %>%

    # Coerce the move object to an ltraj object
    as(., "ltraj") %>%

    # Select the first dataframe in the list
    .[[1]] %>%

    # Select the column for distance in the datafram
    .[["dist"]]

  # Note that the ltraj conversion calculates the distances at the beginning of each step. However, we want them at the end of each step. Let's thus move the vector by one position
  dist <- c(NA, dist[-length(dist)])

  # Again we want the distance before dispersal to remain 0. We can use the vector for days since dispersal for making this possible
  dist[!dispersing] <- 0

  # Ultimately we calculate the cumulative distances and convert them to kilometers
  dist <- round(cumsum(dist) / 1000, 2)

  # Coerce the days vector to a character vector and add some text. Note that we make sure the day is given in 3 digits. This makes sure that the length of the text remains constant so that the textbox is aligned nicely throughout the animation
  days <- paste0("Dispersal-Day: ", sprintf("%03d", days))

  # Coerce the distance vector to a character vector and add some text. As before we make sure the distance is give in three digits
  dist <- paste0("Dispersal-Distance: ", sprintf("%03d", round(dist)), " km")

  # Now we can prepare the animation
  frames <- frames_spatial(
      m
    , tail_colour   = c("yellow")
    , trace_colour  = c("yellow")
    , trace_show    = TRUE
    , path_legend   = FALSE
    , path_size     = 5
    , tail_size     = 2
    , ext           = extent(bbox)
    , map_service   = "mapbox"
    , map_type      = "satellite"
    , map_res       = 1
    , map_token     = "pk.eyJ1IjoiZG9keDkiLCJhIjoiY2p3dnltejJjMGR4YjN5bXp0ZjA2ZXBzMCJ9.4hirgQ-1SfJ2KHI7SR54cQ") %>%

    # Add additional objects to the plot as desired
    add_timestamps(m, type = "label", size = 10) %>%
    add_northarrow(colour = "white", position = "upperright", size = 3) %>%
    add_scalebar(colour = "white", position = "bottomright", distance = 50) %>%
    add_text(days, x = 22.48, y = -20.4, type = "label", size = 8) %>%
    add_text(dist, x = 22.665, y = -20.55, type = "label", size = 8) %>%
    add_labels(x = "Longitude", y = "Latitude") %>%
    add_progress()

  # Look at a frame to make sure everything is correct
  frames[30]

  # Animate the frames
  animate_frames(
      frames
    , out_file = paste0(names[i], "Animation.mov")
    , overwrite = TRUE
    , width = 1980
    , height = 1080
    , end_pause = 2
  )}
