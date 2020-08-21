############################################################
#### Automatic Download of Climate Data in Botswana
############################################################
# Clear R's brain
rm(list = ls())

# Set the working directory
input <- "/home/david/Downloads/Popecol/csv"
output1 <- "/home/david/ownCloud/University/14. FS 19/Masterarbeit/03_Data/00_LargeExtent"
output2 <- "/home/david/ownCloud/University/14. FS 19/Masterarbeit/03_Data/01_SmallExtent"
setwd(input)

# Load packages
library(tidyverse)
library(rvest)
library(lubridate)

# Specify the fixed part of the url
url_fix <- "http://www.sasscalweathernet.org/weatherstat_hourly_we.php?loggerid_crit=67585&date_crit_daily="

# Specify the variable part of the url. In our case this is only the date. We
# want to loop through all dates since the first GPS recording in our data. The
# first GPS fix we have is from 31.12.2011. Unfortunately the weather station
# only recorded data from 30.03.2015 on. Let's download from then
url_var <- seq(from = as.Date("2015-05-30"), to = Sys.Date() - 1, by = "day")
url <- paste0(url_fix, url_var)

# Define the scraping function
scrapethat <- function(x){
  page <- x
  data <- page %>%
    read_html() %>%
    html_nodes("#borderlinie_dott") %>%
    html_text() %>%
    as.data.frame(stringsAsFactors = FALSE)
}

# Run the function on the desired URLs
download <- lapply(url, scrapethat)

# Keep a backup
backup <- download
download <- backup

# We only want to keep the first entry of each dataframe
for (i in 1:length(download)){
  download[[i]] <- download[[i]][1, ]
}

# Put them together into a dataframe
download <- do.call(rbind, download) %>% as.data.frame()

# Split the remaining text into several columns
split <- list()
for (i in 1:nrow(download)){

  # Get the text in the i'th row and turn it into a character row
  row <- as.character(download[i, 1])

  # Replace some of the spacing characters
  row <- gsub("-\r\n\r\n", "\r\n", row)

  # Split the characters according to the remaining spacing characters
  row <- strsplit(row, split = "\r\n")[[1]][c(2:18, 38:445)] %>%
    matrix(., ncol = 17, byrow = TRUE) %>%
    as.data.frame(., stringsAsFactors = FALSE)

  # Use the first row as column header
  names(row) <- t(row[1, ])

  # Remove the first row (since the column header is now set)
  row <- row[-1, ]

  # We can also retrieve the sunrise and sunset times
  row$Sunrise <- substr(download[i, 1], start = 335, stop = 339)
  row$Sunset <- substr(download[i, 1], start = 340, stop = 344)

  # Finally we want to keep track of the dates
  row$Date <- url_var[i]

  # Put the final dataframe into a list
  split[[i]] <- row
}

# For some dates there is no data available. The resulting dataframe is thus a
# bit weird and messes up our structure. Let's remove this data. We can do so by
# checking the column names. If the first column name is not equal to "Hour", we
# can remove the entry
keep <- c()
for (i in 1:length(split)){
  keep[i] <- names(split[[i]][1]) %in% "Hour"
}
split <- split[keep]

# Collapse the list
final <- do.call(rbind, split)

# Coerce the first 17 columns to numeric (last one is the date and thus not numeric)
final[, 1:17] <- sapply(final[, -18], as.numeric)

# Make the column names slightly nicer, i.e. remove the dots
names(final) <- names(final) %>% gsub("\\.", "", .)

# Store the dataframe to file
write.csv(final, "WeatherData.csv")
setwd(output1)
write.csv(final, "00_General_WeatherData_SASSCAL.csv")
setwd(output2)
write.csv(final, "00_General_WeatherData_SASSCAL.csv")
