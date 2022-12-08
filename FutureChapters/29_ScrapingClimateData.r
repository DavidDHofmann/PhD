############################################################
#### Automatic Download of Climate Data in Botswana
############################################################
# Clear R's brain
rm(list = ls())

# Load packages
library(tidyverse)
library(rvest)
library(lubridate)

# Function to scrape data for a desired date
date <- "2021-01-31"
scrapethat <- function(date) {

  # Specify fixed and variable part of the url
  url_fix <- "http://www.sasscalweathernet.org/weatherstat_hourly_we.php?loggerid_crit=67585&date_crit_daily="
  url_var <- ymd(date)
  url <- paste0(url_fix, url_var)

  # Scrape url
  dat <- url %>%
    read_html() %>%
    html_nodes(xpath = '//*[@id="inb_1"]/table[2]') %>%
    html_table()

  # Clean data
  suppressWarnings(
    dat <- dat[[1]] %>%
      setNames(paste0(names(.), slice(., 1))) %>%
      slice(-1) %>%
      mutate(Date = ymd(date)) %>%
      mutate(across(3:18, as.numeric))
  )

  # Return the final dataframe
  return(dat)

}

# Try it
scrapethat("2021-12-01")
