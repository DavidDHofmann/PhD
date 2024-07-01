################################################################################
#### Covariate Table
################################################################################
# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Load required packages
library(tidyverse)    # To wrangle data
library(kableExtra)   # To generate a nice table

# Load covariate table
dat <- read_csv("03_Data/01_RawData/CovariateDescription.csv")
dat <- dat[, -1]

# Some beautifying
names(dat) <- gsub(names(dat), pattern = " ", replacement = "\n")
names(dat) <- linebreak(names(dat), align = "c")
dat$Description <- NULL

# Make a nice table
kbl(dat, booktabs = T, format = "latex", escape = F, align = "lcccc") %>%
  # kable_styling(latex_options = "scale_down") %>%
  pack_rows("(1) Landscape Characteristics", 1, 9, bold = F) %>%
  pack_rows("(2) Climate Descriptors", 10, 11, bold = F) %>%
  pack_rows("(3) Anthropogenic Features", 12, 14, bold = F) %>%
  pack_rows("(4) Light Availability", 15, 16, bold = F) %>%
  row_spec(12:14, background = "#f2f2f2") %>%
  row_spec(4:6, background = "#f2f2f2") %>%
  row_spec(15:16, background = "#f2f2f2") %>%
  # footnote(
  #     general        = "Layers highlighted in gray were combined into a single proxy for human-influence."
  #   , threeparttable = T
  # ) %>%
  writeLines("04_Manuscript/Figures/Covariates.tex")
