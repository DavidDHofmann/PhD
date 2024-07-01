################################################################################
#### Glossary
################################################################################
# Description: Preparing the table for the glossary

# Clear R's brain
rm(list = ls())

# Load required packages
library(tidyverse)    # To wrangle data
library(kableExtra)   # To generate a nice table

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_4"
setwd(wd)

# Load glossary
dat <- read_csv("Glossary.csv") %>%
  arrange(Order) %>%
  select(Term, Definition)

# Make a nice table
kbl(dat, booktabs = T, format = "latex", escape = T, align = "ll", caption = c("Glossary of terms. Terms in the glossary are printed in bold at first occurrence in the main text. Definitions are always given in the context of step selection functions (SSFs)."), label = "Glossary") %>%
  kable_styling(latex_options = c("scale_down", "striped")) %>%
  column_spec(1, bold = T) %>%
  column_spec(2, width = "10 cm") %>%
  writeLines("04_Manuscript/99_Glossary.tex")
