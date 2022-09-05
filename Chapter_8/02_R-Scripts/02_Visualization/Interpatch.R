################################################################################
#### Table of Inter-Patch Connectivity
################################################################################
# Clear R's brain
rm(list = ls())

# Change the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_8"
setwd(wd)

# Load required packages
library(tidyverse)

# Load results on interpatch connectivity
dat <- read_rds("03_Data/03_Results/BootstrappedInterpatchConnectivity.rds")
dat <- subset(dat, FloodLevel != "Mean")

################################################################################
#### Frequency
################################################################################
test <- dat %>%
  select(-c(StepNumber, StepNumberSE)) %>%
  # mutate(Combined = paste(round(dat$Freq), "$\\pm$", round(dat$FreqSE))) %>%
  mutate(Combined = paste(round(dat$Freq), round(dat$FreqSE))) %>%
  select(-c(Freq, FreqSE)) %>%
  pivot_wider(names_from = SourceArea, values_from = Combined) %>%
  arrange(CurrentArea, desc(FloodLevel)) %>%
  select(CurrentArea, FloodLevel, everything())

library(kableExtra)
kbl(test, booktabs = T, format = "latex")

xtable(test) %>% print.xtable(type = "latex", sanitize.text.function = function(x){x})

levels <- unique(dat$FloodLevel)
tables <- list()
for (i in levels) {
  tables[[i]] <- dat %>%
    subset(FloodLevel == i) %>%
    select(-c(StepNumber, StepNumberSE, FloodLevel)) %>%
    nest(Metrics = c(Freq, FreqSE)) %>%
    pivot_wider(names_from = SourceArea, values_from = Metrics)
}

library(xtable)
dataset1 = matrix( c(1,2,3,4, 5, 6 ), 2 , 3)
dataset2 = matrix( c(4,3,5,10, 1, 0), 2, 3)

dataset <- rbind(dataset1,dataset2)  #combine dataset

means <- apply(dataset,1,mean)  #calculate row means
sds <- apply(dataset,1,sd)	#calculate row standard deviation

msd <- paste(round(means,2)," (",round(sds,2),")",sep="") #mean and

rn <- c("Var1","Var2")  #rownames
cn <- c("Dataset1","Dataset2") #column names
tab <- matrix(msd,2,2,dimnames=list(rn,cn))

xtable(tab)

################################################################################
#### Steps
################################################################################

sub <- subset(dat, FloodLevel == "Min")
sub <- select(sub, From, To, Freq)

t1 <- pivot_wider(sub, names_from = From, values_from = Freq)
t2 <- pivot_wider(sub, names_from = From, values_from = Freq)
