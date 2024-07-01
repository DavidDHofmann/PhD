################################################################################
#### Spearman's Rank Correlation for Different Number of "Case" Observaitons
################################################################################
# Description: Assess how Spearman's rank correlation depends on the number
# of random steps per stratum

# Clear R's brain
rm(list = ls())

# Set working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_3"
setwd(wd)

# Load required packages
library(tidyverse)
library(survival)
library(pbmcapply)

# Load custom functions
source("02_R-Scripts/00_Functions.R")

################################################################################
#### Helper Functions
################################################################################
# Function to predict a step-selection score
score <- function(covariates, formula, betas) {
  mat   <- model.matrix(formula, covariates)[, -1, drop = F]
  scori <- exp(as.vector(mat %*% betas))
  return(scori)
}

# Function to compute the rank frequency
rank <- function(dat) {
  ranks <- dat %>%
    arrange(id, desc(score_predicted)) %>%
    group_by(id) %>%
    mutate(rank = 1:n()) %>%
    ungroup() %>%
    subset(case == 1) %>%
    count(rank) %>%
    arrange(rank)
  return(ranks)
}

# Function to simulate some data
simdat <- function(n_case, n_rand, form, betas) {
  n <- n_case * (n_rand + 1)
  df <- tibble(
      id = rep(1:n_case, length.out = n)
    , x1 = rnorm(n)
    , x2 = rnorm(n)
  ) %>% arrange(id)
  df$score <- score(df, form, betas)
  df <- df %>%
    nest(data = -id) %>%
    mutate(data = map(data, function(x) {
        x$case <- 0
        x$case[sample(nrow(x), size = 1, prob = x$score)] <- 1
        return(x)
    })) %>%
    unnest(data) %>%
    arrange(id, desc(case)) %>%
    group_by(id) %>%
    mutate(within_id = 1:n()) %>%
    ungroup()
  return(df)
}

################################################################################
#### Test Functions
################################################################################
# Simulation parameters
n_case <- 100
n_rand <- 5
form   <- ~ x1 + x2
betas  <- c(-0.5, 0.5)

# Simulate
sim <- simdat(n_case, n_rand, form, betas)

# Run conditional logistic regression to try and obtain the betas
mod <- clogit(case ~ x1 + x2 + strata(id), sim)
summary(mod)

# Predict scores based on covariates
sim$score_predicted <- score(sim, ~ x1 + x2, coef(mod))

# Check rank correlation
ran <- rank(sim)
rho <- cor(x = ran$rank, y = ran$n, method = "spearman")

# Visualize
plot(n ~ rank, ran, type = "o", pch = 20, main = paste0("Number Random Steps = ", n_rand))
mtext(paste0("rho = ", round(rho, 2)))

################################################################################
#### Prepare Visualization
################################################################################
# Simulation parameters
n_case <- 100
n_rand <- c(5, 10, 25, 50, 75, 100)
form   <- ~ x1 + x2
betas  <- c(-0.5, 0.5)
reps   <- 1:100

# Simulate data
sim <- simdat(n_case, max(n_rand), form, betas)

# Run conditional logistic regression to try and obtain the betas (let's assume
# we miss the second covariate)
mod <- clogit(case ~ x1 + x2 + strata(id), sim)

# Predict a score
sim$score_predicted <- score(sim, ~ x1 + x2, coef(mod))

# Check how the correlation changes when we consider different amounts of
# random observations
ran <- lapply(n_rand, function(i) {
  ranks <- rank(subset(sim, within_id <= i))
  ranks$nsteps <- i
  return(ranks)
}) %>% do.call(rbind, .)

# Visualize
p1 <- ggplot(ran, aes(x = rank, y = n)) +
  geom_point() +
  facet_wrap(~ nsteps) +
  sm_statCorr(corr_method = "spearman", color = "cornflowerblue") +
  theme_minimal() +
  theme(strip.background =  element_rect(fill = "gray95", color = "transparent")) +
  theme_awesome() +
  xlab("Rank") +
  ylab("Frequency") +
  theme(
      axis.title.y     = element_text(angle = 90)
    , panel.grid.minor = element_blank()
  )

# Store the plot to file
ggsave("04_Manuscript/Figures/RankCorrelationNumberSteps.png"
  , plot   = p1
  , device = png
  , bg     = "white"
  , width  = 5
  , height = 3
  , scale  = 1.8
)
