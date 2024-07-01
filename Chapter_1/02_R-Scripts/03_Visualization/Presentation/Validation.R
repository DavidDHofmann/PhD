################################################################################
#### Plot Movement Model Validation
################################################################################
# Description: Plot of the validation of the movement model

# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/Chapter_1"
setwd(wd)

# Load required packages
library(tidyverse)    # For data wrangling
library(lemon)        # For nice axis labels
library(latex2exp)    # For latex expressions
library(msir)         # For prediction interval around loess
library(ggpubr)       # For putting plots together
library(ggdark)       # To quickly change a bright theme to dark

# Load required data
validation     <- read_rds("03_Data/03_Results/99_ModelValidation.rds")
dat_pref       <- read_rds("03_Data/03_Results/99_ModelValidation(Data).rds")[[1]]
dat_rand       <- read_rds("03_Data/03_Results/99_ModelValidation(Data).rds")[[2]]
validation_con <- read_rds("03_Data/03_Results/99_ConnectivityValidation.rds")

# Put the k-fold cross validation data from observed and random preferences
# together
dat_pref$Group <- "Realized"
dat_rand$Group <- "Random"
dat <- rbind(dat_pref, dat_rand)

# We want to plot this data and add the information from the validation table
# too. Let's prepare a column that indicates the text that we want to plot on
# top of the data. Let's first round the values from the validation table
validation[, c(2:4)] <- round(validation[, c(2:4)], 2)

# Match the Grouping names of the validation table to the groups above
validation$Group <- as.factor(c("Realized", "Random"))

# Create a dataframe which we use to annotate the two facets
text <- data.frame(
    Group = c("Realized", "Random")
  , Text = c(
        paste0(
            "$\\bar{r}_s = "
          , validation$Mean[1]
          , "$, $95\\%-CI = "
          , validation$LCL[1]
          , "$, $"
          , validation$UCL[1], "$"
        )
      , paste0(
            "$\\bar{r}_s = "
          , validation$Mean[2]
          , "$, $95\\%-CI = "
          , validation$LCL[2]
          , "$, $"
          , validation$UCL[2], "$"
        )
    #   "$\\bar{r}_s = -0.45$, $95%-CI = -0.48$, $-0.42$"
    # , "$\\bar{r}_s = 0.04$, $95%-CI = 0.01$, $0.08$"
  )
)

# Reorder the factors in the Group variable
dat$Group <- factor(dat$Group, levels = c("Realized", "Random"))

# Get loess
loess <- dat %>%
  group_by(Group) %>%
  nest() %>%
  mutate(data = lapply(data, function(x){
    l <- loess.sd(x$Frequency ~ x$Rank)
    df <- data.frame(
        Loess = l$y
      , Rank  = l$x
      , Upper = l$upper
      , Lower = l$lower
    )
    return(df)
  })) %>%
  unnest(data)

# Plot the data, once for realized, once for random preferences
p1 <- lapply(unique(loess$Group), function(x) {
  sub_loess <- subset(loess, Group == x)
  sub_text <- subset(text, Group == x)
  sub_data <- subset(dat, Group == x)
  ggplot(sub_loess, aes(x = Rank, y = Loess)) +
    geom_jitter(aes(x = Rank, y = Frequency)
      , data  = sub_data,
      , alpha = 0.1
      , size  = 1
    ) +
    geom_ribbon(aes(ymin = Lower, ymax = Upper)
      , alpha = 0.4
      , fill  = "steelblue"
      , col   = "steelblue"
    ) +
    geom_line(
        linewidth = 1
      , col       = "steelblue"
    ) +
    dark_theme_classic() +
    theme(
        panel.background = element_blank()
      , plot.background  = element_blank()
    ) +
    coord_capped_cart(
        left   = "both"
      , bottom = "both"
      , ylim   = c(0, 50)
    ) +
    geom_text(data = sub_text
      , mapping = aes(
          x = -Inf
        , y = -Inf
        , label = TeX(Text, output = "character")
      )
      , hjust   = -0.05
      , vjust   = -0.5
      , parse   = TRUE
      , size    = 3
    ) +
    ylab("Frequency")
})

# Let's add some more information to the validation of our connectivity map
validation_con <- validation_con %>%
  subset(ModelType == "PSF") %>%
  dplyr::select(Steps, Model) %>%
  unnest(Model) %>%
  setNames(c("Steps", "Coefficient", "SE", "pvalue")) %>%
  mutate(
      LCI_90 = Coefficient - qnorm(1 - (1 - 0.90) / 2) * SE
    , UCI_90 = Coefficient + qnorm(1 - (1 - 0.90) / 2) * SE
    , LCI_95 = Coefficient - qnorm(1 - (1 - 0.95) / 2) * SE
    , UCI_95 = Coefficient + qnorm(1 - (1 - 0.95) / 2) * SE
    , LCI_99 = Coefficient - qnorm(1 - (1 - 0.99) / 2) * SE
    , UCI_99 = Coefficient + qnorm(1 - (1 - 0.99) / 2) * SE
  )

# Add stars indicating the significance
validation_con$Significance <- sapply(1:nrow(validation_con), function(x){
  if (validation_con$pvalue[x] <= 0.01){
    return("***")
  } else if (validation_con$pvalue[x] <= 0.05){
    return("**")
  } else if (validation_con$pvalue[x] <= 0.1){
    return("*")
  }
})

# Prepare dataset for plotting confidence intervals
validation_confs <- validation_con %>%
  dplyr::select(Steps, Coefficient, LCI_90:UCI_99) %>%
  gather(key = confidence_level, value = value, LCI_90:UCI_99) %>%
  separate(col = confidence_level, into = c("Type", "Level"), sep = "_") %>%
  spread(key = Type, value = value) %>%
  mutate(Level = paste0(Level, "%"))

# Plot of path selection model
p1[[3]] <- ggplot(data = validation_con, aes(y = Coefficient, x = as.factor(Steps))) +
  geom_point(shape = 3, size = 2.5, col = "steelblue") +
  geom_errorbar(
      aes(
        ymin = LCI
      , ymax = UCI
      , size = factor(Level)
    )
    , data   = validation_confs
    , width  = 0
    , col    = "steelblue"
    , alpha  = 0.5
  ) +
  geom_text(
      aes(label = Significance, hjust = 0.5, vjust = -0.5, angle = 90)
    , show.legend = F
  ) +
  geom_hline(
      yintercept = 0
    , color      = "darkgrey"
    , lty        = 2
    , lwd        = 0.3
  ) +
  dark_theme_classic() +
  ylim(c(-1, 3)) +
  coord_capped_cart(
      left   = "both"
    , bottom = "both"
  ) +
  labs(
      x = "Number of Simulated Steps"
    # , y = expression(beta*"-Coefficient_{Connectivity}")
    , y = parse(text = TeX("$\\beta-Coefficient_{Connectivity}$"))
  ) +
  scale_size_manual(
      name   = "Confidence Level"
    , values = c(2, 1, 0.3)
  ) +
  theme(
    , legend.position   = "none"
    , legend.text       = element_text(face = 3)
    , legend.title      = element_text(face = 3)
    , panel.background  = element_blank()
    , plot.background   = element_blank()
  ) +
  guides(size = guide_legend(title.position = "top", title.hjust = 0.5, legend.hjust = 4))

# Put plots together
p2 <- ggarrange(p1[[1]], p1[[2]], p1[[3]]
  , ncol          = 3
  , label.x       = -0.10
)

# Store plot
ggsave(plot = p2, "06_Presentation/Validation.png"
  , bg        = "transparent"
  , width     = 10
  , height    = 3
  , scale     = 1.2
)
