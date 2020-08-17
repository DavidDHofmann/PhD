################################################################################
#### Preparation of all Plots for Chapter I
################################################################################
# Clear R's brain
rm(list = ls())

# Set the working directory
wd <- "/home/david/ownCloud/University/15. PhD/00_WildDogs"
setwd(wd)

# Set a seed
set.seed(12345)

# Surpress scientific notation
options(scipen = 999)

# Load required packages
library(rgdal)        # To load and store spatial data
library(raster)       # To manipulate raster data
library(tmap)         # For beautiful maps
library(tidyverse)    # For data wrangling
library(glmmTMB)      # To handle glmm models
library(viridis)      # For nice colors
library(lemon)        # Some ggplot additions
library(cowplot)      # Further ggplot additions
library(imager)       # To import images
library(Cairo)        # To store pdf images (correctly)
library(rgeos)        # To manipulate spatial data
library(grid)         # To overlay ggplots with images
library(ggdark)       # A dark ggplot theme
library(rasterVis)    # For beautiful spatial plots
library(gridExtra)    # To arrange trellis/levelplots
library(xtable)       # To store nice latex tables
library(parallel)     # For parallel computing
library(ggsignif)     # To plot significance stars in ggplot

# Load custom functions
source("Functions.r")

################################################################################
#### Colors
################################################################################
# Create a color ramp
pal1 <- colorRampPalette(colors = c("black", "orange", "yellow"))

# Visualize colors
plot(1:20, 1:20, pty = 20, pch = 20, cex = 15, col = pal1(20))

################################################################################
#### Plot of the Study Area
################################################################################
# Load required data
africa <- "03_Data/02_CleanData/00_General_Africa.shp" %>% readOGR()
kaza <- "03_Data/02_CleanData/00_General_KAZA_KAZA.shp" %>% readOGR()
dogs <- "03_Data/02_CleanData/00_General_WildDogs_IUCN.shp" %>% readOGR()
prot <- "03_Data/02_CleanData/02_LandUseTypes_Protected_PeaceParks(1Class).shp" %>% readOGR()
water <- "03_Data/02_CleanData/03_LandscapeFeatures_MajorWaters_GEOFABRIK.shp" %>% readOGR(.)

# Remove small islands for plotting
africa <- subset(africa, !(ID %in% c(27:41, 689, 690)))

# Create a buffered polygon of africa
africa2 <- gBuffer(africa, width = metersToDegrees(100))

# Get the extent of the KAZA
kaza_ext <- as(extent(kaza), "SpatialPolygons")
crs(kaza_ext) <- CRS("+init=epsg:4326")

# Create labels for countries
labels_countries <- data.frame(
    x = c(16, 24, 11, 38, 34)
  , y = c(-9, -30, -26, -13, -23)
  , Label = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")
)
coordinates(labels_countries) <- c("x", "y")
crs(labels_countries) <- CRS("+init=epsg:4326")

# Create lines pointing towards these countries
l1 <- rbind(c(17, -10), c(19, -12)) %>% spLines()
l2 <- rbind(c(24, -29), c(24, -23)) %>% spLines()
l3 <- rbind(c(13, -25), c(17, -22)) %>% spLines()
l4 <- rbind(c(35, -13), c(30, -14)) %>% spLines()
l5 <- rbind(c(33, -22), c(30, -19)) %>% spLines()
lines_countries <- rbind(l1, l2, l3, l4, l5)
crs(lines_countries) <- crs(labels_countries)

# Prepare a plot of the KAZA. Create labels for countries first.
labels_countries2 <- data.frame(
    x = c(20.39, 23.94, 20.07, 25.69, 28.22)
  , y = c(-15.28, -19.94, -19.39, -15.22, -18.9)
  , Label = c("Angola", "Botswana", "Namibia", "Zambia", "Zimbabwe")
)
coordinates(labels_countries2) <- c("x", "y")
crs(labels_countries2) <- CRS("+init=epsg:4326")

# Create labels for some geographical landmarks
labels_waters <- data.frame(
    x     = c(22.6, 23.7, 27.8, 25.6)
  , y     = c(-19, -18.2, -17, -20.7)
  , Label = c(
    "Okavango\nDelta", "Linyanti\nSwamp", "Lake\nKariba", "Makgadikgadi\nPans"
  )
)
coordinates(labels_waters) <- c("x", "y")
crs(labels_waters) <- CRS("+init=epsg:4326")

# Create labels for some national parks
labels_nationalparks <- data.frame(
    x = c(26.56, 28.61, 21.15, 25.87, 20.38, 23.58, 23.71, 24.51, 20.78)
  , y = c(-19.08, -17.05, -17.26, -14.66, -16.08, -21.4, -19.29, -18.65, -18.81)
  , Label = paste0(c(
      "Hwange", "Matusadona", "Luengue-Luiana", "Kafue", "Mavinga"
    , "Central Kalahari", "Moremi", "Chobe", "Khaudum"
  ), "\nNP")
)
coordinates(labels_nationalparks) <- c("x", "y")
crs(labels_nationalparks) <- CRS("+init=epsg:4326")

# Prepare a map of Africa
p1 <- tm_shape(africa2) +
    tm_polygons(col = "gray70", border.col = "gray70", lwd = 2) +
  # tm_grid(
  #     n.y                 = 5
  #   , n.x                 = 5
  #   , labels.inside.frame = false
  #   , lines               = true
  #   , ticks               = true
  #   , col = "gray12"
  # ) +
  tm_shape(africa) +
    tm_polygons(
        col = "gray40"
      , lwd = 0.7
      , border.col = "gray70"
    ) +
  tm_shape(dogs) +
    tm_polygons(
        col           = "orange"
      , alpha         = 0.8
      , border.alpha  = 0
    ) +
  tm_shape(kaza) +
    tm_borders(
        col = "white"
      , lty = 1
      , lwd = 1
    ) +
  tm_shape(kaza_ext) +
    tm_borders(
        col = pal1(20)[10]
      , lty = 3
      , lwd = 1.5
    ) +
  tm_shape(lines_countries) +
    tm_lines(
        col = "white"
    ) +
  tm_shape(labels_countries) +
    tm_text(
        "Label"
      , size = 0.6
      , col = "white"
    ) +
  tm_layout(
      asp         = 0.8
    , frame       = "black"
    , frame.lwd   = 3
    , legend.show = FALSE
    , bg.color    = "black"
  ) +
  tm_credits("(a)"
    , position  = c("right", "top")
    , size      = 1.5
    , col       = "white"
)

# Prepare a map of Kaza
p2 <- tm_shape(prot) +
    tm_polygons(
      , col           = "gray20"
      , border.col    = "black"
      , border.alpha  = 0.8
      , lwd           = 0.5
      , legend.show   = F
    ) +
  tm_shape(water) +
    tm_polygons(
        col           = "cornflowerblue"
      , border.col    = "cornflowerblue"
      , border.alpha  = 0.6
      , lwd           = 0.2
    ) +
  tm_shape(kaza, is.master = T) +
    tm_borders(
        col = "white"
      , lty = 1
      , lwd = 2
    ) +
  tm_shape(africa) +
    tm_borders(
        col = "gray70"
    ) +
  tm_shape(labels_countries2) +
    tm_text("Label"
      , col       = "gray70"
      , fontface  = 2
      , size      = 1.5
    ) +
  tm_shape(labels_waters) +
    tm_text("Label"
      , fontface  = 3
      , size      = 0.5
      , col       = "white"
    ) +
  tm_shape(labels_nationalparks) +
    tm_text("Label"
      , col       = "gray70"
      , fontface  = 3
      , size      = 0.5
    ) +
  tm_grid(
      n.y                 = 5
    , n.x                 = 5
    , labels.inside.frame = FALSE
    , lines               = FALSE
    , ticks               = TRUE
  ) +
  tm_layout(
    , bg.color = "black"
    , frame                   = "gray20"
    , frame.lwd               = 3
    , asp                     = 1.2
    , legend.outside          = TRUE
    , legend.outside.position = "left"
    , legend.text.color       = "white"
    , legend.title.color      = "black"
    , legend.stack            = "vertical"
    , legend.text.size        = 0.8
    , legend.bg.color         = "black"
  ) +
  tm_scale_bar(
        position  = c("right", "bottom")
      , text.size = 0.5
      , text.col  = "white"
      , width     = 0.125
  ) +
  tm_credits("(b)"
    , position  = c("left", "top")
    , size      = 1.5
    , col       = "white"
  ) +
  tm_compass(
      color.dark  = "white"
    , color.light = "white"
    , text.color  = "white"
    , position    = c("left", "bottom")
)

################################################################################
#### Prepare a common Legend
################################################################################
# We now prepare a legend that the two plots share. However, we have to assign
# to one of the two plots. Let's use p2 for this
p2 <- p2 + tm_add_legend(
    type    = "fill"
  , labels  = c(
      "Wild Dog Populations"
    , "Major Water Areas"
    , "Protected Areas"
  )
  , col = c(
      "orange"
    , "cornflowerblue"
    , "gray20"
  )) + tm_add_legend(
    type    = "line"
  , labels  = c("KAZA-TFCA Borders", "Country Borders")
  , col     = c("white", "gray70")
  , lty     = c(1, 1)
  , lwd     = c(2, 2)
) + tm_layout(legend.frame.lwd = 2, legend.text.size = 1.05)

# Combine plots and store them
CairoPDF("05_Manuscript2/99_StudyArea.pdf", width = 9, height = 5.25, bg = "black")
p2
print(p1, vp = viewport(0.158, 0.341, width = 0.56, height = 0.56))
dev.off()

################################################################################
#### Movement Model AICs
################################################################################
# Load required data
aics <- read_rds("03_Data/03_Results/99_MovementModelAICs.rds")

# Also remove the ModelID column
aics$ModelID <- NULL

# Replace covariates with short code and add cos(ta) and log(sl) as covariatess
aics$Covariates <- aics$Covariates %>%
  gsub(pattern = "\\bWater\\b", replacement = "W") %>%
  gsub(pattern = "\\bShrubs\\b", replacement = "S") %>%
  gsub(pattern = "\\bTrees\\b", replacement = "T") %>%
  gsub(pattern = "\\bProtected\\b", replacement = "P") %>%
  gsub(pattern = "\\bHumansBuff5000\\b", replacement = "HI") %>%
  gsub(pattern = "\\bDistanceToWater\\b", replacement = "DTW") %>%
  gsub(pattern = "\\bDistanceToRoads\\b", replacement = "DTR") %>%
  gsub(pattern = "\\bRoadCrossing\\b", replacement = "RC") %>%
  gsub(pattern = "\\bActivity\\b", replacement = "MA") %>%
  gsub(pattern = ",", replacement = " +") %>%
  paste("cos(ta) + log(sl) +", .)

# We will also add cos(ta) and log(sl) as covariates
aics

# Write the table to a .tex table
print(xtable(aics)
  , floating            = FALSE
  , latex.environments  = NULL
  , booktabs            = TRUE
  , include.rownames    = FALSE
  , type                = "latex"
  , file                = "05_Manuscript2/99_MovementModelAICs.tex"
)

################################################################################
#### Movement Model Results
################################################################################
# Load the movement model
move.mod <- "03_Data/03_Results/99_MovementModel.rds" %>% readRDS()

# Select best model
best <- move.mod$Model[[1]]

# Summary of the best model
summary(best)

# Calculate p-values
coeffs <- getCoeffs(best, pvalue = TRUE)[-1, ] %>%

  # Calculate confidence intervals
  mutate(
      LCI = Coefficient - 1.96 * SE
    , UCI = Coefficient + 1.96 * SE
  )

# Add stars indicating the significance
coeffs$Significance <- sapply(1:nrow(coeffs), function(x){
  if (coeffs$pvalue[x] <= 0.01){
    return("***")
  } else if (coeffs$pvalue[x] <= 0.05){
    return("**")
  } else if (coeffs$pvalue[x] <= 0.1){
    return("*")
  }
})

# Rename covariates
coeffs$Covariate[coeffs$Covariate == "Shrubs"] <- "Shrubs/Grassland"
coeffs$Covariate[coeffs$Covariate == "Trees"] <- "Trees"
coeffs$Covariate[coeffs$Covariate == "Water"] <- "Water"
coeffs$Covariate[coeffs$Covariate == "DistanceToWater"] <- "DistanceToWater"
coeffs$Covariate[coeffs$Covariate == "cos(ta_)"] <- "cos(ta)"
coeffs$Covariate[coeffs$Covariate == "log(sl_)"] <- "log(sl)"
coeffs$Covariate[coeffs$Covariate == "HumansBuff5000"] <- "HumanInfluence"
coeffs$Covariate[coeffs$Covariate == "log(sl_):ActivityMainActivity"] <- "log(sl):MainActivity"
coeffs$Covariate[coeffs$Covariate == "log(sl_):Water"] <- "log(sl):Water"
coeffs$Covariate[coeffs$Covariate == "log(sl_):Trees"] <- "log(sl):Trees"
coeffs$Covariate[coeffs$Covariate == "cos(ta_):DistanceToWater"] <- "cos(ta):DistanceToWater"
coeffs$Covariate[coeffs$Covariate == "cos(ta_):HumansBuff5000"] <- "cos(ta):HumanInfluence"
coeffs$Preference <- ifelse(coeffs$Coefficient > 0, "Preferred", "Avoided")
coeffs$Preference <- factor(coeffs$Preference, levels = c("Preferred", "Avoided"))

# Specify the order in which the coefficients should be plotted
order <- c(
      "cos(ta)"
    , "cos(ta):HumanInfluence"
    , "cos(ta):DistanceToWater"
    , "log(sl)"
    , "log(sl):MainActivity"
    , "log(sl):Water"
    , "log(sl):Trees"
    , "Water"
    , "DistanceToWater"
    , "Trees"
    , "Shrubs/Grassland"
    , "HumanInfluence"
)

# Specify colors of axis labels
labcols <- c("black", "orange")[c(2, 1, 2, 1, 2, 1, 1, 1, 1, 1, 2, 1)]

# Prepare plot with Covariates on the y-axis and the corresponding
# coefficients on the x-axis
p1 <- ggplot(data = coeffs, aes(y = Covariate, x = Coefficient, col = factor(Preference))) +
  geom_point(shape = 1, size = 2) +
  geom_errorbarh(aes(
      xmin = Coefficient - 1.96 * SE
    , xmax = Coefficient + 1.96 * SE
    , height = 0.2)
  ) +
  geom_text(aes(label = Significance, hjust = 0.5, vjust = -0.2), show.legend = F) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey") +
  scale_y_discrete(limits = rev(order)) +
  theme_cowplot() +
  xlim(c(-1, 1)) +
  coord_capped_cart(left = "both", bottom = "both") +
  labs(x = expression(beta*"-Coefficient")) +
  scale_color_manual(values = c("orange", "black")) +
  theme(
      panel.grid.minor = element_line(size = 0.0)
    , panel.grid.major = element_line(size = 0.0)
    , panel.border     = element_blank()
    , axis.line        = element_line()
    , axis.ticks       = element_line(colour = "black")
    , legend.title     = element_blank()
    , axis.text.y      = element_text(color = rev(labcols))
  )

# Load a picture of a wild dog to the plot
dog <- load.image("/home/david/ownCloud/University/15. PhD/99_Various/test.svg")

# Add it to the previous plot
p <- ggdraw() +
  draw_image(dog, x = 0.9, y = 0.7, hjust = 0.5, vjust = 0.5, scale = 0.2) +
  draw_plot(p1)

# Store the plot
CairoPDF("05_Manuscript2/99_MovementModel.pdf", width = 8, height = 3.5)
p
dev.off()

################################################################################
#### Visualize Interactions
################################################################################
# # Prepare interaction graph for interaction between step length and activity
# # phase. Get the data used to train the model
# data <- model.frame(best)
# xVar = "log(sl_)"
# yVar = "MainActivity"
#
# # Coerce the activity phase to a binary variable indicating times of main
# # activity
# data$MainActivity <- data$Activity == "MainActivity"
#
# # Get the beta estimates from the model
# coeffs <- getCoeffs(best)
#
# # Prepare a sequence for which we want to predict the selection score
# seqX <- seq(min(data$`log(sl_)`), max(data$`log(sl_)`), length = 100)
#
# # Write a function to predict values based on x and y
# predictVals <- function(x, y){
#     x * coeffs$Coefficient[coeffs$Covariate == xVar] +
#     x * y * coeffs$Coefficient[coeffs$Covariate == paste0(xVar, ":", yVar)]
# }
#
# # Apply the function to the seqX and seqY vectors
# predicted_Main <- predictVals(seqX, 1)
# predicted_Low <- predictVals(seqX, 0)
#
# # Put values together
# predictions <- rbind(
#     data.frame(Score = predicted_Main, SL = seqX, Activity = "Main")
#   , data.frame(Score = predicted_Low, SL = seqX, Activity = "Low")
# )
#
# # Exponentiate values if desired
# if (exponentiate){
#   predictions$Score <- exp(predictions$Score)
# }
#
# # Normalize values if desired
# normalize <- function(x){
#   (x - min(x)) / (max(x) - min(x))
# }
# if (norm){
#   predictions$Score <- normalize(predictions$Score)
#   predictions$SL <- normalize(predictions$SL)
# }
#
# # Plot
# p0 <- ggplot(predictions, aes(x = SL, y = Score, col = Activity)) +
#   geom_line()

# Prepare rasters for contourplots
r1 <- visInt2(best
  , "log(sl_)"
  , "Water"
  , exponentiate  = F
  , norm          = T
)
r2 <- visInt2(best
  , "cos(ta_)"
  , "DistanceToWater"
  , exponentiate  = F
  , norm          = T
)
r3 <- visInt2(best
  , "cos(ta_)"
  , "HumansBuff5000"
  , exponentiate  = F
  , norm          = T
)
r4 <- visInt2(best
  , "log(sl_)"
  , "Trees"
  , exponentiate  = F
  , norm          = T
)

# Prepare the plots
p1 <- levelplot(r1
  , contour = TRUE
  , xlab    = "log(sl)"
  , ylab    = "Water"
  , margin  = FALSE
  , main    = "(b1)"
  , col.regions = colorRampPalette(c("black", "orange"))(100)
)
p2 <- levelplot(r2
  , contour = TRUE
  , xlab    = "cos(ta)"
  , ylab    = "DistanceToWater"
  , margin  = FALSE
  , main    = "(b2)"
  , col.regions = colorRampPalette(c("black", "orange"))(100)
)
p3 <- levelplot(r3
  , contour = TRUE
  , xlab    = "cos(ta)"
  , ylab    = "Human Influence"
  , margin  = FALSE
  , main    = "(b3)"
  , col.regions = colorRampPalette(c("black", "orange"))(100)
)
p4 <- levelplot(r4
  , contour = TRUE
  , xlab    = "log(sl)"
  , ylab    = "Trees"
  , margin  = FALSE
  , main    = "(b4)"
  , col.regions = colorRampPalette(c("black", "orange"))(100)
)

# Store them
CairoPDF("05_Manuscript2/99_MovementModel(Interactions).pdf"
  , width     = 14
  , height    = 4
)
grid.arrange(p1, p2, p3, p4, nrow = 1)
dev.off()

################################################################################
#### Plot All Movement Models
################################################################################
# Load all movement models
move.mod <- "03_Data/03_Results/99_MovementModel.rds" %>% readRDS()

# Extract their coefficients
coeffs <- move.mod$Model %>%
  lapply(function(x){
    getCoeffs(x)[-1, ]
  }) %>%
  bind_rows(.id = "Model")

# Check all covariates present in the models
covars <- data.frame(Covariate = unique(coeffs$Covariate))

# For each model identify which covariates were not included
coeffs <- coeffs %>%
  group_by(Model) %>%
  nest() %>%
  mutate(data = map(data, function(x){
    joined <- full_join(x, covars)
    joined$Model <- x$Model[1]
    return(joined)
  })) %>%
  unnest() %>%
  ungroup() %>%
  mutate(
      Coefficient = replace_na(Coefficient, 0)
    , Model       = factor(Model, levels = c(1:length(unique(Model))))
  )

# Create a column indicating if the coefficient is present in the model or not
coeffs$InModel <- factor(coeffs$Coefficient != 0, levels = c(TRUE, FALSE))

# Prepare Confidence intervals
coeffs <- coeffs %>%

  # Calculate confidence intervals
  mutate(
      LCI = Coefficient - 1.96 * SE
    , UCI = Coefficient + 1.96 * SE
  )

# Rename covariates
coeffs$Covariate[coeffs$Covariate == "Shrubs"] <- "Shrubs/Grassland"
coeffs$Covariate[coeffs$Covariate == "Trees"] <- "Trees"
coeffs$Covariate[coeffs$Covariate == "Water"] <- "Water"
coeffs$Covariate[coeffs$Covariate == "DistanceToWater"] <- "DistanceToWater"
coeffs$Covariate[coeffs$Covariate == "cos(ta_)"] <- "cos(ta)"
coeffs$Covariate[coeffs$Covariate == "log(sl_)"] <- "log(sl)"
coeffs$Covariate[coeffs$Covariate == "HumansBuff5000"] <- "HumanInfluence"
coeffs$Covariate[coeffs$Covariate == "log(sl_):ActivityMainActivity"] <- "log(sl):MainActivity"
coeffs$Covariate[coeffs$Covariate == "cos(ta_):Water"] <- "cos(ta):Water"
coeffs$Covariate[coeffs$Covariate == "log(sl_):Water"] <- "log(sl):Water"
coeffs$Covariate[coeffs$Covariate == "cos(ta_):Shrubs"] <- "cos(ta):Shrubs/Grassland"
coeffs$Covariate[coeffs$Covariate == "log(sl_):Shrubs"] <- "log(sl):Shrubs/Grassland"
coeffs$Covariate[coeffs$Covariate == "cos(ta_):Trees"] <- "cos(ta):Trees"
coeffs$Covariate[coeffs$Covariate == "log(sl_):Trees"] <- "log(sl):Trees"
coeffs$Covariate[coeffs$Covariate == "cos(ta_):DistanceToWater"] <- "cos(ta):DistanceToWater"
coeffs$Covariate[coeffs$Covariate == "log(sl_):DistanceToWater"] <- "log(sl):DistanceToWater"
coeffs$Covariate[coeffs$Covariate == "cos(ta_):HumansBuff5000"] <- "cos(ta):HumanInfluence"
coeffs$Covariate[coeffs$Covariate == "log(sl_):HumansBuff5000"] <- "log(sl):HumanInfluence"

# Specify the order in which the coefficients should be plotted
order <- c(
      "cos(ta)"
    , "log(sl)"
    , "HumanInfluence"
    , "Shrubs/Grassland"
    , "Trees"
    , "Water"
    , "DistanceToWater"
    , "cos(ta):HumanInfluence"
    , "cos(ta):Shrubs/Grassland"
    , "cos(ta):Trees"
    , "cos(ta):Water"
    , "log(sl):MainActivity"
    , "log(sl):HumanInfluence"
    , "log(sl):Shrubs/Grassland"
    , "log(sl):Trees"
    , "log(sl):Water"
    , "log(sl):DistanceToWater"
)

# Prepare plot with Covariates on the y-axis and the corresponding
# coefficients on the x-axis
p <- ggplot(
      data = coeffs
    , aes(y = Covariate, x = Coefficient, col = InModel, shape = InModel)
  ) +
  geom_point(size = 1.7) +
  geom_errorbarh(aes(
      xmin = Coefficient - 1.96 * SE
    , xmax = Coefficient + 1.96 * SE
    , height = 0.4)
  ) +
  geom_vline(
      xintercept  = 0
    , linetype    = "dashed"
    , color       = "darkgrey"
  ) +
  theme_cowplot() +
  xlim(c(-1, 1)) +
  coord_capped_cart(
      left    = "both"
    , bottom  = "both"
  ) +
  labs(x = expression(beta*"-Coefficient")) +
  scale_color_manual(values = c("black", "orange")) +
  scale_shape_manual(values = c(1, 4)) +
  scale_y_discrete(limits = rev(order)) +
  labs(
      col   = "Covariate in Model:"
    , shape = "Covariate in Model:"
  ) +
  theme(
      panel.grid.minor = element_line(size = 0.0)
    , panel.grid.major = element_line(size = 0.0)
    , panel.border     = element_blank()
    , axis.line        = element_line()
    , axis.ticks       = element_line(colour = "black")
    , legend.position  = "bottom"
    , legend.direction = "horizontal"
    , strip.background = element_rect(fill = "white", colour = "black")
  ) +
  facet_wrap("Model", nrow = 2)

# Store the plot
CairoPDF("05_Manuscript2/99_AllMovementModels.pdf", width = 12, height = 7)
p
dev.off()

################################################################################
#### Plot Source Areas and Source Points
################################################################################
# Load required data
areas1  <- readOGR("03_Data/03_Results/99_SourceAreas.shp")
areas2  <- readOGR("03_Data/03_Results/99_SourceAreas2.shp")
points  <- readOGR("03_Data/03_Results/99_SourcePoints2.shp")
kaza    <- readOGR("03_Data/02_CleanData/00_General_KAZA_KAZA.shp")

# Randomize source points
head(points)
pointsr <- createPoints(areas = areas2, points = points, n = 100, randomize = T)

# Plot of static source points
p1 <- tm_shape(areas1) +
    tm_borders(
        col = "gray70"
      , lwd = 3
    ) +
  tm_shape(areas2) +
    tm_polygons(
        col         = "gray40"
      , border.col  = "black"
    ) +
  tm_shape(points) +
    tm_dots(
        col = "orange"
      , size = 0.1
    ) +
  tm_shape(kaza) +
    tm_borders(
        col = "white"
      , lwd = 1
      , lty = 2
    ) +
  tm_layout(
      frame       = "black"
    , frame.lwd   = 3
    , legend.show = FALSE
    , bg.col      = "black"
  ) +
  tm_scale_bar(
        position    = "left"
      , text.size   = 0.5
      , width       = 0.125
      , text.color  = "white"
  ) +
  tm_compass(
      color.light = "white"
    , color.dark  = "white"
    , text.color  = "white"
  ) +
  tm_credits("(a)"
    , position  = c("left", "top")
    , size      = 1.5
    , col       = "white"
)

# Plot of random source points
p2 <- tm_shape(areas1) +
    tm_borders(
        col = "gray70"
      , lwd = 3
    ) +
  tm_shape(areas2) +
    tm_polygons(
        col         = "gray40"
      , border.col  = "black"
    ) +
  tm_shape(pointsr) +
    tm_dots(
        col = "orange"
      , size = 0.01
    ) +
  tm_shape(kaza) +
    tm_borders(
        col = "white"
      , lwd = 1
      , lty = 2
    ) +
  tm_layout(
      frame       = "black"
    , frame.lwd   = 3
    , legend.show = FALSE
    , bg.col      = "black"
  ) +
  tm_scale_bar(
        position    = "left"
      , text.size   = 0.5
      , width       = 0.125
      , text.color  = "white"
  ) +
  tm_compass(
      color.light = "white"
    , color.dark  = "white"
    , text.color  = "white"
  ) +
  tm_credits("(b)"
    , position  = c("left", "top")
    , size      = 1.5
    , col       = "white"
)

# Put plots together
p <- tmap_arrange(p1, p2, ncol = 2)

# Store the plot
CairoPDF("05_Manuscript2/99_SourcePoints.pdf", width = 12, height = 5.25)
p
dev.off()

################################################################################
#### Heatmaps
################################################################################
# Let's create a table that shows bhattacharyya's affinity between all heatmaps
# and the permeability or the corridor map, as well as the affinity between
# different heatmaps. Load required data
rasterized <- read_rds("03_Data/03_Results/99_RasterizedSimulations.rds")
heatvsheat <- read_rds("03_Data/03_Results/99_HeatmapComparisons.rds")

# Create table that we want to fill
dat <- expand.grid(
      sampling    = c("Static", "Random")
    , comparison  = c("Permeability", "Corridor", "Heatmap")
) %>% arrange(sampling, comparison)
dat[, 3:8] <- NA
dat <- setNames(dat
  , c("Point Sampling", "Comparison", "68", "125", "250", "500", "1000", "2000")
)

heatmaps <- stack("03_Data/03_Results/99_RasterizedSimulations.tif")
# Connect to correct heatmap
for (i in 1:nrow(rasterized)){
  rasterized$heatmap[[i]] <- heatmaps[[i]]
}

# Collapse heatmaps into a single stack
all <- do.call(c, rasterized$heatmap)
all <- stack(all)

# Prepare list indicating the peculiarities of each of the plots
plots <- list()
steps <- rep(c(68, 125, 250, 500, 1000, 2000), 2)
sampling <- c(rep("static", 6), rep("random", 6))

# Prepare all plots
for (i in 1:nlayers(all)){
  plots[[i]] <- tm_shape(all[[i]]) +
    tm_raster(
        palette     = "-Spectral"
      , style       = "cont"
      , legend.show = F
    ) +
    tm_layout(
        title       = paste0("steps: ", steps[i], "\npoint sampling: ", sampling[i])
      , title.color = "white"
    )
}

# Arrange plots
p <- tmap_arrange(
    # plots[[1]]
    plots[[2]]
  # , plots[[3]]
  , plots[[4]]
  # , plots[[5]]
  , plots[[6]]
  # , plots[[7]]
  , plots[[8]]
  # , plots[[9]]
  , plots[[10]]
  # , plots[[11]]
  , plots[[12]]
  , ncol = 3
)

# Store the arranged plot
# CairoPDF("05_Manuscript2/99_Heatmaps.pdf", width = 21, height = 6, bg = "black")
# p
# dev.off()
CairoPDF("05_Manuscript2/99_Heatmaps.pdf", width = 10.6, height = 6, bg = "black")
p
dev.off()

################################################################################
#### Heatmap Metrics
################################################################################
# Load required data
metrics1 <- read_rds("03_Data/03_Results/99_HeatmapMetrics.rds")
metrics2 <- read_rds("03_Data/03_Results/99_HeatmapMetrics2.rds") %>%
  mutate(Correlation = as.numeric(format(round(Correlation, 2), nsmall = 2))) %>%
  mutate(Affinity = as.numeric(format(round(Affinity, 2), nsmall = 2)))

# Summarize bootstrapped affinity metrics
metrics1 <- metrics1 %>%
  dplyr::select(steps, sampling, BhattacharyyaAffinityPerm, BhattacharyyaAffinityCorr) %>%
  gather(key = Map, value = Correlation, 3:4) %>%
  group_by(steps, sampling, Map) %>%
  summarize(Mean = mean(Correlation), SD = sd(Correlation)) %>%
  mutate(
      Mean    = as.numeric(format(round(Mean, 2), nsmall = 2))
    , Metric  = "Affinity"
    , SD      = as.numeric(format(round(SD, 2), nsmall = 2))
    , Map     = substr(Map, start = 22, 25)
  ) %>%
  ungroup()

# Clean heatmap comparisons
metrics2 <- metrics2 %>%
  mutate(
      sampling  = "n.a."
    , Metric    = "Affinity"
    , Map       = "Heat"
    , SD        = NA
    , Mean      = Affinity
  ) %>% select(steps, sampling, Map, Mean, SD, Metric)

# Put metrics together and do some cleaning
metrics <- rbind(metrics1, metrics2) %>%
  subset(Metric == "Affinity") %>%
  mutate(Group = factor(paste("Heat vs.", Map), levels = c("Heat vs. Perm", "Heat vs. Corr", "Heat vs. Heat")))

# Visualize them
ggplot(metrics, aes(x = steps, y = Mean, col = Group)) +
  geom_point(aes(shape = sampling)) +
  geom_line(aes(lty = sampling)) +
  geom_errorbar(aes(ymin = Mean - SD, ymax = Mean + SD)) +
  dark_theme_minimal() +
  scale_color_manual(values = viridis(20)[c(15, 5, 20)]) +
  scale_shape_manual(values = c(2, 4, 6)) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin())

################################################################################
#### Plot of Areas Reached
################################################################################
# Load required data
reached <- read_rds("03_Data/03_Results/99_AreasReached.rds")
areas <- shapefile("03_Data/03_Results/99_SourceAreas2.shp")
sims <- read_rds("03_Data/03_Results/99_DispersalSimulationSub.rds")

# Function to visualize reached areas
visCon <- function(
      n               = 1     # Source point
    , simulations     = NULL  # Simulated Trajectories
    , reached         = NULL  # Tibble indicating areas reached
    , areas           = NULL  # Source area shapefile
    , alpha           = 0.5   # Opacity of tracks
    , steps           = 68    # Number of steps to consider
    , tracks          = F     # Plot tracks?
    , exclude_source  = F     # Source point included?
    , sampling        = c("Static", "Random")
  ){

  # First of all, we need to know the trajectories that reach other areas. These
  # are stored in the "reached" tibble. Let's subset to the repsective entry
  table <- reached[reached$steps == steps & reached$sampling == sampling, ]
  table <- table$AreasReached[[1]]

  # We also want to be able to plot the source area, so lets identify it
  source <- areas[areas$ID == n, ]

  # Get all areas that are reached from this source area
  reach <- subset(areas, ID %in% table$EndPoint[table$StartPoint == n])

  # Maybe the source area should be excluded?
  if (exclude_source){
    reach <- subset(reach, ID != n)
  }

  # Visualize everything
  p <- tm_shape(areas) +
      tm_polygons(col = "gray20", border.col = "black", lwd = 0.2) +
    tm_shape(reach) +
      tm_polygons(col = "orange", border.col = "black", lwd = 0.2) +
    tm_shape(source) +
      tm_polygons(alpha = 1, col = "orange", border.col = "black", lwd = 3) +
    tm_layout(
        bg.color  = "black"
      , frame     = "white"
      , frame.lwd = 3
    )

  # In case tracks should be plotted, add them
  if (tracks){

    # To visualize the trajectories, we need to create them from the simulations
    sims_sub <- simulations[
        simulations$StepNumber <= steps
      & simulations$PointSampling == sampling
      & simulations$StartPoint == n
      , ]
    trajs <- sims2tracks(sims_sub, steps = steps, sampling = sampling)

    # However, we only want to keep the trajectories that are represented in the
    # table.
    trajs <- subset(trajs, ID %in% table$ID)

    # Add them to the plot
    p <- p + tm_shape(trajs) +
        tm_lines(col = darken("red"), alpha = alpha) +
      tm_shape(source) +
        tm_borders(alpha = 1, col = "white", lwd = 0.5)
  }

  # plot
  p
}

# Run a loop to create multiple maps
plots <- list()
for (i in 1:6){
  plots[[i]] <- visCon(
      n           = 52
    , simulations = sims
    , reached     = reached
    , areas       = areas
    , tracks      = T
    , steps       = reached$steps[i]
    , sampling    = reached$sampling[[i]]
    , alpha       = 0.5
  )
}

# We are only going to use three of the plots. Let's do some cosmetics on them
p1 <- plots[[2]] +
  tm_scale_bar(
        position  = c("right", "bottom")
      , text.size = 0.5
      , text.col  = "white"
      , width     = 0.125
  ) +
  tm_credits("(a)"
    , position  = c("left", "top")
    , size      = 1.5
    , col       = "white"
  ) +
  tm_compass(
      color.dark  = "white"
    , color.light = "white"
    , text.color  = "white"
    , position    = c("left", "bottom")
  ) +
  tm_layout(
      title = "Steps: 125"
    , title.position = c("center", "top")
    , title.color = "white"
  )
p2 <- plots[[4]] +
  tm_scale_bar(
        position  = c("right", "bottom")
      , text.size = 0.5
      , text.col  = "white"
      , width     = 0.125
  ) +
  tm_credits("(b)"
    , position  = c("left", "top")
    , size      = 1.5
    , col       = "white"
  ) +
  tm_compass(
      color.dark  = "white"
    , color.light = "white"
    , text.color  = "white"
    , position    = c("left", "bottom")
  ) +
  tm_layout(
      title = "Steps: 500"
    , title.position = c("center", "top")
    , title.color = "white"
  )
p3 <- plots[[6]] +
  tm_scale_bar(
        position  = c("right", "bottom")
      , text.size = 0.5
      , text.col  = "white"
      , width     = 0.125
  ) +
  tm_credits("(c)"
    , position  = c("left", "top")
    , size      = 1.5
    , col       = "white"
  ) +
  tm_compass(
      color.dark  = "white"
    , color.light = "white"
    , text.color  = "white"
    , position    = c("left", "bottom")
  ) +
  tm_layout(
      title = "Steps: 2000"
    , title.position = c("center", "top")
    , title.color = "white"
  )


# Arrange them
p <- tmap_arrange(p1, p2, p3, nrow = 1)

# Store the arranged plot
CairoPDF("05_Manuscript2/99_AreasReached.pdf", width = 11, height = 3.1, bg = "black")
p
dev.off()

################################################################################
#### Plot Areas Reached Over Time
################################################################################
# Unnest the tibble
dat <- reached %>% unnest(UniqueAreasReached) %>% dplyr::select(-AreasReached)

# Prepare list of comparisons
comps <- list(
    c("68", "125")
  , c("250", "500")
  , c("1000", "2000")
)

comps2 <- list(
    c("125", "250")
  , c("500", "1000")
)

# Visualize
dat %>% ggplot(aes(x = factor(steps), y = UniqueAreasReached, col = sampling)) +
  geom_boxplot() +
  geom_signif(comparisons = comps, map_signif_level = T, col = "white") +
  geom_signif(comparisons = comps2, map_signif_level = T, y_position = 35, col = "white") +
  geom_point(
      position  = position_jitterdodge(jitter.width = 0.5, jitter.height = 0.5)
    , alpha   = 0.5
  ) +
  dark_theme_minimal() +
  scale_color_manual(values = viridis(20)[c(15, 20)])

# Plot average unique areas reached over time
dat %>%
  group_by(steps, sampling) %>%
  summarize(AverageUniqueAreasReached = mean(UniqueAreasReached)) %>%
  ggplot(aes(x = steps, y = AverageUniqueAreasReached, col = sampling)) +
    geom_point() +
    geom_line() +
    dark_theme_minimal() +
    scale_color_manual(values = viridis(20)[c(15, 20)])

# Plot median unique areas reached over time
dat %>%
  group_by(steps, sampling) %>%
  summarize(MedianUniqueAreasReached = median(UniqueAreasReached)) %>%
  ggplot(aes(x = steps, y = MedianUniqueAreasReached, col = sampling)) +
    geom_point() +
    geom_line() +
    dark_theme_minimal() +
    scale_color_manual(values = viridis(20)[c(15, 20)])

# Plot average unique areas reached over time by start point
dat %>%
  group_by(steps, sampling) %>%
  ggplot(aes(x = steps, y = UniqueAreasReached, col = sampling)) +
    geom_point() +
    geom_line() +
    facet_wrap("StartPoint") +
    dark_theme_minimal() +
    scale_color_manual(values = viridis(20)[c(15, 20)])
