
library(scales)
library(terra)
library(tidyverse)
library(sf)
library(raster)
hum <- rast("/home/david/ownCloud/University/15. PhD/Chapter_8/03_Data/02_CleanData/HumanInfluence.tif")
# hum <- rast("/home/david/ownCloud/University/15. PhD/Chapter_1/03_Data/02_CleanData/04_AnthropogenicFeatures_HumanInfluenceBuff_FACEBOOK.grd")
# hum <- hum[[1]]
fre <- read_rds("/home/david/ownCloud/University/15. PhD/Chapter_8/03_Data/03_Results/HeatmapsBetweennessGlobal.rds")
fre <- stack(fre$Heatmap[5:6])
fre <- rast(fre)
# hum <- crop(hum, fre)
# hum <- aggregate(hum, fact = 10, fun = "mean")
# hum[hum <= 4] <- 0
plot(hum)
fre <- resample(fre, hum, method = "bilinear")
hum <- hum / minmax(hum)[2]
fre1 <- fre[[1]] / minmax(fre[[1]])[2]
fre2 <- fre[[2]] / minmax(fre[[2]])[2]
test1 <- hum * fre1
test2 <- hum * fre2
par(mfrow = c(2, 1))
png("AlternativeMin.png")
plot(test1)
dev.off()
png("AlternativeMax.png")
plot(test2)
dev.off()
plot(test2)
diff <- test2 - test1
diff <- as.data.frame(diff, xy = T)
ggplot(data = diff, aes(x = x, y = y, fill = HumanInfluence)) +
  geom_raster() +
  scale_fill_gradientn(
      colors  = c("orange", "white", "cornflowerblue")
    , labels  = function(x){format(x, big.mark = "'")}
    , limits  = c(-0.05, 0.05)
    , oob     = squish
    , guide   = guide_colorbar(
        show.limits    = T
      , ticks          = F
      , barheight      = unit(0.2, "cm")
      , barwidth       = unit(6, "cm")
    )
  ) +
  theme(legend.position = "none") +
  coord_equal()
