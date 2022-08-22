################################################################################
#### Helper Functions
################################################################################
# Function to determine the pdf of a mixed von mises distribution
dvonmises <- function(x, kappa, mu) {
  exp(kappa * cos(x - mu)) / (2 * pi * besselI(kappa, nu = 0))
}

# Function to randomly sample from a mixed von mises distribution
rvonmises <- function(n, kappa, mu, by = 0.01) {
  x <- seq(-pi, +pi, by = by)
  probs <- dvonmises(x, kappa = kappa, mu = mu)
  random <- sample(x, size = n, prob = probs, replace = T)
  return(random)
}

# Function that extends a covariate layer and fills the added border with values
# sampled from the layer
extendRaster <- function(x, y){
  extended <- lapply(1:nlayers(x), function(z) {
    layer <- x[[z]]
    na_mask <- is.na(layer)
    vals <- values(layer)
    vals <- na.omit(vals)
    r <- extend(layer, y)
    na_mask <- extend(na_mask, y, value = 0)
    indices <- which(is.na(values(r)))
    values(r)[indices] <- vals[runif(length(indices), min = 1, max = length(vals))]
    r <- mask(r, na_mask, maskvalue = 1, updatevalue = NA)
    return(r)
  })
  return(stack(extended))
}

# Function to interpolate coordinates between two points
interpolatePoints <- function(x1, x2, y1, y2, by = 1) {
  length <- sqrt((x1 - x2) ** 2 + (y1 - y2) ** 2)
  nsegs <- max(ceiling(length / by), 1)
  x <- seq(x1, x2, length.out = nsegs + 1)
  y <- seq(y1, y2, length.out = nsegs + 1)
  return(cbind(x, y))
}
