`%||%` <- function(x, y) if (is.null(x) || length(x) == 0) y else x

wild_cols <- c(
  ink = "#25332f", sand = "#cf9b62", turquoise = "#60bcc4",
  olive = "#9fa669", red = "#9f4a42", cream = "#f5f1e8",
  water = "#a9d7d6", pale = "#e8e5dc", grey = "#73807b"
)

plot_setup <- function() {
  par(mar = c(0, 0, 0, 0), bg = "transparent", xaxs = "i", yaxs = "i")
}

draw_region <- function() {
  plot.new()
  plot.window(xlim = c(0, 100), ylim = c(0, 100), asp = 1)
  rect(0, 0, 100, 100, col = "#edf0e8", border = NA)
  polygon(c(5, 22, 29, 47, 64, 83, 96, 91, 72, 52, 38, 17),
          c(25, 12, 22, 10, 20, 13, 36, 73, 89, 82, 93, 76),
          col = "#e1e5d8", border = "#bdc6b7", lwd = 1.2)
  polygon(c(18, 25, 34, 37, 31, 22, 15), c(55, 73, 79, 66, 52, 43, 47),
          col = wild_cols$water, border = NA)
  polygon(c(26, 32, 42, 48, 43, 34), c(47, 59, 57, 48, 39, 38),
          col = "#c5dfd5", border = NA)
  for (x in seq(12, 92, by = 13)) abline(v = x, col = "#d7dcd2", lwd = .55)
  for (y in seq(15, 90, by = 15)) abline(h = y, col = "#d7dcd2", lwd = .55)
  text(c(22, 74, 81), c(86, 85, 25), c("OKAVANGO", "CHOBE", "MAKGADIKGADI"), col = "#89938d", cex = .65, font = 2)
}

draw_landscape <- function(mode = "tracks") {
  plot_setup()
  set.seed(54)
  draw_region()
  patches <- data.frame(x = c(24, 37, 69, 82, 60, 44), y = c(61, 45, 73, 31, 24, 79), r = c(5, 4, 5, 6, 4, 3))
  symbols(patches$x, patches$y, circles = patches$r, inches = FALSE,
          add = TRUE, bg = adjustcolor(wild_cols$olive, .33), fg = NA)

  if (mode == "observations") {
    for (i in seq_len(13)) {
      n <- sample(15:38, 1)
      x <- cumsum(rnorm(n, .9, 2.8)) + sample(10:48, 1)
      y <- cumsum(rnorm(n, 0, 2.4)) + sample(30:80, 1)
      lines(x, y, col = adjustcolor(wild_cols$sand, .35), lwd = 1)
      points(x, y, pch = 16, cex = .22, col = adjustcolor(wild_cols$ink, .62))
    }
  } else if (mode == "steps") {
    start <- c(43, 48)
    ends <- cbind(start[1] + sin(seq(-pi, pi, length.out = 26)) * runif(26, 10, 30),
                  start[2] + cos(seq(-pi, pi, length.out = 26)) * runif(26, 10, 30))
    segments(start[1], start[2], ends[, 1], ends[, 2], col = adjustcolor(wild_cols$turquoise, .32), lwd = .8)
    segments(start[1], start[2], 57, 63, col = wild_cols$red, lwd = 3)
    points(start[1], start[2], pch = 21, bg = wild_cols$ink, col = "white", cex = 1.2)
    points(ends[, 1], ends[, 2], pch = 16, col = adjustcolor(wild_cols$turquoise, .7), cex = .55)
    points(57, 63, pch = 21, bg = wild_cols$red, col = "white", cex = 1.15)
    legend("bottomleft", inset = .03, legend = c("Observed", "Available"), col = c(wild_cols$red, wild_cols$turquoise), lwd = c(3, 1), bty = "n", cex = .72)
  } else {
    n_paths <- if (mode == "simulation") 42 else 16
    starts <- matrix(c(24, 61, 37, 45), ncol = 2, byrow = TRUE)
    for (i in seq_len(n_paths)) {
      s <- starts[sample(1:2, 1), ]
      target <- patches[sample(3:nrow(patches), 1), 1:2]
      n <- 35
      t <- seq(0, 1, length.out = n)
      x <- s[1] * (1 - t) + target$x * t + cumsum(rnorm(n, 0, .8))
      y <- s[2] * (1 - t) + target$y * t + cumsum(rnorm(n, 0, .75))
      lines(x, y, col = adjustcolor(if (mode == "simulation") wild_cols$sand else wild_cols$ink, if (mode == "simulation") .22 else .27), lwd = if (mode == "simulation") 1 else 1.3)
    }
    points(starts[, 1], starts[, 2], pch = 21, bg = wild_cols$red, col = "white", cex = 1.4, lwd = 1.2)
    points(patches[3:6, 1], patches[3:6, 2], pch = 21, bg = wild_cols$olive, col = "white", cex = 1.15)
  }
  box(col = "#d5dbd2")
}

draw_sampling_plot <- function() {
  par(mar = c(3.6, 4, .5, .5), bg = "transparent")
  set.seed(11)
  x <- pmax(.2, c(rnorm(1650, 4, .32), rnorm(260, 8, .65), rexp(180, .25)))
  hist(x, breaks = seq(0, 18, .5), col = wild_cols$sand, border = "white",
       xlab = "Hours between successive fixes", ylab = "Number of fixes", main = "", axes = FALSE)
  axis(1, at = seq(0, 18, 2), col = NA, col.axis = wild_cols$grey, cex.axis = .8)
  axis(2, at = pretty(c(0, 500)), col = NA, col.axis = wild_cols$grey, las = 1, cex.axis = .8)
  abline(v = 4, col = wild_cols$red, lwd = 2)
}

draw_gps_overview <- function() {
  old <- par(no.readonly = TRUE)
  on.exit(par(old), add = TRUE)
  layout(matrix(1:3, nrow = 1), widths = c(1.15, 1, 1))
  set.seed(31)

  # Number of fixes per individual, matching the chapter's gpsOverview design.
  fixes <- sort(round(rlnorm(28, log(5200), .36)))
  par(mar = c(4.2, 1.2, 2.2, .7), bg = "transparent")
  barplot(fixes, horiz = TRUE, axes = FALSE, border = NA,
          col = colorRampPalette(c(wild_cols$sand, wild_cols$turquoise))(28),
          xlab = "Number of fixes")
  axis(1, at = pretty(c(0, max(fixes))), labels = format(pretty(c(0, max(fixes))), big.mark = ","),
       col = NA, col.axis = wild_cols$grey, cex.axis = .72)
  title("Fixes by individual", adj = 0, cex.main = .9, col.main = wild_cols$ink)

  # Fix collection times by individual and hour.
  par(mar = c(4.2, 1.2, 2.2, .7), bg = "transparent")
  plot(0, 0, type = "n", xlim = c(0, 23), ylim = c(.5, 28.5), axes = FALSE,
       xlab = "Hour of day", ylab = "")
  for (id in 1:28) {
    hrs <- c(3, 7, 15, 19, 23) + rnorm(5, 0, .13)
    points(hrs, rep(id, 5), pch = 16, cex = .45,
           col = adjustcolor(if (id %% 2) wild_cols$turquoise else wild_cols$olive, .7))
  }
  axis(1, at = c(0, 4, 8, 12, 16, 20, 24), col = NA, col.axis = wild_cols$grey, cex.axis = .72)
  title("Fix timing", adj = 0, cex.main = .9, col.main = wild_cols$ink)

  # Successive-fix time lags by individual.
  par(mar = c(4.2, 1.2, 2.2, .7), bg = "transparent")
  lag_data <- lapply(1:28, function(i) pmax(.4, rnorm(55, 4 + runif(1, -.3, .4), .42)))
  boxplot(lag_data, horizontal = TRUE, axes = FALSE, outline = FALSE,
          col = adjustcolor(wild_cols$sand, .65), border = wild_cols$sand,
          xlab = "Successive-fix lag (hours)")
  axis(1, at = seq(2, 7, 1), col = NA, col.axis = wild_cols$grey, cex.axis = .72)
  abline(v = 4, col = wild_cols$red, lwd = 1.3)
  title("Time lag", adj = 0, cex.main = .9, col.main = wild_cols$ink)
}

draw_step_distribution <- function(n_random = 25) {
  par(mar = c(3.6, 4, .5, .5), bg = "transparent")
  set.seed(100 + n_random)
  x <- rgamma(2400, shape = 1.55, scale = 2.3)
  hist(x, breaks = seq(0, 18, .6), probability = TRUE, col = adjustcolor(wild_cols$turquoise, .6), border = "white", xlab = "Step length (km)", ylab = "Density", main = "", axes = FALSE)
  curve(dgamma(x, 1.55, scale = 2.3), from = 0, to = 18, add = TRUE, col = wild_cols$red, lwd = 2.4)
  axis(1, at = seq(0, 18, 3), col = NA, col.axis = wild_cols$grey, cex.axis = .8)
  axis(2, col = NA, col.axis = wild_cols$grey, las = 1, cex.axis = .8)
}

draw_covariate <- function(name) {
  plot_setup()
  set.seed(match(name, c("Seasonal water", "Distance to water", "Tree cover", "Human influence")) %||% 1)
  z <- outer(seq(-2, 2, length.out = 85), seq(-2, 2, length.out = 85), function(x, y) sin(x * 1.4) + cos(y * 1.8) + exp(-(x^2 + y^2)))
  z <- z + matrix(rnorm(length(z), 0, .15), nrow(z))
  pal <- colorRampPalette(c("#f3efe5", "#d6ddbf", "#9fb89a", "#5f8d86", "#335e62"))(80)
  image(seq(0, 100, length.out = nrow(z)), seq(0, 100, length.out = ncol(z)), z, col = pal, axes = FALSE, useRaster = TRUE)
  contour(seq(0, 100, length.out = nrow(z)), seq(0, 100, length.out = ncol(z)), z, add = TRUE, drawlabels = FALSE, col = adjustcolor("white", .25), lwd = .6)
  polygon(c(5, 22, 29, 47, 64, 83, 96, 91, 72, 52, 38, 17), c(25, 12, 22, 10, 20, 13, 36, 73, 89, 82, 93, 76), border = "white", col = NA, lwd = 1.5)
  text(6, 6, name, adj = 0, col = "white", font = 2, cex = .9)
  box(col = "#d5dbd2")
}

draw_coefficients <- function(terms = NULL) {
  par(mar = c(3.2, 11.5, .7, 1), bg = "transparent")
  defaults <- c("DistanceToWater_Along", "Water_End", "TreeCover_End", "Humans_Along")
  terms <- head(terms %||% defaults, 5)
  labs <- c(terms, "Turning persistence", "Log step length", "Step length")
  labs <- gsub("_", " · ", labs)
  set.seed(44 + length(labs))
  est <- c(seq(-.58, .62, length.out = length(terms)), .43, -.24, .11)
  est <- est + rnorm(length(est), 0, .06)
  se <- runif(length(est), .05, .13)
  y <- seq_along(labs)
  plot(c(-1, 1.1), c(.5, length(labs) + .5), type = "n", axes = FALSE, xlab = "Selection coefficient (β)", ylab = "")
  abline(v = 0, col = "#bcc4bf", lwd = 1)
  segments(est - 1.96 * se, y, est + 1.96 * se, y, col = wild_cols$grey, lwd = 1.4)
  points(est, y, pch = 21, bg = ifelse(est > 0, wild_cols$turquoise, wild_cols$red), col = "white", cex = 1.25)
  axis(1, col = NA, col.axis = wild_cols$grey, cex.axis = .8)
  axis(2, at = y, labels = labs, las = 1, col = NA, col.axis = wild_cols$ink, cex.axis = .78)
}

draw_validation <- function() {
  par(mar = c(3.6, 4, .6, .6), bg = "transparent")
  ranks <- 1:10
  set.seed(2)
  vals <- round(c(3, 5, 7, 9, 13, 17, 23, 31, 42, 58) + rnorm(10, 0, 1.5))
  plot(ranks, vals, type = "b", pch = 21, bg = wild_cols$olive, col = wild_cols$olive, lwd = 2, axes = FALSE, xlab = "Used-step rank", ylab = "Frequency (%)", ylim = c(0, 65))
  polygon(c(ranks, rev(ranks)), c(pmax(vals - 4, 0), rev(vals + 4)), col = adjustcolor(wild_cols$olive, .16), border = NA)
  lines(ranks, vals, col = wild_cols$olive, lwd = 2)
  points(ranks, vals, pch = 21, bg = wild_cols$olive, col = "white", cex = .9)
  axis(1, at = ranks, col = NA, col.axis = wild_cols$grey, cex.axis = .8)
  axis(2, col = NA, col.axis = wild_cols$grey, las = 1, cex.axis = .8)
  text(2.2, 57, "Spearman ρ = 0.94", col = wild_cols$ink, font = 2, cex = .85)
}

draw_connectivity <- function(type) {
  plot_setup()
  set.seed(81)
  n <- 105
  x <- seq(-2.8, 2.8, length.out = n)
  y <- seq(-2.8, 2.8, length.out = n)
  if (type == "Betweenness") {
    z <- outer(x, y, function(a, b) exp(-((b - .45 * sin(a * 1.6))^2) / .13) * exp(-a^2 / 4))
    pal <- colorRampPalette(c("#f0eee7", "#dbcf9c", "#cf9b62", "#9f4a42", "#5a2b2b"))(100)
  } else if (type == "Interpatch connectivity") {
    z <- outer(x, y, function(a, b) exp(-((a + 1)^2 + (b - 1)^2) / .4) + exp(-((a - 1.2)^2 + (b + .8)^2) / .5) + .45 * exp(-(b - .3 * a)^2 / .25))
    pal <- colorRampPalette(c("#f0eee7", "#d7dfca", "#9fa669", "#58716b", "#25332f"))(100)
  } else {
    z <- outer(x, y, function(a, b) exp(-((a + .8)^2 + (b - .7)^2) / .8) + .8 * exp(-((a - 1.1)^2 + (b + .9)^2) / .65) + .45 * exp(-(b - .2 * a)^2 / .8))
    pal <- colorRampPalette(c("#f0eee7", "#d8dfc8", "#9fc2ad", "#60bcc4", "#335e62"))(100)
  }
  image(seq(0, 100, length.out = n), seq(0, 100, length.out = n), z, col = pal, axes = FALSE, useRaster = TRUE)
  contour(seq(0, 100, length.out = n), seq(0, 100, length.out = n), z, levels = pretty(range(z), 7), add = TRUE, drawlabels = FALSE, col = adjustcolor("white", .25), lwd = .6)
  points(c(22, 35), c(64, 46), pch = 21, bg = wild_cols$red, col = "white", cex = 1.25)
  points(c(70, 83, 59), c(73, 31, 24), pch = 21, bg = wild_cols$olive, col = "white", cex = 1)
  box(col = "#d5dbd2")
}
