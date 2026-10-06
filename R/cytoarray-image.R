# Image registration and spot measurement.

.as_gray <- function(img) {
  if (EBImage::colorMode(img) == EBImage::Color) EBImage::channel(img, "gray") else img
}

.image_matrix <- function(img) {
  mat <- EBImage::imageData(img)
  if (length(dim(mat)) == 3) mat <- mat[, , 1]
  mat
}

#' Read and background-flatten a membrane image
#' @param brush opening-disc diameter in pixels.
#' @param invert_if_mean_above invert scans brighter than this mean intensity.
#' @export
preprocess_image <- function(img_path, brush = 51, invert_if_mean_above = 0.5) {
  img <- .as_gray(EBImage::readImage(img_path))
  if (mean(img) > invert_if_mean_above) img <- 1 - img
  bg <- EBImage::opening(img, EBImage::makeBrush(brush, shape = "disc"))
  pmax(img - bg, 0)
}

#' Detect the six corner reference spots.
#' `min_area_fraction` is relative to the largest detected blob. Corner fractions
#' are fractions of the detected array width or height, measured from each edge.
#' @export
detect_ref_spots <- function(img,
                             threshold_percentiles = c(99.9, 99.5, 99.0, 98.0),
                             threshold_scale = 0.7,
                             min_area_fraction = 0.25,
                             array_edge_quantiles = c(0.02, 0.98),
                             corner_x_fraction = 0.30,
                             corner_y_fraction = 0.18,
                             pair_vertical_weight = 0.2) {
  mat <- .image_matrix(img)
  take2 <- function(df, sx, sy) {
    if (nrow(df) < 2) return(df[0, ])
    df <- df[df$area >= min_area_fraction * max(df$area), , drop = FALSE]
    if (nrow(df) < 2) return(df[0, ])
    anchor_index <- which.max(sx * df$cx + sy * df$cy)
    anchor <- df[anchor_index, ]
    rest <- df[-anchor_index, , drop = FALSE]
    mate <- which.min(abs(rest$cx - anchor$cx) + pair_vertical_weight * abs(rest$cy - anchor$cy))
    out <- rbind(anchor, rest[mate, ]); out[order(out$cy), ]
  }
  tl <- tr <- br <- data.frame()
  for (pct in threshold_percentiles) {
    threshold <- stats::quantile(as.vector(mat), pct / 100)
    lab <- EBImage::bwlabel(EBImage::Image(mat > threshold * threshold_scale, dim = dim(mat)))
    feat <- EBImage::computeFeatures.moment(lab, mat)
    shp <- EBImage::computeFeatures.shape(lab, mat)
    if (is.null(feat)) next
    blobs <- data.frame(cx = feat[, "m.cx"], cy = feat[, "m.cy"], area = shp[, "s.area"])
    if (nrow(blobs) < 6) next
    x0 <- stats::quantile(blobs$cx, array_edge_quantiles[1])
    x1 <- stats::quantile(blobs$cx, array_edge_quantiles[2])
    y0 <- stats::quantile(blobs$cy, array_edge_quantiles[1])
    y1 <- stats::quantile(blobs$cy, array_edge_quantiles[2])
    bw <- x1 - x0; bh <- y1 - y0
    tl <- take2(blobs[blobs$cx < x0 + corner_x_fraction * bw & blobs$cy < y0 + corner_y_fraction * bh, ], -1, -1)
    tr <- take2(blobs[blobs$cx > x1 - corner_x_fraction * bw & blobs$cy < y0 + corner_y_fraction * bh, ], 1, -1)
    br <- take2(blobs[blobs$cx > x1 - corner_x_fraction * bw & blobs$cy > y1 - corner_y_fraction * bh, ], 1, 1)
    if (nrow(tl) >= 2 && nrow(tr) >= 2 && nrow(br) >= 2) break
  }
  if (nrow(tl) < 2 || nrow(tr) < 2 || nrow(br) < 2)
    stop("Could not auto-detect reference spots; use ref_mode='manual'.")
  data.frame(mx = c(tl$cx, tr$cx, br$cx), my = c(tl$cy, tr$cy, br$cy),
             zone = c("TL1", "TL2", "TR1", "TR2", "BR1", "BR2"), stringsAsFactors = FALSE)
}

#' Fit the template-to-membrane affine map.
#' @export
compute_affine <- function(ref_template, ref_membrane) {
  G <- cbind(ref_template$tx, ref_template$ty, 1)
  P <- cbind(ref_membrane$mx, ref_membrane$my)
  M <- qr.solve(G, P)
  attr(M, "error") <- mean(sqrt(rowSums((G %*% M - P)^2)))
  M
}

#' Map template spots onto the membrane.
#' @export
map_spots <- function(template_spots, M, uniform_radius = NULL, min_radius_px = 4) {
  mapped <- cbind(template_spots$tx, template_spots$ty, 1) %*% M
  s <- template_spots; s$mx <- mapped[, 1]; s$my <- mapped[, 2]
  scale <- sqrt(M[1, 1]^2 + M[2, 1]^2)
  if (is.null(uniform_radius)) {
    typical_radius <- as.numeric(names(sort(table(template_spots$tr), decreasing = TRUE))[1])
    uniform_radius <- max(min_radius_px, round(typical_radius * scale))
  }
  s$mr <- uniform_radius; s
}

#' Refine spot positions by bounded local-intensity centroid.
#' @export
snap_spots <- function(img, spots, win = 1.4, max_shift = 0.8,
                       min_snr = 1.3, min_window_px = 3, mad_epsilon = 1e-9) {
  mat <- .image_matrix(img)
  w <- dim(mat)[1]; h <- dim(mat)[2]
  r <- spots$mr[1]; wr <- round(win * r)
  spots$snapped <- FALSE
  for (i in seq_len(nrow(spots))) {
    cx <- round(spots$mx[i]); cy <- round(spots$my[i])
    if (is.na(cx) || is.na(cy) || cx < 1 || cx > w || cy < 1 || cy > h) next
    xs <- max(1, cx - wr):min(w, cx + wr)
    ys <- max(1, cy - wr):min(h, cy + wr)
    if (length(xs) < min_window_px || length(ys) < min_window_px) next
    sub <- mat[xs, ys]
    bg <- stats::median(sub)
    sd_ <- stats::mad(sub) + mad_epsilon
    if ((max(sub) - bg) / sd_ < min_snr) next
    wgt <- pmax(sub - (bg + min_snr * sd_), 0)
    if (sum(wgt) <= 0) next
    gx <- sum(outer(xs, rep(1, length(ys))) * wgt) / sum(wgt)
    gy <- sum(outer(rep(1, length(xs)), ys) * wgt) / sum(wgt)
    if (sqrt((gx - spots$mx[i])^2 + (gy - spots$my[i])^2) <= max_shift * r) {
      spots$mx[i] <- gx
      spots$my[i] <- gy
      spots$snapped[i] <- TRUE
    }
  }
  spots
}

#' Measure circular apertures with optional local-annulus background.
#' `sat_frac` is a processed-image maximum-intensity proxy, not device saturation.
#' @export
measure_spots <- function(img, spots, bg_annulus = TRUE,
                          annulus_outer_radius = 1.8, annulus_gap_px = 1,
                          min_aperture_fraction = 0.5, saturation_fraction = 0.999) {
  mat <- .image_matrix(img)
  w <- dim(mat)[1]; h <- dim(mat)[2]; r <- spots$mr[1]
  inner <- expand.grid(dx = -r:r, dy = -r:r)
  inner <- inner[inner$dx^2 + inner$dy^2 <= r^2, ]
  outer_radius <- round(r * annulus_outer_radius)
  ring <- expand.grid(dx = -outer_radius:outer_radius, dy = -outer_radius:outer_radius)
  ring <- ring[ring$dx^2 + ring$dy^2 > (r + annulus_gap_px)^2 & ring$dx^2 + ring$dy^2 <= outer_radius^2, ]
  image_max <- max(mat)
  vals <- mapply(function(mx, my) {
    cx <- round(mx); cy <- round(my)
    if (is.na(cx) || cx < 1 || cx > w || cy < 1 || cy > h)
      return(c(mean = NA, total = NA, peak = NA, bg = NA, sat = NA))
    ix <- cx + inner$dx; iy <- cy + inner$dy
    ok <- ix >= 1 & ix <= w & iy >= 1 & iy <= h
    if (sum(ok) < min_aperture_fraction * length(ok))
      return(c(mean = NA, total = NA, peak = NA, bg = NA, sat = NA))
    pv <- mat[cbind(ix[ok], iy[ok])]
    bgv <- 0
    if (bg_annulus) {
      rx <- cx + ring$dx; ry <- cy + ring$dy
      rok <- rx >= 1 & rx <= w & ry >= 1 & ry <= h
      if (any(rok)) bgv <- stats::median(mat[cbind(rx[rok], ry[rok])])
    }
    c(mean = mean(pv) - bgv,
      total = sum(pmax(pv - bgv, 0)),
      peak = max(pv), bg = bgv,
      sat = mean(pv >= image_max * saturation_fraction))
  }, spots$mx, spots$my)
  spots$mean_int <- pmax(vals["mean", ], 0)
  spots$integ_int <- vals["total", ]
  spots$max_int <- vals["peak", ]
  spots$local_bg <- vals["bg", ]
  spots$sat_frac <- vals["sat", ]
  spots
}

#' Draw a verification overlay.
#' @export
draw_overlay <- function(img, spots, save_path = NULL, original_path = NULL,
                         high_intensity_quantile = 0.75) {
  disp <- if (!is.null(original_path)) {
    .as_gray(EBImage::readImage(original_path))
  } else 1 - img
  if (!is.null(save_path)) {
    grDevices::png(save_path, width = dim(disp)[1] * 2,
                   height = dim(disp)[2] * 2, res = 150)
  }
  EBImage::display(disp, method = "raster")
  theta <- seq(0, 2 * pi, length.out = 50)
  r <- spots$mr[1]
  vint <- if (is.null(spots$mean_int)) rep(NA_real_, nrow(spots)) else spots$mean_int
  snp <- if (is.null(spots$snapped)) rep(NA, nrow(spots)) else spots$snapped
  high_intensity <- stats::quantile(vint, high_intensity_quantile, na.rm = TRUE)
  invisible(mapply(function(mx, my, v, s) {
    col <- if (is.na(v)) "red" else if (isTRUE(s)) "green" else if (!is.na(high_intensity) && v > high_intensity) "orange" else "cyan"
    graphics::lines(mx + r * cos(theta), my + r * sin(theta), col = col, lwd = 0.8)
  }, spots$mx, spots$my, vint, snp))
  if (!is.null(save_path)) grDevices::dev.off()
}


