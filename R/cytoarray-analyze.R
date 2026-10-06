# Per-membrane driver.

#' Quantify one membrane
#' @param brush opening-disc diameter in pixels for image background flattening.
#' @param invert_if_mean_above invert scans brighter than this mean intensity.
#' @param snap_max_shift maximum centroid move in spot radii.
#' @param snap_min_snr minimum local signal-to-noise ratio required for a snap.
#' @param annulus_outer_radius outer annulus radius in spot radii.
#' @param registration_warning_px empirical registration-error warning threshold
#'   in pixels; it is not a validated acceptance criterion.
#' @param qc_saturation_fraction fraction of pixels at the processed-image maximum
#'   that flags a spot in `n_saturated`; this is a proxy, not device saturation.
#' @return list(results = spots table, M = affine, qc = one-row QC data.frame)
#' @export
analyze_membrane <- function(img_path, layout,
                             ref_mode = c("auto", "manual"),
                             refine = c("both", "auto", "manual", "none"),
                             uniform_radius = NULL, reuse_M = NULL,
                             verify_png = TRUE, brush = 51,
                             invert_if_mean_above = 0.5,
                             snap_max_shift = 0.4, snap_min_snr = 1.3,
                             annulus_outer_radius = 1.8,
                             registration_warning_px = 30,
                             qc_saturation_fraction = 0.05) {
  ref_mode <- match.arg(ref_mode); refine <- match.arg(refine)
  if (!interactive() && (ref_mode == "manual" || refine %in% c("manual", "both")))
    stop("Manual reference or refinement requires an interactive R session; choose ref_mode='auto' and refine='auto' or 'none' for scripted runs.")

  img <- preprocess_image(img_path, brush = brush, invert_if_mean_above = invert_if_mean_above)
  M <- if (!is.null(reuse_M)) reuse_M else if (ref_mode == "manual") {
    EBImage::display(img, method = "raster")
    message("Click the 6 corner reference spots IN ORDER: top-left upper, top-left lower, top-right upper, top-right lower, bottom-right upper, bottom-right lower")
    pts <- graphics::locator(6, type = "p", col = "red", pch = 3, cex = 2)
    if (is.null(pts) || length(pts$x) < 6) stop("Manual reference needs 6 clicks; got ", length(pts$x), ". Re-run.")
    compute_affine(layout$ref, data.frame(mx = pts$x, my = pts$y))
  } else compute_affine(layout$ref, detect_ref_spots(img))
  reg_err <- attr(M, "error"); if (is.null(reg_err)) reg_err <- NA

  spots <- map_spots(layout$spots, M, uniform_radius)
  if (refine %in% c("auto", "both")) spots <- snap_spots(img, spots, max_shift = snap_max_shift, min_snr = snap_min_snr)
  if (refine %in% c("manual", "both")) {
    spots <- measure_spots(img, spots, annulus_outer_radius = annulus_outer_radius)
    spots <- .manual_adjust(img, spots)
  }
  spots <- measure_spots(img, spots, annulus_outer_radius = annulus_outer_radius)
  spots$file <- basename(img_path)
  if (verify_png) draw_overlay(img, spots, sub("\\.(tif|tiff|png|jpg)$", "_verify.png", img_path, ignore.case = TRUE), original_path = img_path)

  qc <- data.frame(
    file = basename(img_path), reg_error_px = round(reg_err, 1),
    n_measured = sum(!is.na(spots$mean_int)),
    n_snapped = if (is.null(spots$snapped)) NA else sum(spots$snapped, na.rm = TRUE),
    n_saturated = sum(spots$sat_frac > qc_saturation_fraction, na.rm = TRUE),
    ref_cv = round(.cv(spots$mean_int[spots$position %in% layout$ref_pos]), 2),
    stringsAsFactors = FALSE
  )
  if (!is.na(reg_err) && reg_err > registration_warning_px)
    warning(sprintf("%s: high registration error (%.1f px) — check reference spots / use ref_mode='manual'.", basename(img_path), reg_err))
  list(results = spots, M = M, qc = qc)
}

.manual_adjust <- function(img, spots) {
  message("Manual adjust: click a circle, then its correct centre; ESC/right-click to stop.")
  r <- spots$mr[1]; theta <- seq(0, 2 * pi, length.out = 50)
  repeat {
    draw_overlay(img, spots)
    c1 <- graphics::locator(1); if (is.null(c1)) break
    k <- which.min((spots$mx - c1$x)^2 + (spots$my - c1$y)^2)
    graphics::lines(spots$mx[k] + r * cos(theta), spots$my[k] + r * sin(theta), col = "magenta", lwd = 2)
    c2 <- graphics::locator(1); if (is.null(c2)) break
    spots$mx[k] <- c2$x; spots$my[k] <- c2$y
    if (!is.null(spots$snapped)) spots$snapped[k] <- TRUE
  }
  spots
}

.cv <- function(x) { x <- x[is.finite(x)]; if (length(x) < 2 || mean(x) == 0) return(NA); stats::sd(x) / mean(x) }
