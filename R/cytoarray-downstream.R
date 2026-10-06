# Normalisation and differential analysis.

#' Normalise quantified spots to an analyte x sample matrix
#'
#' Duplicate spots are averaged. Shared dual-gene spots stay as one composite
#' measurement by default. All between-membrane scaling is off by default.
#' @param wide data.frame with `position` and one column per membrane
#' @param layout layout from [load_layout()]
#' @param neg_bg subtract negative-control background per membrane
#' @param ref_norm scale membranes by reference spots
#' @param median_norm scale membranes to a common analyte median
#' @param total_norm scale membranes to a common total signal
#' @param split_shared keep the legacy split into two identical gene-name rows;
#'   FALSE retains one row per printed analyte (111 for ARY028).
#' @return list with analyte matrix, sample names, and normalisation factors
#' @export
normalize_array <- function(wide, layout, neg_bg = FALSE, ref_norm = FALSE,
                            median_norm = FALSE, total_norm = FALSE,
                            split_shared = FALSE) {
  samples <- setdiff(names(wide), c("position", "spot_id", "tx", "ty", "tr"))
  signal <- as.matrix(wide[, samples, drop = FALSE])

  bg <- stats::setNames(rep(0, length(samples)), samples)
  if (neg_bg) {
    bg <- colMeans(signal[match(layout$neg_pos, wide$position), , drop = FALSE], na.rm = TRUE)
    bg[!is.finite(bg)] <- 0
  }
  signal <- sweep(signal, 2, bg, "-")
  signal[signal < 0] <- 0

  nf <- stats::setNames(rep(1, length(samples)), samples)
  if (ref_norm) {
    refm <- colMeans(signal[wide$position %in% layout$ref_pos, , drop = FALSE], na.rm = TRUE)
    refm[!is.finite(refm) | refm <= 0] <- NA
    nf <- mean(refm, na.rm = TRUE) / refm
    nf[!is.finite(nf)] <- 1
    signal <- sweep(signal, 2, nf, "*")
  }

  am <- layout$analyte_map[
    !layout$analyte_map$analyte %in% c("Reference", "Negative Control"), ]
  v1 <- signal[match(am$pos1, wide$position), , drop = FALSE]
  v2 <- signal[match(am$pos2, wide$position), , drop = FALSE]
  n_obs <- (!is.na(v1)) + (!is.na(v2))
  v1[is.na(v1)] <- 0
  v2[is.na(v2)] <- 0
  mat <- (v1 + v2) / n_obs
  rownames(mat) <- am$analyte

  if (split_shared) for (sp in layout$splits) {
    k <- which(rownames(mat) == sp[1])
    if (!length(k)) next
    # Same shared measurement under both printed names; not independent proteins.
    shared <- mat[k, , drop = FALSE]
    r1 <- r2 <- shared
    rownames(r1) <- rep(sp[2], nrow(shared))
    rownames(r2) <- rep(sp[3], nrow(shared))
    mat <- rbind(mat[-k, , drop = FALSE], r1, r2)
  }
  mat[!is.finite(mat)] <- NA

  mf <- NULL
  if (median_norm) {
    md <- apply(mat, 2, stats::median, na.rm = TRUE)
    mf <- mean(md) / md
    mat <- sweep(mat, 2, mf, "*")
  }
  tf <- NULL
  if (total_norm) {
    cs <- colSums(mat, na.rm = TRUE)
    tf <- mean(cs) / cs
    mat <- sweep(mat, 2, tf, "*")
  }
  list(mat = mat, samples = samples,
       factors = list(bg = bg, ref = nf, median = mf, total = tf))
}

#' Differential analysis across groups with limma
#'
#' @param mat analyte x sample matrix from [normalize_array()]
#' @param group group labels, ordered or named by `colnames(mat)`
#' @param control control group label
#' @param offset optional pseudo-count added before log2
#' @param offset_quantile quantile of positive intensities used when `offset` is
#'   NULL; default 0.02 preserves the existing analysis
#' @param trend estimate an intensity-dependent prior variance trend; TRUE by default.
#' @param robust use robust empirical Bayes moderation; TRUE by default.
#' @details Each matrix column must represent an independent biological sample.
#'   Each group requires at least two such samples. Duplicate spots or repeated
#'   scans of a membrane are technical repeats and must not be separate columns.
#'   The function cannot infer sample provenance from intensities or labels.
#'   A single pooled membrane per group supports descriptive fold changes only.
#' @return data.frame with analytes and per-contrast logFC, P, and FDR; its
#'   `offset` attribute records the pseudo-count used
#' @import statmod
#' @export
diff_array <- function(mat, group, control, offset = NULL, offset_quantile = 0.02,
                       trend = TRUE, robust = TRUE) {
  mat <- as.matrix(mat)
  if (is.null(colnames(mat))) stop("mat needs sample column names")
  if (!is.null(names(group))) {
    if (anyDuplicated(names(group)) || !setequal(names(group), colnames(mat)))
      stop("group names must match matrix samples")
    group <- group[colnames(mat)]
  } else if (length(group) != ncol(mat)) {
    stop("group length must equal ncol(mat)")
  }
  if (anyNA(group)) stop("group contains NA")
  group <- factor(group)
  if (nlevels(group) < 2L) stop("diff_array needs a control and at least one treatment group")
  if (any(table(group) < 2L))
    stop("Each group needs at least two independent biological samples. Duplicate spots and repeated scans are technical repeats; use descriptive fold changes for single or pooled membranes.")
  if (!control %in% levels(group)) stop("control is not present in group")
  group <- stats::relevel(group, ref = control)

  if (is.null(offset)) {
    positive <- mat[is.finite(mat) & mat > 0]
    if (!length(positive)) stop("no positive values available to estimate offset")
    if (length(offset_quantile) != 1L || offset_quantile < 0 || offset_quantile > 1)
      stop("offset_quantile must be between 0 and 1")
    offset <- as.numeric(stats::quantile(positive, offset_quantile, na.rm = TRUE))
  }
  if (length(offset) != 1L || !is.finite(offset) || offset <= 0)
    stop("offset must be positive and finite")

  design <- stats::model.matrix(~group)
  if (robust && !requireNamespace("statmod", quietly = TRUE))
    stop("robust = TRUE requires statmod; install it before running this comparison.")
  fit <- limma::eBayes(limma::lmFit(log2(mat + offset), design),
                      trend = trend, robust = robust)
  out <- data.frame(analyte = rownames(mat), stringsAsFactors = FALSE)
  for (co in setdiff(levels(group), control)) {
    tt <- limma::topTable(fit, coef = paste0("group", co), number = Inf, sort.by = "none")
    out[[paste0(co, ".logFC")]] <- tt$logFC
    out[[paste0(co, ".P")]] <- tt$P.Value
    out[[paste0(co, ".FDR")]] <- tt$adj.P.Val
  }
  out <- out[order(-abs(out[[grep("logFC", names(out))[1]]])), ]
  attr(out, "offset") <- offset
  out
}

