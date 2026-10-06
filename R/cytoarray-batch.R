# Folder-in / folder-out workflow: quantify a folder, keep the run in the
# environment, then re-do individual membranes by hand as needed.

#' Quantify every membrane in a folder
#'
#' Keeps the result in a variable (a `cytoarray_run`): per-membrane spots/affine/QC
#' plus the layout and paths, so single membranes can be re-done with [remanual()].
#' Overlays + intensity tables are written to `output`.
#'
#' @param input,output input image folder and output folder
#' @param template,product array layout (per-spot template CSV + built-in product)
#' @param refine,ref_mode passed to [analyze_membrane()] ("auto" / "manual" / ...)
#' @param pattern image filename pattern
#' @param include optional explicit image basenames to include
#' @return a `cytoarray_run` object (print it to see per-membrane QC)
#' @export
run_arrays <- function(input, output, template = NULL,
                       product = "ARY028", refine = "auto", ref_mode = "auto",
                       pattern = "\\.(tif|tiff)$", include = NULL) {
  layout <- load_layout(template, product)
  dir.create(output, showWarnings = FALSE, recursive = TRUE)
  fl <- list.files(input, pattern, full.names = TRUE, ignore.case = TRUE)
  fl <- fl[!grepl("_verify", basename(fl), ignore.case = TRUE)]
  if (!is.null(include)) {
    missing <- setdiff(include, basename(fl))
    if (length(missing)) stop("include samples not found: ", paste(missing, collapse = ", "))
    fl <- fl[basename(fl) %in% include]
  }
  stopifnot(length(fl) > 0)
  message(length(fl), " membranes in ", input)
  mem <- lapply(fl, .quant_one, layout = layout, output = output,
                ref_mode = ref_mode, refine = refine)
  names(mem) <- vapply(mem, `[[`, "", "id")
  incomplete <- names(mem)[vapply(mem, function(x) {
    nrow(x$results) != nrow(layout$spots) || anyNA(x$results$mean_int)
  }, logical(1))]
  if (length(incomplete))
    warning("Incomplete spot measurements in: ", paste(incomplete, collapse = ", "),
            "; inspect QC before downstream analysis.", call. = FALSE)
  arr <- structure(list(membranes = mem),
                   layout = layout, input = input, output = output,
                   class = "cytoarray_run")
  write_tables(arr); arr
}

#' Re-do chosen membranes by hand (default: manual reference + manual fine-tune)
#'
#' @param arr a `cytoarray_run` from [run_arrays()]
#' @param ids membrane id(s) to redo (the names shown when you print `arr`)
#' @param ref_mode,refine how to redo (default fully manual)
#' @return the updated `cytoarray_run` (tables are rewritten)
#' @export
remanual <- function(arr, ids, ref_mode = "manual", refine = "both") {
  lay <- attr(arr, "layout"); out <- attr(arr, "output")
  miss <- setdiff(ids, names(arr$membranes))
  if (length(miss)) stop("unknown membranes: ", paste(miss, collapse = ", "))
  for (id in ids)
    arr$membranes[[id]] <- .quant_one(arr$membranes[[id]]$path, lay, out, ref_mode, refine)
  write_tables(arr); arr
}

#' (Re)write intensity tables + QC from a run
#' @return the run, invisibly; writes intensity_spots_raw / intensity_analyte_normalized / QC_per_membrane
#' @export
write_tables <- function(arr, output = attr(arr, "output")) {
  dir.create(output, showWarnings = FALSE, recursive = TRUE)
  lay <- attr(arr, "layout")
  wide <- .wide_from_run(arr)
  raw <- wide
  pos2an <- stats::setNames(c(lay$analyte_map$analyte, lay$analyte_map$analyte),
                            c(lay$analyte_map$pos1, lay$analyte_map$pos2))
  raw$analyte <- pos2an[raw$position]
  raw <- raw[, c("spot_id", "position", "analyte", setdiff(names(raw), c("spot_id", "position", "analyte"))), drop = FALSE]
  utils::write.csv(raw, file.path(output, "intensity_spots_raw.csv"), row.names = FALSE)
  m <- normalize_array(wide, lay, total_norm = FALSE)$mat
  utils::write.csv(data.frame(analyte = rownames(m), m, check.names = FALSE),
                   file.path(output, "intensity_analyte_normalized.csv"), row.names = FALSE)
  utils::write.csv(qc_table(arr), file.path(output, "QC_per_membrane.csv"), row.names = FALSE)
  invisible(arr)
}

#' Per-membrane QC table from a run
#' @export
qc_table <- function(arr) do.call(rbind, lapply(arr$membranes, `[[`, "qc"))

#' Normalised analyte x membrane matrix from a run (for downstream stats)
#' @export
intensity_matrix <- function(arr) {
  lay <- attr(arr, "layout")
  wide <- .wide_from_run(arr)
  normalize_array(wide, lay, total_norm = FALSE)$mat
}

#' @export
print.cytoarray_run <- function(x, ...) {
  cat("cytoarray_run:", length(x$membranes), "membranes @", attr(x, "input"), "\n")
  print(qc_table(x)[, c("file", "reg_error_px", "n_measured", "n_snapped")], row.names = FALSE)
  cat("\n redo one by hand:  arr <- remanual(arr, \"<id>\")\n")
  invisible(x)
}

# internal: quantify one membrane -> a run entry
.quant_one <- function(path, layout, output, ref_mode, refine) {
  id <- sub("\\.(tif|tiff|png|jpg)$", "", basename(path), ignore.case = TRUE)
  m  <- analyze_membrane(path, layout, ref_mode = ref_mode, refine = refine, verify_png = FALSE)
  draw_overlay(NULL, m$results,
               file.path(output, paste0(id, "_verify.png")), original_path = path)
  list(path = path, id = id, results = m$results, M = m$M, qc = m$qc)
}

.wide_from_run <- function(arr) {
  mem <- arr$membranes
  data.frame(spot_id = mem[[1]]$results$spot_id,
             position = mem[[1]]$results$position,
             do.call(cbind, lapply(mem, function(e) e$results$mean_int)),
             check.names = FALSE)
}
