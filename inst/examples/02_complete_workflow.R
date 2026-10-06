# Open in RStudio. Edit the settings below, then run the script interactively.
# Sample IDs are image file names without .tif/.tiff, e.g. ctrl_1.tif -> ctrl_1.
input <- "path/to/membrane_images"
output <- "path/to/array_results"
group <- c(ctrl_1 = "ctrl", ctrl_2 = "ctrl", ctrl_3 = "ctrl",
           treatA_1 = "treatA", treatA_2 = "treatA", treatA_3 = "treatA",
           treatB_1 = "treatB", treatB_2 = "treatB", treatB_3 = "treatB")
control <- "ctrl"

library(cytoarray)
if (packageVersion("cytoarray") < "0.0.4") stop("Install cytoarray 0.0.4 or later before using this demo.")
# Quantify all TIFFs. Click TL upper/lower, TR upper/lower, BR upper/lower.
# Adjust an off-centre spot by clicking its circle, then its centre.
arr <- run_arrays(input, output, ref_mode = "manual", refine = "both")
print(qc_table(arr))
saveRDS(arr, file.path(output, "quantification.rds"))

# Inspect every verify.png and the QC table before continuing.
# To redo one: arr <- remanual(arr, ids = "ctrl_1"), then saveRDS() again.
# You can instead run the group-comparison section later after loading:
# arr <- readRDS(file.path(output, "quantification.rds"))
readline("Review the overlays and QC. Press Enter to compare groups, or Esc to stop: ")

# GROUP COMPARISONS: run this section only after accepting the measurements.
M <- intensity_matrix(arr)
# Supply the original intensities: diff_array() applies its log2 transform itself.
# Fit all groups once; compare treatA and treatB with the same control.
res <- diff_array(M, group, control = control, offset_quantile = 0.02,
                  trend = TRUE, robust = TRUE)
write.csv(res, file.path(output, "group_comparisons_all.csv"), row.names = FALSE)

# Save a separate FDR < 0.05 table for each treatment versus control.
for (treatment in setdiff(unique(group), control)) {
  columns <- c("analyte", paste0(treatment, c(".logFC", ".P", ".FDR")))
  fdr <- res[[paste0(treatment, ".FDR")]]
  selected <- res[!is.na(fdr) & fdr < 0.05, columns, drop = FALSE]
  write.csv(selected, file.path(output, paste0(treatment, "_FDR_lt_0.05.csv")),
            row.names = FALSE)
}
# Positive logFC means higher signal in the treatment than in the control.
writeLines(c(paste("control:", control), "offset_quantile: 0.02",
             paste("offset:", format(attr(res, "offset"), digits = 17)),
             "trend: TRUE", "robust: TRUE", "FDR: BH per contrast"),
           file.path(output, "comparison_settings.txt"))
write.csv(data.frame(sample = names(group), group = unname(group)),
          file.path(output, "sample_groups.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output, "sessionInfo.txt"))
