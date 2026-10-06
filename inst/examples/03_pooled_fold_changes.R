# Descriptive changes for one pooled membrane per group: no P values or FDR.
# Use a saved run from this pooled experiment, separate from replicated experiments.
object_file <- "path/to/pooled_results/quantification.rds"
output <- "path/to/pooled_results"
control_sample <- "ctrl_pool"
treatment_samples <- c("treatA_pool", "treatB_pool")

library(cytoarray)
if (packageVersion("cytoarray") < "0.0.4") stop("Install cytoarray 0.0.4 or later before using this demo.")
arr <- readRDS(object_file)
M <- intensity_matrix(arr)

# POOLED FOLD CHANGES: keep this experiment's control and treatment membranes.
M <- M[, c(control_sample, treatment_samples), drop = FALSE]
positive <- M[is.finite(M) & M > 0]
if (!length(positive)) stop("No positive intensities available to estimate the offset.")
offset <- as.numeric(quantile(positive, 0.02, na.rm = TRUE))
result <- data.frame(analyte = rownames(M))
for (sample in treatment_samples) {
  logFC <- log2((M[, sample] + offset) / (M[, control_sample] + offset))
  result[[paste0(sample, ".log2FC")]] <- logFC
  result[[paste0(sample, ".fold_change")]] <- 2^logFC
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)
write.csv(result, file.path(output, "pooled_descriptive_fold_changes.csv"), row.names = FALSE)
writeLines(c(paste("control_sample:", control_sample), "offset_quantile: 0.02",
             paste("offset:", format(offset, digits = 17)), "Descriptive only: no P values or FDR"),
           file.path(output, "pooled_settings.txt"))
