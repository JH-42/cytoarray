# Open in RStudio, edit these two folders, then click Source.
input <- "path/to/membrane_images"
output <- "path/to/array_results"

library(cytoarray)
# Click six references: top-left upper/lower, top-right upper/lower,
# then bottom-right upper/lower.
# To move a spot, click its circle, then its centre. Esc/right-click ends adjustment.
arr <- run_arrays(input, output, ref_mode = "manual", refine = "both")
# For one image, add include = "ctrl_1.tif" to run_arrays() above.
print(qc_table(arr))
# Check the verification images, then save the run.
saveRDS(arr, file.path(output, "quantification.rds"))
# Redo one: arr <- remanual(arr, ids = "ctrl_1"), then save again.
