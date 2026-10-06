# Getting started

Install cytoarray using the [README](README.md), then follow these steps in RStudio.

## 1. Prepare two folders

Put one membrane in each `.tif` or `.tiff` file, for example `ctrl_1.tif`.
Choose a new folder for the results. Your original images are not changed.
The ARY028 layout is included, so you do not need to supply a template.

Open [01_quantify_arrays.R](inst/examples/01_quantify_arrays.R), change `input`
and `output`, then click **Source**:

```r
input <- "path/to/membrane_images"
output <- "path/to/array_results"
library(cytoarray)
arr <- run_arrays(input, output, ref_mode = "manual", refine = "both")
print(qc_table(arr))
saveRDS(arr, file.path(output, "quantification.rds"))
```

The script processes all TIFFs in `input`. For only one membrane, use:

```r
arr <- run_arrays(input, output, include = "ctrl_1.tif",
                  ref_mode = "manual", refine = "both")
```

`input` is still the folder. `include` is the complete file name in that folder.
To select several images, use `include = c("ctrl_1.tif", "treat_1.tif")`.

## 2. Click the six reference spots

For each membrane, click the centres in this order:

1. Top-left: upper, then lower.
2. Top-right: upper, then lower.
3. Bottom-right: upper, then lower.

<img src="man/figures/reference_points.png" width="350" alt="Reference spots numbered in click order"/>

## 3. Adjust and inspect the circles

With `refine = "both"`, automatic refinement runs first, then manual adjustment.
Click an off-centre circle, then its correct centre. End adjustment with Esc or
right-click, depending on your R graphics device, to continue to the next membrane.

<img src="man/figures/example_overlay.png" width="350" alt="Manual spot-adjustment view using saved final coordinates"/>

The view shows saved final spot coordinates. During manual adjustment, the image
background is processed for clearer viewing; saved verification images use the
original grayscale image. Green means adjusted, red a missing measurement, orange
a strong unadjusted signal, and cyan another measured position.

| Setting | What it does |
| --- | --- |
| `include = NULL` | Process all TIFFs |
| `include = "ctrl_1.tif"` | Process one named TIFF |
| `ref_mode = "manual"` | Click the six reference spots |
| `ref_mode = "auto"` | Find the reference spots automatically |
| `refine = "both"` | Automatic adjustment, then manual adjustment |
| `refine = "auto"` | Automatic adjustment only |

The function defaults are auto/auto. The starter script explicitly uses manual/both,
which needs an interactive R session. Use auto/auto for unattended scripts and
inspect the output afterwards. Other settings are described in `?run_arrays`.

## 4. Find and save your results

| File in the output folder | Contents |
| --- | --- |
| `intensity_spots_raw.csv` | Intensity of each printed spot |
| `intensity_analyte_normalized.csv` | Target-by-membrane intensity table |
| `QC_per_membrane.csv` | QC summary for each membrane |
| `<sample>_verify.png` | Spot positions drawn on the membrane |
| `quantification.rds` | Saved run object, created by the starter script |

Check every verification image and the QC table before downstream analysis.
To redo one membrane, use its file name without the extension:

```r
arr <- remanual(arr, ids = "ctrl_1")
saveRDS(arr, file.path(output, "quantification.rds"))
```

This updates its overlay and the output tables. Save again after any corrections.
To reopen the saved object later, use `arr <- readRDS("path/to/quantification.rds")`.

## 5. Compare groups

Open [02_complete_workflow.R](inst/examples/02_complete_workflow.R). At the top,
edit `input`, `output`, the named `group` vector and `control`. For example,
`ctrl_1 = "ctrl"` assigns `ctrl_1.tif` to the control group. The example has three
independent samples each for ctrl, treatA and treatB; use your own sample IDs.

The script quantifies the images, saves the object, then pauses for you to inspect
the overlays and QC. Press Enter to continue only after checking the measurements.
To compare groups later, load the saved object with the commented `readRDS()` line
and run only the group-comparison section; do not quantify the images again.

```r
M <- intensity_matrix(arr)
res <- diff_array(M, group, control = control)
```

Supply the intensity matrix directly, without taking log2 first. One fit uses all
groups to compare each treatment against the chosen control. The default ARY028
matrix has 111 rows. The script writes:

- `group_comparisons_all.csv`: all targets, with treatment `.logFC`, `.P` and `.FDR` columns.
- `<treatment>_FDR_lt_0.05.csv`: targets with FDR < 0.05 for that comparison.
- `sample_groups.csv`, `comparison_settings.txt` and `sessionInfo.txt`: settings to keep.

Positive logFC means higher signal in the treatment than in the control. Empty
filtered tables mean no targets passed that FDR cutoff. Use independent biological
samples, not duplicate spots or repeated scans as replicates; one pooled membrane
per group supports descriptive fold changes only. See `?diff_array` for options.

For **one pooled membrane per group**, use the optional
[pooled fold-change script](inst/examples/03_pooled_fold_changes.R) on that experiment's
saved object. Edit the object path, output folder and control/treatment sample IDs.
It exports descriptive log2 fold changes and fold changes, without P values or FDR.
Pathway enrichment uses fgsea separately; it is not included in cytoarray.

For help and license, see the [README](README.md#help).
