# cytoarray 0.0.5

- Keep a text-only colour legend and the complete workflow examples.
- Update release documentation; runtime functions are unchanged.

# cytoarray 0.0.4

- Keep one row per shared printed measurement by default (111 ARY028 analytes).
  `normalize_array(split_shared = TRUE)` requests the legacy split format.
- Default `diff_array()` to `trend = TRUE, robust = TRUE`, matching manuscript
  Fig.6D MC38 comparisons; preserve positive 2% offset and per-contrast BH FDR.
- Refuse inference when any group has fewer than two sample columns. Columns must
  be independent biological samples; duplicate spots and repeated scans are technical.
- Provide English README, starter guide and editable folder workflow, with dependency references.
- Leave pixel preprocessing, registration, refinement and spot measurement unchanged.
  Exact manuscript replay uses saved final coordinates, not new auto-detection.

# cytoarray 0.0.3

- Remove result heatmap generation: `plot_array_heatmap()`, its namespace export,
  help page, and the `pheatmap` dependency.
- Retain image QC overlays, registration, manual correction, spot measurement,
  matrix/table export and group differential analysis. Their function arguments,
  defaults and bodies are unchanged from the source snapshot used for this release.
- Replace the local project README with portable installation and usage guidance.
- Restore the user-specified author identity (Junhao Wang) and GPL-3 metadata.
  Do not carry the source snapshot's conflicting MIT `LICENSE` file into this package.
- Use version 0.0.3 because two code states already carried version 0.0.2:
  the installed package lacks the newer source's incomplete-spot warning.
  This candidate retains that warning and is not yet published or globally installed.
- Shared-measurement splitting remains the legacy matrix-export behavior;
  the README explains the statistical implications and provides an explicit ARY028 adapter.
