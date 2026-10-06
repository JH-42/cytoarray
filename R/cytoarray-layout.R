# Array layout metadata. A "layout" bundles everything product-specific so the
# image/analysis code stays generic: corner reference template coords, the
# reference & negative control positions, and the position->analyte map
# (with the two duplicate spot positions per analyte).

#' Built-in layout for the R&D Systems Proteome Profiler Mouse XL array (ARY028)
#' @return a list with `ref`, `ref_pos`, `neg_pos`, `analyte_map`, `splits`
#' @export
ary028_layout <- function() {
  ref <- data.frame(
    tx   = c(80, 80, 636, 638, 628, 628),
    ty   = c(128, 188, 130, 188, 1470, 1526),
    zone = c("TL1", "TL2", "TR1", "TR2", "BR1", "BR2"),
    stringsAsFactors = FALSE
  )
  pos1 <- c("A1","A3","A5","A7","A9","A11","A13","A15","A17","A19","A21","A23",
            "B3","B5","B7","B9","B11","B13","B15","B17","B19","B21",
            "C3","C5","C7","C9","C11","C13","C15","C17","C19","C21",
            "D1","D3","D5","D7","D9","D11","D13","D15","D17","D19","D21","D23",
            "E1","E3","E5","E7","E9","E11","E13","E15","E17","E19","E21","E23",
            "F1","F3","F5","F7","F9","F11","F13","F15","F17","F19","F21","F23",
            "G1","G3","G5","G7","G9","G11","G13","G15","G17","G19","G21","G23",
            "H1","H3","H5","H7","H9","H11","H13","H15","H17","H19","H21","H23",
            "I1","I3","I5","I7","I9","I11","I13","I15","I17","I19","I21",
            "J1","J3","J5","J7","J9","J11","J13","J15","J17","J19","J21","J23")
  pos2 <- vapply(pos1, function(p) {
    L <- sub("[0-9]+$", "", p); N <- as.integer(sub("^[A-J]", "", p)); paste0(L, N + 1)
  }, character(1))
  analyte <- c("Reference","Adipoq","Areg","Angpt1","Angpt2","Angptl3","Tnfsf13b","Cd93","Ccl2",
               "Ccl3/Ccl4","Ccl5","Reference","Ccl6","Ccl11","Ccl12","Ccl17","Ccl19","Ccl20",
               "Ccl21a","Ccl22","Cd14","Cd40","Cd160","Rarres2","Chi3l1","F3","C5/Hc","Cfd","Crp",
               "Cx3cl1","Cxcl1","Cxcl2","Cxcl9","Cxcl10","Cxcl11","Cxcl13","Cxcl16","Cst3","Dkk1",
               "Dpp4","Egf","Eng","Col18a1","Ahsg","Fgf1","Fgf21","Flt3l","Gas6","Csf3","Gdf15",
               "Csf2","Hgf","Icam1","Ifng","Igfbp1","Igfbp2","Igfbp3","Igfbp5","Igfbp6","Il1a",
               "Il1b","Il1rn","Il2","Il3","Il4","Il5","Il6","Il7","Il10","Il11","Il12b","Il13",
               "Il15","Il17a","Il22","Il23a","Il27","Il28a/Il28b","Il33","Ldlr","Lep","Lif","Lcn2",
               "Cxcl5","Csf1","Mmp2","Mmp3","Mmp9","Mpo","Spp1","Tnfrsf11b","Tymp","Pdgfb","Apcs",
               "Ptx3","Postn","Dlk1","Prl2c2","Pcsk9","Ager","Rbp4","Reg3g","Retn","Reference",
               "Sele","Selp","Serpine1","Serpinf1","Thpo","Havcr1","Tnf","Vcam1","Vegfa","Wisp1",
               "Negative Control")
  list(
    product     = "ARY028",
    ref         = ref,
    ref_pos     = c("A1","A2","A23","A24","J1","J2"),
    neg_pos     = c("J23","J24"),
    analyte_map = data.frame(pos1 = pos1, pos2 = pos2, analyte = analyte, stringsAsFactors = FALSE),
    splits      = list(c("Ccl3/Ccl4","Ccl3","Ccl4"),
                       c("Il28a/Il28b","Il28a","Il28b"),
                       c("C5/Hc","C5","Hc"))
  )
}

#' Load a layout: per-spot template coords + product metadata
#'
#' @param template_csv path to a CSV with columns tx, ty, tr, position, spot_id
#'   (one row per printed spot; same coordinate frame as the layout's `ref`).
#'   Defaults to the bundled ARY028 template (works out of the box once installed).
#' @param product one of the built-in products (currently "ARY028"), or a list
#'   returned by [ary028_layout()] for a custom array.
#' @param default_radius_px fallback spot radius when template has no `tr` column.
#' @return a list = product metadata + `$spots` (the per-spot template table)
#' @export
load_layout <- function(template_csv = NULL, product = "ARY028", default_radius_px = 12) {
  meta <- if (is.list(product)) product
          else switch(product, ARY028 = ary028_layout(),
                      stop("Unknown product: ", product))
  if (is.null(template_csv))
    template_csv <- system.file("extdata", "ARY028_template.csv", package = "cytoarray")
  if (!nzchar(template_csv) || !file.exists(template_csv))
    stop("template_csv not found; pass a path, or install the package for the bundled ARY028 template.")
  spots <- utils::read.csv(template_csv, stringsAsFactors = FALSE)
  need <- c("tx", "ty", "position", "spot_id")
  if (!all(need %in% names(spots)))
    stop("template_csv must contain columns: ", paste(need, collapse = ", "))
  if (is.null(spots$tr)) spots$tr <- default_radius_px
  meta$spots <- spots
  meta
}
