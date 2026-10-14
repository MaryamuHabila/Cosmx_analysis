
WORKS 
## =============================================================================
## 03_QC_low_counts_per_cell.R
## =============================================================================
## PURPOSE
##   Vignette step 4.1, first filter only (total counts per cell), 6000-plex panel:
##     threshold = min(100, 10th percentile of nCount_RNA)
##     flag      = nCount_RNA < threshold
##   Approach A: one threshold for the whole slide (vignette as written).
##   Approach B: the same rule applied separately within each tissue section
##               (6 patients: AAA2, AAA38, AAA40, AAA45, AAA46, AAA48).
##   Cells are only FLAGGED here. Removal is vignette step 4.4.
##
## RUN AFTER 02_seurat_inspect_and_split_sections.R (needs metadata with section_id).
## HOW TO RUN: one chunk at a time (select, Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 15: settings, input, output folders
## =============================================================================
suppressPackageStartupMessages(library(data.table))   # tables and grouped summaries

output_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"

## Reuse metadata from script 02 if it is in the session; otherwise reload the saved copy.
if (!exists("metadata")) metadata <- readRDS(file.path(output_dir, "metadata_unfiltered.RDS"))
stopifnot("section_id" %in% colnames(metadata))       # section labels from Chunk 9
stopifnot(identical(rownames(metadata), metadata$cell_id))
section_names <- sort(unique(metadata$section_id))    # AAA2 ... AAA48

## Vignette recommendations, 6000plex: 50 (permissive) or 100 (conservative);
## all panels: no more than 10% of cells (conservative).
count_min <- 100    # count threshold
q_cap     <- 0.10   # never flag more than the lowest 10% of cells

## Output folders
dir_qc  <- file.path(output_dir, "QC")
dir_all <- file.path(dir_qc, "QC_low_counts_per_cell_all_cells")
dir_sec <- file.path(dir_qc, "QC_low_counts_per_tissue_section")
dir.create(dir_all, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_sec, showWarnings = FALSE, recursive = TRUE)

## Cells with 0 counts: log2(0) = -Inf, so hist() silently omits them from the
## plot, but they ARE counted in the flags below.
sum(metadata$nCount_RNA == 0)

## Clear any leftover graphics state (e.g. an old par(mfrow)) that can cause
## "figure margins too large".
graphics.off()


## =============================================================================
## CHUNK 16: helper functions
## =============================================================================
## Threshold rule from the vignette: min(count threshold, 10th percentile).
qc_threshold <- function(x) min(count_min, unname(quantile(x, q_cap)))

## One summary row for a group of cells.
##   x   = nCount_RNA of those cells
##   thr = threshold applied
##   q10 = 10th percentile used in the threshold rule
qc_row <- function(x, thr, q10, group) {
  thr <- unname(thr); q10 <- unname(q10)
  n   <- length(x)                           # total cells
  nf  <- sum(x < thr)                        # flagged cells
  data.table(
    group              = group,
    n_cells            = n,
    q10_counts         = q10,                                  # 10th percentile value
    threshold          = thr,                                  # min(100, q10)
    threshold_source   = ifelse(count_min <= q10, paste0(count_min, " counts"), "10th percentile"),
    n_below_100_counts = sum(x < count_min),                   # cells < 100 counts
    pct_below_100      = round(100 * sum(x < count_min) / n, 2),
    n_below_q10        = sum(x < q10),                         # cells < 10th percentile
    pct_below_q10      = round(100 * sum(x < q10) / n, 2),
    n_flagged          = nf,
    pct_flagged        = round(100 * nf / n, 2),
    n_retained         = n - nf,
    flagged_label      = paste0(format(nf, big.mark = ","), " of ",
                                format(n, big.mark = ","), " cells flagged")
  )
}

## Vignette histogram, plus the flagged label on top and the <100 / <q10 counts.
##   r = one row from qc_row()
qc_hist <- function(x, r, title) {
  hist(log2(x), breaks = 100, xlab = "Log2 counts per cell", ylab = "N cells", main = "")  # vignette
  abline(v = log2(r$threshold), col = "red")                                             # vignette
  legend("topleft", legend = paste0(round(100 * r$n_flagged / r$n_cells, 1), "% rejected"))  # vignette
  legend("topright", bty = "n", cex = 0.8, legend = c(
    paste0("< ", count_min, " counts: ", format(r$n_below_100_counts, big.mark = ","), " cells"),
    paste0("< 10th pct (", round(r$q10_counts, 1), " counts): ",
           format(r$n_below_q10, big.mark = ","), " cells"),
    paste0("threshold = min(", count_min, ", 10th pct) = ", round(r$threshold, 1))))
  mtext(title,           side = 3, line = 2.0, font = 2, cex = 0.9)   # line 1 above plot
  mtext(r$flagged_label, side = 3, line = 0.7, cex = 0.85)            # "X of N cells flagged"
}

## Draw to a fixed-size PNG (never hits pane-size limits), then show the same
## plot in a separate window not constrained by the RStudio Plots pane.
show_and_save <- function(draw_fun, file, width, height) {
  png(file, width = width, height = height, units = "in", res = 300)
  draw_fun()
  dev.off()
  dev.new(width = width, height = height, noRStudioGD = TRUE)
  draw_fun()
  message("Saved: ", file)
}


## =============================================================================
## CHUNK 17: APPROACH A - all cells, one threshold for the whole slide
## =============================================================================
x_all   <- metadata$nCount_RNA          # total counts of every cell
q10_all <- quantile(x_all, q_cap)       # 10th percentile of the whole slide
thr_all <- qc_threshold(x_all)          # min(100, q10_all)

flag_all <- x_all < thr_all             # TRUE = flagged
names(flag_all) <- metadata$cell_id     # keep cell IDs with the flags

## Row 1 = whole slide; rows 2-7 = how that SAME global threshold hits each section.
qc_all_summary <- rbind(
  qc_row(x_all, thr_all, q10_all, "All cells"),
  rbindlist(lapply(section_names, function(sec) {
    qc_row(x_all[metadata$section_id == sec], thr_all, q10_all, sec)
  }))
)
stopifnot(sum(flag_all) == qc_all_summary[group == "All cells", n_flagged])   # cross-check
print(qc_all_summary)

## Plot
draw_all <- function() {
  op <- par(mar = c(4.5, 4.5, 4, 1)); on.exit(par(op))   # room for 2 title lines
  qc_hist(x_all, qc_all_summary[group == "All cells"],
          "All cells - one threshold for the whole slide")
}
show_and_save(draw_all, file.path(dir_all, "QC_low_counts_all_cells_hist.png"), 9, 6)


## =============================================================================
## CHUNK 18: APPROACH B - each tissue section (patient) with its own threshold
## =============================================================================
q10_sec  <- thr_sec <- setNames(numeric(length(section_names)), section_names)
flag_sec <- setNames(logical(nrow(metadata)), metadata$cell_id)   # all FALSE to start

for (sec in section_names) {
  idx <- metadata$section_id == sec          # cells of this section
  x   <- metadata$nCount_RNA[idx]
  q10_sec[sec]  <- quantile(x, q_cap)        # this section's 10th percentile
  thr_sec[sec]  <- qc_threshold(x)           # min(100, this section's q10)
  flag_sec[idx] <- x < thr_sec[sec]
}

## One row per section, then a total row (sum across sections).
qc_sec_summary <- rbindlist(lapply(section_names, function(sec) {
  qc_row(metadata$nCount_RNA[metadata$section_id == sec], thr_sec[sec], q10_sec[sec], sec)
}))
tot <- qc_sec_summary[, .(n_cells = sum(n_cells), n_below_100_counts = sum(n_below_100_counts),
                          n_below_q10 = sum(n_below_q10), n_flagged = sum(n_flagged))]
tot[, `:=`(group = "All sections (sum)",
           pct_below_100 = round(100 * n_below_100_counts / n_cells, 2),
           pct_below_q10 = round(100 * n_below_q10 / n_cells, 2),
           pct_flagged   = round(100 * n_flagged / n_cells, 2),
           n_retained    = n_cells - n_flagged,
           flagged_label = paste0(format(n_flagged, big.mark = ","), " of ",
                                  format(n_cells, big.mark = ","), " cells flagged"))]
qc_sec_summary <- rbind(qc_sec_summary, tot, use.names = TRUE, fill = TRUE)  # threshold = NA for total
stopifnot(sum(flag_sec) == qc_sec_summary[group == "All sections (sum)", n_flagged])  # cross-check
print(qc_sec_summary)

## Plot: 2 x 3 grid, one panel per section, each with its own threshold.
draw_sec <- function() {
  op <- par(mfrow = c(2, 3), mar = c(4.5, 4.5, 4, 1), oma = c(0, 0, 2.5, 0)); on.exit(par(op))
  for (sec in section_names) {
    qc_hist(metadata$nCount_RNA[metadata$section_id == sec],
            qc_sec_summary[group == sec], sec)
  }
  s <- qc_sec_summary[group == "All sections (sum)"]
  mtext(paste0("Per-section thresholds: ", s$flagged_label, " (", s$pct_flagged, "%)"),
        side = 3, outer = TRUE, line = 0.8, font = 2)
}
show_and_save(draw_sec, file.path(dir_sec, "QC_low_counts_per_section_grid.png"), 15, 9)


## =============================================================================
## CHUNK 19: compare A vs B, print, save
## =============================================================================
## Per section: flagged under A, under B, both, and only one approach.
comparison <- data.table(section_id = metadata$section_id, A = flag_all, B = flag_sec)[, .(
  n_cells               = .N,
  flagged_A_global      = sum(A),
  flagged_B_per_section = sum(B),
  flagged_both          = sum(A & B),
  only_A                = sum(A & !B),
  only_B                = sum(!A & B)
), by = section_id][order(section_id)]
comparison <- rbind(comparison,
                    comparison[, lapply(.SD, sum), .SDcols = -"section_id"][, section_id := "All cells"],
                    use.names = TRUE)
print(comparison)

## Console summary
message("Approach A (whole slide):  ", qc_all_summary[group == "All cells", flagged_label],
        " (", qc_all_summary[group == "All cells", pct_flagged], "%)")
message("Approach B (per section):  ", qc_sec_summary[group == "All sections (sum)", flagged_label],
        " (", qc_sec_summary[group == "All sections (sum)", pct_flagged], "%)")

## Save tables and flags (flags are named by cell_id; reload with readRDS()).
fwrite(qc_all_summary, file.path(dir_all, "QC_low_counts_all_cells_summary.csv"))
saveRDS(flag_all,      file.path(dir_all, "flag_low_counts_all_cells.RDS"))
fwrite(qc_sec_summary, file.path(dir_sec, "QC_low_counts_per_section_summary.csv"))
saveRDS(flag_sec,      file.path(dir_sec, "flag_low_counts_per_section.RDS"))
fwrite(comparison,     file.path(dir_qc,  "QC_low_counts_A_vs_B_comparison.csv"))
message("QC outputs saved to ", dir_qc)




graphics.off()                                          # close any plot device left open by an error
dir_qc                                                  # full path of the QC folder
list.files(dir_qc, recursive = TRUE, full.names = TRUE) # every file saved inside it


















## =============================================================================
## 04_QC_cell_area.R
## =============================================================================
## PURPOSE
##   Vignette step 4.1 (second part): flag cells with excessively large areas
##   (likely segmentation errors), then add them to the low-count flags:
##     area_threshold <- 30000                       # pixels, as in the vignette
##     flag <- flag | (metadata$Area > area_threshold)
##
##   Reported two ways, matching script 03:
##     A) all cells      -> combined with flag_all     (one low-count threshold)
##     B) per section    -> combined with flag_section (low-count threshold per section)
##
## NOTE: the area threshold is a FIXED number (30,000 px), not a percentile.
##   So the SAME cells are area-flagged in A and B; only the tallies per
##   section and the combined totals (which include different low-count flags) differ.
##
## Run after 03_QC_low_counts_per_cell.R. Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 20: inputs, settings, folders
## =============================================================================
library(data.table)

## If starting a fresh session, reload metadata and the low-count flags:
# output_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"
# metadata   <- readRDS(file.path(output_dir, "metadata_unfiltered.RDS"))
# flag_all   <- readRDS(file.path(output_dir, "QC", "QC_low_counts_per_cell_all_cells",  "flag_low_counts_all_cells.RDS"))
# flag_sec   <- readRDS(file.path(output_dir, "QC", "QC_low_counts_per_tissue_section", "flag_low_counts_per_section.RDS"))

section_names <- sort(unique(metadata$section_id))   # the 6 sections

flag_section <- flag_sec   # script 03 named the per-section low-count flags "flag_sec"

## Low-count flags must be in the same cell order as metadata
stopifnot(identical(names(flag_all),     metadata$cell_id))
stopifnot(identical(names(flag_section), metadata$cell_id))

## Vignette uses metadata$Area, in PIXELS. Check the column exists.
stopifnot("Area" %in% colnames(metadata))
stopifnot(!anyNA(metadata$Area))

area_threshold <- 30000                   # vignette value, pixels
px_um          <- 0.120280945             # µm per pixel (same value as script 02, Chunk 6)
area_threshold_um2 <- area_threshold * px_um^2   # same threshold in µm², for reporting only
area_threshold_um2                                # ~434 µm²

## Sanity check: Area (px) x pixel area should reproduce AtoMx's Area.um2
if ("Area.um2" %in% colnames(metadata)) {
  print(all.equal(metadata$Area * px_um^2, metadata$Area.um2, tolerance = 1e-3))  # TRUE = consistent
}

## Output folders: QC/Area/<approach>
dir_area         <- file.path(output_dir, "QC", "Area")
dir_area_all     <- file.path(dir_area, "Area_per_cell_all_cells")
dir_area_section <- file.path(dir_area, "Area_per_tissue_section")
dir.create(dir_area_all,     showWarnings = FALSE, recursive = TRUE)
dir.create(dir_area_section, showWarnings = FALSE, recursive = TRUE)

## Area flag (identical for both approaches, because the threshold is fixed)
flag_area <- setNames(metadata$Area > area_threshold, metadata$cell_id)

graphics.off()   # clear any plot device left open (avoids "figure margins too large")


## =============================================================================
## CHUNK 21: helper functions
## =============================================================================
## fmt(): numbers with thousands separators, e.g. 123,709
fmt <- function(n) format(n, big.mark = ",")

## area_row(): one summary row for a group of cells.
##   area_flag = area flags of the group; low_flag = low-count flags of the group
area_row <- function(area_flag, low_flag, label) {
  n        <- length(area_flag)
  combined <- low_flag | area_flag                    # vignette: flag <- flag | (Area > threshold)
  data.table(
    group                   = label,
    n_cells_total           = n,
    area_threshold_px       = area_threshold,
    area_threshold_um2      = round(area_threshold_um2, 1),
    n_area_flagged          = sum(area_flag),         # cells > 30,000 px
    pct_area_flagged        = round(100 * mean(area_flag), 2),
    n_lowcount_flagged      = sum(low_flag),          # from script 03
    n_area_and_lowcount     = sum(area_flag & low_flag),   # already flagged by low counts
    n_area_only_new         = sum(area_flag & !low_flag),  # newly removed by the area step
    n_flagged_combined      = sum(combined),          # total removed so far
    pct_flagged_combined    = round(100 * mean(combined), 2),  # vignette: mean(flag)
    n_kept                  = n - sum(combined),
    area_flagged_label      = paste0(fmt(sum(area_flag)), " of ", fmt(n), " cells flagged (area)"),
    combined_flagged_label  = paste0(fmt(sum(combined)),  " of ", fmt(n), " cells flagged (low counts + area)")
  )
}

## area_hist(): vignette histogram + red line + "% rejected" legend (top right).
area_hist <- function(area, title, cex_leg = 1) {
  hist(area, breaks = 100, main = title, xlab = "Cell area (pixels)")
  abline(v = area_threshold, col = "red")
  legend("topright", cex = cex_leg, legend = c(
    paste0(round(mean(area > area_threshold) * 100, 2), "% rejected"),   # vignette text
    paste0(fmt(sum(area > area_threshold)), " cells > ", fmt(area_threshold), " px")))
}

## removed_table(): one row per area-flagged cell, with its area.
removed_table <- function(keep_rows, low_flag) {
  d <- as.data.table(metadata[keep_rows, c("cell_id", "section_id", "fov", "Area", "nCount_RNA")])
  if ("Area.um2" %in% colnames(metadata)) d[, Area.um2 := metadata$Area.um2[keep_rows]]
  d[, also_lowcount_flagged := low_flag[keep_rows]]   # TRUE = was already removed by low counts
  d[order(section_id, -Area)]                          # largest first within each section
}

## show_and_save(): same as script 03. Save PNG at a fixed size, then show the
## plot in a separate window not limited by the RStudio Plots pane.
show_and_save <- function(draw_fun, file, width, height) {
  png(file, width = width, height = height, units = "in", res = 300)
  draw_fun()
  dev.off()
  dev.new(width = width, height = height, noRStudioGD = TRUE)
  draw_fun()
  message("Saved: ", file)
}


## =============================================================================
## CHUNK 22: APPROACH A - all cells
## =============================================================================
area_all_summary <- rbind(
  area_row(flag_area, flag_all, "All cells"),
  rbindlist(lapply(section_names, function(sec) {
    idx <- metadata$section_id == sec
    area_row(flag_area[idx], flag_all[idx], sec)      # same thresholds, tallied per section
  }))
)

## Combined flag for approach A (carried forward to FOV QC)
flag_all_combined <- flag_all | flag_area
stopifnot(sum(flag_all_combined) == area_all_summary[group == "All cells", n_flagged_combined])

## Plot 1: area only (cells removed because Area > 30,000 px)
draw_area_all <- function() {
  op <- par(mar = c(4.5, 4.5, 4.5, 1)); on.exit(par(op))   # top margin fits the 2-line title
  r <- area_all_summary[group == "All cells"]
  area_hist(metadata$Area, paste0("All cells\n", r$area_flagged_label))
}
show_and_save(draw_area_all, file.path(dir_area_all, "hist_cell_area_all_cells.png"), 9, 6)

## Plot 2: low counts + area combined (total removed so far)
draw_combined_all <- function() {
  op <- par(mar = c(4.5, 4.5, 4.5, 1)); on.exit(par(op))
  r <- area_all_summary[group == "All cells"]
  hist(metadata$Area, breaks = 100,
       main = paste0("All cells\n", r$combined_flagged_label),   # low counts + area, of total cells
       xlab = "Cell area (pixels)")
  abline(v = area_threshold, col = "red")                        # area threshold
  legend("topright", legend = c(
    paste0(r$pct_flagged_combined, "% rejected"),                # combined % removed
    paste0("low counts: ", fmt(r$n_lowcount_flagged)),           # from script 03
    paste0("area: ",       fmt(r$n_area_flagged)),               # Area > 30,000 px
    paste0("in both: ",    fmt(r$n_area_and_lowcount))))         # counted once in the total
}
show_and_save(draw_combined_all, file.path(dir_area_all, "hist_lowcounts_plus_area_all_cells.png"), 9, 6)

## Table of every area-flagged cell
area_removed_all <- removed_table(which(flag_area), flag_all)

## Save
fwrite(area_all_summary, file.path(dir_area_all, "Area_QC_all_cells_summary.csv"))
fwrite(area_removed_all, file.path(dir_area_all, "Area_QC_all_cells_removed_cells.csv"))
saveRDS(flag_all_combined, file.path(dir_area_all, "flag_lowcounts_plus_area_all_cells.RDS"))

## Print
print(area_all_summary)
print(area_removed_all)


## =============================================================================
## CHUNK 23: APPROACH B - per tissue section
## =============================================================================
rows <- lapply(section_names, function(sec) {
  idx <- metadata$section_id == sec
  area_row(flag_area[idx], flag_section[idx], sec)    # low-count flags from per-section thresholds
})
area_section_summary <- rbindlist(rows)

## Combined row: sums over the 6 sections
tot <- area_section_summary[, .(
  group               = "All sections combined",
  n_cells_total       = sum(n_cells_total),
  area_threshold_px   = area_threshold,
  area_threshold_um2  = round(area_threshold_um2, 1),
  n_area_flagged      = sum(n_area_flagged),
  n_lowcount_flagged  = sum(n_lowcount_flagged),
  n_area_and_lowcount = sum(n_area_and_lowcount),
  n_area_only_new     = sum(n_area_only_new),
  n_flagged_combined  = sum(n_flagged_combined),
  n_kept              = sum(n_kept))]
tot[, `:=`(pct_area_flagged       = round(100 * n_area_flagged     / n_cells_total, 2),
           pct_flagged_combined   = round(100 * n_flagged_combined / n_cells_total, 2),
           area_flagged_label     = paste0(fmt(n_area_flagged),     " of ", fmt(n_cells_total), " cells flagged (area)"),
           combined_flagged_label = paste0(fmt(n_flagged_combined), " of ", fmt(n_cells_total), " cells flagged (low counts + area)"))]
area_section_summary <- rbindlist(list(area_section_summary, tot), use.names = TRUE, fill = TRUE)

## Combined flag for approach B
flag_section_combined <- flag_section | flag_area
stopifnot(sum(flag_section_combined) == tot$n_flagged_combined)

## Grid plot 1: 2 x 3, area only (cells removed because Area > 30,000 px)
draw_area_grid <- function() {
  op <- par(mfrow = c(2, 3),               # 2 rows x 3 columns
            mar = c(4, 4, 4.5, 0.5),       # panel margins: bottom, left, top, right
            oma = c(0.5, 0.5, 2.5, 0.5),   # outer margin; top space for the overall title
            cex = 0.7)                     # smaller text so all 6 panels fit
  on.exit(par(op))                         # restore settings afterwards
  for (sec in section_names) {
    r <- area_section_summary[group == sec]
    area_hist(metadata$Area[metadata$section_id == sec],
              paste0(sec, "\n", r$area_flagged_label),       # e.g. "AAA2 / 120 of 20,000 cells flagged (area)"
              cex_leg = 0.8)
  }
  mtext(paste0("Area > ", fmt(area_threshold), " px - all sections: ", tot$area_flagged_label),
        outer = TRUE, font = 2, cex = 1)
}
show_and_save(draw_area_grid, file.path(dir_area_section, "hist_cell_area_per_section_grid.png"), 12, 8)

## Grid plot 2: 2 x 3, low counts + area combined (total removed so far)
draw_combined_grid <- function() {
  op <- par(mfrow = c(2, 3), mar = c(4, 4, 4.5, 0.5), oma = c(0.5, 0.5, 2.5, 0.5), cex = 0.7)
  on.exit(par(op))
  for (sec in section_names) {
    r <- area_section_summary[group == sec]
    hist(metadata$Area[metadata$section_id == sec], breaks = 100,
         main = paste0(sec, "\n", r$combined_flagged_label),  # low counts + area, of total cells
         xlab = "Cell area (pixels)")
    abline(v = area_threshold, col = "red")                   # area threshold
    legend("topright", cex = 0.8, legend = c(
      paste0(r$pct_flagged_combined, "% rejected"),           # combined % removed
      paste0("low counts: ", fmt(r$n_lowcount_flagged)),      # from script 03
      paste0("area: ",       fmt(r$n_area_flagged)),          # Area > 30,000 px
      paste0("in both: ",    fmt(r$n_area_and_lowcount))))    # counted once in the total
  }
  mtext(paste0("Low counts + area - all sections: ", tot$combined_flagged_label),
        outer = TRUE, font = 2, cex = 1)
}
show_and_save(draw_combined_grid, file.path(dir_area_section, "hist_lowcounts_plus_area_per_section_grid.png"), 12, 8)

## Table of every area-flagged cell (also_lowcount_flagged uses per-section flags here)
area_removed_section <- removed_table(which(flag_area), flag_section)

## Save
fwrite(area_section_summary, file.path(dir_area_section, "Area_QC_per_section_summary.csv"))
fwrite(area_removed_section, file.path(dir_area_section, "Area_QC_per_section_removed_cells.csv"))
saveRDS(flag_section_combined, file.path(dir_area_section, "flag_lowcounts_plus_area_per_section.RDS"))

## Print
print(area_section_summary)
print(area_removed_section)
area_removed_section[, .N, by = section_id]   # area-flagged cells per section


















## =============================================================================
## 05_1_QC_FOV.R
## =============================================================================






## =============================================================================
## 05_1_QC_FOV.R (Used)   (alternative to 05_QC_FOV.R - run ONE of the two, not both)
## =============================================================================
## PURPOSE
##   FOV QC exactly as described in the CosMx Scratch Space post
##   "FOV QC from single-cell gene expression in spatial dataset" (Danaher, 2024),
##   with no changes to the method:
##     - 7x7 grid per FOV (squares with < 10 cells dropped)
##     - 10 most similar squares in OTHER FOVs, max 1 per FOV
##     - signal loss: FOV flagged if > 75% of its squares lost > 60% of expected counts
##     - bias: FOV flagged if >= 2 bits of one reporter cycle lose > 60% of expected signal
##   ">60% loss" is the function default max_totalcounts_loss = max_prop_loss = 0.6
##   (i.e. the "< 40% of expected" cut-off). It is not an extra step; it is how the
##   post defines a "low" square, so it is kept unchanged.
##
##   DIFFERENCE FROM 05_QC_FOV.R: runFOVQC() is run SEPARATELY for each tissue
##   section, so comparator squares come only from the same section
##   (post: "We advise this approach be applied separately to each slide or tissue").
##   NB: the runFOVQC "tissue" argument does NOT do this - it only renames FOVs,
##   comparators are still searched across all tissues. Hence the loop below.
##
##   The 6 per-section results are merged into ONE object, "fovqcresult", with the
##   same structure and the same variable/file names as 05_QC_FOV.R, so script 06
##   onwards runs unchanged. Outputs go to QC/05_1_QC_FOV.
##
## Run after script 02 (needs counts, xy, metadata with section_id) in the session.
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 24: packages, settings, folders
## =============================================================================
## Needed once: install.packages(c("FNN", "scales", "circlize", "ggplot2"))
##              BiocManager::install("ComplexHeatmap")
library(data.table)       # tables
library(Matrix)           # sparse counts
library(ggplot2)          # bar plot (vignette 2)
library(ComplexHeatmap)   # bias heatmap (vignette 2)
library(circlize)         # colorRamp2() for the heatmap colours (vignette 2)

output_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"

## If starting a fresh session, reload the unfiltered objects saved by script 02:
# counts   <- readRDS(file.path(output_dir, "counts_unfiltered.RDS"))
# xy       <- readRDS(file.path(output_dir, "xy_unfiltered.RDS"))
# metadata <- readRDS(file.path(output_dir, "metadata_unfiltered.RDS"))

## FOV QC runs on all cells (before any cells are removed).
stopifnot(identical(rownames(counts), metadata$cell_id))
stopifnot(identical(rownames(xy),     metadata$cell_id))
section_names <- sort(unique(metadata$section_id))

## Output folder (same variable name as script 05, new folder)
dir_fov <- file.path(output_dir, "QC", "05_1_QC_FOV")
dir.create(dir_fov, showWarnings = FALSE, recursive = TRUE)

graphics.off()   # clear any open plot device


## =============================================================================
## CHUNK 25: FOV QC tool and barcode map (vignette)
## =============================================================================
## Source the FOV QC functions (runFOVQC, mapFlaggedFOVs, FOVSignalLossSpatialPlot, ...)
source("https://raw.githubusercontent.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/Main/_code/FOV%20QC/FOV%20QC%20utils.R")

## Thresholds = the function defaults (the post's ">60% loss"). Read them from the
## function itself so they are exactly what is used; kept as variables, as in script 05.
max_totalcounts_loss <- formals(runFOVQC)$max_totalcounts_loss   # 0.6
max_prop_loss        <- formals(runFOVQC)$max_prop_loss          # 0.6
stopifnot(max_totalcounts_loss == 0.6, max_prop_loss == 0.6)     # stop if the tool's defaults ever change

## Gene -> barcode map for every panel
allbarcodes <- readRDS(url("https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/raw/Main/_code/FOV%20QC/barcodes_by_panel.RDS"))
names(allbarcodes)                    # Hs_IO, Hs_UCC, Hs_6k, Mm_Neuro, Mm_UCC, Hs_WTX
barcodemap <- allbarcodes$Hs_6k       # <------- 6k panel
head(barcodemap)

## Genes shared by the panel map and the counts matrix (runFOVQC uses only these)
length(intersect(barcodemap$gene, colnames(counts)))


## =============================================================================
## CHUNK 26: run FOV QC per tissue section (vignette call, defaults)
## =============================================================================
## Row indices of the cells in each section, in section_names order.
idx_by_section <- split(seq_len(nrow(counts)), metadata$section_id)[section_names]

## FOVs per section. Each square needs 10 OTHER FOVs in its section for 10
## comparators; with <= 10 FOVs it gets fewer (the tool averages what it finds).
n_fov_per_section <- sapply(idx_by_section, function(i) length(unique(metadata$fov[i])))
print(n_fov_per_section)
if (any(n_fov_per_section <= 10)) {
  warning("Section(s) with <= 10 FOVs (fewer than 10 comparators possible): ",
          paste(names(n_fov_per_section)[n_fov_per_section <= 10], collapse = ", "))
}

## Run runFOVQC on one section at a time: counts, xy and fov of that section only,
## so comparator squares can only come from FOVs of the same section.
fovqcresult_by_section <- list()
for (sec in section_names) {
  i <- idx_by_section[[sec]]                                   # cells of this section
  message("---- FOV QC: ", sec, " (", length(i), " cells, ", n_fov_per_section[[sec]], " FOVs)")
  fovqcresult_by_section[[sec]] <- runFOVQC(counts     = counts[i, , drop = FALSE],
                                            xy         = xy[i, , drop = FALSE],
                                            fov        = metadata$fov[i],
                                            barcodemap = barcodemap)   # defaults: 0.6 / 0.6
}
saveRDS(fovqcresult_by_section, file.path(dir_fov, "fovqcresult_by_section.RDS"))   # per-section results


## =============================================================================
## CHUNK 26b: merge the 6 section results into one "fovqcresult" (same structure
##            as a single runFOVQC() output, so all later code and plots work)
## =============================================================================
## FOV IDs are unique across the slide (checked in chunk 27), so grid IDs
## ("fov_x_y") are unique too and the pieces can simply be stacked.
res_list   <- unname(fovqcresult_by_section)                                      # drop list names so they are not pasted onto element names
cell_order <- order(unlist(idx_by_section, use.names = FALSE))                    # puts stacked cells back in original order
stack_vec  <- function(get) do.call(c, lapply(res_list, get))                     # join vectors, keep element names
stack_mat  <- function(get) do.call(rbind, lapply(res_list, get))                 # join matrices by row

fovqcresult <- list(
  flaggedfovs                   = stack_vec(function(r) r$flaggedfovs),
  flaggedfovs_fortotalcounts    = stack_vec(function(r) r$flaggedfovs_fortotalcounts),
  flaggedfovs_forbias           = stack_vec(function(r) r$flaggedfovs_forbias),
  flagged_fov_x_gene            = stack_mat(function(r) r$flagged_fov_x_gene),    # NULL if no FOV flagged for bias
  flags_per_fov_x_reportercycle = stack_mat(function(r) r$flags_per_fov_x_reportercycle),
  fovstats = list(flag      = stack_mat(function(r) r$fovstats$flag),
                  bias      = stack_mat(function(r) r$fovstats$bias),
                  p         = stack_mat(function(r) r$fovstats$p),
                  propagree = stack_mat(function(r) r$fovstats$propagree)),
  resid                         = stack_mat(function(r) r$resid),
  totalcountsresids             = stack_vec(function(r) setNames(as.vector(r$totalcountsresids),
                                                                 names(r$totalcountsresids))),   # 1-d array -> named vector
  gridinfo = list(gridid  = stack_vec(function(r) r$gridinfo$gridid)[cell_order],  # per cell, original order
                  gridfov = stack_vec(function(r) r$gridinfo$gridfov)),            # per grid square
  xy                            = stack_mat(function(r) r$xy)[cell_order, , drop = FALSE],
  fov                           = stack_vec(function(r) r$fov)[cell_order]
)
## Empty flag lists come back as NULL from c(); make them character(0) like runFOVQC.
for (nm in c("flaggedfovs", "flaggedfovs_fortotalcounts", "flaggedfovs_forbias")) {
  if (is.null(fovqcresult[[nm]])) fovqcresult[[nm]] <- character(0)
}

## Checks: merged object lines up with the original cells
stopifnot(identical(fovqcresult$fov, as.character(metadata$fov)))
stopifnot(identical(rownames(fovqcresult$xy), rownames(xy)))
stopifnot(!anyDuplicated(rownames(fovqcresult$resid)))

saveRDS(fovqcresult, file.path(dir_fov, "fovqcresult.RDS"))   # same file name as script 05

## FOVs flagged for any reason, for loss of signal (low counts), for bias
fovqcresult$flaggedfovs
fovqcresult$flaggedfovs_fortotalcounts
fovqcresult$flaggedfovs_forbias


## =============================================================================
## CHUNK 27: summary - FOVs flagged of total, whole slide and per section
## =============================================================================
fmt <- function(n) format(n, big.mark = ",")

## One row per FOV: its section and whether it was flagged, and why.
## (runFOVQC returns FOV IDs as text, so compare as text.)
fov_table <- unique(data.table(fov = as.character(metadata$fov), section_id = metadata$section_id))
stopifnot(!anyDuplicated(fov_table$fov))                     # each FOV in one section only
fov_table[, flagged_low_counts := fov %in% fovqcresult$flaggedfovs_fortotalcounts]
fov_table[, flagged_bias       := fov %in% fovqcresult$flaggedfovs_forbias]
fov_table[, flagged_any        := fov %in% fovqcresult$flaggedfovs]
fov_table[, reason := fifelse(flagged_low_counts & flagged_bias, "low counts + bias",
                              fifelse(flagged_low_counts, "low counts",
                                      fifelse(flagged_bias, "bias", "not flagged")))]
setorder(fov_table, section_id, fov)

## Per section: flagged FOVs of total FOVs
fov_qc_summary <- fov_table[, .(
  n_fov_total         = .N,
  n_fov_flagged_any   = sum(flagged_any),
  n_fov_low_counts    = sum(flagged_low_counts),
  n_fov_bias          = sum(flagged_bias),
  flagged_fov_ids     = paste(sort(as.integer(fov[flagged_any])), collapse = ", ")
), by = section_id]

## Whole-slide row
fov_qc_summary <- rbind(fov_qc_summary, fov_table[, .(
  section_id        = "All sections",
  n_fov_total       = .N,
  n_fov_flagged_any = sum(flagged_any),
  n_fov_low_counts  = sum(flagged_low_counts),
  n_fov_bias        = sum(flagged_bias),
  flagged_fov_ids   = paste(sort(as.integer(fov[flagged_any])), collapse = ", "))])
fov_qc_summary[, pct_fov_flagged := round(100 * n_fov_flagged_any / n_fov_total, 1)]
fov_qc_summary[, flagged_label   := paste0(n_fov_flagged_any, " of ", n_fov_total, " FOVs flagged")]

print(fov_qc_summary)
print(fov_table[flagged_any == TRUE])   # every flagged FOV, its section and reason

fwrite(fov_qc_summary,               file.path(dir_fov, "FOV_QC_summary_per_section.csv"))
fwrite(fov_table[flagged_any == TRUE], file.path(dir_fov, "FOV_QC_flagged_FOVs.csv"))
fwrite(fov_table,                    file.path(dir_fov, "FOV_QC_all_FOVs.csv"))

tot <- fov_qc_summary[section_id == "All sections"]


## =============================================================================
## CHUNK 28: plotting helpers
## =============================================================================
## Save PNG at a fixed size, then show in a separate window (avoids
## "figure margins too large" in the RStudio Plots pane).
show_and_save <- function(draw_fun, file, width, height) {
  png(file, width = width, height = height, units = "in", res = 300)
  draw_fun()
  dev.off()
  dev.new(width = width, height = height, noRStudioGD = TRUE)
  draw_fun()
  message("Saved: ", file)
}

## Slide-shaped figure size (longest side 11 in).
rx <- diff(range(xy[, 1])); ry <- diff(range(xy[, 2]))
map_w <- rx * 11 / max(rx, ry) + 1
map_h <- ry * 11 / max(rx, ry) + 1.5

## Summary text on the spatial plots: overall line under the title,
## "AAAx: n of N FOVs flagged" above each section.
add_fov_labels <- function(type = c("any", "low_counts")) {
  type <- match.arg(type)
  col_n <- if (type == "any") "n_fov_flagged_any" else "n_fov_low_counts"
  what  <- if (type == "any") "FOVs flagged (any reason)" else "FOVs flagged (low counts)"
  mtext(paste0("All sections (QC run per section): ", tot[[col_n]], " of ", tot$n_fov_total, " ", what),
        side = 3, line = 0.2, cex = 0.8)
  for (sec in section_names) {
    inds <- metadata$section_id == sec
    r    <- fov_qc_summary[section_id == sec]
    text(mean(range(xy[inds, 1])), max(xy[inds, 2]), pos = 3, xpd = NA, cex = 0.7, font = 2,
         labels = paste0(sec, ": ", r[[col_n]], " of ", r$n_fov_total, " FOVs"))
  }
}


## =============================================================================
## CHUNK 29: PLOT 1 - map of flagged FOVs, all 6 sections on one slide (mapFlaggedFOVs)
## =============================================================================
## Blue = all FOVs; red with FOV number = flagged for low counts OR bias.
draw_map <- function() {
  op <- par(mar = c(4, 4, 4, 1)); on.exit(par(op))
  mapFlaggedFOVs(fovqcresult)
  add_fov_labels("any")
}
show_and_save(draw_map, file.path(dir_fov, "01_flagged_FOVs_map.png"), map_w, map_h)


## =============================================================================
## CHUNK 30: PLOT 2 - log2FC of total counts in space, all 6 sections (FOVSignalLossSpatialPlot)
## =============================================================================
## Cells coloured by log2 fold-change in total counts vs the 10 similar squares in
## other FOVs OF THE SAME SECTION; yellow boxes = FOVs flagged for low counts.
draw_signal_loss <- function() {
  op <- par(mar = c(4, 4, 4, 1)); on.exit(par(op))
  FOVSignalLossSpatialPlot(fovqcresult)
  add_fov_labels("low_counts")
}
show_and_save(draw_signal_loss, file.path(dir_fov, "02_signal_loss_low_counts.png"), map_w, map_h)


## =============================================================================
## CHUNK 31: PLOT 3 - reporter bias heatmap (vignette 2)
## =============================================================================
## Vignette 2 function, unchanged except for three guards marked "## guard":
##   drop = FALSE (works with 1 FOV), na.rm = TRUE, and colour breaks kept increasing.
FOVEffectsHeatmap <- function(res, fovs = NULL,
                              cluster_rows = FALSE,
                              cluster_columns = TRUE) {
  if (is.null(fovs)) {
    fovs <- rownames(res$fovstats$bias)
  }
  transposed_matrix <- t(res$fovstats$bias[fovs, , drop = FALSE])   ## guard: drop = FALSE
  lo <- min(transposed_matrix, na.rm = TRUE)                         ## guard: na.rm
  hi <- max(transposed_matrix, na.rm = TRUE)
  if (lo >= 0) lo <- -1e-6                                           ## guard: breaks must increase
  if (hi <= 0) hi <-  1e-6
  col_fun = colorRamp2(
    breaks = c(lo, 0, hi),
    colors = c("magenta", "white", "cyan")
  )
  ht_obj <- ComplexHeatmap::Heatmap(
    matrix = transposed_matrix,
    name = "FOV bias",
    col = col_fun,
    cluster_rows = cluster_rows,
    cluster_columns = cluster_columns,
    show_row_dend = TRUE,
    show_column_dend = TRUE,
    row_names_gp = gpar(fontsize = 7),
    column_names_gp = gpar(fontsize = 8)
  )
  return(ht_obj)
}

## Vignette 2 plots the flagged FOVs only (columns = FOVs, rows = reporter cycle x colour).
if (length(fovqcresult$flaggedfovs) > 0) {
  ht_title <- paste0("Reporter bias, flagged FOVs: ", tot$flagged_label,
                     " (bias: ", tot$n_fov_bias, ", low counts: ", tot$n_fov_low_counts, ")")
  draw_heatmap <- function() {
    draw(FOVEffectsHeatmap(fovqcresult, fovs = fovqcresult$flaggedfovs,
                           cluster_columns = length(fovqcresult$flaggedfovs) > 1),
         column_title = ht_title)
  }
  show_and_save(draw_heatmap, file.path(dir_fov, "03_reporter_bias_heatmap.png"), 10, 15)  # vignette size
} else {
  message("No flagged FOVs: heatmap not drawn.")
}

## Which reporter cycles were flagged in which FOVs (the flags themselves;
## heatmap colours show bias values, not flags). >= 0.5 = cycle flagged.

flags_rc <- fovqcresult$flags_per_fov_x_reportercycle                        # FOVs x reporter cycles (all 6 sections stacked)
flags_rc_flagged <- flags_rc[rowSums(flags_rc >= 0.5) > 0, , drop = FALSE]   # keep FOVs with >= 1 flagged cycle

if (nrow(flags_rc_flagged) > 0) {                                            # guard: only if a cycle was flagged
  print(flags_rc_flagged)
  fwrite(as.data.table(flags_rc_flagged, keep.rownames = "fov"),             # row names (FOV IDs) -> column "fov"
         file.path(dir_fov, "FOV_QC_flagged_reportercycles.csv"))
} else {
  message("No reporter cycles flagged in any FOV (0 FOVs flagged for bias): file not written.")
}
## =============================================================================
## CHUNK 32: PLOT 4 - genes flagged in multiple FOVs (vignette 2 bar plot)
## =============================================================================
## Genes come only from reporter cycles flagged for bias in flagged FOVs.
if (!is.null(fovqcresult$flagged_fov_x_gene)) {
  flagged_genes <- fovqcresult$flagged_fov_x_gene[, 2]
  flagged_genes <- flagged_genes[!is.na(flagged_genes)]
  flagged_genes <- as.data.frame(rev(sort(table(flagged_genes))))
  colnames(flagged_genes)[1] <- "target"
  
  ## Per-gene table for every panel gene: number of FOVs it was flagged in (0 = none).
  ## (Vignette 2 stores this in adata$var$fov_qc_flagged.)
  df <- data.frame(target = colnames(counts))
  df <- merge(df, flagged_genes, by = "target", all.x = TRUE, all.y = FALSE)
  df$Freq[is.na(df$Freq)] <- 0L
  colnames(df)[2] <- "fov_qc_flagged"
  fwrite(df, file.path(dir_fov, "FOV_QC_genes_n_flagged_FOVs.csv"))
  
  p <- ggplot(flagged_genes, aes(x = factor(Freq), fill = factor(Freq))) +
    geom_bar() +
    geom_text(stat = "count", aes(label = after_stat(count)), vjust = -0.5) +
    labs(x = "Number of FOVs", y = "Number of Flagged Genes",
         title = paste0(nrow(flagged_genes), " genes flagged"),
         subtitle = paste0("in ", tot$n_fov_bias, " of ", tot$n_fov_total, " FOVs flagged for bias")) +
    theme_minimal() +
    theme(legend.position = "none")
  print(p)
  ggsave(p, filename = file.path(dir_fov, "04_genes_flagged_in_n_FOVs_barplot.png"),
         width = 5, height = 5, dpi = 150)
} else {
  message("No FOVs flagged for bias: no genes flagged, bar plot not drawn.")
}


## =============================================================================
## CHUNK 33: FOVs flagged, on screen
## =============================================================================
message("Whole slide (QC run per section): ", tot$flagged_label, " (", tot$pct_fov_flagged,
        "%) - low counts: ", tot$n_fov_low_counts, ", bias: ", tot$n_fov_bias)
for (sec in section_names) {
  r <- fov_qc_summary[section_id == sec]
  message(sec, ": ", r$flagged_label, " (", r$pct_fov_flagged, "%) - low counts: ",
          r$n_fov_low_counts, ", bias: ", r$n_fov_bias,
          if (nchar(r$flagged_fov_ids) > 0) paste0(" - FOVs ", r$flagged_fov_ids) else "")
}
message("All FOV QC outputs saved to ", dir_fov)




















## =============================================================================
## 06_QC_ncount_nfeature.R
## =============================================================================
## PURPOSE
##   Vignette 2, section 3.2.1 (Counts and Features): permissive, symmetrical
##   threshold of 2.5 SD around the mean, in log10 space, for
##     nCount_RNA   (transcripts per cell)
##     nFeature_RNA (genes per cell)
##   A cell PASSES if all four hold:
##     features_min <= nFeature_RNA <= features_max  and
##     counts_min   <= nCount_RNA   <= counts_max
##   otherwise it is "filtered".
##
##   A) all cells: one set of thresholds for the whole slide
##   B) per tissue section: thresholds computed within each section (patient)
##   Not combined with the low-count / area / FOV flags here.
##
## Run after script 02 (needs metadata with section_id). Chunk by chunk.
## =============================================================================


## =============================================================================
## CHUNK 34: packages, settings, folders
## =============================================================================
## Needed once: install.packages(c("ggExtra", "scales"))
library(data.table)   # tables
library(ggplot2)      # plots
library(ggExtra)      # ggMarginal(): density curves on the joint plot (vignette)
library(scales)       # log10 axis labels (vignette)

output_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"
# metadata <- readRDS(file.path(output_dir, "metadata_unfiltered.RDS"))   # fresh session only

section_names <- sort(unique(metadata$section_id))
n_sd <- 2.5                                  # vignette: 2.5 SD around the mean in log10 space

## Plot colours (display only; thresholds and counts are unchanged)
qc_cols <- c("pass" = "#1F4E9E", "filtered" = "#E8590C")   # strong blue, strong orange

dir_nc     <- file.path(output_dir, "QC", "ncount_and_nfeature")
dir_nc_all <- file.path(dir_nc, "all_cells")
dir_nc_sec <- file.path(dir_nc, "per_tissue_section")
dir.create(dir_nc_all, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_nc_sec, showWarnings = FALSE, recursive = TRUE)

## Cells with 0 counts or 0 genes: log10(0) = -Inf would make the mean -Inf,
## so thresholds are computed on cells with > 0; zero cells are then
## "filtered" by the minimum thresholds.
sum(metadata$nCount_RNA == 0); sum(metadata$nFeature_RNA == 0)


## =============================================================================
## CHUNK 35: helper functions
## =============================================================================
fmt <- function(n) format(n, big.mark = ",")

## Thresholds (vignette formulas), pass/filtered per cell, one summary row.
##   md = metadata rows of one group; label = group name
nc_nf_qc <- function(md, label) {
  pos <- md$nCount_RNA > 0 & md$nFeature_RNA > 0          # exclude log10(0)
  lc  <- log10(md$nCount_RNA[pos])
  lf  <- log10(md$nFeature_RNA[pos])
  counts_min   <- round(10^(mean(lc) - n_sd * sd(lc)))   # vignette
  features_min <- round(10^(mean(lf) - n_sd * sd(lf)))
  counts_max   <- round(10^(mean(lc) + n_sd * sd(lc)))
  features_max <- round(10^(mean(lf) + n_sd * sd(lf)))
  
  status <- ifelse(md$nFeature_RNA >= features_min & md$nFeature_RNA <= features_max &
                     md$nCount_RNA   >= counts_min   & md$nCount_RNA   <= counts_max,
                   "pass", "filtered")                  # vignette
  n  <- nrow(md); nf <- sum(status == "filtered")
  
  row <- data.table(
    group                = label,
    n_cells              = n,
    min_nCount_RNA       = min(md$nCount_RNA),           # data range
    max_nCount_RNA       = max(md$nCount_RNA),
    min_nFeature_RNA     = min(md$nFeature_RNA),
    max_nFeature_RNA     = max(md$nFeature_RNA),
    counts_min_threshold   = counts_min,                 # chosen thresholds
    counts_max_threshold   = counts_max,
    features_min_threshold = features_min,
    features_max_threshold = features_max,
    n_below_counts_min   = sum(md$nCount_RNA   < counts_min),
    n_above_counts_max   = sum(md$nCount_RNA   > counts_max),
    n_below_features_min = sum(md$nFeature_RNA < features_min),
    n_above_features_max = sum(md$nFeature_RNA > features_max),
    n_filtered           = nf,                           # outside any threshold
    pct_filtered         = round(100 * nf / n, 2),
    n_pass               = n - nf,
    filtered_label       = paste0(fmt(nf), " of ", fmt(n), " cells filtered")
  )
  list(status = status, row = row)
}

## The four vignette plots for one group, saved into out_dir.
##   md = metadata rows; status = pass/filtered; r = summary row
nc_nf_plots <- function(md, status, r, out_dir, tag) {
  meta <- data.frame(nCount_RNA = md$nCount_RNA, nFeature_RNA = md$nFeature_RNA,
                     x_slide_mm = md$x_slide_mm, y_slide_mm = md$y_slide_mm,
                     counts_filter = status)
  meta <- meta[order(meta$counts_filter == "filtered"), ]   # draw filtered cells last, on top
  thr_txt <- paste0("counts ", r$counts_min_threshold, "-", fmt(r$counts_max_threshold),
                    "; features ", r$features_min_threshold, "-", fmt(r$features_max_threshold))
  ttl <- paste0(tag, ": ", r$filtered_label, " (", r$pct_filtered, "%)")
  
  ## Joint distribution (vignette)
  p <- ggplot(data = meta, aes(x = nFeature_RNA, y = nCount_RNA, color = counts_filter)) +
    geom_point(alpha = 0.3, size = 0.3, pch = 16) +          # vignette: alpha 0.002, pch "." (too faint)
    scale_color_manual(values = qc_cols) +
    geom_vline(xintercept = c(r$features_min_threshold, r$features_max_threshold), lty = "dotted") +
    geom_hline(yintercept = c(r$counts_min_threshold, r$counts_max_threshold), lty = "dotted") +
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = trans_format("log10", math_format(10^.x))) +
    theme_minimal() + coord_fixed() + theme(legend.position = "none") +
    labs(title = ttl, subtitle = thr_txt)
  p <- ggMarginal(p, type = "density")
  ggsave(file.path(out_dir, paste0("01_qc_joint_distribution_", tag, ".png")), p, width = 8, height = 8)
  
  ## Transcripts per cell (vignette)
  p1 <- ggplot(data = meta, aes(x = nCount_RNA)) + geom_density(alpha = 0.01) + theme_minimal() +
    geom_vline(xintercept = c(r$counts_min_threshold, r$counts_max_threshold), lty = "dotted") +
    labs(title = ttl, subtitle = paste0("nCount_RNA thresholds: ", r$counts_min_threshold, " - ",
                                        fmt(r$counts_max_threshold), " | below: ", fmt(r$n_below_counts_min),
                                        ", above: ", fmt(r$n_above_counts_max)))
  ggsave(file.path(out_dir, paste0("02_qc_dist_counts_", tag, ".png")), p1, width = 8, height = 8)
  
  ## Features per cell (vignette)
  p2 <- ggplot(data = meta, aes(x = nFeature_RNA)) + geom_density(alpha = 0.01) + theme_minimal() +
    geom_vline(xintercept = c(r$features_min_threshold, r$features_max_threshold), lty = "dotted") +
    labs(title = ttl, subtitle = paste0("nFeature_RNA thresholds: ", r$features_min_threshold, " - ",
                                        fmt(r$features_max_threshold), " | below: ", fmt(r$n_below_features_min),
                                        ", above: ", fmt(r$n_above_features_max)))
  ggsave(file.path(out_dir, paste0("03_qc_dist_features_", tag, ".png")), p2, width = 8, height = 8)
  
  ## Filtered cells in space (vignette)
  p3 <- ggplot(data = meta) +
    geom_point(aes(x = x_slide_mm, y = y_slide_mm, color = counts_filter), size = 0.05, alpha = 0.4) +  # vignette: 0.001 / 0.02
    coord_fixed() + theme_bw() +
    scale_color_manual(values = qc_cols) +
    facet_wrap(~counts_filter) +
    guides(color = guide_legend(override.aes = list(size = 7, alpha = 1))) +
    labs(title = ttl)
  ggsave(file.path(out_dir, paste0("04_qc_xy_flagged_", tag, ".png")), p3, width = 8, height = 8)
  
  message("Saved 4 plots for ", tag, " in ", out_dir)
}


## =============================================================================
## CHUNK 36: APPROACH A - all cells, one set of thresholds
## =============================================================================
res_all <- nc_nf_qc(metadata, "All cells")
nc_nf_all_summary <- res_all$row
print(nc_nf_all_summary)

nc_nf_plots(metadata, res_all$status, nc_nf_all_summary, dir_nc_all, "all_cells")

fwrite(nc_nf_all_summary, file.path(dir_nc_all, "ncount_nfeature_all_cells_summary.csv"))
filter_nc_nf_all <- setNames(res_all$status == "filtered", metadata$cell_id)   # TRUE = filtered
saveRDS(filter_nc_nf_all, file.path(dir_nc_all, "flag_ncount_nfeature_all_cells.RDS"))


## =============================================================================
## CHUNK 37: APPROACH B - each tissue section, its own thresholds
## =============================================================================
status_sec <- setNames(character(nrow(metadata)), metadata$cell_id)
rows <- list()
for (sec in section_names) {
  idx <- metadata$section_id == sec
  res <- nc_nf_qc(metadata[idx, ], sec)
  status_sec[idx] <- res$status
  rows[[sec]] <- res$row
  d <- file.path(dir_nc_sec, sec); dir.create(d, showWarnings = FALSE)
  nc_nf_plots(metadata[idx, ], res$status, res$row, d, sec)            # 4 plots per section
}
nc_nf_sec_summary <- rbindlist(rows)

## Total row (sum over sections; thresholds differ per section so left NA)
tot <- nc_nf_sec_summary[, .(n_cells = sum(n_cells),
                             min_nCount_RNA = min(min_nCount_RNA), max_nCount_RNA = max(max_nCount_RNA),
                             min_nFeature_RNA = min(min_nFeature_RNA), max_nFeature_RNA = max(max_nFeature_RNA),
                             n_below_counts_min = sum(n_below_counts_min), n_above_counts_max = sum(n_above_counts_max),
                             n_below_features_min = sum(n_below_features_min), n_above_features_max = sum(n_above_features_max),
                             n_filtered = sum(n_filtered))]
tot[, `:=`(group = "All sections (sum)", pct_filtered = round(100 * n_filtered / n_cells, 2),
           n_pass = n_cells - n_filtered,
           filtered_label = paste0(fmt(n_filtered), " of ", fmt(n_cells), " cells filtered"))]
nc_nf_sec_summary <- rbind(nc_nf_sec_summary, tot, use.names = TRUE, fill = TRUE)
print(nc_nf_sec_summary)

## Grid: joint distribution per section with its own threshold lines
md_plot <- data.frame(section_id = metadata$section_id, nCount_RNA = metadata$nCount_RNA,
                      nFeature_RNA = metadata$nFeature_RNA, counts_filter = status_sec)
md_plot <- md_plot[order(md_plot$counts_filter == "filtered"), ]   # filtered drawn on top
thr_lines <- nc_nf_sec_summary[group %in% section_names]
thr_lines[, section_id := paste0(group, ": ", filtered_label)]
md_plot$section_id <- thr_lines$section_id[match(md_plot$section_id, thr_lines$group)]  # panel titles with counts
p_grid <- ggplot(md_plot, aes(nFeature_RNA, nCount_RNA, color = counts_filter)) +
  geom_point(alpha = 0.3, size = 0.3, pch = 16) +
  scale_color_manual(values = qc_cols) +
  geom_vline(data = thr_lines, aes(xintercept = features_min_threshold), lty = "dotted") +
  geom_vline(data = thr_lines, aes(xintercept = features_max_threshold), lty = "dotted") +
  geom_hline(data = thr_lines, aes(yintercept = counts_min_threshold),   lty = "dotted") +
  geom_hline(data = thr_lines, aes(yintercept = counts_max_threshold),   lty = "dotted") +
  scale_x_log10(labels = trans_format("log10", math_format(10^.x))) +
  scale_y_log10(labels = trans_format("log10", math_format(10^.x))) +
  facet_wrap(~section_id, nrow = 2) + theme_minimal() + theme(legend.position = "none") +
  labs(title = paste0("Per-section thresholds: ", tot$filtered_label, " (", tot$pct_filtered, "%)"))
ggsave(file.path(dir_nc_sec, "01_qc_joint_distribution_per_section_grid.png"), p_grid, width = 12, height = 8)

fwrite(nc_nf_sec_summary, file.path(dir_nc_sec, "ncount_nfeature_per_section_summary.csv"))
filter_nc_nf_sec <- setNames(status_sec == "filtered", metadata$cell_id)   # TRUE = filtered
saveRDS(filter_nc_nf_sec, file.path(dir_nc_sec, "flag_ncount_nfeature_per_section.RDS"))


## =============================================================================
## CHUNK 38: results on screen
## =============================================================================
show_cols <- c("group", "n_cells", "min_nCount_RNA", "max_nCount_RNA", "min_nFeature_RNA", "max_nFeature_RNA",
               "counts_min_threshold", "counts_max_threshold", "features_min_threshold", "features_max_threshold",
               "n_filtered", "pct_filtered")
print(rbind(nc_nf_all_summary[, ..show_cols], nc_nf_sec_summary[, ..show_cols]))

message("A, all cells:          ", nc_nf_all_summary$filtered_label, " (", nc_nf_all_summary$pct_filtered, "%)")
message("B, per tissue section: ", tot$filtered_label, " (", tot$pct_filtered, "%)")
message("Outputs saved to ", dir_nc)











## =============================================================================
## 07_QC_FOV_borders.R
## =============================================================================
## PURPOSE
##   Vignette 2, section 3.2.2 (Optional) Cells near FOV borders:
##     SplitRatioToLocal = 0      -> cell not on an FOV border (kept)
##     SplitRatioToLocal > 0      -> border cell; value = its area / mean cell area in the FOV
##     cutoff 0.5 (log2 = -1)     -> border cells smaller than half the average cell are removed
##   remove = SplitRatioToLocal > 0 & SplitRatioToLocal < 0.5
##
##   The ratio is already relative to each FOV (AtoMx computes it), and the
##   cutoff is fixed, so the SAME cells are removed whether you look at all
##   cells or at each section. What differs is only how the counts are tallied:
##     all_cells/            : whole slide, and per FOV
##     per_tissue_section/   : per section, and per FOV within each section
##   Not combined with the other QC flags here.
##
## Run after script 02 (needs metadata with section_id). Chunk by chunk.
## =============================================================================


## =============================================================================
## CHUNK 39: packages, settings, folders
## =============================================================================
## Needed once: install.packages(c("patchwork", "dplyr", "scales"))
library(data.table)   # tables
library(dplyr)        # %>%, filter, mutate, count (vignette)
library(ggplot2)      # plots
library(patchwork)    # inset_element(): pie chart inside the histogram (vignette)
library(scales)

output_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"
# metadata <- readRDS(file.path(output_dir, "metadata_unfiltered.RDS"))   # fresh session only

stopifnot("SplitRatioToLocal" %in% colnames(metadata))   # AtoMx column used by the vignette
stopifnot(!anyNA(metadata$SplitRatioToLocal))
section_names <- sort(unique(metadata$section_id))

SplitRatioToLocalThreshold <- 0.5                 # vignette cutoff (log2 = -1)
fb_cols <- c("removed" = "#E8590C", "kept" = "grey80")   # plot colours only

dir_fb     <- file.path(output_dir, "QC", "FOV_Borders")
dir_fb_all <- file.path(dir_fb, "all_cells")
dir_fb_sec <- file.path(dir_fb, "per_tissue_section")
dir.create(dir_fb_all, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_fb_sec, showWarnings = FALSE, recursive = TRUE)

## Per-cell flags
meta <- as.data.frame(metadata)
meta$border_cell    <- meta$SplitRatioToLocal > 0                                   # touches FOV edge
meta$border_removed <- meta$border_cell & meta$SplitRatioToLocal < SplitRatioToLocalThreshold
meta$border_status  <- ifelse(meta$border_removed, "removed", "kept")

fmt <- function(n) format(n, big.mark = ",")


## =============================================================================
## CHUNK 40: helper functions
## =============================================================================
## Summary for any grouping (grp = NULL whole slide, "section_id", or c("section_id","fov"))
fb_summary <- function(d, grp = NULL) {
  as.data.table(d)[, .(
    n_cells             = .N,
    n_border_cells      = sum(border_cell),                      # SplitRatioToLocal > 0
    pct_border_cells    = round(100 * mean(border_cell), 2),
    n_removed           = sum(border_removed),                   # border & ratio < 0.5
    pct_removed         = round(100 * mean(border_removed), 2),
    n_kept              = .N - sum(border_removed),
    removed_label       = paste0(fmt(sum(border_removed)), " of ", fmt(.N), " cells removed (border, ratio < 0.5)")
  ), by = grp]
}

## Vignette figure: histogram of log2(SplitRatioToLocal) for border cells + pie inset
fb_vignette_plot <- function(d, title) {
  p_detail <- ggplot(
    data = d %>% filter(SplitRatioToLocal != 0),
    aes(x = log2(SplitRatioToLocal))) +
    geom_histogram(bins = 100, fill = "steelblue") +
    labs(x = expression(log[2]("SplitRatioToLocal")),
         y = "Number of FOV border cells",
         title = title,
         subtitle = paste0(fmt(sum(d$border_removed)), " of ", fmt(sum(d$border_cell)),
                           " border cells removed (", fmt(nrow(d)), " cells total)")) +
    theme_minimal() +
    geom_vline(xintercept = log2(SplitRatioToLocalThreshold), color = "darkorange", lty = "dashed")
  
  pie_data <- d %>%
    mutate(Category = if_else(SplitRatioToLocal > 0, "Border cell", "Non-border cell")) %>%
    count(Category) %>%
    mutate(Percentage = (n / sum(n)) * 100,
           Label = paste0(Category, " (", round(Percentage, 1), "%)"))
  
  p_pie <- ggplot(pie_data, aes(x = 2, y = n, fill = Category)) +
    geom_bar(stat = "identity", width = 1, color = "white") +
    coord_polar("y", start = 0) +
    scale_fill_manual(values = c("Border cell" = "steelblue", "Non-border cell" = "gray90")) +
    geom_text(aes(label = Label), position = position_stack(vjust = 0.5), size = 3, fontface = "bold") +
    theme_void() + xlim(0.5, 2.5) + theme(legend.position = "none")
  
  p_detail + inset_element(p_pie, left = 0.6, bottom = 0.6, right = 1.0, top = 1.0)
}

## Spatial plot: removed border cells on top, per-FOV removed count at each FOV centre
fb_spatial_plot <- function(d, fov_tab, title, subtitle, label_size = 2) {
  d <- d[order(d$border_status == "removed"), ]            # removed drawn last
  ggplot(d, aes(x_slide_mm, y_slide_mm, colour = border_status)) +
    geom_point(size = 0.05, alpha = 0.5) +
    scale_colour_manual(values = fb_cols) +
    geom_text(data = fov_tab, aes(fov_x, fov_y, label = n_removed), inherit.aes = FALSE,
              size = label_size, colour = "black") +
    coord_fixed() + theme_minimal() +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    labs(title = title, subtitle = subtitle, colour = NULL,
         caption = "Numbers = border cells removed in that FOV", x = "x (mm)", y = "y (mm)")
}

## FOV centres (for labels)
fov_centres <- as.data.table(meta)[, .(fov_x = mean(range(x_slide_mm)),
                                       fov_y = mean(range(y_slide_mm))), by = .(section_id, fov)]


## =============================================================================
## CHUNK 41: ALL CELLS - whole slide and per FOV
## =============================================================================
fb_all_summary <- fb_summary(meta)[, group := "All cells"]
setcolorder(fb_all_summary, "group")
fb_all_per_fov <- fb_summary(meta, grp = c("section_id", "fov"))[order(fov)]
fb_all_per_fov <- merge(fb_all_per_fov, fov_centres, by = c("section_id", "fov"))

print(fb_all_summary)
print(fb_all_per_fov[n_removed > 0, .(fov, section_id, n_cells, n_border_cells, n_removed, pct_removed)])

## Vignette histogram + pie
p <- fb_vignette_plot(meta, "All cells")
ggsave(file.path(dir_fb_all, "01_splitratiotolocal_all_cells.png"), p, width = 6, height = 6)
print(p)

## Spatial: whole slide
p_sp <- fb_spatial_plot(meta, fb_all_per_fov, "All cells - FOV border cells removed",
                        paste0(fb_all_summary$removed_label, " (", fb_all_summary$pct_removed, "%)"),
                        label_size = 1.8)
ggsave(file.path(dir_fb_all, "02_spatial_border_cells_removed_all_cells.png"), p_sp, width = 12, height = 12, dpi = 300)
print(p_sp)

## Save
fwrite(fb_all_summary, file.path(dir_fb_all, "FOV_borders_all_cells_summary.csv"))
fwrite(fb_all_per_fov[, !c("fov_x", "fov_y")], file.path(dir_fb_all, "FOV_borders_all_cells_per_FOV.csv"))
flag_border_removed <- setNames(meta$border_removed, meta$cell_id)          # TRUE = removed
saveRDS(flag_border_removed, file.path(dir_fb_all, "flag_FOV_border_cells.RDS"))


## =============================================================================
## CHUNK 42: PER TISSUE SECTION - per section and per FOV within each section
## =============================================================================
fb_sec_summary <- fb_summary(meta, grp = "section_id")[order(section_id)]
fb_sec_summary <- rbind(fb_sec_summary, fb_summary(meta)[, section_id := "All sections"], use.names = TRUE)
print(fb_sec_summary)
fwrite(fb_sec_summary, file.path(dir_fb_sec, "FOV_borders_per_section_summary.csv"))

for (sec in section_names) {
  d_sec   <- meta[meta$section_id == sec, ]
  fov_sec <- fb_all_per_fov[section_id == sec][order(fov)]     # same cells, same counts
  r       <- fb_sec_summary[section_id == sec]
  d_out   <- file.path(dir_fb_sec, sec); dir.create(d_out, showWarnings = FALSE)
  
  ## Per-FOV table for this section
  fwrite(fov_sec[, !c("fov_x", "fov_y")], file.path(d_out, paste0("FOV_borders_", sec, "_per_FOV.csv")))
  
  ## Vignette histogram + pie
  ggsave(file.path(d_out, paste0("01_splitratiotolocal_", sec, ".png")),
         fb_vignette_plot(d_sec, sec), width = 6, height = 6)
  
  ## Spatial with per-FOV counts
  ggsave(file.path(d_out, paste0("02_spatial_border_cells_removed_", sec, ".png")),
         fb_spatial_plot(d_sec, fov_sec, paste0(sec, " - FOV border cells removed"),
                         paste0(r$removed_label, " (", r$pct_removed, "%)"), label_size = 3),
         width = 9, height = 9, dpi = 300)
  
  message(sec, ": ", r$removed_label, " (", r$pct_removed, "%)")
}

## Grid: all 6 sections, removed border cells on top
d_grid <- meta
lab <- setNames(paste0(fb_sec_summary$section_id, ": ", fmt(fb_sec_summary$n_removed), " of ",
                       fmt(fb_sec_summary$n_cells), " removed"), fb_sec_summary$section_id)
d_grid$panel <- lab[d_grid$section_id]
d_grid <- d_grid[order(d_grid$border_status == "removed"), ]
tot <- fb_sec_summary[section_id == "All sections"]
p_grid <- ggplot(d_grid, aes(x_slide_mm, y_slide_mm, colour = border_status)) +
  geom_point(size = 0.05, alpha = 0.5) +
  scale_colour_manual(values = fb_cols) +
  facet_wrap(~panel, nrow = 2, scales = "free") +
  theme_minimal() +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  labs(title = paste0("FOV border cells removed - all sections: ", tot$removed_label,
                      " (", tot$pct_removed, "%)"), colour = NULL, x = "x (mm)", y = "y (mm)")
ggsave(file.path(dir_fb_sec, "02_spatial_border_cells_removed_per_section_grid.png"), p_grid,
       width = 14, height = 9, dpi = 300)
print(p_grid)


## =============================================================================
## CHUNK 43: results on screen
## =============================================================================
message("All cells: ", fb_all_summary$removed_label, " (", fb_all_summary$pct_removed, "%); ",
        fmt(fb_all_summary$n_border_cells), " border cells in total")
for (sec in section_names) {
  r <- fb_sec_summary[section_id == sec]
  message(sec, ": ", r$removed_label, " (", r$pct_removed, "%) across ",
          nrow(fb_all_per_fov[section_id == sec & n_removed > 0]), " of ",
          nrow(fb_all_per_fov[section_id == sec]), " FOVs")
}
message("Outputs saved to ", dir_fb)


















## =============================================================================
## 08_QC_SBR.R
## =============================================================================
## PURPOSE
##   Vignette 2, section 3.2.3: spatially smoothed Signal to Background Ratio (SBR)
##     S   = mean target counts per cell   (nCount_RNA / number of targets)
##     B   = mean negative-probe counts    (nCount_negprobes / number of negprobes)
##     both smoothed over space with liana's Gaussian kernel, then
##     SBR = smoothed S / (smoothed B + 1e-9)
##   Vignette decision: low-SBR cells (log2(SBR) < 0) that are spatially
##   concentrated are flagged for removal.
##
##   Reproduces li.ut.spatial_neighbors(bandwidth = 0.010, cutoff = 0.080,
##   kernel = 'gaussian', set_diag = TRUE, standardize = TRUE) in R
##   (liana source, v1.10.0):
##     - 100 nearest neighbours (liana default max_neighbours = 100) + the cell itself
##     - weight = exp(-d^2 / (2 * bandwidth^2)); coordinates in mm (as in the vignette),
##       so bandwidth 0.010 mm = 10 µm
##     - cutoff applies to the WEIGHT: weights <= 0.08 set to 0
##       (0.08 corresponds to a distance of ~22.5 µm)
##     - each row scaled to sum to 1
##
##   Smoothing the mean is exact: mean over genes of (W %*% X) = W %*% (mean over
##   genes of X), so no smoothed 123k x 6k matrix is needed.
##
##   Thresholds do not depend on the panel: S and B are per-probe MEANS, so the
##   6k vs WTX difference in number of targets is divided out. log2(SBR) < 0
##   means target signal is below background.
##
##   Neighbours are within ~22.5 µm and sections lie > 1 mm apart, so each cell's
##   SBR is the same whether computed on all cells or per section; per section
##   only changes the tallies and plots.
##
## Run after script 02 (needs metadata with section_id, negmean; counts;
## negprobes_counts). Chunk by chunk.
## =============================================================================

## =============================================================================
## 08_QC_SBR.R
## =============================================================================
## PURPOSE
##   Vignette 2, section 3.2.3: spatially smoothed Signal to Background Ratio (SBR)
##     S   = mean target counts per cell   (nCount_RNA / number of targets)
##     B   = mean negative-probe counts    (nCount_negprobes / number of negprobes)
##     both smoothed over space with liana's Gaussian kernel, then
##     SBR = smoothed S / (smoothed B + 1e-9)
##   Vignette decision: low-SBR cells (log2(SBR) < 0) that are spatially
##   concentrated are flagged for removal.
##
##   Reproduces li.ut.spatial_neighbors(bandwidth = 0.010, cutoff = 0.080,
##   kernel = 'gaussian', set_diag = TRUE, standardize = TRUE) in R
##   (liana source, v1.10.0):
##     - 100 nearest neighbours (liana default max_neighbours = 100) + the cell itself
##     - weight = exp(-d^2 / (2 * bandwidth^2)); coordinates in mm (as in the vignette),
##       so bandwidth 0.010 mm = 10 µm
##     - cutoff applies to the WEIGHT: weights <= 0.08 set to 0
##       (0.08 corresponds to a distance of ~22.5 µm)
##     - each row scaled to sum to 1
##
##   Smoothing the mean is exact: mean over genes of (W %*% X) = W %*% (mean over
##   genes of X), so no smoothed 123k x 6k matrix is needed.
##
##   Thresholds do not depend on the panel: S and B are per-probe MEANS, so the
##   6k vs WTX difference in number of targets is divided out. log2(SBR) < 0
##   means target signal is below background.
##
##   Neighbours are within ~22.5 µm and sections lie > 1 mm apart, so each cell's
##   SBR is the same whether computed on all cells or per section; per section
##   only changes the tallies and plots.
##
## Run after script 02 (needs metadata with section_id, negmean; counts;
## negprobes_counts). Chunk by chunk.
## =============================================================================


## =============================================================================
## CHUNK 44: packages, settings, folders
## =============================================================================
library(data.table)
library(Matrix)     # sparse weight matrix
library(dbscan)     # kNN() nearest neighbours
library(ggplot2)

output_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"
## Fresh session only:
# metadata         <- readRDS(file.path(output_dir, "metadata_unfiltered.RDS"))
# counts           <- readRDS(file.path(output_dir, "counts_unfiltered.RDS"))
# negprobes_counts <- readRDS(file.path(output_dir, "negprobes_unfiltered.RDS"))

section_names <- sort(unique(metadata$section_id))

## Vignette parameters
bandwidth      <- 0.010   # mm (10 µm)
weight_cutoff  <- 0.080   # kernel weight below which a neighbour is ignored
max_neighbours <- 100     # liana default
epsilon        <- 1e-9
sbr_threshold  <- 0       # log2(SBR) < 0 flagged (vignette)

dir_sbr     <- file.path(output_dir, "QC", "SBR")
dir_sbr_all <- file.path(dir_sbr, "all_cells")
dir_sbr_sec <- file.path(dir_sbr, "per_tissue_section")
dir.create(dir_sbr_all, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_sbr_sec, showWarnings = FALSE, recursive = TRUE)

fmt <- function(n) format(n, big.mark = ",")


## =============================================================================
## CHUNK 45: per-cell mean signal and background (not yet smoothed)
## =============================================================================
n_targets   <- ncol(counts)             # targets in the panel (vignette: columns of adata.X)
n_negprobes <- ncol(negprobes_counts)   # negative probes (vignette: adata.obsm['counts_neg'])
n_targets; n_negprobes

mean_target <- metadata$nCount_RNA       / n_targets     # = rowMeans(counts)
mean_neg    <- metadata$nCount_negprobes / n_negprobes   # = rowMeans(negprobes_counts) = negmean


## =============================================================================
## CHUNK 46: Gaussian connectivity matrix (liana spatial_neighbors, in R)
## =============================================================================
xy_mm <- as.matrix(metadata[, c("x_slide_mm", "y_slide_mm")])
nn    <- dbscan::kNN(xy_mm, k = max_neighbours)            # 100 nearest, self excluded

n <- nrow(xy_mm)
i <- c(seq_len(n), rep(seq_len(n), max_neighbours))       # rows: self, then neighbours
j <- c(seq_len(n), as.vector(nn$id))                       # columns
d <- c(rep(0, n),  as.vector(nn$dist))                     # distance 0 to self (set_diag = TRUE)
w <- exp(-(d^2) / (2 * bandwidth^2))                       # Gaussian kernel
keep <- w > weight_cutoff                                  # cutoff on weight
W <- sparseMatrix(i = i[keep], j = j[keep], x = w[keep], dims = c(n, n))
W <- Diagonal(x = 1 / rowSums(W)) %*% W                    # standardize: rows sum to 1

summary(rowSums(W > 0))    # neighbours used per cell (incl. itself)


## =============================================================================
## CHUNK 47: smoothed SBR
## =============================================================================
smoothed_mean_target <- as.vector(W %*% mean_target)
smoothed_mean_neg    <- as.vector(W %*% mean_neg)
metadata$SBR      <- smoothed_mean_target / (smoothed_mean_neg + epsilon)
metadata$log2_SBR <- log2(metadata$SBR)
metadata$low_SBR  <- metadata$log2_SBR < sbr_threshold     # TRUE = flagged

summary(metadata$log2_SBR)


## =============================================================================
## CHUNK 48: helpers - bins (vignette), summary row, spatial plot (vignette)
## =============================================================================
sbr_breaks <- c(-Inf, -1, 0, 1, 2, 3, Inf)
sbr_bins   <- c("<= -1", "-1 to 0", "0 to 1", "1 to 2", "2 to 3", "> 3")
bin_colors <- c("<= -1" = "#D73027", "-1 to 0" = "#F46D43", "0 to 1" = "#D9EF8B",
                "1 to 2" = "#A6D96A", "2 to 3" = "#1A9850", "> 3" = "#006837")   # vignette

## include.lowest so a cell with SBR = 0 (log2 = -Inf) falls in "<= -1"
sbr_cut <- function(x) cut(x, breaks = sbr_breaks, labels = sbr_bins, include.lowest = TRUE)

sbr_row <- function(md, label) {
  b <- table(sbr_cut(md$log2_SBR))
  r <- data.table(group = label, n_cells = nrow(md),
                  median_log2_SBR = round(median(md$log2_SBR), 2),
                  n_low_SBR = sum(md$low_SBR),
                  pct_low_SBR = round(100 * mean(md$low_SBR), 2))
  r[, paste0("n_", sbr_bins) := as.list(as.integer(b))]
  r[, low_SBR_label := paste0(fmt(n_low_SBR), " of ", fmt(n_cells), " cells log2(SBR) < 0")]
  r
}

## Vignette plot; legend labels carry each bin's % of cells in this group
sbr_spatial <- function(md, title, subtitle) {
  pct  <- round(100 * table(sbr_cut(md$log2_SBR)) / nrow(md), 1)
  labs_bin <- paste0(sbr_bins, " (", pct, "%)")
  cols <- setNames(bin_colors, labs_bin)
  md$bin <- factor(labs_bin[as.integer(sbr_cut(md$log2_SBR))], levels = labs_bin)
  ggplot(md) +
    geom_point(aes(x = x_slide_mm, y = y_slide_mm, color = bin), size = 0.001, alpha = 1) +
    scale_color_manual(values = cols, drop = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 7, alpha = 1))) +
    labs(color = "log2(SBR)", title = title, subtitle = subtitle) +
    coord_fixed() + theme_bw()
}


## =============================================================================
## CHUNK 49: ALL CELLS
## =============================================================================
sbr_all_summary <- sbr_row(metadata, "All cells")
print(sbr_all_summary)

p <- sbr_spatial(metadata, "Spatially smoothed SBR - all cells",
                 paste0(sbr_all_summary$low_SBR_label, " (", sbr_all_summary$pct_low_SBR, "%)"))
ggsave(file.path(dir_sbr_all, "01_qc_xy_sbr_all_cells.png"), p, width = 10, height = 10)

## log2(SBR) distribution with the threshold
p_h <- ggplot(metadata, aes(log2_SBR)) + geom_histogram(bins = 100, fill = "steelblue") +
  geom_vline(xintercept = sbr_threshold, colour = "#D73027", lty = "dashed") + theme_minimal() +
  labs(x = "log2(SBR)", y = "Number of cells", title = "All cells",
       subtitle = sbr_all_summary$low_SBR_label)
ggsave(file.path(dir_sbr_all, "02_log2_SBR_histogram_all_cells.png"), p_h, width = 7, height = 5)

fwrite(sbr_all_summary, file.path(dir_sbr_all, "SBR_all_cells_summary.csv"))
fwrite(data.table(cell_id = metadata$cell_id, section_id = metadata$section_id, fov = metadata$fov,
                  SBR = metadata$SBR, log2_SBR = metadata$log2_SBR, low_SBR = metadata$low_SBR),
       file.path(dir_sbr_all, "SBR_per_cell.csv"))
saveRDS(setNames(metadata$low_SBR, metadata$cell_id), file.path(dir_sbr_all, "flag_low_SBR.RDS"))


## =============================================================================
## CHUNK 50: PER TISSUE SECTION
## =============================================================================
sbr_sec_summary <- rbindlist(lapply(section_names, function(sec)
  sbr_row(metadata[metadata$section_id == sec, ], sec)))
sbr_sec_summary <- rbind(sbr_sec_summary, sbr_row(metadata, "All sections"))
print(sbr_sec_summary)
fwrite(sbr_sec_summary, file.path(dir_sbr_sec, "SBR_per_section_summary.csv"))

for (sec in section_names) {
  md <- metadata[metadata$section_id == sec, ]
  r  <- sbr_sec_summary[group == sec]
  ggsave(file.path(dir_sbr_sec, paste0("01_qc_xy_sbr_", sec, ".png")),
         sbr_spatial(md, paste0("Spatially smoothed SBR - ", sec),
                     paste0(r$low_SBR_label, " (", r$pct_low_SBR, "%)")),
         width = 9, height = 9)
}

## Grid: all 6 sections, one colour scale, panel titles with counts
md_g <- metadata
md_g$bin <- sbr_cut(md_g$log2_SBR)
lab <- setNames(paste0(sbr_sec_summary$group, ": ", sbr_sec_summary$low_SBR_label), sbr_sec_summary$group)
md_g$panel <- lab[md_g$section_id]
tot <- sbr_sec_summary[group == "All sections"]
p_grid <- ggplot(md_g) +
  geom_point(aes(x_slide_mm, y_slide_mm, color = bin), size = 0.001) +
  scale_color_manual(values = bin_colors, drop = FALSE) +
  facet_wrap(~panel, nrow = 2, scales = "free") + theme_bw() +
  guides(color = guide_legend(override.aes = list(size = 7, alpha = 1))) +
  labs(color = "log2(SBR)", title = paste0("Spatially smoothed SBR - ", tot$low_SBR_label,
                                           " (", tot$pct_low_SBR, "%)"))
ggsave(file.path(dir_sbr_sec, "02_qc_xy_sbr_per_section_grid.png"), p_grid, width = 14, height = 9)

## Whole slide, cells in their original positions, each section labelled
## above it with its own count of low-SBR cells.
sec_lab <- as.data.table(metadata)[, .(x = mean(range(x_slide_mm)), y = max(y_slide_mm)), by = section_id]
sec_lab <- merge(sec_lab, sbr_sec_summary[, .(section_id = group, n_low_SBR, n_cells, pct_low_SBR)],
                 by = "section_id")
sec_lab[, label := paste0(section_id, ": ", fmt(n_low_SBR), " of ", fmt(n_cells),
                          " log2(SBR) < 0 (", pct_low_SBR, "%)")]

p_slide <- sbr_spatial(metadata, "Spatially smoothed SBR - per tissue section",
                       paste0("All sections: ", tot$low_SBR_label, " (", tot$pct_low_SBR, "%)")) +
  geom_label(data = sec_lab, aes(x = x, y = y, label = label), vjust = -0.3,
             size = 3, fontface = "bold", fill = "white", label.size = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.08)))   # room for the top labels
ggsave(file.path(dir_sbr_sec, "03_qc_xy_sbr_per_section_on_slide.png"), p_slide, width = 12, height = 12)
print(p_slide)


## =============================================================================
## CHUNK 51: results on screen
## =============================================================================
message("All cells: ", sbr_all_summary$low_SBR_label, " (", sbr_all_summary$pct_low_SBR, "%)")
for (sec in section_names) {
  r <- sbr_sec_summary[group == sec]
  message(sec, ": ", r$low_SBR_label, " (", r$pct_low_SBR, "%)")
}
message("Outputs saved to ", dir_sbr)











































## =============================================================================
## 09_QC_Proseg.R
## =============================================================================



## =============================================================================
## 09_QC_Proseg.R
## =============================================================================
## PURPOSE
##   Segmentation QC with ProSeg, run SEPARATELY on each of the 6 tissue
##   sections (six ProSeg runs). Continues from scripts 02-08.
##   Results are saved as files only (QC/QC_Proseg/); metadata is NOT changed.
##
## WHAT PROSEG DOES
##   AtoMx assigned each transcript to a cell from image-based outlines.
##   ProSeg re-segments the same transcripts from their positions and gene
##   identities (probabilistic segmentation), starting from AtoMx's assignment.
##   If AtoMx segmented a cell well, ProSeg keeps its transcripts together.
##
## QC METRIC: purity, per AtoMx cell
##   purity = transcripts in the ProSeg cell that received most of them
##            / all of the cell's transcripts that ProSeg assigned to a cell
##   flagged       = purity < 0.5 (fewer than half the transcripts stay together)
##   not evaluated = < 3 evaluable transcripts (purity NA; NOT flagged)
##
## ProSeg options checked against the ProSeg source (v3.2.0). Transcript output
## columns used: transcript_id, observed_x, observed_y, gene, assignment, background.
##
## OUTPUT: QC/QC_Proseg/
##   01_prepared_transcripts/   ProSeg input, one file per section
##   02_proseg_raw_output/<section>/   ProSeg's own output + run log
##   03_results/                CSVs and plots (cells and transcripts)
##
## REQUIREMENT: ProSeg installed (Terminal: cargo install proseg).
## Run in a fresh R session (Session > Restart R), chunk by chunk.
## =============================================================================


## =============================================================================
## CHUNK 1: packages, paths, settings
## =============================================================================
suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

output_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"
slide_dir  <- "/Users/maryamuhabilausman/Desktop/AAA_CosMx_Analysis/Downloaded_data/flatFiles/AC20260629_AAA3"

## Section files written by the earlier version of this script are reused
dir_ps  <- file.path(output_dir, "QC", "QC_Proseg", "per_tissue_section")
dir_in  <- file.path(dir_ps, "01_prepared_transcripts")
dir_raw <- file.path(dir_ps, "02_proseg_raw_output")
dir_res <- file.path(dir_ps, "03_results")
dir_sec_plots <- file.path(dir_res, "per_section_plots")
dir_tx_out    <- file.path(dir_res, "flagged_cell_transcripts")
for (d in c(dir_in, dir_raw, dir_sec_plots, dir_tx_out)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

pixel_size_um        <- 0.12028   # µm per CosMx pixel (ProSeg's CosMx value)
purity_threshold     <- 0.5       # flag cells with purity below this
min_transcripts_eval <- 3         # fewer evaluable transcripts -> not evaluated
coord_tolerance_um   <- 0.05      # allowed rounding difference in the linking check
excluded_genes_regex <- "^(FalseCode|SystemControl|NegPrb|Negative)"   # control probes (ProSeg's CosMx default)
max_transcript_nucleus_distance <- NULL   # NULL = ProSeg default (60 µm)
n_threads    <- NULL              # NULL = all CPU cores; e.g. 6 keeps the laptop usable
rerun_proseg <- FALSE             # FALSE = reuse ProSeg output that already exists

cargo_bin <- file.path(Sys.getenv("HOME"), ".cargo", "bin")
if (dir.exists(cargo_bin) && !grepl(cargo_bin, Sys.getenv("PATH"), fixed = TRUE))
  Sys.setenv(PATH = paste(cargo_bin, Sys.getenv("PATH"), sep = ":"))
stopifnot("proseg not found - install with: cargo install proseg" = nzchar(Sys.which("proseg")))

fmt <- function(n) format(n, big.mark = ",")


## =============================================================================
## CHUNK 2: cell metadata (read only, not modified)
## =============================================================================
metadata <- readRDS(file.path(output_dir, "metadata_unfiltered.RDS"))
stopifnot(identical(rownames(metadata), metadata$cell_id))
section_names <- sort(unique(metadata$section_id))

fov_to_section <- unique(data.table(fov = metadata$fov, section_id = metadata$section_id))
stopifnot("an FOV lies in two sections" = !anyDuplicated(fov_to_section$fov))

## Cell names are "c_<slide>_<FOV>_<cell>"; the prefix "c_<slide>" is rebuilt per transcript
cell_prefix <- unique(sub("_[^_]+_[^_]+$", "", metadata$cell_id))
stopifnot("cell_id prefix not unique" = length(cell_prefix) == 1)


## =============================================================================
## CHUNK 3: prepare one ProSeg input file per section
##          (skipped if all 6 files already exist)
## =============================================================================
in_file  <- function(sec) file.path(dir_in, paste0(sec, "_tx_prepared.csv.gz"))
lk_file  <- function(sec) file.path(dir_in, paste0(sec, "_row_id_lookup.RDS"))
prepared <- all(file.exists(sapply(section_names, in_file)), file.exists(sapply(section_names, lk_file)))

if (!prepared) {
  txfile <- list.files(slide_dir, pattern = "tx_file.*\\.gz$", full.names = TRUE)[1]
  stopifnot("tx_file not found" = !is.na(txfile))
  cols_needed <- c("fov", "cell_ID", "x_global_px", "y_global_px", "target", "z", "CellComp")
  cols_read   <- intersect(cols_needed, names(fread(txfile, nrows = 0)))
  tx <- fread(txfile, select = cols_read)
  stopifnot(all(c("fov", "cell_ID", "x_global_px", "y_global_px", "target") %in% colnames(tx)))
  
  ## Control probes are left in (removing rows would copy the 73M-row table);
  ## ProSeg drops them itself with --excluded-genes.
  tx[fov_to_section, on = "fov", section_id := i.section_id]
  if (anyNA(tx$section_id)) tx <- tx[!is.na(section_id)]
  
  ## Cell name per transcript; cell_ID = 0 = outside every cell
  tx[, seurat_cell := ""]
  tx[cell_ID != 0, seurat_cell := paste0(cell_prefix, "_", fov, "_", cell_ID)]
  pct_matched <- 100 * mean(unique(tx$seurat_cell[tx$seurat_cell != ""]) %in% metadata$cell_id)
  message("Transcript cell names found in metadata: ", round(pct_matched, 2), "%")
  stopifnot(pct_matched > 95)
  
  keep_cols <- intersect(c("fov", "seurat_cell", "cell_ID", "x_global_px", "y_global_px",
                           "target", "z", "CellComp"), colnames(tx))
  for (sec in section_names) {
    d <- tx[section_id == sec, ..keep_cols]
    d[, row_id := seq_len(.N)]                     # unique ID carried through ProSeg
    saveRDS(d[, .(row_id, seurat_cell, x_global_px, y_global_px)], lk_file(sec))
    fwrite(d, in_file(sec), compress = "gzip")
    message(sec, ": ", fmt(nrow(d)), " transcripts written")
  }
  rm(tx, d); invisible(gc())
}

## Columns present in the prepared files (decides the optional ProSeg options)
in_cols <- names(fread(in_file(section_names[1]), nrows = 0))
has_z           <- "z" %in% in_cols
has_compartment <- "CellComp" %in% in_cols


## =============================================================================
## CHUNK 4: run ProSeg on each section
## =============================================================================
##   --transcript-id-column row_id : carry row_id to the output (for linking)
##   --coordinate-scale 0.12028    : pixels -> µm
##   --cell-assignment-column      : AtoMx's cell per transcript, the starting point
build_proseg_args <- function(tx_path, out_path) {
  a <- c(tx_path,
         "--output-transcript-metadata", file.path(out_path, "transcript-metadata.csv"),
         "--output-cell-metadata",       file.path(out_path, "cell-metadata.csv"),
         "--output-spatialdata",         file.path(out_path, "spatialdata.zarr"),
         "--overwrite", "--exclude-spatialdata-transcripts",
         "--transcript-id-column", "row_id",
         "--gene-column", "target", "--x-column", "x_global_px", "--y-column", "y_global_px",
         "--fov-column", "fov",
         "--cell-id-column", "seurat_cell", "--cell-id-unassigned", "",
         "--cell-assignment-column", "cell_ID", "--cell-assignment-unassigned", "0",
         "--non-unique-cell-ids",
         "--excluded-genes", excluded_genes_regex,
         "--coordinate-scale", as.character(pixel_size_um))
  a <- if (has_z) c(a, "--z-column", "z") else c(a, "--z-column", "x_global_px", "--ignore-z-coord")
  if (has_compartment) a <- c(a, "--compartment-column", "CellComp", "--compartment-nuclear", "Nuclear")
  if (!is.null(max_transcript_nucleus_distance)) a <- c(a, "--max-transcript-nucleus-distance", max_transcript_nucleus_distance)
  if (!is.null(n_threads)) a <- c(a, "--nthreads", n_threads)
  a
}
quote_for_shell <- function(x) ifelse(x == "" | grepl("[^A-Za-z0-9_./=-]", x), shQuote(x), x)

for (sec in section_names) {
  out_path <- file.path(dir_raw, sec); dir.create(out_path, showWarnings = FALSE)
  tm_path  <- file.path(out_path, "transcript-metadata.csv")
  if (file.exists(tm_path) && !rerun_proseg) { message(sec, ": ProSeg output exists - reused"); next }
  message(sec, ": running ProSeg (can take a long time)...")
  log_lines <- system2("proseg", quote_for_shell(build_proseg_args(in_file(sec), out_path)),
                       stdout = TRUE, stderr = TRUE)
  writeLines(log_lines, file.path(out_path, "proseg_run_log.txt"))
  if (!file.exists(tm_path)) warning(sec, ": ProSeg wrote no output - see proseg_run_log.txt")
  else message(sec, ": done")
}


## =============================================================================
## CHUNK 5: link ProSeg output to cells; purity per cell; transcripts of flagged cells
## =============================================================================
score_section <- function(sec) {
  tm_path <- file.path(dir_raw, sec, "transcript-metadata.csv")
  if (!file.exists(tm_path)) { warning(sec, ": no ProSeg output - skipped"); return(NULL) }
  ps <- fread(tm_path, select = c("transcript_id", "observed_x", "observed_y", "gene", "assignment", "background"))
  lookup <- readRDS(lk_file(sec))
  
  ## row_id is 1..N, so ProSeg's transcript_id = lookup row number. Link by index,
  ## then confirm by coordinates (ProSeg position = pixel position x 0.12028).
  stopifnot(identical(lookup$row_id, seq_len(nrow(lookup))),
            all(ps$transcript_id >= 1 & ps$transcript_id <= nrow(lookup)))
  idx <- as.integer(ps$transcript_id)
  max_dev <- max(pmax(abs(ps$observed_x - lookup$x_global_px[idx] * pixel_size_um),
                      abs(ps$observed_y - lookup$y_global_px[idx] * pixel_size_um)))
  if (anyDuplicated(idx) || max_dev > coord_tolerance_um) {
    warning(sec, ": linking check failed (max deviation ", signif(max_dev, 3), " µm) - skipped"); return(NULL)
  }
  ps[, seurat_cell := lookup$seurat_cell[idx]]
  if (!is.logical(ps$background)) ps[, background := tolower(as.character(background)) %in% c("true", "1")]
  
  ## Evaluable transcript: in an AtoMx cell, not background, assigned to a ProSeg cell
  ev <- ps[seurat_cell != "" & !background & !is.na(assignment)]
  counts_pc <- ev[, .N, by = .(seurat_cell, assignment)]
  p <- counts_pc[, .(n_transcripts_evaluated = sum(N), n_in_top_proseg_cell = max(N),
                     top_proseg_cell = assignment[which.max(N)],
                     n_proseg_cells_spanned = .N), by = seurat_cell]
  p[, purity := n_in_top_proseg_cell / n_transcripts_evaluated]
  p[n_transcripts_evaluated < min_transcripts_eval, purity := NA_real_]
  
  ## Every transcript of the flagged cells, with what ProSeg did to it
  flagged_ids <- p[!is.na(purity) & purity < purity_threshold, seurat_cell]
  tx_f <- ps[seurat_cell %chin% flagged_ids]
  tx_f[p, on = "seurat_cell", top_proseg_cell := i.top_proseg_cell]
  tx_f[, transcript_fate := fifelse(background, "background (noise)",
                                    fifelse(is.na(assignment), "not assigned by ProSeg",
                                            fifelse(assignment == top_proseg_cell, "kept in main ProSeg cell",
                                                    "moved to another ProSeg cell")))]
  tx_f <- tx_f[, .(cell_id = seurat_cell, section_id = sec, transcript_id, gene,
                   x_um = observed_x, y_um = observed_y, proseg_cell = assignment,
                   main_proseg_cell_of_this_cell = top_proseg_cell, background, transcript_fate)]
  
  report <- data.table(section_id = sec, transcripts_in_input = nrow(lookup),
                       transcripts_kept_by_proseg = nrow(ps),
                       transcripts_evaluable = nrow(ev),
                       max_coordinate_deviation_um = signif(max_dev, 3))
  message(sec, ": linked and verified; ", fmt(nrow(p)), " cells scored, ",
          fmt(length(flagged_ids)), " flagged")
  list(p = p, tx_f = tx_f, report = report)
}

sc <- setNames(lapply(section_names, score_section), section_names)
p_all       <- rbindlist(lapply(sc, `[[`, "p"))
link_report <- rbindlist(lapply(sc, `[[`, "report"))

## One row per cell in metadata (same order); not evaluated -> not flagged
i <- match(metadata$cell_id, p_all$seurat_cell)
res <- data.table(cell_id = metadata$cell_id, section_id = metadata$section_id, fov = metadata$fov,
                  x_slide_mm = metadata$x_slide_mm, y_slide_mm = metadata$y_slide_mm,
                  nCount_RNA = metadata$nCount_RNA,
                  n_transcripts_evaluated = p_all$n_transcripts_evaluated[i],
                  n_in_main_proseg_cell   = p_all$n_in_top_proseg_cell[i],
                  n_proseg_cells_spanned  = p_all$n_proseg_cells_spanned[i],
                  purity                  = p_all$purity[i])
res[, evaluated := !is.na(purity)]
res[, flagged := evaluated & purity < purity_threshold]
res[, status := factor(fifelse(!evaluated, "not evaluated",
                               fifelse(flagged, "flagged (purity < 0.5)", "pass")),
                       levels = c("pass", "not evaluated", "flagged (purity < 0.5)"))]


## =============================================================================
## CHUNK 6: save results (CSV)
## =============================================================================
summary_sec <- res[, .(n_cells = .N, n_evaluated = sum(evaluated), n_not_evaluated = sum(!evaluated),
                       median_purity = round(median(purity, na.rm = TRUE), 3),
                       n_flagged = sum(flagged),
                       pct_flagged_of_all_cells = round(100 * mean(flagged), 2),
                       pct_flagged_of_evaluated = round(100 * sum(flagged) / sum(evaluated), 2)),
                   by = section_id][order(section_id)]
summary_sec <- rbind(summary_sec, res[, .(section_id = "All sections", n_cells = .N,
                                          n_evaluated = sum(evaluated), n_not_evaluated = sum(!evaluated),
                                          median_purity = round(median(purity, na.rm = TRUE), 3), n_flagged = sum(flagged),
                                          pct_flagged_of_all_cells = round(100 * mean(flagged), 2),
                                          pct_flagged_of_evaluated = round(100 * sum(flagged) / sum(evaluated), 2))])
summary_sec[, flagged_label := paste0(fmt(n_flagged), " of ", fmt(n_cells), " cells flagged")]

## Transcripts of flagged cells: what happened to them, per section
tx_fate <- rbindlist(lapply(sc, function(x) x$tx_f[, .N, by = .(section_id, transcript_fate)]))
tx_fate_wide <- dcast(tx_fate, section_id ~ transcript_fate, value.var = "N", fill = 0)

fwrite(summary_sec,  file.path(dir_res, "01_ProSeg_flagged_cells_summary_per_section.csv"))
fwrite(res,          file.path(dir_res, "02_ProSeg_purity_all_cells.csv"))
fwrite(res[flagged == TRUE], file.path(dir_res, "03_ProSeg_flagged_cells.csv"))
fwrite(tx_fate_wide, file.path(dir_res, "04_ProSeg_flagged_cells_transcript_fate_per_section.csv"))
fwrite(link_report,  file.path(dir_res, "05_ProSeg_linking_report_per_section.csv"))
for (sec in names(sc)) if (!is.null(sc[[sec]]))
  fwrite(sc[[sec]]$tx_f, file.path(dir_tx_out, paste0(sec, "_flagged_cell_transcripts.csv.gz")), compress = "gzip")
saveRDS(setNames(res$flagged, res$cell_id), file.path(dir_res, "flag_ProSeg_low_purity.RDS"))

print(summary_sec); print(tx_fate_wide); print(link_report)


## =============================================================================
## CHUNK 7: plots
## =============================================================================
st_cols   <- c("pass" = "grey80", "not evaluated" = "grey40", "flagged (purity < 0.5)" = "red")
fate_cols <- c("kept in main ProSeg cell" = "grey60", "moved to another ProSeg cell" = "red",
               "background (noise)" = "blue", "not assigned by ProSeg" = "orange")
dp  <- res[order(status)]                                   # flagged drawn on top
tot <- summary_sec[section_id == "All sections"]
lab <- merge(res[, .(x = mean(range(x_slide_mm)), y = max(y_slide_mm)), by = section_id],
             summary_sec[, .(section_id, flagged_label, pct_flagged_of_evaluated)], by = "section_id")
lab[, label := paste0(section_id, ": ", flagged_label, " (", pct_flagged_of_evaluated, "% of evaluated)")]

## 1. Whole slide: 6 sections in position, each run separately, counts above each
p_slide <- ggplot(dp, aes(x_slide_mm, y_slide_mm, colour = status)) +
  geom_point(size = 0.05) + scale_colour_manual(values = st_cols, drop = FALSE) +
  geom_label(data = lab, aes(x, y, label = label), inherit.aes = FALSE, vjust = -0.3,
             size = 3, fontface = "bold", fill = "white", label.size = 0.2) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.08))) +
  coord_fixed() + theme_bw() + guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(title = paste0("ProSeg low purity, per tissue section: ", tot$flagged_label), 
       x = "x (mm)", y = "y (mm)", colour = NULL)
ggsave(file.path(dir_res, "01_ProSeg_flagged_cells_on_slide.png"), p_slide, width = 12, height = 12, dpi = 300)
print(p_slide)

## 2. Grid, one panel per section
pl <- setNames(paste0(summary_sec$section_id, ": ", summary_sec$flagged_label), summary_sec$section_id)
dp[, panel := pl[section_id]]
p_grid <- ggplot(dp, aes(x_slide_mm, y_slide_mm, colour = status)) +
  geom_point(size = 0.05) + scale_colour_manual(values = st_cols, drop = FALSE) +
  facet_wrap(~panel, nrow = 2, scales = "free") + theme_bw() +
  guides(colour = guide_legend(override.aes = list(size = 3))) +
  labs(title = paste0("ProSeg low purity, per tissue section: ", tot$flagged_label),
       x = "x (mm)", y = "y (mm)", colour = NULL)
ggsave(file.path(dir_res, "02_ProSeg_flagged_cells_per_section_grid.png"), p_grid, width = 16, height = 10, dpi = 300)

## 3. Purity histograms
p_hist <- ggplot(res[evaluated == TRUE], aes(purity)) +
  geom_histogram(bins = 50, fill = "steelblue") +
  geom_vline(xintercept = purity_threshold, colour = "red", linetype = "dashed") +
  facet_wrap(~section_id, nrow = 2, scales = "free_y") + theme_classic() +
  labs(title = "ProSeg purity of evaluated cells", x = "ProSeg purity", y = "N cells")
ggsave(file.path(dir_res, "03_ProSeg_purity_histograms.png"), p_hist, width = 14, height = 8, dpi = 300)

## 4. Per section: flagged cells, and the transcripts of flagged cells by fate
for (sec in section_names) {
  r <- summary_sec[section_id == sec]
  ggsave(file.path(dir_sec_plots, paste0(sec, "_01_ProSeg_flagged_cells.png")),
         ggplot(dp[section_id == sec], aes(x_slide_mm, y_slide_mm, colour = status)) +
           geom_point(size = 0.1) + scale_colour_manual(values = st_cols, drop = FALSE) +
           coord_fixed() + theme_bw() + guides(colour = guide_legend(override.aes = list(size = 3))) +
           labs(title = paste0(sec, ": ", r$flagged_label),
                subtitle = paste0(fmt(r$n_evaluated), " evaluated, median purity ", r$median_purity),
                x = "x (mm)", y = "y (mm)", colour = NULL),
         width = 8, height = 8, dpi = 300)
  
  tf <- if (!is.null(sc[[sec]])) sc[[sec]]$tx_f else NULL
  if (!is.null(tf) && nrow(tf) > 0) {
    tf <- tf[order(transcript_fate == "moved to another ProSeg cell")]   # moved drawn on top
    n_f <- tf[, .N, by = transcript_fate]
    sub <- paste(paste0(n_f$transcript_fate, ": ", fmt(n_f$N)), collapse = " | ")
    ggsave(file.path(dir_sec_plots, paste0(sec, "_02_ProSeg_flagged_cell_transcripts.png")),
           ggplot(tf, aes(x_um / 1000, y_um / 1000, colour = transcript_fate)) +
             geom_point(size = 0.02) + scale_colour_manual(values = fate_cols, drop = FALSE) +
             coord_fixed() + theme_bw() + guides(colour = guide_legend(override.aes = list(size = 3))) +
             labs(title = paste0(sec, ": transcripts of ", fmt(r$n_flagged), " flagged cells"),
                  subtitle = sub, x = "x (mm)", y = "y (mm)", colour = NULL),
           width = 9, height = 9, dpi = 300)
  }
}


## =============================================================================
## CHUNK 8: results on screen
## =============================================================================
for (sec in c(section_names, "All sections")) {
  r <- summary_sec[section_id == sec]
  message(sec, ": ", r$flagged_label, " (", r$pct_flagged_of_evaluated, "% of ",
          fmt(r$n_evaluated), " evaluated)")
}
message("ProSeg results saved to ", dir_res)















## =============================================================================
## 10_QC_descriptive_plots.R
## =============================================================================
## PURPOSE
##   Vignette 1, section 4.3 Descriptive QC plots (no filtering, no action):
##     1. total counts over space
##     2. spatially smoothed, totalcount-normalised background over space
##        (mean of nCount_negprobes / nCount_RNA over each cell's 50 nearest neighbours)
##   Drawn three ways:
##     all_cells/            whole slide (vignette as written)
##     per_tissue_section/   one plot per section; background smoothed within the
##                           section and coloured on that section's own scale
##                           + one plot of all sections on the slide, labelled per section
##   Vignette code used as written; changes are marked "## added".
##
## Run after script 02 (needs metadata with section_id, xy). Chunk by chunk.
## =============================================================================

##install.packages('InSituCor')
## =============================================================================
## CHUNK 52: packages, inputs, folders
## =============================================================================
## Needed once:
##install.packages(c("scales", "spatstat.geom", "devtools"))
##devtools::install_github("https://github.com/Nanostring-Biostats/InSituCor")
##renv::install("Nanostring-Biostats/InSituCor")   # installs from GitHub into your renv project library

remotes::install_github("Nanostring-Biostats/InSituCor", upgrade = "never")
library(scales)   # viridis_pal() (vignette)
library(data.table)
library(InSituCor)


output_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"

## Load metadata and xy from disk if they are not already in this R session
if (!exists("metadata")) metadata <- readRDS(file.path(output_dir, "metadata_unfiltered.RDS"))
if (!exists("xy"))       xy       <- readRDS(file.path(output_dir, "xy_unfiltered.RDS"))

stopifnot(identical(rownames(xy), metadata$cell_id),   # same cells, same order
          "section_id" %in% colnames(metadata))        # section labels were saved


dir_desc     <- file.path(output_dir, "QC", "Descriptive_plots")
dir_desc_all <- file.path(dir_desc, "all_cells")
dir_desc_sec <- file.path(dir_desc, "per_tissue_section")
dir.create(dir_desc_all, showWarnings = FALSE, recursive = TRUE)
dir.create(dir_desc_sec, showWarnings = FALSE, recursive = TRUE)
fmt <- function(n) format(n, big.mark = ",")

## Save PNG at the shape of the cells plotted, then show in a separate window   ## added
show_and_save <- function(draw_fun, file, pts) {
  rx <- diff(range(pts[, 1])); ry <- diff(range(pts[, 2]))
  w <- 10 * rx / max(rx, ry) + 0.5; h <- 10 * ry / max(rx, ry) + 0.8
  png(file, width = w, height = h, units = "in", res = 300); draw_fun(); dev.off()
  dev.new(width = w, height = h, noRStudioGD = TRUE); draw_fun()
  message("Saved: ", file)
}

## Section label positions (top centre of each section)                      ## added
sec_pos <- as.data.table(metadata)[, .(x = mean(range(x_slide_mm)), y = max(y_slide_mm)), by = section_id]


## =============================================================================
## CHUNK 53: colour functions (vignette formulas)
## =============================================================================
## Plot 1: fixed scale, 10 to 5000 counts, log2
counts_col <- function(n) {
  viridis_pal(option = "B")(101)[
    round(1 + 100 * (log2(pmin(pmax(n, 10), 5000)) - log2(10)) / log2(5000))]
}
counts_legend <- function() {
  legendvals <- c(10, 100, 1000, 5000)
  legend("bottomright", pch = 16,
         col = c("white", viridis_pal(option = "B")(101)[round(1 + 100 * (log2(legendvals) - log2(10)) / log2(5000))]),
         legend = c("Total counts:", legendvals))
}

## Plot 2: scale relative to the 99.95th percentile of the values plotted
bg_col <- function(v) {
  viridis_pal(option = "B")(101)[round(1 + 100 * pmin(v / quantile(v, 0.9995), 1))]
}
bg_legend <- function(v) {
  legendvals <- signif(median(v) * c(0.5,1,2,5,10), 2)
  legend("bottomright", pch = 16,
         col = c("white", viridis_pal(option = "B")(101)[signif(1 + 100 * pmin(legendvals / quantile(v, 0.9995), 1), 2)]),
         legend = c("Spatially smoothed totalcounts-normalized background:", paste0(c(0.5,1,2,5,10), "x median")))
}


## =============================================================================
## CHUNK 54: smoothed background (vignette), whole slide and within each section
## =============================================================================
## Background rate per cell. A cell with 0 total counts gives x/0 (Inf/NaN),
## which would spread NaN to its 50 neighbours; such cells are set to 0.      ## added
bg_rate <- metadata$nCount_negprobes / metadata$nCount_RNA
message("Cells with nCount_RNA = 0 (rate set to 0): ", sum(!is.finite(bg_rate)))
bg_rate[!is.finite(bg_rate)] <- 0

## Whole slide: vignette uses subset = metadata$slide_ID (one slide here)
slide_subset <- if ("slide_ID" %in% colnames(metadata)) metadata$slide_ID else rep(1, nrow(metadata))
neighbors <- InSituCor:::nearestNeighborGraph(xy[, 1], xy[, 2], N = 50, subset = slide_subset)
smoothedgbrate <- as.vector(InSituCor:::neighbor_mean(x = bg_rate, neighbors = neighbors))

## Per section: neighbours restricted to the same section                    ## added
neighbors_sec <- InSituCor:::nearestNeighborGraph(xy[, 1], xy[, 2], N = 50, subset = metadata$section_id)
smoothedgbrate_sec <- as.vector(InSituCor:::neighbor_mean(x = bg_rate, neighbors = neighbors_sec))

## Summary per section (printed and saved)                                   ## added
desc_summary <- data.table(section_id = metadata$section_id, nCount_RNA = metadata$nCount_RNA,
                           bg = smoothedgbrate_sec)[, .(
                             n_cells                    = .N,
                             median_nCount_RNA          = median(nCount_RNA),
                             median_smoothed_background = signif(median(bg), 3),
                             p99.95_smoothed_background = signif(quantile(bg, 0.9995), 3)), by = section_id][order(section_id)]
desc_summary <- rbind(desc_summary, data.table(section_id = "All cells", n_cells = nrow(metadata),
                                               median_nCount_RNA = median(metadata$nCount_RNA),
                                               median_smoothed_background = signif(median(smoothedgbrate), 3),
                                               p99.95_smoothed_background = signif(quantile(smoothedgbrate, 0.9995), 3)))
print(desc_summary)
fwrite(desc_summary, file.path(dir_desc, "descriptive_summary_per_section.csv"))


## =============================================================================
## CHUNK 55: ALL CELLS (vignette as written)
## =============================================================================
draw_counts_all <- function() {
  par(mar = c(0,0,0,0))
  plot(xy, pch = 16, cex = 0.1, asp = 1, col = counts_col(metadata$nCount_RNA))
  counts_legend()
}
show_and_save(draw_counts_all, file.path(dir_desc_all, "01_total_counts_over_space_all_cells.png"), xy)

draw_bg_all <- function() {
  par(mar = c(0,0,0,0))
  plot(xy, pch = 16, cex = 0.1, asp = 1, col = bg_col(smoothedgbrate))
  bg_legend(smoothedgbrate)
}
show_and_save(draw_bg_all, file.path(dir_desc_all, "02_smoothed_background_over_space_all_cells.png"), xy)


## =============================================================================
## CHUNK 56: PER TISSUE SECTION - one plot per section
## =============================================================================
for (sec in section_names) {
  idx <- metadata$section_id == sec
  xy_s <- xy[idx, , drop = FALSE]
  n_s  <- metadata$nCount_RNA[idx]
  bg_s <- smoothedgbrate_sec[idx]
  d_out <- file.path(dir_desc_sec, sec); dir.create(d_out, showWarnings = FALSE)
  
  ## 1. total counts (same fixed scale as the whole slide)
  show_and_save(function() {
    par(mar = c(0, 0, 2, 0))                                                ## added: title room
    plot(xy_s, pch = 16, cex = 0.1, asp = 1, col = counts_col(n_s))
    title(paste0(sec, ": ", fmt(sum(idx)), " cells, median ", round(median(n_s)), " counts"), cex.main = 0.9)
    counts_legend()
  }, file.path(d_out, paste0("01_total_counts_over_space_", sec, ".png")), xy_s)
  
  ## 2. smoothed background (scale relative to this section)
  show_and_save(function() {
    par(mar = c(0, 0, 2, 0))
    plot(xy_s, pch = 16, cex = 0.1, asp = 1, col = bg_col(bg_s))
    title(paste0(sec, ": median smoothed background ", signif(median(bg_s), 3)), cex.main = 0.9)
    bg_legend(bg_s)
  }, file.path(d_out, paste0("02_smoothed_background_over_space_", sec, ".png")), xy_s)
}


## =============================================================================
## CHUNK 57: PER TISSUE SECTION - all sections on one plot (slide positions)
## =============================================================================
## Same colour scale for every section, so sections can be compared;
## each section labelled with its own values.                                 ## added
add_sec_labels <- function(lab) {
  for (k in seq_len(nrow(sec_pos))) {
    text(sec_pos$x[k], sec_pos$y[k], pos = 3, xpd = NA, cex = 0.7, font = 2,
         labels = lab(sec_pos$section_id[k]))
  }
}

show_and_save(function() {
  par(mar = c(0, 0, 2, 0))
  plot(xy, pch = 16, cex = 0.1, asp = 1, col = counts_col(metadata$nCount_RNA))
  title("Total counts per tissue section", cex.main = 0.9)
  add_sec_labels(function(s) {
    r <- desc_summary[section_id == s]
    paste0(s, ": ", fmt(r$n_cells), " cells, median ", round(r$median_nCount_RNA), " counts")
  })
  counts_legend()
}, file.path(dir_desc_sec, "01_total_counts_over_space_per_section_on_slide.png"), xy)

show_and_save(function() {
  par(mar = c(0, 0, 2, 0))
  plot(xy, pch = 16, cex = 0.1, asp = 1, col = bg_col(smoothedgbrate_sec))
  title("Smoothed background per tissue section (neighbours within section, common scale)", cex.main = 0.9)
  add_sec_labels(function(s) paste0(s, ": median ", desc_summary[section_id == s, median_smoothed_background]))
  bg_legend(smoothedgbrate_sec)
}, file.path(dir_desc_sec, "02_smoothed_background_over_space_per_section_on_slide.png"), xy)

message("Descriptive plots saved to ", dir_desc)













## =============================================================================
## 10_QC_combined_flags.R
## =============================================================================
## PURPOSE
##   Combine the cells flagged by the QC scripts, two ways:
##
##   PART A - PER TISSUE SECTION, QCs 03 to 09 (7 QCs)
##     A cell is removed if ANY of the 7 QCs flags it; each cell is counted ONCE.
##     Outputs: CSV of every cell with the QCs that flagged it, counts per
##     section, slide plot (red = removed), per-section plots, UpSet plot of cells,
##     and the same at TRANSCRIPT level (transcripts in removed cells, in red,
##     plus an UpSet plot of transcripts).
##
##   PART B - ALL CELLS analysis, QCs 03 to 08 (6 QCs; ProSeg was run per section only)
##     Same rule; slide plot with the number of removed cells, CSVs, UpSet plot of cells.
##
##   Flags used (all saved by earlier scripts):
##     QC                     PART A (per section)                    PART B (all cells)
##     03 Low counts          flag_low_counts_per_section.RDS         flag_low_counts_all_cells.RDS
##     04 Large area          Area > 30,000 px (fixed threshold: identical in A and B)
##     05 FOV QC              fovqcresult.RDS (run once on the whole slide: identical in A and B)
##     06 nCount/nFeature     flag_ncount_nfeature_per_section.RDS    flag_ncount_nfeature_all_cells.RDS
##     07 FOV border          flag_FOV_border_cells.RDS (fixed threshold: identical in A and B)
##     08 Low SBR             flag_low_SBR.RDS (fixed threshold: identical in A and B)
##     09 ProSeg low purity   flag_ProSeg_low_purity.RDS (per section) -
##   Nothing is removed from the data here; the final flags are saved for removal.
##
## Run in a fresh R session after scripts 02-09. Chunk by chunk.
## =============================================================================


## =============================================================================
## CHUNK 1: packages, paths, folders
## =============================================================================
## Needed once: install.packages("UpSetR")
library(data.table)
library(UpSetR)       # upset(), fromList(), fromExpression()

output_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"
qc_dir     <- file.path(output_dir, "QC")

metadata <- readRDS(file.path(output_dir, "metadata_unfiltered.RDS"))
xy       <- readRDS(file.path(output_dir, "xy_unfiltered.RDS"))            # x_mm, y_mm
stopifnot(identical(rownames(metadata), metadata$cell_id),
          identical(rownames(xy), metadata$cell_id))
section_names <- sort(unique(metadata$section_id))

dir_A <- file.path(qc_dir, "independent_cells_per_tissue_section_03_to_09_QCs")
dir_B <- file.path(qc_dir, "independent_cells_all_cells_03_to_08_QCs")
for (d in c(file.path(dir_A, "cells", "per_tissue_section"), file.path(dir_A, "transcripts", "per_tissue_section"),
            file.path(dir_B, "per_tissue_section")))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

## ProSeg input files (per section) hold every transcript with its cell and position
dir_ps_in <- file.path(qc_dir, "QC_Proseg", "per_tissue_section", "01_prepared_transcripts")
control_regex <- "^(FalseCode|SystemControl|NegPrb|Negative)"   # control probes: not RNA
px_to_mm      <- 0.120280945 / 1000                             # same conversion as script 02

fmt <- function(n) format(n, big.mark = ",")


## =============================================================================
## CHUNK 2: load every flag
## =============================================================================
rd <- function(...) {                                           # read one saved flag, check cell order
  f <- readRDS(file.path(qc_dir, ...))
  stopifnot(identical(names(f), metadata$cell_id))
  unname(f)
}
fovqcresult <- readRDS(file.path(qc_dir, "FOV", "fovqcresult.RDS"))
length(fovqcresult$flaggedfovs)                                 # flagged FOVs

area_flag   <- metadata$Area > 30000
fov_flag    <- is.element(metadata$fov, fovqcresult$flaggedfovs)
border_flag <- rd("FOV_Borders", "all_cells", "flag_FOV_border_cells.RDS")
sbr_flag    <- rd("SBR", "all_cells", "flag_low_SBR.RDS")

flags_A <- list(                                                 # PART A: per tissue section
  `03_Low_counts`        = rd("QC_low_counts_per_tissue_section", "flag_low_counts_per_section.RDS"),
  `04_Large_area`        = area_flag,
  `05_FOV_QC`            = fov_flag,
  `06_nCount_nFeature`   = rd("ncount_and_nfeature", "per_tissue_section", "flag_ncount_nfeature_per_section.RDS"),
  `07_FOV_border`        = border_flag,
  `08_Low_SBR`           = sbr_flag,
  `09_ProSeg_low_purity` = rd("QC_Proseg", "per_tissue_section", "03_results", "flag_ProSeg_low_purity.RDS"))

flags_B <- list(                                                 # PART B: all cells
  `03_Low_counts`      = rd("QC_low_counts_per_cell_all_cells", "flag_low_counts_all_cells.RDS"),
  `04_Large_area`      = area_flag,
  `05_FOV_QC`          = fov_flag,
  `06_nCount_nFeature` = rd("ncount_and_nfeature", "all_cells", "flag_ncount_nfeature_all_cells.RDS"),
  `07_FOV_border`      = border_flag,
  `08_Low_SBR`         = sbr_flag)


## =============================================================================
## CHUNK 3: helper functions
## =============================================================================
## Per-cell table: one TRUE/FALSE column per QC, number of QCs, the QCs listed
## by name, and removed = flagged by at least one QC (each cell once).
build_cell_table <- function(fl) {
  d <- data.table(cell_id = metadata$cell_id, section_id = metadata$section_id, fov = metadata$fov,
                  x_slide_mm = metadata$x_slide_mm, y_slide_mm = metadata$y_slide_mm,
                  nCount_RNA = metadata$nCount_RNA)
  for (q in names(fl)) set(d, j = q, value = fl[[q]])
  d[, n_QCs := rowSums(.SD), .SDcols = names(fl)]
  removed_by <- rep("", nrow(d))                        # e.g. "03_Low_counts; 08_Low_SBR"
  for (q in names(fl)) removed_by <- ifelse(fl[[q]], ifelse(removed_by == "", q, paste(removed_by, q, sep = "; ")), removed_by)
  d[, QCs_removed_by := removed_by]
  d[, removed := n_QCs > 0]
  d
}

## Per-section summary: cells, cells flagged by each QC, independent removed, kept
summarise_cells <- function(d, qcs) {
  one <- function(x, label) {
    r <- data.table(section_id = label, n_cells = nrow(x))
    for (q in qcs) set(r, j = paste0("n_", q), value = sum(x[[q]]))
    r[, `:=`(total_flag_events = sum(x$n_QCs),                     # a cell counted once per QC
             n_removed_independent = sum(x$removed),              # each cell once
             pct_removed = round(100 * mean(x$removed), 2),
             n_kept = sum(!x$removed),
             n_removed_by_1_QC = sum(x$n_QCs == 1),
             n_removed_by_2_or_more_QCs = sum(x$n_QCs >= 2))]
    r[, removed_label := paste0(fmt(n_removed_independent), " of ", fmt(n_cells), " cells removed")]
    r
  }
  rbind(rbindlist(lapply(section_names, function(s) one(d[section_id == s], s))), one(d, "All sections"))
}

## Save PNG at slide shape, then show in a separate window
show_and_save <- function(draw_fun, file, pts, extra_h = 0.8) {
  rx <- diff(range(pts[, 1])); ry <- diff(range(pts[, 2]))
  w <- 10 * rx / max(rx, ry) + 0.5; h <- 10 * ry / max(rx, ry) + extra_h
  png(file, width = w, height = h, units = "in", res = 300); draw_fun(); dev.off()
  dev.new(width = w, height = h, noRStudioGD = TRUE); draw_fun()
  message("Saved: ", file)
}

sec_pos <- as.data.table(metadata)[, .(x = mean(range(x_slide_mm)), y = max(y_slide_mm)), by = section_id]

## Whole slide: grey80 = kept, red = removed (drawn on top); counts above each section
draw_cells_slide <- function(d, s, title) {
  o <- order(d$removed)
  par(mar = c(0, 0, 2, 0))
  plot(xy[o, ], pch = 16, cex = 0.1, asp = 1, col = c("grey80", "red")[1 + d$removed[o]])
  tot <- s[section_id == "All sections"]
  title(paste0(title, ": ", tot$removed_label, " (", tot$pct_removed, "%)"), cex.main = 0.9)
  for (k in seq_len(nrow(sec_pos))) {
    r <- s[section_id == sec_pos$section_id[k]]
    text(sec_pos$x[k], sec_pos$y[k], pos = 3, xpd = NA, cex = 0.7, font = 2,
         labels = paste0(r$section_id, ": ", r$removed_label, " (", r$pct_removed, "%)"))
  }
}

## 2 x 3 grid, one panel per section
draw_cells_grid <- function(d, s, title) {
  op <- par(mfrow = c(2, 3), mar = c(0.5, 0.5, 2.5, 0.5), oma = c(0, 0, 2.5, 0)); on.exit(par(op))
  for (sec in section_names) {
    idx <- which(d$section_id == sec); o <- idx[order(d$removed[idx])]
    plot(xy[o, ], pch = 16, cex = 0.1, asp = 1, axes = FALSE, xlab = "", ylab = "",
         col = c("grey80", "red")[1 + d$removed[o]])
    r <- s[section_id == sec]
    title(paste0(sec, ": ", r$removed_label, " (", r$pct_removed, "%)"), cex.main = 1)
  }
  tot <- s[section_id == "All sections"]
  mtext(paste0(title, ": ", tot$removed_label, " (", tot$pct_removed, "%)"), outer = TRUE, font = 2)
}

## One plot per section
save_cells_per_section <- function(d, s, out_dir, tag) {
  for (sec in section_names) {
    idx <- which(d$section_id == sec); o <- idx[order(d$removed[idx])]
    r <- s[section_id == sec]
    png(file.path(out_dir, paste0(sec, "_removed_cells_", tag, ".png")), width = 8, height = 8, units = "in", res = 300)
    par(mar = c(0, 0, 2, 0))
    plot(xy[o, ], pch = 16, cex = 0.15, asp = 1, col = c("grey80", "red")[1 + d$removed[o]])
    title(paste0(sec, " (", tag, "): ", r$removed_label, " (", r$pct_removed, "%)"), cex.main = 0.9)
    dev.off()
  }
}


## =============================================================================
## CHUNK 4: PART A - cells, per tissue section, QCs 03 to 09
## =============================================================================
qcs_A   <- names(flags_A)
cells_A <- build_cell_table(flags_A)
sum_A   <- summarise_cells(cells_A, qcs_A)
print(sum_A)

fwrite(sum_A,                      file.path(dir_A, "cells", "01_removed_cells_summary_per_section.csv"))
fwrite(cells_A,                    file.path(dir_A, "cells", "02_all_cells_with_QCs_removed_by.csv"))
fwrite(cells_A[removed == TRUE],   file.path(dir_A, "cells", "03_removed_cells_with_QCs_removed_by.csv"))
saveRDS(setNames(cells_A$removed, cells_A$cell_id), file.path(dir_A, "final_QC_flag_per_tissue_section_03_to_09.RDS"))

## Plots: whole slide, grid, one per section
show_and_save(function() draw_cells_slide(cells_A, sum_A, "Removed cells, per tissue section QCs 03-09"),
              file.path(dir_A, "cells", "04_removed_cells_on_slide.png"), xy)
png(file.path(dir_A, "cells", "05_removed_cells_per_section_grid.png"), width = 15, height = 10, units = "in", res = 300)
draw_cells_grid(cells_A, sum_A, "Removed cells, per tissue section QCs 03-09"); dev.off()
save_cells_per_section(cells_A, sum_A, file.path(dir_A, "cells", "per_tissue_section"), "QC03_09_per_section")

## UpSet plot of cells (each cell in every QC set that flags it)
cell_sets <- lapply(qcs_A, function(q) cells_A$cell_id[cells_A[[q]]]); names(cell_sets) <- qcs_A
cell_sets <- Filter(function(x) length(x) > 0, cell_sets)
png(file.path(dir_A, "cells", "06_upset_removed_cells.png"), width = 14, height = 8, units = "in", res = 300)
print(upset(fromList(cell_sets), nintersects = NA, order.by = "freq", nsets = length(cell_sets),
            text.scale = 1.4, mainbar.y.label = "Cells in intersection", sets.x.label = "Cells per QC"))
dev.off()


## =============================================================================
## CHUNK 5: PART A - transcripts in removed cells (red), per tissue section QCs 03 to 09
## =============================================================================
## Transcripts are read section by section from the ProSeg input files (script 09),
## which hold every transcript with its cell and position. Control probes and
## transcripts outside any cell are left out (cell QC does not apply to them).
## Each transcript takes the QC flags of its cell.
md_counts <- data.table(cell_id = metadata$cell_id, nCount_RNA = metadata$nCount_RNA)
cell_key <- cells_A[, .(cell_id, removed, key = gsub("; ", "&", QCs_removed_by))]  # UpSet key, e.g. "03_Low_counts&08_Low_SBR"

## Whole-slide PNG is drawn while sections are read, so all transcripts are never in memory at once
rng_x <- range(xy[, 1]); rng_y <- range(xy[, 2])
w_s <- 10 * diff(rng_x) / max(diff(rng_x), diff(rng_y)) + 0.5
h_s <- 10 * diff(rng_y) / max(diff(rng_x), diff(rng_y)) + 0.8
png(file.path(dir_A, "transcripts", "01_transcripts_in_removed_cells_on_slide.png"),
    width = w_s, height = h_s, units = "in", res = 300)
dev_slide <- dev.cur()
par(mar = c(0, 0, 2, 0))
plot(NA, xlim = rng_x, ylim = rng_y, asp = 1, axes = FALSE, xlab = "", ylab = "")

tx_sum <- list(); tx_keys <- list(); tx_check <- list()
for (sec in section_names) {
  f <- file.path(dir_ps_in, paste0(sec, "_tx_prepared.csv.gz"))
  stopifnot("ProSeg input file missing - run script 09 first" = file.exists(f))
  tx <- fread(f, select = c("seurat_cell", "x_global_px", "y_global_px", "target"))
  ctrl <- grep(control_regex, unique(tx$target), value = TRUE)
  tx <- tx[seurat_cell != "" & !target %chin% ctrl, .(seurat_cell, x = x_global_px * px_to_mm, y = y_global_px * px_to_mm)]
  tx[cell_key, on = c(seurat_cell = "cell_id"), `:=`(removed = i.removed, key = i.key)]
  stopifnot("transcripts with unknown cell" = !anyNA(tx$removed))
  
  ## Check: transcripts per cell should equal the cell's nCount_RNA
  n_per_cell <- tx[, .N, by = seurat_cell]
  n_per_cell[md_counts, on = c(seurat_cell = "cell_id"), nCount_RNA := i.nCount_RNA]
  tx_check[[sec]] <- data.table(section_id = sec, pct_cells_matching_nCount_RNA =
                                  round(100 * mean(n_per_cell$N == n_per_cell$nCount_RNA), 2))
  
  tx_sum[[sec]]  <- data.table(section_id = sec, n_transcripts = nrow(tx),
                               n_transcripts_removed = sum(tx$removed),
                               pct_transcripts_removed = round(100 * mean(tx$removed), 2),
                               n_transcripts_kept = sum(!tx$removed))
  tx_keys[[sec]] <- tx[removed == TRUE, .N, by = key]
  
  ## Per-section plot: grey = transcripts in kept cells, red = in removed cells (on top)
  r <- tx_sum[[sec]]
  png(file.path(dir_A, "transcripts", "per_tissue_section", paste0(sec, "_transcripts_in_removed_cells.png")),
      width = 9, height = 9, units = "in", res = 300)
  par(mar = c(0, 0, 2, 0))
  plot(tx$x[!tx$removed], tx$y[!tx$removed], pch = ".", col = "grey80", asp = 1, axes = FALSE, xlab = "", ylab = "",
       xlim = range(tx$x), ylim = range(tx$y))
  points(tx$x[tx$removed], tx$y[tx$removed], pch = ".", col = "red")
  title(paste0(sec, ": ", fmt(r$n_transcripts_removed), " of ", fmt(r$n_transcripts),
               " transcripts in removed cells (", r$pct_transcripts_removed, "%)"), cex.main = 0.9)
  dev.off()
  
  ## Add this section to the whole-slide plot
  dev.set(dev_slide)
  points(tx$x[!tx$removed], tx$y[!tx$removed], pch = ".", col = "grey80")
  points(tx$x[tx$removed],  tx$y[tx$removed],  pch = ".", col = "red")
  message(sec, ": ", fmt(r$n_transcripts_removed), " of ", fmt(r$n_transcripts), " transcripts in removed cells")
  rm(tx, n_per_cell); invisible(gc())
}
tx_summary <- rbindlist(tx_sum)
tx_summary <- rbind(tx_summary, tx_summary[, .(section_id = "All sections", n_transcripts = sum(n_transcripts),
                                               n_transcripts_removed = sum(n_transcripts_removed),
                                               pct_transcripts_removed = round(100 * sum(n_transcripts_removed) / sum(n_transcripts), 2),
                                               n_transcripts_kept = sum(n_transcripts_kept))])
tx_summary[, removed_label := paste0(fmt(n_transcripts_removed), " of ", fmt(n_transcripts), " transcripts")]

## Finish the whole-slide transcript plot: title and per-section counts
dev.set(dev_slide)
tt <- tx_summary[section_id == "All sections"]
title(paste0("Transcripts in removed cells, per tissue section QCs 03-09: ", tt$removed_label,
             " (", tt$pct_transcripts_removed, "%)"), cex.main = 0.9)
for (k in seq_len(nrow(sec_pos))) {
  r <- tx_summary[section_id == sec_pos$section_id[k]]
  text(sec_pos$x[k], sec_pos$y[k], pos = 3, xpd = NA, cex = 0.7, font = 2,
       labels = paste0(r$section_id, ": ", r$removed_label, " (", r$pct_transcripts_removed, "%)"))
}
dev.off()

print(tx_summary); print(rbindlist(tx_check))
fwrite(tx_summary,          file.path(dir_A, "transcripts", "02_transcripts_in_removed_cells_summary_per_section.csv"))
fwrite(rbindlist(tx_check), file.path(dir_A, "transcripts", "03_check_transcripts_per_cell_vs_nCount_RNA.csv"))

## UpSet plot of transcripts: intersection sizes = transcripts in cells with that QC combination
key_counts <- rbindlist(tx_keys)[, .(N = sum(N)), by = key]
fwrite(key_counts[order(-N)], file.path(dir_A, "transcripts", "04_transcripts_per_QC_combination.csv"))
png(file.path(dir_A, "transcripts", "05_upset_transcripts_in_removed_cells.png"), width = 14, height = 8, units = "in", res = 300)
print(upset(fromExpression(setNames(key_counts$N, key_counts$key)), nintersects = NA, order.by = "freq",
            nsets = uniqueN(unlist(strsplit(key_counts$key, "&"))), text.scale = 1.4,
            mainbar.y.label = "Transcripts in intersection", sets.x.label = "Transcripts per QC"))
dev.off()


## =============================================================================
## CHUNK 6: PART B - cells, all-cells analysis, QCs 03 to 08
## =============================================================================
qcs_B   <- names(flags_B)
cells_B <- build_cell_table(flags_B)
sum_B   <- summarise_cells(cells_B, qcs_B)
print(sum_B)

fwrite(sum_B,                    file.path(dir_B, "01_removed_cells_summary_per_section.csv"))
fwrite(cells_B,                  file.path(dir_B, "02_all_cells_with_QCs_removed_by.csv"))
fwrite(cells_B[removed == TRUE], file.path(dir_B, "03_removed_cells_with_QCs_removed_by.csv"))
saveRDS(setNames(cells_B$removed, cells_B$cell_id), file.path(dir_B, "final_QC_flag_all_cells_03_to_08.RDS"))

show_and_save(function() draw_cells_slide(cells_B, sum_B, "Removed cells, all-cells QCs 03-08"),
              file.path(dir_B, "04_removed_cells_on_slide.png"), xy)
png(file.path(dir_B, "05_removed_cells_per_section_grid.png"), width = 15, height = 10, units = "in", res = 300)
draw_cells_grid(cells_B, sum_B, "Removed cells, all-cells QCs 03-08"); dev.off()
save_cells_per_section(cells_B, sum_B, file.path(dir_B, "per_tissue_section"), "QC03_08_all_cells")

## UpSet plot of cells (each cell in every QC set that flags it)
cell_sets_B <- lapply(qcs_B, function(q) cells_B$cell_id[cells_B[[q]]]); names(cell_sets_B) <- qcs_B
cell_sets_B <- Filter(function(x) length(x) > 0, cell_sets_B)
png(file.path(dir_B, "06_upset_removed_cells.png"), width = 14, height = 8, units = "in", res = 300)
print(upset(fromList(cell_sets_B), nintersects = NA, order.by = "freq", nsets = length(cell_sets_B),
            text.scale = 1.4, mainbar.y.label = "Cells in intersection", sets.x.label = "Cells per QC"))
dev.off()


## =============================================================================
## CHUNK 7: results on screen
## =============================================================================
message("\nPART A - per tissue section, QCs 03-09 (cells):")
for (sec in c(section_names, "All sections")) {
  r <- sum_A[section_id == sec]; message("  ", sec, ": ", r$removed_label, " (", r$pct_removed, "%)")
}
message("\nPART A - transcripts in removed cells:")
for (sec in c(section_names, "All sections")) {
  r <- tx_summary[section_id == sec]; message("  ", sec, ": ", r$removed_label, " (", r$pct_transcripts_removed, "%)")
}
message("\nPART B - all-cells analysis, QCs 03-08 (cells):")
for (sec in c(section_names, "All sections")) {
  r <- sum_B[section_id == sec]; message("  ", sec, ": ", r$removed_label, " (", r$pct_removed, "%)")
}
message("\nPART A saved to ", dir_A, "\nPART B saved to ", dir_B)








## =============================================================================
## 11_add_condition_and_filter_seurat.R same as 11 above but i get to choose 03, ..., 09
## =============================================================================
## PURPOSE
##   1. Original Seurat object + 2 new metadata columns, saved under a NEW name
##      (seurat_with_sections.RDS itself is never overwritten or updated):
##        condition = "AAA" for every cell (TAA added later)
##        sample_id = AAA section of each cell (AAA2, AAA38, AAA40, AAA45,
##                    AAA46, AAA48), copied from section_id used in scripts 02-10
##   2. A new, QC-filtered Seurat object: the same counts and the same metadata
##      rows, for the cells KEPT after removing the cells flagged by the QCs YOU
##      CHOOSE (any of 03, 04, 05, 06, 07, 08, 09; per-tissue-section flags from
##      script 10, Part A). A cell is removed if ANY chosen QC flags it.
##      Each choice is saved under its own name, e.g. ..._removed_QC03_04_06.RDS,
##      so different choices never overwrite each other.
##
## WHY NOT subset()
##   The object was saved with an older SeuratObject. subset() and
##   UpdateSeuratObject() both fail on its spatial part (cell outlines), with
##   "protection stack overflow". So, as in the earlier filtering script, the
##   kept cells are selected directly from the counts and the metadata, and a new
##   Seurat object is built from them. Nothing is recalculated: counts and
##   metadata values are copied exactly (checked). Cell centroids are added back
##   as the spatial part; cell positions are also in metadata (x_slide_mm,
##   y_slide_mm), which is what every plot so far has used.
##
## Run in a fresh R session, chunk by chunk.
## =============================================================================


## =============================================================================
## CHUNK 1: packages and paths
## =============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(data.table)
})

output_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"
seurat_in  <- file.path(output_dir, "seurat_with_sections.RDS")                     # script 02 (read only)
qc_dir_A   <- file.path(output_dir, "QC", "independent_cells_per_tissue_section_03_to_09_QCs")
flag_table <- file.path(qc_dir_A, "cells", "02_all_cells_with_QCs_removed_by.csv")     # script 10: TRUE/FALSE per QC per cell
flag_csv   <- file.path(qc_dir_A, "cells", "01_removed_cells_summary_per_section.csv")  # script 10: counts (for the 03-09 check)

## >>>>> CHOOSE WHICH QCs TO REMOVE CELLS BY <<<<<
##   "03" Low counts        "04" Large area        "05" FOV QC
##   "06" nCount/nFeature   "07" FOV border        "08" Low SBR
##   "09" ProSeg low purity
qcs_to_remove <- c("03", "04", "05", "07")        # e.g. later: c("03", "04", "05", "09"); all: c("03","04","05","06","07","08","09")

qc_columns <- c(`03` = "03_Low_counts", `04` = "04_Large_area", `05` = "05_FOV_QC",
                `06` = "06_nCount_nFeature", `07` = "07_FOV_border", `08` = "08_Low_SBR",
                `09` = "09_ProSeg_low_purity")             # column names in script 10's per-cell table
stopifnot("unknown QC code in qcs_to_remove" = all(qcs_to_remove %in% names(qc_columns)))
qcs_to_remove <- sort(unique(qcs_to_remove))
qc_tag <- paste0("QC", paste(qcs_to_remove, collapse = "_"))                        # e.g. "QC03_04_06"
seurat_out_all      <- file.path(output_dir, "seurat_with_sections_condition.RDS")          # all cells + 2 columns
dir_filtered        <- file.path(output_dir, "QC_filtered_objects", qc_tag)                 # one folder per choice
seurat_out_filtered <- file.path(dir_filtered, paste0("seurat_QC_filtered_removed_", qc_tag, ".RDS"))
dir.create(dir_filtered, showWarnings = FALSE, recursive = TRUE)
fmt <- function(n) format(n, big.mark = ",")


## =============================================================================
## CHUNK 2: load the original object, add the 2 columns, save under a new name
## =============================================================================
obj <- readRDS(seurat_in)
metadata_qc <- readRDS(file.path(output_dir, "metadata_unfiltered.RDS"))   # metadata used by scripts 03-10

## Same cells, same order, same section labels as every QC script used
stopifnot("cells differ from metadata_unfiltered" = identical(colnames(obj), metadata_qc$cell_id),
          "section labels differ from metadata_unfiltered" =
            identical(unname(as.character(obj$section_id)), as.character(metadata_qc$section_id)))

obj$condition <- "AAA"                                   # all 6 sections are AAA
obj$sample_id <- unname(as.character(obj$section_id))    # which AAA section each cell is from
stopifnot(!anyNA(obj$condition), !anyNA(obj$sample_id),
          setequal(unique(obj$sample_id), c("AAA2", "AAA38", "AAA40", "AAA45", "AAA46", "AAA48")))
table(obj$condition, obj$sample_id)                      # cells per condition x section

if (!file.exists(seurat_out_all)) {                     # written once; same for every QC choice
  saveRDS(obj, seurat_out_all)                           # original + 2 columns, nothing else changed
  message("Saved (all cells + condition, sample_id): ", seurat_out_all)
}


## =============================================================================
## CHUNK 3: which cells to remove (the chosen per-tissue-section QCs)
## =============================================================================
## Script 10's per-cell table: one TRUE/FALSE column per QC (per-section results)
qc_tab <- fread(flag_table, select = c("cell_id", "section_id", unname(qc_columns)))
stopifnot(setequal(qc_tab$cell_id, colnames(obj)))
qc_tab <- qc_tab[match(colnames(obj), cell_id)]          # align rows to the Seurat cells
stopifnot(identical(qc_tab$cell_id, colnames(obj)),
          identical(qc_tab$section_id, unname(obj$sample_id)))

chosen  <- unname(qc_columns[qcs_to_remove])
flagmat <- as.matrix(qc_tab[, ..chosen])                 # cells x chosen QCs
removed <- rowSums(flagmat) > 0                          # removed if ANY chosen QC flags it (each cell once)
keep_cells <- colnames(obj)[!removed]                    # kept cells, original order

## Per section: cells flagged by each chosen QC, removed (once), kept
check <- data.table(section_id = qc_tab$section_id, removed = removed, flagmat)[
  , c(list(n_cells = .N), lapply(.SD, sum), list(n_removed = sum(removed), n_kept = sum(!removed))),
  by = section_id, .SDcols = chosen][order(section_id)]
print(check)

## Each chosen QC's count per section must equal script 10's summary
s10 <- fread(flag_csv)[section_id != "All sections"]
for (q in chosen) {
  s10_q <- setNames(s10[[paste0("n_", q)]], s10$section_id)[check$section_id]
  stopifnot("QC counts differ from script 10" = all(check[[q]] == s10_q))
}
## With all 7 QCs chosen, removed cells must equal script 10's total
if (length(chosen) == length(qc_columns))
  stopifnot(all(check$n_removed == setNames(s10$n_removed_independent, s10$section_id)[check$section_id]))


## =============================================================================
## CHUNK 4: build the filtered object from the kept cells' counts and metadata
## =============================================================================
get_counts <- function(o, assay) {                       # same accessor as script 02
  tryCatch(LayerData(o, assay = assay, layer = "counts"),
           error = function(e) GetAssayData(o, assay = assay, slot = "counts"))
}
md_keep <- obj[[]][keep_cells, ]                         # metadata rows of kept cells (all columns)

## RNA counts -> new object; metadata copied in
obj_filtered <- CreateSeuratObject(counts = get_counts(obj, "RNA")[, keep_cells],
                                   assay = "RNA", meta.data = md_keep)

## Control assays (negative probes, false codes), same kept cells
for (a in setdiff(Assays(obj), "RNA")) {
  obj_filtered[[a]] <- CreateAssayObject(counts = get_counts(obj, a)[, keep_cells])
}

## Spatial part: cell centroids of the kept cells (same coordinates as the original)
img_name  <- Images(obj)[1]                                              # "AC20260629_AAA3"
centroids <- GetTissueCoordinates(obj[[img_name]], which = "centroids")
cent_keep <- centroids[match(keep_cells, centroids$cell), c("x", "y", "cell")]
stopifnot(!anyNA(cent_keep$cell))
obj_filtered[[img_name]] <- CreateFOV(coords = CreateCentroids(cent_keep), type = "centroids",
                                      assay = "RNA", key = "fov_")

## Metadata: put back the original rows exactly (CreateSeuratObject may re-type
## columns such as orig.ident); same cells, same order
stopifnot(identical(rownames(md_keep), colnames(obj_filtered)))
obj_filtered@meta.data <- md_keep


## =============================================================================
## CHUNK 5: checks - the filtered object must match the original for every kept cell
## =============================================================================
same_counts <- function(A, B) identical(dim(A), dim(B)) && identical(dimnames(A), dimnames(B)) &&
  max(abs(A - B)) == 0                                   # every count identical
stopifnot(
  "cells differ"       = identical(colnames(obj_filtered), keep_cells),
  "genes differ"       = identical(rownames(obj_filtered), rownames(obj)),
  "assays differ"      = setequal(Assays(obj_filtered), Assays(obj)),
  "RNA counts differ"  = same_counts(get_counts(obj_filtered, "RNA"), get_counts(obj, "RNA")[, keep_cells]),
  "metadata differs"   = identical(obj_filtered[[]], md_keep),
  "nCount_RNA differs" = all(Matrix::colSums(get_counts(obj_filtered, "RNA")) == md_keep$nCount_RNA),
  "centroids differ"   = isTRUE(all.equal(GetTissueCoordinates(obj_filtered[[img_name]], which = "centroids")[, c("x", "y")],
                                          cent_keep[, c("x", "y")], check.attributes = FALSE)))
for (a in setdiff(Assays(obj), "RNA"))
  stopifnot(max(abs(get_counts(obj_filtered, a) - get_counts(obj, a)[, keep_cells])) == 0)

after <- as.data.table(table(sample_id = obj_filtered$sample_id))
setnames(after, "N", "n_cells_after_QC")
after <- merge(after, check[, .(sample_id = section_id, n_cells_before_QC = n_cells, n_removed)], by = "sample_id")
after[, QCs_removed := paste(chosen, collapse = "; ")]
stopifnot(all(after$n_cells_after_QC == after$n_cells_before_QC - after$n_removed))
print(after)
obj_filtered


## =============================================================================
## CHUNK 6: save
## =============================================================================
saveRDS(obj_filtered, seurat_out_filtered)
fwrite(merge(after, check[, !c("n_cells", "n_removed", "n_kept")], by.x = "sample_id", by.y = "section_id"),
       file.path(dir_filtered, paste0("cells_before_and_after_", qc_tag, "_per_section.csv")))
message("Saved (removed cells flagged by ", paste(chosen, collapse = ", "), "): ", seurat_out_filtered)
message("Cells: ", fmt(ncol(obj)), " before, ", fmt(sum(removed)), " removed, ", fmt(ncol(obj_filtered)), " kept")

## Downstream examples (not run):
# obj_filtered <- readRDS(seurat_out_filtered)
# sub <- subset(obj_filtered, subset = sample_id %in% c("AAA48", "AAA38"))   # chosen AAA sections
# Idents(obj_filtered) <- "sample_id"                                         # compare sections









## =============================================================================
## 11b_QC_flags_for_chosen_QCs.R
## =============================================================================
## PURPOSE
##   Script 10 (Part A) outputs, for ONLY the QCs chosen in script 11
##   (qcs_to_remove). A cell is removed if ANY chosen QC flags it; each cell is
##   counted ONCE.
##   Outputs, one folder per choice, inside script 11's folder for that choice:
##     QC_filtered_objects/<QCtag>/QC_flags/
##       cells/        01 summary per section, 02 all cells with the QCs that flagged
##                     them, 03 removed cells only (CSVs); 04 slide plot, 05 2 x 3 grid,
##                     per_tissue_section/ one plot per section, 06 UpSet plot of cells
##       transcripts/  01 slide plot, per_tissue_section/ plots, 02 summary CSV,
##                     03 transcripts-per-cell vs nCount_RNA check, 04 transcripts per
##                     QC combination, 05 UpSet plot of transcripts
##       final_QC_flag_per_tissue_section_<QCtag>.RDS   (TRUE = removed, named by cell_id)
##
##   Flags are read from script 10's per-cell table (02_all_cells_with_QCs_removed_by.csv),
##   the same table script 11 filters on, so both scripts remove the same cells.
##   Helper functions and plots are copied unchanged from script 10.
##   Nothing is removed from the data here.
##
## Run after scripts 09-10, with the SAME qcs_to_remove as script 11.
## Fresh R session, chunk by chunk.
## =============================================================================


## =============================================================================
## CHUNK 1: packages, paths, QC choice, folders
## =============================================================================
library(data.table)
library(UpSetR)       # upset(), fromList(), fromExpression()

output_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"
qc_dir     <- file.path(output_dir, "QC")
qc_dir_A   <- file.path(qc_dir, "independent_cells_per_tissue_section_03_to_09_QCs")   # script 10, Part A
flag_table <- file.path(qc_dir_A, "cells", "02_all_cells_with_QCs_removed_by.csv")     # TRUE/FALSE per QC per cell

## >>>>> CHOOSE WHICH QCs (must match qcs_to_remove in script 11) <<<<<
##   "03" Low counts        "04" Large area        "05" FOV QC
##   "06" nCount/nFeature   "07" FOV border        "08" Low SBR
##   "09" ProSeg low purity
qcs_to_remove <- c("03", "04", "05", "07")

qc_columns <- c(`03` = "03_Low_counts", `04` = "04_Large_area", `05` = "05_FOV_QC",
                `06` = "06_nCount_nFeature", `07` = "07_FOV_border", `08` = "08_Low_SBR",
                `09` = "09_ProSeg_low_purity")                       # column names in script 10's per-cell table
stopifnot("unknown QC code in qcs_to_remove" = all(qcs_to_remove %in% names(qc_columns)))
qcs_to_remove <- sort(unique(qcs_to_remove))
qc_tag   <- paste0("QC", paste(qcs_to_remove, collapse = "_"))       # e.g. "QC03_04_05_07" (same as script 11)
qc_label <- paste0("QCs ", paste(qcs_to_remove, collapse = ", "))    # plot titles, e.g. "QCs 03, 04, 05, 07"
chosen   <- unname(qc_columns[qcs_to_remove])                        # chosen column names

dir_filtered <- file.path(output_dir, "QC_filtered_objects", qc_tag) # script 11's folder for this choice
dir_C <- file.path(dir_filtered, "QC_flags")                         # this script's outputs
for (d in c(file.path(dir_C, "cells", "per_tissue_section"), file.path(dir_C, "transcripts", "per_tissue_section")))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

metadata <- readRDS(file.path(output_dir, "metadata_unfiltered.RDS"))
xy       <- readRDS(file.path(output_dir, "xy_unfiltered.RDS"))      # x_mm, y_mm
stopifnot(identical(rownames(metadata), metadata$cell_id),
          identical(rownames(xy), metadata$cell_id))
section_names <- sort(unique(metadata$section_id))

## ProSeg input files (per section) hold every transcript with its cell and position
dir_ps_in <- file.path(qc_dir, "QC_Proseg", "per_tissue_section", "01_prepared_transcripts")
control_regex <- "^(FalseCode|SystemControl|NegPrb|Negative)"       # control probes: not RNA
px_to_mm      <- 0.120280945 / 1000                                 # same conversion as script 02

fmt <- function(n) format(n, big.mark = ",")


## =============================================================================
## CHUNK 2: load the chosen flags from script 10's per-cell table
## =============================================================================
qc_tab <- fread(flag_table, select = c("cell_id", chosen))           # cell_id + chosen QC columns only
stopifnot(setequal(qc_tab$cell_id, metadata$cell_id))
qc_tab <- qc_tab[match(metadata$cell_id, cell_id)]                  # align rows to metadata order
stopifnot(identical(qc_tab$cell_id, metadata$cell_id))

flags_C <- lapply(setNames(chosen, chosen), function(q) qc_tab[[q]]) # list: one TRUE/FALSE vector per chosen QC
stopifnot("flags must be TRUE/FALSE with no NA" =
            all(vapply(flags_C, is.logical, logical(1))) && !anyNA(unlist(flags_C)))


## =============================================================================
## CHUNK 3: helper functions (unchanged from script 10)
## =============================================================================
## Per-cell table: one TRUE/FALSE column per QC, number of QCs, the QCs listed
## by name, and removed = flagged by at least one QC (each cell once).
build_cell_table <- function(fl) {
  d <- data.table(cell_id = metadata$cell_id, section_id = metadata$section_id, fov = metadata$fov,
                  x_slide_mm = metadata$x_slide_mm, y_slide_mm = metadata$y_slide_mm,
                  nCount_RNA = metadata$nCount_RNA)
  for (q in names(fl)) set(d, j = q, value = fl[[q]])
  d[, n_QCs := rowSums(.SD), .SDcols = names(fl)]
  removed_by <- rep("", nrow(d))                        # e.g. "03_Low_counts; 07_FOV_border"
  for (q in names(fl)) removed_by <- ifelse(fl[[q]], ifelse(removed_by == "", q, paste(removed_by, q, sep = "; ")), removed_by)
  d[, QCs_removed_by := removed_by]
  d[, removed := n_QCs > 0]
  d
}

## Per-section summary: cells, cells flagged by each QC, independent removed, kept
summarise_cells <- function(d, qcs) {
  one <- function(x, label) {
    r <- data.table(section_id = label, n_cells = nrow(x))
    for (q in qcs) set(r, j = paste0("n_", q), value = sum(x[[q]]))
    r[, `:=`(total_flag_events = sum(x$n_QCs),                     # a cell counted once per QC
             n_removed_independent = sum(x$removed),              # each cell once
             pct_removed = round(100 * mean(x$removed), 2),
             n_kept = sum(!x$removed),
             n_removed_by_1_QC = sum(x$n_QCs == 1),
             n_removed_by_2_or_more_QCs = sum(x$n_QCs >= 2))]
    r[, removed_label := paste0(fmt(n_removed_independent), " of ", fmt(n_cells), " cells removed")]
    r
  }
  rbind(rbindlist(lapply(section_names, function(s) one(d[section_id == s], s))), one(d, "All sections"))
}

## Save PNG at slide shape, then show in a separate window
show_and_save <- function(draw_fun, file, pts, extra_h = 0.8) {
  rx <- diff(range(pts[, 1])); ry <- diff(range(pts[, 2]))
  w <- 10 * rx / max(rx, ry) + 0.5; h <- 10 * ry / max(rx, ry) + extra_h
  png(file, width = w, height = h, units = "in", res = 300); draw_fun(); dev.off()
  dev.new(width = w, height = h, noRStudioGD = TRUE); draw_fun()
  message("Saved: ", file)
}

sec_pos <- as.data.table(metadata)[, .(x = mean(range(x_slide_mm)), y = max(y_slide_mm)), by = section_id]

## Whole slide: grey80 = kept, red = removed (drawn on top); counts above each section
draw_cells_slide <- function(d, s, title) {
  o <- order(d$removed)
  par(mar = c(0, 0, 2, 0))
  plot(xy[o, ], pch = 16, cex = 0.1, asp = 1, col = c("grey80", "red")[1 + d$removed[o]])
  tot <- s[section_id == "All sections"]
  title(paste0(title, ": ", tot$removed_label, " (", tot$pct_removed, "%)"), cex.main = 0.9)
  for (k in seq_len(nrow(sec_pos))) {
    r <- s[section_id == sec_pos$section_id[k]]
    text(sec_pos$x[k], sec_pos$y[k], pos = 3, xpd = NA, cex = 0.7, font = 2,
         labels = paste0(r$section_id, ": ", r$removed_label, " (", r$pct_removed, "%)"))
  }
}

## 2 x 3 grid, one panel per section
draw_cells_grid <- function(d, s, title) {
  op <- par(mfrow = c(2, 3), mar = c(0.5, 0.5, 2.5, 0.5), oma = c(0, 0, 2.5, 0)); on.exit(par(op))
  for (sec in section_names) {
    idx <- which(d$section_id == sec); o <- idx[order(d$removed[idx])]
    plot(xy[o, ], pch = 16, cex = 0.1, asp = 1, axes = FALSE, xlab = "", ylab = "",
         col = c("grey80", "red")[1 + d$removed[o]])
    r <- s[section_id == sec]
    title(paste0(sec, ": ", r$removed_label, " (", r$pct_removed, "%)"), cex.main = 1)
  }
  tot <- s[section_id == "All sections"]
  mtext(paste0(title, ": ", tot$removed_label, " (", tot$pct_removed, "%)"), outer = TRUE, font = 2)
}

## One plot per section
save_cells_per_section <- function(d, s, out_dir, tag) {
  for (sec in section_names) {
    idx <- which(d$section_id == sec); o <- idx[order(d$removed[idx])]
    r <- s[section_id == sec]
    png(file.path(out_dir, paste0(sec, "_removed_cells_", tag, ".png")), width = 8, height = 8, units = "in", res = 300)
    par(mar = c(0, 0, 2, 0))
    plot(xy[o, ], pch = 16, cex = 0.15, asp = 1, col = c("grey80", "red")[1 + d$removed[o]])
    title(paste0(sec, " (", tag, "): ", r$removed_label, " (", r$pct_removed, "%)"), cex.main = 0.9)
    dev.off()
  }
}


## =============================================================================
## CHUNK 4: cells flagged by the chosen QCs, per tissue section
## =============================================================================
cells_C <- build_cell_table(flags_C)
sum_C   <- summarise_cells(cells_C, chosen)
print(sum_C)

fwrite(sum_C,                    file.path(dir_C, "cells", "01_removed_cells_summary_per_section.csv"))
fwrite(cells_C,                  file.path(dir_C, "cells", "02_all_cells_with_QCs_removed_by.csv"))
fwrite(cells_C[removed == TRUE], file.path(dir_C, "cells", "03_removed_cells_with_QCs_removed_by.csv"))
saveRDS(setNames(cells_C$removed, cells_C$cell_id),
        file.path(dir_C, paste0("final_QC_flag_per_tissue_section_", qc_tag, ".RDS")))

## Check 1: removed per section must equal script 11's numbers for the same choice (if script 11 has run)
s11_csv <- file.path(dir_filtered, paste0("cells_before_and_after_", qc_tag, "_per_section.csv"))
if (file.exists(s11_csv)) {
  s11 <- fread(s11_csv)
  cmp <- merge(sum_C[section_id != "All sections", .(sample_id = section_id, n_removed_independent)],
               s11[, .(sample_id, n_removed)], by = "sample_id")
  stopifnot("removed cells differ from script 11" = all(cmp$n_removed_independent == cmp$n_removed))
  message("Check passed: removed cells per section match script 11.")
} else message("Script 11 not run yet for ", qc_tag, ": comparison skipped.")

## Check 2: with all 7 QCs chosen, the flag must equal script 10's final flag cell for cell
if (length(chosen) == length(qc_columns)) {
  f10 <- readRDS(file.path(qc_dir_A, "final_QC_flag_per_tissue_section_03_to_09.RDS"))
  stopifnot("flag differs from script 10" = identical(unname(f10[cells_C$cell_id]), cells_C$removed))
}

## Plots: whole slide, grid, one per section
show_and_save(function() draw_cells_slide(cells_C, sum_C, paste0("Removed cells, per tissue section ", qc_label)),
              file.path(dir_C, "cells", "04_removed_cells_on_slide.png"), xy)
png(file.path(dir_C, "cells", "05_removed_cells_per_section_grid.png"), width = 15, height = 10, units = "in", res = 300)
draw_cells_grid(cells_C, sum_C, paste0("Removed cells, per tissue section ", qc_label)); dev.off()
save_cells_per_section(cells_C, sum_C, file.path(dir_C, "cells", "per_tissue_section"), paste0(qc_tag, "_per_section"))

## UpSet plot of cells (each cell in every QC set that flags it)
## UpSetR needs at least 2 non-empty sets, so it is skipped for a single-QC choice.
cell_sets <- lapply(chosen, function(q) cells_C$cell_id[cells_C[[q]]]); names(cell_sets) <- chosen
cell_sets <- Filter(function(x) length(x) > 0, cell_sets)
if (length(cell_sets) >= 2) {
  png(file.path(dir_C, "cells", "06_upset_removed_cells.png"), width = 14, height = 8, units = "in", res = 300)
  print(upset(fromList(cell_sets), nintersects = NA, order.by = "freq", nsets = length(cell_sets),
              text.scale = 1.4, mainbar.y.label = "Cells in intersection", sets.x.label = "Cells per QC"))
  dev.off()
} else message("UpSet plot of cells skipped: fewer than 2 QCs flag any cell.")


## =============================================================================
## CHUNK 5: transcripts in removed cells (red), chosen QCs
## =============================================================================
## Transcripts are read section by section from the ProSeg input files (script 09),
## which hold every transcript with its cell and position. Control probes and
## transcripts outside any cell are left out (cell QC does not apply to them).
## Each transcript takes the QC flags of its cell.
md_counts <- data.table(cell_id = metadata$cell_id, nCount_RNA = metadata$nCount_RNA)
cell_key <- cells_C[, .(cell_id, removed, key = gsub("; ", "&", QCs_removed_by))]  # UpSet key, e.g. "03_Low_counts&07_FOV_border"

## Whole-slide PNG is drawn while sections are read, so all transcripts are never in memory at once
rng_x <- range(xy[, 1]); rng_y <- range(xy[, 2])
w_s <- 10 * diff(rng_x) / max(diff(rng_x), diff(rng_y)) + 0.5
h_s <- 10 * diff(rng_y) / max(diff(rng_x), diff(rng_y)) + 0.8
png(file.path(dir_C, "transcripts", "01_transcripts_in_removed_cells_on_slide.png"),
    width = w_s, height = h_s, units = "in", res = 300)
dev_slide <- dev.cur()
par(mar = c(0, 0, 2, 0))
plot(NA, xlim = rng_x, ylim = rng_y, asp = 1, axes = FALSE, xlab = "", ylab = "")

tx_sum <- list(); tx_keys <- list(); tx_check <- list()
for (sec in section_names) {
  f <- file.path(dir_ps_in, paste0(sec, "_tx_prepared.csv.gz"))
  stopifnot("ProSeg input file missing - run script 09 first" = file.exists(f))
  tx <- fread(f, select = c("seurat_cell", "x_global_px", "y_global_px", "target"))
  ctrl <- grep(control_regex, unique(tx$target), value = TRUE)
  tx <- tx[seurat_cell != "" & !target %chin% ctrl, .(seurat_cell, x = x_global_px * px_to_mm, y = y_global_px * px_to_mm)]
  tx[cell_key, on = c(seurat_cell = "cell_id"), `:=`(removed = i.removed, key = i.key)]
  stopifnot("transcripts with unknown cell" = !anyNA(tx$removed))
  
  ## Check: transcripts per cell should equal the cell's nCount_RNA
  n_per_cell <- tx[, .N, by = seurat_cell]
  n_per_cell[md_counts, on = c(seurat_cell = "cell_id"), nCount_RNA := i.nCount_RNA]
  tx_check[[sec]] <- data.table(section_id = sec, pct_cells_matching_nCount_RNA =
                                  round(100 * mean(n_per_cell$N == n_per_cell$nCount_RNA), 2))
  
  tx_sum[[sec]]  <- data.table(section_id = sec, n_transcripts = nrow(tx),
                               n_transcripts_removed = sum(tx$removed),
                               pct_transcripts_removed = round(100 * mean(tx$removed), 2),
                               n_transcripts_kept = sum(!tx$removed))
  tx_keys[[sec]] <- tx[removed == TRUE, .N, by = key]
  
  ## Per-section plot: grey = transcripts in kept cells, red = in removed cells (on top)
  r <- tx_sum[[sec]]
  png(file.path(dir_C, "transcripts", "per_tissue_section", paste0(sec, "_transcripts_in_removed_cells_", qc_tag, ".png")),
      width = 9, height = 9, units = "in", res = 300)
  par(mar = c(0, 0, 2, 0))
  plot(tx$x[!tx$removed], tx$y[!tx$removed], pch = ".", col = "grey80", asp = 1, axes = FALSE, xlab = "", ylab = "",
       xlim = range(tx$x), ylim = range(tx$y))
  points(tx$x[tx$removed], tx$y[tx$removed], pch = ".", col = "red")
  title(paste0(sec, " (", qc_label, "): ", fmt(r$n_transcripts_removed), " of ", fmt(r$n_transcripts),
               " transcripts in removed cells (", r$pct_transcripts_removed, "%)"), cex.main = 0.9)
  dev.off()
  
  ## Add this section to the whole-slide plot
  dev.set(dev_slide)
  points(tx$x[!tx$removed], tx$y[!tx$removed], pch = ".", col = "grey80")
  points(tx$x[tx$removed],  tx$y[tx$removed],  pch = ".", col = "red")
  message(sec, ": ", fmt(r$n_transcripts_removed), " of ", fmt(r$n_transcripts), " transcripts in removed cells")
  rm(tx, n_per_cell); invisible(gc())
}
tx_summary <- rbindlist(tx_sum)
tx_summary <- rbind(tx_summary, tx_summary[, .(section_id = "All sections", n_transcripts = sum(n_transcripts),
                                               n_transcripts_removed = sum(n_transcripts_removed),
                                               pct_transcripts_removed = round(100 * sum(n_transcripts_removed) / sum(n_transcripts), 2),
                                               n_transcripts_kept = sum(n_transcripts_kept))])
tx_summary[, removed_label := paste0(fmt(n_transcripts_removed), " of ", fmt(n_transcripts), " transcripts")]

## Finish the whole-slide transcript plot: title and per-section counts
dev.set(dev_slide)
tt <- tx_summary[section_id == "All sections"]
title(paste0("Transcripts in removed cells, per tissue section ", qc_label, ": ", tt$removed_label,
             " (", tt$pct_transcripts_removed, "%)"), cex.main = 0.9)
for (k in seq_len(nrow(sec_pos))) {
  r <- tx_summary[section_id == sec_pos$section_id[k]]
  text(sec_pos$x[k], sec_pos$y[k], pos = 3, xpd = NA, cex = 0.7, font = 2,
       labels = paste0(r$section_id, ": ", r$removed_label, " (", r$pct_transcripts_removed, "%)"))
}
dev.off()

print(tx_summary); print(rbindlist(tx_check))
fwrite(tx_summary,          file.path(dir_C, "transcripts", "02_transcripts_in_removed_cells_summary_per_section.csv"))
fwrite(rbindlist(tx_check), file.path(dir_C, "transcripts", "03_check_transcripts_per_cell_vs_nCount_RNA.csv"))

## UpSet plot of transcripts: intersection sizes = transcripts in cells with that QC combination
key_counts <- rbindlist(tx_keys)[, .(N = sum(N)), by = key]
fwrite(key_counts[order(-N)], file.path(dir_C, "transcripts", "04_transcripts_per_QC_combination.csv"))
n_tx_sets <- uniqueN(unlist(strsplit(key_counts$key, "&")))
if (n_tx_sets >= 2) {
  png(file.path(dir_C, "transcripts", "05_upset_transcripts_in_removed_cells.png"), width = 14, height = 8, units = "in", res = 300)
  print(upset(fromExpression(setNames(key_counts$N, key_counts$key)), nintersects = NA, order.by = "freq",
              nsets = n_tx_sets, text.scale = 1.4,
              mainbar.y.label = "Transcripts in intersection", sets.x.label = "Transcripts per QC"))
  dev.off()
} else message("UpSet plot of transcripts skipped: fewer than 2 QCs flag any cell.")


## =============================================================================
## CHUNK 6: results on screen
## =============================================================================
message("\n", qc_label, " - cells removed per tissue section:")
for (sec in c(section_names, "All sections")) {
  r <- sum_C[section_id == sec]; message("  ", sec, ": ", r$removed_label, " (", r$pct_removed, "%)")
}
message("\n", qc_label, " - transcripts in removed cells:")
for (sec in c(section_names, "All sections")) {
  r <- tx_summary[section_id == sec]; message("  ", sec, ": ", r$removed_label, " (", r$pct_transcripts_removed, "%)")
}
message("\nSaved to ", dir_C)






















