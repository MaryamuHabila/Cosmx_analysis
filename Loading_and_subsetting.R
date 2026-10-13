## Rule used throughout: one column    -> obj$column
##                       several columns -> obj[[c("col1", "col2")]]
##                       all metadata    -> obj[[]]
##                       a count matrix  -> LayerData(obj, assay = , layer = "counts")
## New per-cell values are written straight into obj: obj$new_column <- ...
##loads the toolkits;

suppressPackageStartupMessages({
  library(Seurat)        
  library(SeuratObject)  
  library(Matrix)        
  library(data.table)    
  library(dbscan)        
  library(ggplot2)  
  library(R.utils)  
})

## =============================================================================
## file locations and output folders
## =============================================================================
## all inputs and outputs
seurat_file <- "/Users/maryamuhabilausman/Desktop/AAA_CosMx_Analysis/Downloaded_data/seuratObject_AC20260629_AAA3.RDS"  # AtoMx Seurat export
slide_dir   <- "/Users/maryamuhabilausman/Desktop/AAA_CosMx_Analysis/Downloaded_data/flatFiles/AC20260629_AAA3"      # AtoMx flat files (FOV positions)
output_dir  <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA"                                                # all results go here

dir_plots    <- file.path(output_dir, "spatial_plots")   # builds the path text ".../AAA/spatial_plots"
dir_sections <- file.path(output_dir, "per_section")     # parent folder for one sub-folder per section
dir.create(dir_plots,    recursive = TRUE, showWarnings = FALSE)  # creates the folder; recursive = parents too; silent if it exists
dir.create(dir_sections, recursive = TRUE, showWarnings = FALSE)


## =============================================================================
##load and inspect the Seurat object
## =============================================================================
##confirms the input exists and is a Seurat object, then shows what it holds
## before anything is computed.
stopifnot(file.exists(seurat_file))   # stop now if the path is wrong
obj <- readRDS(seurat_file)           # open the .RDS file and store the restored object as obj
stopifnot(inherits(obj, "Seurat"))    # every Seurat function below needs this class

## ---- 2a: the object as a whole ----
obj                    # 6526 features (6183 genes + 343 control probes), 123,709 cells, 3 assays, 1 FOV object
dim(obj)               # genes x cells of the default assay (6183 x 123,709)
Assays(obj)            # "RNA" (genes), "negprobes" and "falsecode" (background controls)
DefaultAssay(obj)      # "RNA": used whenever no assay is named
Layers(obj[["RNA"]])   # "counts" (raw) and "data"; this script uses counts only
Images(obj)            # "AC20260629_AAA3": spatial object holding centroids, polygons, transcripts
head(Cells(obj))       # cell names (= colnames(obj))
head(Features(obj))    # gene names (= rownames(obj))

## ---- 2b: per-cell metadata (one row per cell, one column per AtoMx variable) ----
dim(obj[[]])           # 123,709 cells x 66 columns
colnames(obj[[]])      # every metadata column name
head(obj[[]])          # first 6 cells. REDUNDANT: obj[[]][1:10, ] shows the same for 10 cells
str(obj[[]])           # type of each column (int, num, chr, logi)
summary(obj[[]])       # min/median/max per column; long, but one-off
# View(obj[[]])        # spreadsheet view; slow for 123k rows, so run only when needed
table(obj$fov)         # cells per FOV (267 FOVs)
summary(obj$Area.um2)  # cell area in µm²


## =============================================================================
##  structural validation
## =============================================================================
## stops the run if an assay or metadata column the workflow needs is
## missing, or if AtoMx cell IDs disagree with Seurat's cell names.
required_assays <- c("RNA", "negprobes", "falsecode")
missing_assays  <- setdiff(required_assays, Assays(obj))            # names required but not present
if (length(missing_assays)) stop("Missing assay(s): ", paste(missing_assays, collapse = ", "))

required_metadata <- c("cell_id", "fov", "x_slide_mm", "y_slide_mm", "Area.um2",
                       "nCount_RNA", "nFeature_RNA", "nCount_negprobes", "nCount_falsecode")
missing_metadata  <- setdiff(required_metadata, colnames(obj[[]]))  # columns required but not present
if (length(missing_metadata)) stop("Missing metadata column(s): ", paste(missing_metadata, collapse = ", "))

stopifnot(!anyDuplicated(obj$cell_id))                          # every cell ID is unique
stopifnot(identical(Cells(obj), unname(obj$cell_id)))           # cell_id = Seurat cell names, same order
stopifnot(identical(unname(obj$cell_id), unname(obj$cell_ID)))  # AtoMx's duplicate ID column agrees; only cell_id is used
uniqueN(obj$fov)                                                # EXPECT 267 FOVs


## =============================================================================
## count matrices (cells x features) and x/y coordinates
## =============================================================================
##Seurat stores counts as features x cells; t() transposes to cells x features, the orientation used by the NanoString CosMx workflow.
## Negative probes target sequences absent from the human genome; false codes are barcodes assigned to no probe. Both measure background, not expression.
counts           <- t(LayerData(obj, assay = "RNA",       layer = "counts"))  # cells x genes
negprobes_counts <- t(LayerData(obj, assay = "negprobes", layer = "counts"))  # cells x negative probes
falsecode_counts <- t(LayerData(obj, assay = "falsecode", layer = "counts"))  # cells x false codes
## ALTERNATIVE, same matrix: t(obj[["RNA"]]$counts)

dim(counts);           counts[1:5, 1:5]            # "." = zero in a sparse matrix
dim(negprobes_counts); negprobes_counts[1:5, 1:5]
dim(falsecode_counts); falsecode_counts[1:5, 1:5]
class(counts)                                      # dgCMatrix = compressed sparse matrix
1 - nnzero(counts) / prod(dim(counts))             # fraction of zeros: why sparse storage is used

## x_slide_mm / y_slide_mm are AtoMx's whole-slide coordinates, already in mm.
## (From pixels it would be px * 0.120280945 µm/px / 1000 = mm; not needed here.)
xy <- as.matrix(obj[[c("x_slide_mm", "y_slide_mm")]])  # two columns straight from obj; row names = cell names
colnames(xy) <- c("x_mm", "y_mm")
stopifnot(all(is.finite(xy)))                          # no NA/NaN/Inf before clustering or area maths
head(xy)
apply(xy, 2, range)                                    # slide extent in mm (y min ≈ -0.51 is valid)


## =============================================================================
##  integrity checks
## =============================================================================
## DEFENCE: every matrix lists the same cells in the same order as obj, and
## summing each cell's counts reproduces the totals AtoMx recorded.
stopifnot(identical(rownames(counts),           Cells(obj)))
stopifnot(identical(rownames(negprobes_counts), Cells(obj)))
stopifnot(identical(rownames(falsecode_counts), Cells(obj)))
stopifnot(identical(rownames(xy),               Cells(obj)))
stopifnot(all(rowSums(counts)           == obj$nCount_RNA))        # per-cell RNA total
stopifnot(all(rowSums(negprobes_counts) == obj$nCount_negprobes))  # per-cell negative-probe total
stopifnot(all(rowSums(falsecode_counts) == obj$nCount_falsecode))  # per-cell false-code total
message("Integrity checks passed: identities and per-cell totals match.")


## =============================================================================
## CHUNK 6: mean background per cell
## =============================================================================
## a control SUM scales with how many control probes the panel has; dividing by that number gives a per-probe MEAN, comparable across panels.
n_negprobes <- ncol(negprobes_counts)   # number of negative probes in this panel
n_falsecode <- ncol(falsecode_counts)   # number of false codes in this panel
n_negprobes; n_falsecode                # together they should equal 343

obj$negmean       <- obj$nCount_negprobes / n_negprobes   # new column written straight into obj
obj$falsecodemean <- obj$nCount_falsecode / n_falsecode
summary(obj$negmean)
summary(obj$falsecodemean)


## =============================================================================
## spatial components: centroids, polygons, transcripts, FOV positions
## =============================================================================
## the spatial data live in an FOV-class object inside obj.
img_name <- Images(obj)[1]   # name of the spatial object ("AC20260629_AAA3")
fov_obj  <- obj[[img_name]]  # the spatial object itself
names(fov_obj)               # EXPECT "centroids" "segmentation" "molecules"

## ---- 7a: centroids = one centre point per cell ----
centroids_df <- GetTissueCoordinates(fov_obj, which = "centroids")   # columns x, y, cell
head(centroids_df)
summary(centroids_df[, c("x", "y")])                                 # values < ~50, so mm
## Are the Seurat centroids the same numbers as AtoMx's x_slide_mm / y_slide_mm?
## TRUE justifies using the metadata coordinates (xy) for everything below.
isTRUE(all.equal(centroids_df$x, xy[centroids_df$cell, "x_mm"], check.attributes = FALSE)) &&
  isTRUE(all.equal(centroids_df$y, xy[centroids_df$cell, "y_mm"], check.attributes = FALSE))

## ---- 7b: polygons = segmentation outline of each cell (one row per vertex) ----
polygons_df <- NULL
if ("segmentation" %in% names(fov_obj)) {
  polygons_df <- GetTissueCoordinates(fov_obj, which = "segmentation")  # columns x, y, cell
  print(head(polygons_df, 10))
  print(nrow(polygons_df))                  # total vertices
  print(length(unique(polygons_df$cell)))   # cells with an outline; EXPECT 123,709
  print(mean(table(polygons_df$cell)))      # mean vertices per cell
}

## ---- 7c: transcripts = x, y of each detected RNA molecule ----
## The full table has hundreds of millions of rows, so load only the chosen gene(s).
genes_to_show <- intersect("AC024901.1", rownames(obj))   # keep genes that are in the panel
tx_df <- NULL
if ("molecules" %in% names(fov_obj)) {
  tx_df <- GetTissueCoordinates(fov_obj[["molecules"]], features = genes_to_show)  # columns x, y, molecule
  print(head(tx_df))
  print(table(tx_df$molecule))              # transcripts per chosen gene
}

## ---- 7d: FOV positions from the AtoMx flat file (for cross-checking only) ----
fovfile <- list.files(slide_dir, pattern = "fov_positions.*", full.names = TRUE)[1]
fov_positions_flat <- NULL
if (!is.na(fovfile)) {
  fov_positions_flat <- fread(fovfile)
  print(head(fov_positions_flat))           # one row per FOV and its slide position
}


## =============================================================================
## CHUNK 8: section labels (the one label AtoMx does not provide)
## =============================================================================
## DEFENCE: sections are defined from position only (no expression). DBSCAN
## chains cells < 1 mm apart (eps = 1 mm); minPts = 1 means every cell is
## assigned and none is noise. Cluster numbers are then mapped to tissue IDs by
## matching the plot to the AtoMx slide thumbnail.
if ("section_id" %in% colnames(obj[[]])) {
  message("section_id already present - reused, no clustering run.")
} else {
  clust_raw <- dbscan(xy, eps = 1, minPts = 1)$cluster      # one cluster number per cell, in Cells(obj) order
  obj$section_id_placeholder <- paste0("section", clust_raw)
  print(table(obj$section_id_placeholder))                  # cells per cluster
  if (length(unique(clust_raw)) != 6) warning("Expected 6 clusters - adjust eps.")
  
  ## Cluster centres: where each cluster sits on the slide.
  print(aggregate(cbind(x_slide_mm, y_slide_mm) ~ section_id_placeholder,
                  data = obj[[c("x_slide_mm", "y_slide_mm", "section_id_placeholder")]], FUN = mean))
  
  ## Plot to match clusters to the AtoMx slide thumbnail.
  p_placeholder <- ggplot(obj[[c("x_slide_mm", "y_slide_mm", "section_id_placeholder")]],
                          aes(x_slide_mm, y_slide_mm, colour = section_id_placeholder)) +
    geom_point(size = 0.05, alpha = 0.3) + coord_fixed() + theme_minimal() +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(title = "dbscan clusters - match to AtoMx thumbnail", x = "x (mm)", y = "y (mm)")
  ggsave(file.path(dir_plots, "placeholder_overview_plot.png"), p_placeholder, width = 8, height = 8, dpi = 300)
  print(p_placeholder)
  
  ## dbscan numbers clusters in the order it meets cells; it cannot know tissue IDs.
  ## >>> CONFIRM this mapping against the plot and the AtoMx thumbnail. <<<
  section_name_mapping <- c(section1 = "AAA38", section2 = "AAA45", section3 = "AAA46",
                            section4 = "AAA40", section5 = "AAA48", section6 = "AAA2")
  stopifnot(all(unique(obj$section_id_placeholder) %in% names(section_name_mapping)))
  obj$section_id <- unname(section_name_mapping[obj$section_id_placeholder])
}
table(obj$section_id)   # cells per named section


## =============================================================================
## CHUNK 9: tallies per FOV and per section
## =============================================================================
## DEFENCE: every number is a tally of existing AtoMx columns, except tissue
## area, which belongs to a whole section and must be computed.
## md holds only the columns the summaries need, taken straight from obj.
md <- as.data.table(obj[[c("cell_id", "section_id", "fov", "x_slide_mm", "y_slide_mm",
                           "Area.um2", "nCount_RNA", "nFeature_RNA", "negmean", "falsecodemean")]])

## Sanity check: an FOV should lie in only one section.
fov_multi <- md[, .(n_sections = uniqueN(section_id)), by = fov][n_sections > 1]
if (nrow(fov_multi) > 0) { warning("FOVs spanning >1 section:"); print(fov_multi) }

## ---- 9a: per-FOV table (one row per FOV; .N = cells in that FOV) ----
fov_summary <- md[, .(
  n_cells              = .N,
  total_transcripts    = sum(nCount_RNA),
  median_tx_per_cell   = median(nCount_RNA),
  median_cell_area_um2 = median(Area.um2),
  fov_x_min_mm = min(x_slide_mm), fov_x_max_mm = max(x_slide_mm),
  fov_y_min_mm = min(y_slide_mm), fov_y_max_mm = max(y_slide_mm),
  fov_x_centre_mm = mean(range(x_slide_mm)), fov_y_centre_mm = mean(range(y_slide_mm))
), by = .(section_id, fov)][order(section_id, fov)]
fov_summary

## Cells per FOV, one column per section (0 = FOV not in that section).
cells_per_fov_per_section <- dcast(fov_summary, fov ~ section_id, value.var = "n_cells", fill = 0)

## ---- 9b: tissue-section area ----
## Occupied-bin area: lay a 50 µm grid; count squares holding >= 1 cell centre; x square area.
##   Follows the true outline, including a ring-shaped aortic wall.
## Convex-hull area: smallest convex polygon around all cells (shoelace formula).
##   Overestimates rings because it includes the lumen.
bin_mm <- 0.05
occupied_area_mm2 <- function(x, y, b) nrow(unique(data.table(floor(x / b), floor(y / b)))) * b^2
hull_area_mm2 <- function(x, y) {
  h <- chull(x, y); xs <- x[h]; ys <- y[h]                       # hull vertices in order
  0.5 * abs(sum(xs * c(ys[-1], ys[1]) - c(xs[-1], xs[1]) * ys))  # shoelace formula
}

## ---- 9c: per-section table ----
section_summary <- md[, .(
  n_fov                 = uniqueN(fov),
  n_cells               = .N,
  total_transcripts     = sum(nCount_RNA),
  median_tx_per_cell    = median(nCount_RNA),
  median_genes_per_cell = median(nFeature_RNA),
  median_cell_area_um2  = median(Area.um2),
  summed_cell_area_mm2  = sum(Area.um2) / 1e6,                     # µm² -> mm²
  tissue_area_bins_mm2  = occupied_area_mm2(x_slide_mm, y_slide_mm, bin_mm),
  tissue_area_hull_mm2  = hull_area_mm2(x_slide_mm, y_slide_mm),
  mean_negmean          = mean(negmean),
  mean_falsecodemean    = mean(falsecodemean)
), by = section_id][order(section_id)]
section_summary


## =============================================================================
## CHUNK 10: per-section objects (subsets only; nothing recomputed)
## =============================================================================
## DEFENCE: each section's tables are row subsets of the whole-slide objects;
## its Seurat object is subset(obj), which keeps assays and spatial data linked.
## Transcripts carry no cell label, so each is given the section of its nearest
## cell centre (an approximation).
if (!is.null(tx_df) && nrow(tx_df) > 0) {
  nn <- kNN(xy, k = 1, query = as.matrix(tx_df[, c("x", "y")]))$id[, 1]  # row of nearest cell
  tx_df$section_id <- unname(obj$section_id)[nn]
}

section_names <- sort(unique(obj$section_id))
metadata_by_section <- expression_by_section <- negprobes_by_section <- list()
falsecode_by_section <- xy_by_section <- polygons_by_section <- list()
transcripts_by_section <- seurat_by_section <- list()

for (sec in section_names) {
  cells_sec <- Cells(obj)[obj$section_id == sec]   # cell names of this section
  
  ## Row subsets of the whole-slide objects.
  metadata_by_section[[sec]]   <- obj[[]][cells_sec, ]          # metadata rows, straight from obj
  expression_by_section[[sec]] <- counts[cells_sec, ]
  negprobes_by_section[[sec]]  <- negprobes_counts[cells_sec, ]
  falsecode_by_section[[sec]]  <- falsecode_counts[cells_sec, ]
  xy_by_section[[sec]]         <- xy[cells_sec, ]
  
  ## assign() also creates named copies (AAA2_metadata ...) so View(AAA2_metadata) works.
  ## OPTIONAL convenience: the lists above already hold the same data.
  assign(paste0(sec, "_metadata"),          metadata_by_section[[sec]],     envir = .GlobalEnv)
  assign(paste0(sec, "_expression_matrix"), expression_by_section[[sec]],   envir = .GlobalEnv)
  assign(paste0(sec, "_negprobes"),         negprobes_by_section[[sec]],    envir = .GlobalEnv)
  assign(paste0(sec, "_falsecode"),         falsecode_by_section[[sec]],    envir = .GlobalEnv)
  assign(paste0(sec, "_xy"),                xy_by_section[[sec]],           envir = .GlobalEnv)
  assign(paste0(sec, "_fov_summary"),       fov_summary[section_id == sec], envir = .GlobalEnv)
  
  if (!is.null(polygons_df)) {
    polygons_by_section[[sec]] <- polygons_df[polygons_df$cell %in% cells_sec, ]
    assign(paste0(sec, "_polygons"), polygons_by_section[[sec]], envir = .GlobalEnv)
  }
  if (!is.null(tx_df) && nrow(tx_df) > 0) {
    transcripts_by_section[[sec]] <- tx_df[tx_df$section_id == sec, ]
    assign(paste0(sec, "_transcripts"), transcripts_by_section[[sec]], envir = .GlobalEnv)
  }
  
  ## Seurat object of this section only (assays + spatial data kept together).
  seurat_by_section[[sec]] <- tryCatch(subset(obj, cells = cells_sec),
                                       error = function(e) { message(sec, ": subset failed - ", e$message); NULL })
  assign(paste0(sec, "_seurat"), seurat_by_section[[sec]], envir = .GlobalEnv)
  ## SEURAT ALTERNATIVE, same result: subset(obj, subset = section_id == sec)
  
  ## Check the subsets still line up.
  stopifnot(identical(rownames(expression_by_section[[sec]]), rownames(metadata_by_section[[sec]])))
  message(sec, ": ", length(cells_sec), " cells, ", uniqueN(metadata_by_section[[sec]]$fov), " FOVs")
}
ls(pattern = "^AAA")   # every per-section object just created


## =============================================================================
## CHUNK 11: plots for visual verification
## =============================================================================
## DEFENCE: confirms the section labels and FOV layout visually against the slide.
## 11a: whole slide, sections coloured and labelled at their centre.
centroids_sec <- md[, .(x = mean(x_slide_mm), y = mean(y_slide_mm)), by = section_id]
p_overview <- ggplot(md, aes(x_slide_mm, y_slide_mm, colour = section_id)) +
  geom_point(size = 0.05, alpha = 0.3) +
  geom_text(data = centroids_sec, aes(x, y, label = section_id), colour = "black", size = 5, fontface = "bold") +
  coord_fixed() + theme_minimal() +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "AAA CosMx - 6 tissue sections", x = "x (mm)", y = "y (mm)")
ggsave(file.path(dir_plots, "section_overview_plot.png"), p_overview, width = 8, height = 8, dpi = 300)
print(p_overview)
## SEURAT ALTERNATIVE: ImageDimPlot(obj, fov = img_name, group.by = "section_id")
## (SpatialDimPlot is for Visium, not CosMx.) Axis orientation may differ from this plot.

## 11b: each section, cells coloured by FOV, FOV numbers at FOV centres.
for (sec in section_names) {
  d <- md[section_id == sec]; fl <- fov_summary[section_id == sec]
  p <- ggplot(d, aes(x_slide_mm, y_slide_mm, colour = factor(fov))) +
    geom_point(size = 0.1, alpha = 0.5, show.legend = FALSE) +
    geom_text(data = fl, aes(fov_x_centre_mm, fov_y_centre_mm, label = fov), colour = "black", size = 2.5) +
    coord_fixed() + theme_minimal() +
    labs(title = paste0(sec, " - ", nrow(d), " cells, ", nrow(fl), " FOVs"), x = "x (mm)", y = "y (mm)")
  ggsave(file.path(dir_plots, paste0(sec, "_spatial_by_fov.png")), p, width = 7, height = 7, dpi = 300)
}

## 11c: cell polygons of one FOV, filled by transcripts per cell.
## group = cell draws each cell's vertices as one closed shape.
if (!is.null(polygons_df)) {
  example_fov <- fov_summary[section_id == section_names[1]]$fov[1]
  poly_fov <- polygons_df[polygons_df$cell %in% md[fov == example_fov]$cell_id, ]
  poly_fov$nCount_RNA <- obj[[]][poly_fov$cell, "nCount_RNA"]   # each vertex gets its cell's count
  p_poly <- ggplot(poly_fov, aes(x, y, group = cell, fill = nCount_RNA)) +
    geom_polygon(colour = "grey30", linewidth = 0.1) +
    coord_fixed() + theme_minimal() + scale_fill_viridis_c() +
    labs(title = paste0("Cell polygons - FOV ", example_fov, " (", section_names[1], ")"))
  ggsave(file.path(dir_plots, paste0("polygons_FOV", example_fov, ".png")), p_poly, width = 7, height = 7, dpi = 300)
}

## 11d: transcript locations of the chosen genes.
if (!is.null(tx_df) && nrow(tx_df) > 0) {
  p_tx <- ggplot(tx_df, aes(x, y, colour = molecule)) +
    geom_point(size = 0.01, alpha = 0.3) + coord_fixed() + theme_minimal() +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(title = "Transcript locations", x = "x (mm)", y = "y (mm)")
  ggsave(file.path(dir_plots, "transcripts_selected_genes.png"), p_tx, width = 8, height = 8, dpi = 300)
}


## =============================================================================
## CHUNK 12: save
## =============================================================================
## DEFENCE: everything is saved so later steps reload with readRDS() instead of
## re-running this script.
saveRDS(counts,           file.path(output_dir, "counts_unfiltered.RDS"))
saveRDS(negprobes_counts, file.path(output_dir, "negprobes_unfiltered.RDS"))
saveRDS(falsecode_counts, file.path(output_dir, "falsecode_unfiltered.RDS"))
saveRDS(obj[[]],          file.path(output_dir, "metadata_unfiltered.RDS"))   # metadata straight from obj
saveRDS(xy,               file.path(output_dir, "xy_unfiltered.RDS"))
saveRDS(obj,              file.path(output_dir, "seurat_with_sections.RDS"))
if (!is.null(polygons_df)) saveRDS(polygons_df, file.path(output_dir, "polygons_unfiltered.RDS"))
fwrite(fov_summary,               file.path(output_dir, "fov_summary.csv"))
fwrite(section_summary,           file.path(output_dir, "section_summary.csv"))
fwrite(cells_per_fov_per_section, file.path(output_dir, "cells_per_fov_per_section.csv"))

## Per-section files, one folder per section.
for (sec in section_names) {
  d <- file.path(dir_sections, sec); dir.create(d, showWarnings = FALSE)
  saveRDS(metadata_by_section[[sec]],   file.path(d, paste0(sec, "_metadata.RDS")))
  saveRDS(expression_by_section[[sec]], file.path(d, paste0(sec, "_expression_matrix.RDS")))
  saveRDS(negprobes_by_section[[sec]],  file.path(d, paste0(sec, "_negprobes.RDS")))
  saveRDS(falsecode_by_section[[sec]],  file.path(d, paste0(sec, "_falsecode.RDS")))
  saveRDS(xy_by_section[[sec]],         file.path(d, paste0(sec, "_xy.RDS")))
  fwrite(fov_summary[section_id == sec], file.path(d, paste0(sec, "_fov_summary.csv")))
  if (!is.null(polygons_df))              saveRDS(polygons_by_section[[sec]],    file.path(d, paste0(sec, "_polygons.RDS")))
  if (length(transcripts_by_section))     saveRDS(transcripts_by_section[[sec]], file.path(d, paste0(sec, "_transcripts.RDS")))
  if (!is.null(seurat_by_section[[sec]])) saveRDS(seurat_by_section[[sec]],      file.path(d, paste0(sec, "_seurat.RDS")))
}
message("All outputs saved to ", output_dir)


## =============================================================================
## CHUNK 13: how to VIEW things
## =============================================================================
## Whole slide
# View(obj[[]])                           # all cells, all metadata
# View(section_summary)                   # one row per section: FOVs, cells, transcripts, areas, background
# View(fov_summary)                       # one row per FOV: cells, transcripts, position
# View(cells_per_fov_per_section)
# View(as.matrix(counts[1:50, 1:50]))     # corner of the whole-slide expression matrix
# head(xy); head(polygons_df); head(tx_df)
#
## One section (swap AAA2 for AAA38, AAA40, AAA45, AAA46, AAA48)
# View(AAA2_metadata)
# View(AAA2_fov_summary)
# dim(AAA2_expression_matrix)
# View(as.matrix(AAA2_expression_matrix[1:50, 1:50]))
# View(as.matrix(AAA2_negprobes[1:50, ]))
# View(as.matrix(AAA2_falsecode[1:50, ]))
# head(AAA2_xy); head(AAA2_polygons); head(AAA2_transcripts)
# AAA2_seurat
#
## Useful one-liners
# table(AAA2_metadata$fov)                                       # cells per FOV
# summary(AAA2_metadata$nCount_RNA)                              # transcripts per cell
# summary(AAA2_metadata$Area.um2)                                # cell areas (µm²)
# sort(colSums(AAA2_expression_matrix), decreasing = TRUE)[1:20] # top 20 genes
## =============================================================================