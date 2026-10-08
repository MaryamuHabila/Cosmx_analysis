

## =============================================================================
## CHUNK 1: Paths and output folders
## =============================================================================
## PURPOSE
##   Start from the AtoMx Seurat object, view everything it contains, then split
##   it into the 6 tissue sections (AAA2, AAA38, AAA40, AAA45, AAA46, AAA48).
##
## GUIDING PRINCIPLE: reuse, do not recalculate
##   AtoMx already recorded most per-cell information in the metadata:
##     cell ID ................ cell_id      (already unique: "c_<slide>_<FOV>_<cell>")
##     FOV of each cell ....... fov
##     x, y position (mm) ..... x_slide_mm, y_slide_mm   (already in mm; no px conversion needed)
##     cell area (µm²) ........ Area.um2
##     transcripts per cell ... nCount_RNA
##     genes per cell ......... nFeature_RNA
##     negative-probe sum ..... nCount_negprobes
##     false-code sum ......... nCount_falsecode
##   So every per-section table below is a SUBSET of these existing rows/columns.
##
##   Only three things are not in the AtoMx output and must be derived:
##     1. section_id    : AtoMx does not know which tissue piece a cell is in.
##     2. negmean / falsecodemean : sum / number of control probes.
##     3. tissue-section area : a property of a whole section, not of a cell.
##   Counts of cells/FOVs/transcripts per FOV or per section are simple
##   tallies of existing columns, not new measurements.
##
## HOW TO RUN: select a chunk in RStudio, press Cmd+Enter, read the output,
## then move on.
## =============================================================================

## A package is a collection of functions; library() loads it into this session.
## install.packages() is needed only once per computer.
# install.packages(c("Seurat", "dbscan", "data.table", "ggplot2"))

suppressPackageStartupMessages({
  library(Seurat)        # plotting and subsetting of Seurat objects
  library(SeuratObject)  # accessors: LayerData(), GetTissueCoordinates()
  library(Matrix)        # sparse matrices: store only non-zero values
  library(data.table)    # fast grouped summaries (per FOV, per section)
  library(dbscan)        # spatial clustering, used once to find the 6 sections
  library(ggplot2)       # plots
})

## Paths tell R where to read the input and where to write results.
## The Seurat object I want to analyse is located at this particular place on my computer."
##Tell R where input is
##↓
##Tell R where output should go
##↓
##Create folders for output


seurat_file <- "/Users/maryamuhabilausman/Desktop/AAA_CosMx_Analysis/Downloaded_data/seuratObject_AC20260629_AAA3.RDS"

##  This is where my CosMx slide-related files are
slide_dir   <- "/Users/maryamuhabilausman/Desktop/AAA_CosMx_Analysis/Downloaded_data/flatFiles/AC20260629_AAA3"  # only for the optional FOV-position file

##This is the folder where I want my analysis results to go
output_dir  <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"


##Take my output folder and make a path to a subfolder called spatial_plots
dir_plots    <- file.path(output_dir, "spatial_plots")   # figures
dir_sections <- file.path(output_dir, "per_section")     # one sub-folder per section
dir.create(dir_plots,    showWarnings = FALSE, recursive = TRUE)  # recursive = also create parent folders
dir.create(dir_sections, showWarnings = FALSE, recursive = TRUE)


##Alternatively Use seurat
##seurat_file <- "/Users/maryamuhabilausman/Desktop/AAA_CosMx_Analysis/Downloaded_data/seuratObject_AC20260629_AAA3.RDS"

##output_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"

##dir_plots <- file.path(output_dir, "spatial_plots")
##dir_sections <- file.path(output_dir, "per_section")

##dir.create(dir_plots, showWarnings = FALSE, recursive = TRUE)
##dir.create(dir_sections, showWarnings = FALSE, recursive = TRUE)












## =============================================================================
## CHUNK 2: Opening the Seurat object
## =============================================================================
## A Seurat object is one container holding assays (count matrices),
## meta.data (one row per cell) and images (spatial coordinates).
##"Take the .RDS file stored at seurat_file, open it, and store whatever is inside it as obj

obj <- readRDS(seurat_file)   # readRDS() restores an R object saved on disk

obj                  # Gene expression data, metadata, spatial information, images/FOV information, assays (RNA, negprobes and falsecode) 
# Show me what this object looks like.
#A summary of the Seurat object
dim(obj)     # What are the dimensions of this object? genes x cells of the active assay (RNA)
dim(obj[["RNA"]]$counts)  #interested in the raw counts matrix and your object contains it there
dim(obj[["RNA"]])  # Give me the dimensions of the RNA assay.
Assays(obj)          #What assays are stored inside my Seurat object, RNA = real genes; negprobes and falsecode = background controls
DefaultAssay(obj)    # Which assay has Seurat currently designated as the default?, assay used when none is named
Images(obj)          # "What spatial image/FOV objects are stored in my Seurat object?", name of the spatial object (centroids, polygons, transcripts)
head(colnames(obj))  # What are the names of the columns, cell names
head(rownames(obj))  # Only show me the first few 
## 6183 genes x 123,709 cells; assays RNA, negprobes, falsecode.

##Alternatively Use seurat
##obj <- readRDS("/Users/maryamuhabilausman/Desktop/AAA_CosMx_Analysis/Downloaded_data/seuratObject_AC20260629_AAA3.RDS")





## =============================================================================
## CHUNK 3: 
## =============================================================================
## Get expression/count matrices, 
## Get cell metadata
## Get spatial information
## Make an x/y coordinate table
## Check that everything matches
## Calculate background measures
## Identify the 6 tissue sections
## Calculate statistics for FOVs and sections
## Split everything into individual sections
## Make plots to visually check things
## Save everything
## Show you how to inspect the results

##pulling things out of Seurat and turning them into separate R matrices/data frames/lists.
## A mixture of Seurat, ordinary R, data.table, ggplot2, and DBSCAN.

## o = which Seurat object?
## assay = which assay?
  
get_counts <- function(o, assay) {
  tryCatch(LayerData(o, assay = assay, layer = "counts"),
           error = function(e) GetAssayData(o, assay = assay, slot = "counts"))
}

## "Take obj and get the counts from the RNA assay."
## normally is genes x cells t() transposes to cells x features, the orientation used by the NanoString workflow.
##Get the raw RNA counts from Seurat and transpose them so cells are rows and genes are columns.

counts           <- Matrix::t(get_counts(obj, "RNA"))        # cells x genes
negprobes_counts <- Matrix::t(get_counts(obj, "negprobes"))  # cells x negative probes
falsecode_counts <- Matrix::t(get_counts(obj, "falsecode"))  # cells x false codes


## Alternative maybe use LayerData(obj, assay = "RNA", layer = "counts")

## i should now see cells x genes
dim(counts);           counts[1:5, 1:5]           # "." = zero in a sparse matrix
dim(negprobes_counts); negprobes_counts[1:5, 1:5]
dim(falsecode_counts); falsecode_counts[1:5, 1:5]
class(counts)   ##a sparse matrix.                                  # dgCMatrix = compressed sparse matrix

## Sparsity: fraction of zeros; explains why a sparse format is used.
1 - Matrix::nnzero(counts) / prod(dim(counts))

##
## Negative probes target sequences absent from the human genome.
## False codes are barcodes not assigned to any probe.
## Both measure background: counts that are noise, not real expression.





## Expression = what the cell contains.

## Metadata = information describing the cell.

## "Take the Seurat object's metadata and store it as a separate object called metadata."
metadata <- obj[[]]

## How many cells and metadata variables do I have?
dim(metadata)          # cells x variables
## Show me every metadata variable.
colnames(metadata)     # every variable AtoMx recorded per cell
head(metadata, 5)      # first 5 cells

##Tell me what type each column contains
str(metadata)          # data type of each column
# View(metadata)       # spreadsheet view



## =============================================================================
## CHUNK 5: CENTROIDS, CELL POLYGONS, TRANSCRIPT LOCATIONS, FOV POSITIONS
## =============================================================================
## Spatial data live in an FOV-class object inside obj.
## It tells  which spatial image/FOV objects are stored in the Seurat object.
img_name <- Images(obj)[1]

##"Take that spatial object out of the Seurat object and store it as fov_obj."
##"What components are inside this spatial object?"
fov_obj  <- obj[[img_name]]
names(fov_obj)   # EXPECT: "centroids" "segmentation" "molecules"

## ---- 5a: centroids = one centre point per cell ----
## Give me the coordinates of the cell centres
centroids_df <- GetTissueCoordinates(fov_obj, which = "centroids")   # columns x, y, cell
head(centroids_df)
summary(centroids_df[, c("x", "y")])   # values < ~50 → units are mm

## "Are the coordinates in the spatial object exactly the same as the x/y coordinates already stored in my metadata?"
##Confirm centroids are the same numbers as x_slide_mm / y_slide_mm.
## If TRUE, the metadata columns can be used as the xy coordinates directly.
## Seurat spatial coordinates =  AtoMx metadata coordinates
isTRUE(all.equal(centroids_df$x, metadata[centroids_df$cell, "x_slide_mm"])) &&
  isTRUE(all.equal(centroids_df$y, metadata[centroids_df$cell, "y_slide_mm"]))


## ---- 5b: cell polygons = segmentation outline of each cell ----
## "Does fov_obj contain something called segmentation?"

## Each cell has several rows: one per outline vertex, in drawing order.
## This extracts the coordinates defining each cell's outline.
## How many polygon vertices are there in total?
## "How many unique cells have segmentation outlines?"
## "On average, how many vertices does each cell polygon have?"
polygons_df <- NULL
if ("segmentation" %in% names(fov_obj)) {
  polygons_df <- GetTissueCoordinates(fov_obj, which = "segmentation")  # columns x, y, cell
  print(head(polygons_df, 10))
  print(nrow(polygons_df))                  # total vertices
  print(length(unique(polygons_df$cell)))   # cells with an outline; expect 123,709
  print(mean(table(polygons_df$cell)))      # mean vertices per cell
}

## ---- 5c: transcript locations = x, y of each detected RNA molecule ----
## The full table has hundreds of millions of rows, so load only a few genes.
genes_to_show <- intersect(c("AC024901.1"), rownames(obj))  # keep genes present in the panel
tx_df <- NULL
if ("molecules" %in% names(fov_obj)) {
  tx_df <- GetTissueCoordinates(fov_obj[["molecules"]], features = genes_to_show)  # columns x, y, molecule
  print(head(tx_df))
  print(table(tx_df$molecule))              # transcripts per chosen gene
}

## ---- 5d: FOV positions ----
## Two sources:
## (i) the AtoMx fov_positions flat file
## (ii) the extent of each FOV's cells, computed in Chunk 10 (fov_summary).
fovfile <- list.files(slide_dir, pattern = "fov_positions.*", full.names = TRUE)[1]
fov_positions_flat <- NULL
if (!is.na(fovfile)) {
  fov_positions_flat <- fread(fovfile)
  print(head(fov_positions_flat))           # one row per FOV and its slide position
}


## =============================================================================
## x_slide_mm / y_slide_mm are AtoMx's global coordinates, already in mm.
## (From pixels it would be: px * 120.280945 nm/px / 1e6 = mm; not needed here.)
xy <- as.matrix(metadata[, c("x_slide_mm", "y_slide_mm")])
colnames(xy) <- c("x_mm", "y_mm")
## "Use each cell's ID as the row name."
rownames(xy) <- metadata$cell_id   # same cell names as counts and metadata
head(xy)
apply(xy, 2, range)                # slide extent in mm


## =============================================================================
## CHUNK 7: integrity checks
## =============================================================================
## stopifnot() stops the script if a statement is FALSE, so errors cannot pass unnoticed.

## 1) cell_id is already unique across the slide (it encodes slide, FOV and cell).
stopifnot(!anyDuplicated(metadata$cell_id))

## 2) counts, control matrices, metadata and xy list the same cells in the same order.
stopifnot(identical(rownames(metadata),         metadata$cell_id))
stopifnot(identical(rownames(counts),           metadata$cell_id))
stopifnot(identical(rownames(negprobes_counts), metadata$cell_id))
stopifnot(identical(rownames(falsecode_counts), metadata$cell_id))
stopifnot(identical(rownames(xy),               metadata$cell_id))

## 3) per cell, summing counts across features reproduces AtoMx's recorded totals.
stopifnot(max(abs(Matrix::rowSums(counts)           - metadata$nCount_RNA))       == 0)
stopifnot(max(abs(Matrix::rowSums(negprobes_counts) - metadata$nCount_negprobes)) == 0)
stopifnot(max(abs(Matrix::rowSums(falsecode_counts) - metadata$nCount_falsecode)) == 0)
message("Integrity checks passed: identities and per-cell totals match.")


## =============================================================================
## CHUNK 8: mean background per cell (negprobes and falsecode)
## =============================================================================
## The SUM of control counts depends on how many control probes the panel has.
## Dividing by that number gives a per-probe MEAN, comparable across panels.
n_negprobes <- ncol(negprobes_counts)   # number of negative probes in this panel
n_falsecode <- ncol(falsecode_counts)   # number of false codes in this panel
n_negprobes; n_falsecode

metadata$negmean       <- metadata$nCount_negprobes / n_negprobes
metadata$falsecodemean <- metadata$nCount_falsecode / n_falsecode
summary(metadata$negmean)
summary(metadata$falsecodemean)


## =============================================================================
## CHUNK 9: section label (the one thing AtoMx does not provide)
## =============================================================================
## If section_id already exists (e.g. you load seurat_with_sections.RDS), reuse it.
## Otherwise derive it once with dbscan on xy (position only, no gene expression):
##   eps = 1 mm    : cells closer than 1 mm are chained into the same piece;
##                   separate tissue pieces lie further apart than that.
##   minPts = 1    : every cell is assigned a cluster; none labelled noise.
if ("section_id" %in% colnames(metadata)) {
  message("section_id already present - reused, no clustering run.")
} else {
  clust_raw <- dbscan::dbscan(xy, eps = 1, minPts = 1)$cluster
  metadata$section_id_placeholder <- paste0("section", clust_raw)
  print(table(metadata$section_id_placeholder))           # cells per cluster
  if (length(unique(clust_raw)) != 6) warning("Expected 6 clusters - adjust eps.")
  
  ## Cluster centres: where each cluster sits on the slide.
  print(aggregate(cbind(x_slide_mm, y_slide_mm) ~ section_id_placeholder, data = metadata, FUN = mean))
  
  ## Plot to match clusters against the AtoMx slide thumbnail.
  p_placeholder <- ggplot(metadata, aes(x_slide_mm, y_slide_mm, colour = section_id_placeholder)) +
    geom_point(size = 0.05, alpha = 0.3) + coord_fixed() + theme_minimal() +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(title = "dbscan clusters - match to AtoMx thumbnail", x = "x (mm)", y = "y (mm)")
  ggsave(file.path(dir_plots, "placeholder_overview_plot.png"), p_placeholder, width = 8, height = 8, dpi = 300)
  print(p_placeholder)
  
  ## dbscan numbers clusters in the order it meets cells; it cannot know tissue IDs.
  ## >>> CONFIRM this mapping against the plot and AtoMx thumbnail before trusting it. <<<
  section_name_mapping <- c(section1 = "AAA38", section2 = "AAA45", section3 = "AAA46",
                            section4 = "AAA40", section5 = "AAA48", section6 = "AAA2")
  stopifnot(all(unique(metadata$section_id_placeholder) %in% names(section_name_mapping)))
  metadata$section_id <- unname(section_name_mapping[metadata$section_id_placeholder])
}
table(metadata$section_id)   # cells per named section

## Copy the new columns into the Seurat object so Seurat functions can use them.
obj$section_id    <- metadata$section_id
obj$negmean       <- metadata$negmean
obj$falsecodemean <- metadata$falsecodemean


## =============================================================================
## CHUNK 10: tallies per FOV and per section (all from existing columns)
## =============================================================================
md <- as.data.table(metadata)   # data.table allows fast "by group" summaries

## Sanity check: an FOV should lie in only one section.
fov_multi <- md[, .(n_sections = uniqueN(section_id)), by = fov][n_sections > 1]
if (nrow(fov_multi) > 0) { warning("FOVs spanning >1 section:"); print(fov_multi) }

## ---- 10a: per-FOV table ----
## by = .(section_id, fov) → one row per FOV; .N = number of cells in that FOV.
## FOV position = extent and centre of its cells' coordinates.
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

## ---- 10b: tissue-section area (a section-level property, so it must be computed) ----
## Occupied-bin area: lay a 50 µm grid; count squares containing ≥1 cell; multiply by square area.
##   Follows the true outline, including a ring-shaped aortic wall.
## Convex-hull area: area of the smallest convex polygon enclosing all cells.
##   Overestimates rings/curves because it includes the lumen.
bin_mm <- 0.05
occupied_area_mm2 <- function(x, y, b) nrow(unique(data.table(floor(x / b), floor(y / b)))) * b^2
hull_area_mm2 <- function(x, y) {
  h <- chull(x, y); xs <- x[h]; ys <- y[h]                       # hull vertices
  0.5 * abs(sum(xs * c(ys[-1], ys[1]) - c(xs[-1], xs[1]) * ys))  # shoelace formula
}

## ---- 10c: per-section table ----
section_summary <- md[, .(
  n_fov                 = uniqueN(fov),
  n_cells               = .N,
  total_transcripts     = sum(nCount_RNA),
  median_tx_per_cell    = median(nCount_RNA),
  median_genes_per_cell = median(nFeature_RNA),
  median_cell_area_um2  = median(Area.um2),
  summed_cell_area_mm2  = sum(Area.um2) / 1e6,                     # µm² → mm²
  tissue_area_bins_mm2  = occupied_area_mm2(x_slide_mm, y_slide_mm, bin_mm),
  tissue_area_hull_mm2  = hull_area_mm2(x_slide_mm, y_slide_mm),
  mean_negmean          = mean(negmean),
  mean_falsecodemean    = mean(falsecodemean)
), by = section_id][order(section_id)]
section_summary


## =============================================================================
## CHUNK 11: per-section objects with clear names (subsets, nothing recomputed)
## =============================================================================
## Transcripts carry no cell label in the Seurat object, so each chosen-gene
## transcript is given the section of its nearest cell centre.
if (!is.null(tx_df) && nrow(tx_df) > 0) {
  nn <- dbscan::kNN(xy, k = 1, query = as.matrix(tx_df[, c("x", "y")]))$id[, 1]
  tx_df$section_id <- metadata$section_id[nn]
}

section_names <- sort(unique(metadata$section_id))
metadata_by_section <- expression_by_section <- negprobes_by_section <- list()
falsecode_by_section <- xy_by_section <- polygons_by_section <- list()
transcripts_by_section <- seurat_by_section <- list()

for (sec in section_names) {
  cells_sec <- metadata$cell_id[metadata$section_id == sec]   # cells of this section
  
  ## Row subsets of the combined objects.
  metadata_by_section[[sec]]   <- metadata[cells_sec, ]
  expression_by_section[[sec]] <- counts[cells_sec, ]
  negprobes_by_section[[sec]]  <- negprobes_counts[cells_sec, ]
  falsecode_by_section[[sec]]  <- falsecode_counts[cells_sec, ]
  xy_by_section[[sec]]         <- xy[cells_sec, ]
  
  ## assign() creates a named object, e.g. AAA2_metadata, so View(AAA2_metadata) works.
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
  
  ## Check the subsets still line up.
  stopifnot(identical(rownames(expression_by_section[[sec]]), rownames(metadata_by_section[[sec]])))
  message(sec, ": ", length(cells_sec), " cells, ", uniqueN(metadata_by_section[[sec]]$fov), " FOVs")
}
ls(pattern = "^AAA")   # every per-section object just created


## =============================================================================
## CHUNK 12: plots for visual verification
## =============================================================================
## 12a: whole slide, sections coloured and labelled at their centre.
centroids_sec <- md[, .(x = mean(x_slide_mm), y = mean(y_slide_mm)), by = section_id]
p_overview <- ggplot(md, aes(x_slide_mm, y_slide_mm, colour = section_id)) +
  geom_point(size = 0.05, alpha = 0.3) +
  geom_text(data = centroids_sec, aes(x, y, label = section_id), colour = "black", size = 5, fontface = "bold") +
  coord_fixed() + theme_minimal() +
  guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
  labs(title = "AAA CosMx - 6 tissue sections", x = "x (mm)", y = "y (mm)")
ggsave(file.path(dir_plots, "section_overview_plot.png"), p_overview, width = 8, height = 8, dpi = 300)
print(p_overview)

## 12b: each section, cells coloured by FOV, FOV numbers at FOV centres.
for (sec in section_names) {
  d <- md[section_id == sec]; fl <- fov_summary[section_id == sec]
  p <- ggplot(d, aes(x_slide_mm, y_slide_mm, colour = factor(fov))) +
    geom_point(size = 0.1, alpha = 0.5, show.legend = FALSE) +
    geom_text(data = fl, aes(fov_x_centre_mm, fov_y_centre_mm, label = fov), colour = "black", size = 2.5) +
    coord_fixed() + theme_minimal() +
    labs(title = paste0(sec, " - ", nrow(d), " cells, ", nrow(fl), " FOVs"), x = "x (mm)", y = "y (mm)")
  ggsave(file.path(dir_plots, paste0(sec, "_spatial_by_fov.png")), p, width = 7, height = 7, dpi = 300)
}

## 12c: cell polygons of one FOV, filled by transcripts per cell.
## group = cell tells ggplot to draw each cell's vertices as one closed shape.
if (!is.null(polygons_df)) {
  example_fov <- fov_summary[section_id == section_names[1]]$fov[1]
  poly_fov <- polygons_df[polygons_df$cell %in% md[fov == example_fov]$cell_id, ]
  poly_fov$nCount_RNA <- metadata[poly_fov$cell, "nCount_RNA"]
  p_poly <- ggplot(poly_fov, aes(x, y, group = cell, fill = nCount_RNA)) +
    geom_polygon(colour = "grey30", linewidth = 0.1) +
    coord_fixed() + theme_minimal() + scale_fill_viridis_c() +
    labs(title = paste0("Cell polygons - FOV ", example_fov, " (", section_names[1], ")"))
  ggsave(file.path(dir_plots, paste0("polygons_FOV", example_fov, ".png")), p_poly, width = 7, height = 7, dpi = 300)
}

## 12d: transcript locations of the chosen genes.
if (!is.null(tx_df) && nrow(tx_df) > 0) {
  p_tx <- ggplot(tx_df, aes(x, y, colour = molecule)) +
    geom_point(size = 0.01, alpha = 0.3) + coord_fixed() + theme_minimal() +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(title = "Transcript locations", x = "x (mm)", y = "y (mm)")
  ggsave(file.path(dir_plots, "transcripts_selected_genes.png"), p_tx, width = 8, height = 8, dpi = 300)
}


## =============================================================================
## CHUNK 13: save
## =============================================================================
## Combined objects: reload later with readRDS() instead of repeating this script.
saveRDS(counts,           file.path(output_dir, "counts_unfiltered.RDS"))
saveRDS(negprobes_counts, file.path(output_dir, "negprobes_unfiltered.RDS"))
saveRDS(falsecode_counts, file.path(output_dir, "falsecode_unfiltered.RDS"))
saveRDS(metadata,         file.path(output_dir, "metadata_unfiltered.RDS"))
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
  if (!is.null(polygons_df))                    saveRDS(polygons_by_section[[sec]],    file.path(d, paste0(sec, "_polygons.RDS")))
  if (length(transcripts_by_section))           saveRDS(transcripts_by_section[[sec]], file.path(d, paste0(sec, "_transcripts.RDS")))
  if (!is.null(seurat_by_section[[sec]]))       saveRDS(seurat_by_section[[sec]],      file.path(d, paste0(sec, "_seurat.RDS")))
}
message("All outputs saved to ", output_dir)


## =============================================================================
## CHUNK 14: how to VIEW things
## =============================================================================
## Whole slide
# View(metadata)                          # all cells, all metadata
# View(section_summary)                   # one row per section: FOVs, cells, transcripts, areas, background
# View(fov_summary)                       # one row per FOV: cells, transcripts, position
# View(cells_per_fov_per_section)
# View(as.matrix(counts[1:50, 1:50]))     # corner of the combined expression matrix
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
# table(AAA2_metadata$fov)                       # cells per FOV
# summary(AAA2_metadata$nCount_RNA)              # transcripts per cell
# summary(AAA2_metadata$Area.um2)                # cell areas (µm²)
# sort(Matrix::colSums(AAA2_expression_matrix), decreasing = TRUE)[1:20]   # top 20 genes
## =============================================================================























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
## 05_QC_FOV.R
## =============================================================================
## PURPOSE
##   Vignette step 4.2 (FOV QC), default thresholds (0.6 / 0.6):
##     runFOVQC(counts, xy, fov, barcodemap)       # vignette 1 call
##   Two checks, both inside runFOVQC():
##     1. low total counts  -> fovqcresult$flaggedfovs_fortotalcounts
##     2. reporter bias     -> fovqcresult$flaggedfovs_forbias
##   Plots:
##     vignette 1: mapFlaggedFOVs(), FOVSignalLossSpatialPlot()
##     vignette 2: bias heatmap (ComplexHeatmap), genes-flagged bar plot
##   Summary: FOVs flagged of total FOVs, whole slide and per tissue section.
##   FOV flags are NOT combined with the low-count/area flags here.
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

## Vignette 1 runs FOV QC on all cells (before any cells are removed).
stopifnot(identical(rownames(counts), metadata$cell_id))
stopifnot(identical(rownames(xy),     metadata$cell_id))
section_names <- sort(unique(metadata$section_id))

## Thresholds: function defaults, written out so they are visible.
max_totalcounts_loss <- 0.6   # low counts: grid square flagged if < 40% of expected total counts
max_prop_loss        <- 0.6   # reporter bias: bit flagged if < 40% of expected signal

## Output folder
dir_fov <- file.path(output_dir, "QC", "FOV")
dir.create(dir_fov, showWarnings = FALSE, recursive = TRUE)

graphics.off()   # clear any open plot device


## =============================================================================
## CHUNK 25: FOV QC tool and barcode map (vignette)
## =============================================================================
## Source the FOV QC functions (runFOVQC, mapFlaggedFOVs, FOVSignalLossSpatialPlot, ...)
source("https://raw.githubusercontent.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/Main/_code/FOV%20QC/FOV%20QC%20utils.R")

## Gene -> barcode map for every panel
allbarcodes <- readRDS(url("https://github.com/Nanostring-Biostats/CosMx-Analysis-Scratch-Space/raw/Main/_code/FOV%20QC/barcodes_by_panel.RDS"))
names(allbarcodes)                    # Hs_IO, Hs_UCC, Hs_6k, Mm_Neuro, Mm_UCC, Hs_WTX
barcodemap <- allbarcodes$Hs_6k       # <------- 6k panel (vignettes use Hs_WTX)
head(barcodemap)

## Genes shared by the panel map and the counts matrix (runFOVQC uses only these)
length(intersect(barcodemap$gene, colnames(counts)))


## =============================================================================
## CHUNK 26: run FOV QC (vignette 1 call, default thresholds)
## =============================================================================
fovqcresult <- runFOVQC(counts = counts, xy = xy, fov = metadata$fov, barcodemap = barcodemap,
                        max_prop_loss = max_prop_loss,
                        max_totalcounts_loss = max_totalcounts_loss)
saveRDS(fovqcresult, file.path(dir_fov, "fovqcresult.RDS"))

## FOVs flagged for any reason, for loss of signal (low counts), for bias (vignette)
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
  mtext(paste0("All sections: ", tot[[col_n]], " of ", tot$n_fov_total, " ", what),
        side = 3, line = 0.2, cex = 0.8)
  for (sec in section_names) {
    inds <- metadata$section_id == sec
    r    <- fov_qc_summary[section_id == sec]
    text(mean(range(xy[inds, 1])), max(xy[inds, 2]), pos = 3, xpd = NA, cex = 0.7, font = 2,
         labels = paste0(sec, ": ", r[[col_n]], " of ", r$n_fov_total, " FOVs"))
  }
}


## =============================================================================
## CHUNK 29: PLOT 1 - map of flagged FOVs (vignette 1 and 2: mapFlaggedFOVs)
## =============================================================================
## Blue = all FOVs; red with FOV number = flagged for low counts OR bias.
draw_map <- function() {
  op <- par(mar = c(4, 4, 4, 1)); on.exit(par(op))
  mapFlaggedFOVs(fovqcresult)
  add_fov_labels("any")
}
show_and_save(draw_map, file.path(dir_fov, "01_flagged_FOVs_map.png"), map_w, map_h)


## =============================================================================
## CHUNK 30: PLOT 2 - signal loss (low counts) (vignette 1 and 2: FOVSignalLossSpatialPlot)
## =============================================================================
## Cells coloured by log2 fold-change in total counts vs similar regions in
## other FOVs; yellow boxes = FOVs flagged for low counts.
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
flags_rc <- fovqcresult$flags_per_fov_x_reportercycle
flags_rc_flagged <- flags_rc[rowSums(flags_rc >= 0.5) > 0, , drop = FALSE]
print(flags_rc_flagged)
fwrite(as.data.table(flags_rc_flagged, keep.rownames = "fov"),
       file.path(dir_fov, "FOV_QC_flagged_reportercycles.csv"))


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
  ## vignette's scale_x_discrete(breaks = 1:13) was specific to its data, so not used
  print(p)
  ggsave(p, filename = file.path(dir_fov, "04_genes_flagged_in_n_FOVs_barplot.png"),
         width = 5, height = 5, dpi = 150)                          # vignette: 4 x 5 in, 150 dpi
} else {
  message("No FOVs flagged for bias: no genes flagged, bar plot not drawn.")
}


## =============================================================================
## CHUNK 33: FOVs flagged, on screen
## =============================================================================
message("Whole slide: ", tot$flagged_label, " (", tot$pct_fov_flagged, "%) - low counts: ",
        tot$n_fov_low_counts, ", bias: ", tot$n_fov_bias)
for (sec in section_names) {
  r <- fov_qc_summary[section_id == sec]
  message(sec, ": ", r$flagged_label, " (", r$pct_fov_flagged, "%) - low counts: ",
          r$n_fov_low_counts, ", bias: ", r$n_fov_bias,
          if (nchar(r$flagged_fov_ids) > 0) paste0(" - FOVs ", r$flagged_fov_ids) else "")
}
message("All FOV QC outputs saved to ", dir_fov)














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
flags_rc <- fovqcresult$flags_per_fov_x_reportercycle
flags_rc_flagged <- flags_rc[rowSums(flags_rc >= 0.5) > 0, , drop = FALSE]
print(flags_rc_flagged)
fwrite(as.data.table(flags_rc_flagged, keep.rownames = "fov"),
       file.path(dir_fov, "FOV_QC_flagged_reportercycles.csv"))


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
## 06_ncount_nfeature.R
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
  thr_txt <- paste0("counts ", r$counts_min_threshold, "-", fmt(r$counts_max_threshold),
                    "; features ", r$features_min_threshold, "-", fmt(r$features_max_threshold))
  ttl <- paste0(tag, ": ", r$filtered_label, " (", r$pct_filtered, "%)")
  
  ## Joint distribution (vignette)
  p <- ggplot(data = meta, aes(x = nFeature_RNA, y = nCount_RNA, color = counts_filter)) +
    geom_point(alpha = 0.002, pch = ".") +
    scale_color_manual(values = c("pass" = "dodgerblue", "filtered" = "darkorange")) +
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
    geom_point(aes(x = x_slide_mm, y = y_slide_mm, color = counts_filter), size = 0.001, alpha = 0.02) +
    coord_fixed() + theme_bw() +
    scale_color_manual(values = c("pass" = "dodgerblue", "filtered" = "darkorange")) +
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
thr_lines <- nc_nf_sec_summary[group %in% section_names]
thr_lines[, section_id := paste0(group, ": ", filtered_label)]
md_plot$section_id <- thr_lines$section_id[match(md_plot$section_id, thr_lines$group)]  # panel titles with counts
p_grid <- ggplot(md_plot, aes(nFeature_RNA, nCount_RNA, color = counts_filter)) +
  geom_point(alpha = 0.01, pch = ".") +
  scale_color_manual(values = c("pass" = "dodgerblue", "filtered" = "darkorange")) +
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
## 09_QC_Proseg.R including run on ALL Cells
## =============================================================================
## PURPOSE
##   Segmentation QC with ProSeg, continuing from scripts 02-08.
##   Run two ways, each saved in its own folder under QC/QC_Proseg/:
##     all_cells/            ONE ProSeg run on the transcripts of all 6 sections
##     per_tissue_section/   SIX ProSeg runs, one per section
##   Results are saved as files only; metadata is NOT changed.
##
## WHAT PROSEG DOES
##   AtoMx assigned each transcript to a cell from image-based outlines.
##   ProSeg re-segments the same transcripts from their positions and gene
##   identities (probabilistic segmentation), starting from AtoMx's assignment.
##   If AtoMx segmented a cell well, ProSeg keeps its transcripts together.
##
## QC METRIC: purity, per AtoMx cell (same rule as the earlier ProSeg script)
##   purity = transcripts in the ProSeg cell that received most of them
##            / all of the cell's transcripts that ProSeg assigned to a cell
##   flagged       = purity < 0.5 (fewer than half the transcripts stay together)
##   not evaluated = < 3 evaluable transcripts (purity NA; NOT flagged)
##
## ProSeg options checked against the ProSeg source (v3.2.0): all flags used
## below exist; transcript-metadata output columns used here are
## transcript_id, observed_x, observed_y, assignment, background.
##
## REQUIREMENT: ProSeg installed (Terminal: cargo install proseg).
## The all-cells run holds every transcript in memory at once, so it needs
## more RAM and time than any single-section run.
## =============================================================================


## =============================================================================
## CHUNK 1: packages, paths, settings
## =============================================================================
suppressPackageStartupMessages({
  library(data.table)   # fast reading/writing of the large transcript file
  library(ggplot2)
})

output_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"
slide_dir  <- "/Users/maryamuhabilausman/Desktop/AAA_CosMx_Analysis/Downloaded_data/flatFiles/AC20260629_AAA3"

dir_ps <- file.path(output_dir, "QC", "QC_Proseg")
approaches <- c("all_cells", "per_tissue_section")
for (a in approaches) for (s in c("01_prepared_transcripts", "02_proseg_raw_output", "03_results"))
  dir.create(file.path(dir_ps, a, s), showWarnings = FALSE, recursive = TRUE)

pixel_size_um        <- 0.12028   # µm per CosMx pixel (ProSeg's CosMx value)
purity_threshold     <- 0.5       # flag cells with purity below this
min_transcripts_eval <- 3         # fewer evaluable transcripts -> not evaluated
coord_tolerance_um   <- 0.05      # allowed rounding difference in the linking check
excluded_genes_regex <- "^(FalseCode|SystemControl|NegPrb|Negative)"   # control probes (ProSeg's CosMx default)
max_transcript_nucleus_distance <- NULL   # NULL = ProSeg default (60 µm)
n_threads    <- NULL              # NULL = all CPU cores; e.g. 6 keeps the laptop usable
rerun_proseg <- FALSE             # FALSE = reuse ProSeg output that already exists

## R must find the proseg program installed by cargo
cargo_bin <- file.path(Sys.getenv("HOME"), ".cargo", "bin")
if (dir.exists(cargo_bin) && !grepl(cargo_bin, Sys.getenv("PATH"), fixed = TRUE))
  Sys.setenv(PATH = paste(cargo_bin, Sys.getenv("PATH"), sep = ":"))
stopifnot("proseg not found - install with: cargo install proseg" = nzchar(Sys.which("proseg")))

fmt <- function(n) format(n, big.mark = ",")


## =============================================================================
## CHUNK 2: free memory, then cell metadata (read only, not modified)
## =============================================================================
## The transcript file has ~73 million rows. R's default limit is 16 GB, and the
## objects left from scripts 02-10 (counts, Seurat objects, per-section copies)
## already use much of it. Everything they hold is saved on disk, so remove them
## from memory and keep only what this script needs.
rm(list = setdiff(ls(), c("output_dir", "slide_dir", "dir_ps", "approaches", "pixel_size_um",
                          "purity_threshold", "min_transcripts_eval", "coord_tolerance_um",
                          "excluded_genes_regex", "max_transcript_nucleus_distance",
                          "n_threads", "rerun_proseg", "fmt")))
invisible(gc())

metadata <- readRDS(file.path(output_dir, "metadata_unfiltered.RDS"))
stopifnot(identical(rownames(metadata), metadata$cell_id))
section_names <- sort(unique(metadata$section_id))

## FOV -> section: every transcript inherits its FOV's section
fov_to_section <- unique(data.table(fov = metadata$fov, section_id = metadata$section_id))
stopifnot("an FOV lies in two sections" = !anyDuplicated(fov_to_section$fov))

## Cell names are "c_<slide>_<FOV>_<cell>"; the prefix "c_<slide>" is rebuilt per transcript
cell_prefix <- unique(sub("_[^_]+_[^_]+$", "", metadata$cell_id))
stopifnot("cell_id prefix not unique" = length(cell_prefix) == 1)
cell_prefix


## =============================================================================
## CHUNK 3: load transcripts and link each to its cell
## =============================================================================
## One row per detected transcript. Only the needed columns are read (memory).
txfile <- list.files(slide_dir, pattern = "tx_file.*\\.gz$", full.names = TRUE)[1]
stopifnot("tx_file not found" = !is.na(txfile))

cols_needed <- c("fov", "cell_ID", "x_global_px", "y_global_px", "target", "z", "CellComp")
cols_read   <- intersect(cols_needed, names(fread(txfile, nrows = 0)))
tx <- fread(txfile, select = cols_read)
nrow(tx); colnames(tx)
stopifnot(all(c("fov", "cell_ID", "x_global_px", "y_global_px", "target") %in% colnames(tx)))

## Control probes are NOT removed here: removing rows makes a full copy of the
## 73-million-row table (this is what hit the 16 GB limit). ProSeg drops them
## itself (--excluded-genes), so they never reach its output or the purity score.
control_targets <- grep(excluded_genes_regex, unique(tx$target), value = TRUE)
n_control <- tx[target %chin% control_targets, .N]
message("Control-probe transcripts (dropped by ProSeg): ", fmt(n_control))

## Section per transcript (join by reference, no copy)
tx[fov_to_section, on = "fov", section_id := i.section_id]
if (anyNA(tx$section_id)) tx <- tx[!is.na(section_id)]
table(tx$section_id)

## Cell name per transcript; cell_ID = 0 = outside every cell
tx[, seurat_cell := ""]
tx[cell_ID != 0, seurat_cell := paste0(cell_prefix, "_", fov, "_", cell_ID)]
assigned_names <- unique(tx$seurat_cell[tx$seurat_cell != ""])
pct_matched <- 100 * mean(assigned_names %in% metadata$cell_id)
message("Transcript cell names found in metadata: ", round(pct_matched, 2), "%")
stopifnot("transcript cell names do not match metadata cell names" = pct_matched > 95)
rm(assigned_names); gc()

## Optional columns ProSeg can use
z_col <- if ("z" %in% colnames(tx)) "z" else NA
has_compartment <- "CellComp" %in% colnames(tx)


## =============================================================================
## CHUNK 4: write ProSeg input files (all cells: 1 file; per section: 6 files)
## =============================================================================
## row_id = unique number per transcript that ProSeg carries to its output
## (--transcript-id-column), so each output transcript links back to its cell.
keep_cols <- c("fov", "seurat_cell", "cell_ID", "x_global_px", "y_global_px", "target")
if (!is.na(z_col))   keep_cols <- c(keep_cols, z_col)
if (has_compartment) keep_cols <- c(keep_cols, "CellComp")

write_input <- function(d, dir_in, name) {
  d[, row_id := seq_len(.N)]
  saveRDS(d[, .(row_id, seurat_cell, x_global_px, y_global_px)],
          file.path(dir_in, paste0(name, "_row_id_lookup.RDS")))
  fwrite(d, file.path(dir_in, paste0(name, "_tx_prepared.csv.gz")), compress = "gzip")
  n_rna <- d[!target %chin% control_targets, .N]           # transcripts ProSeg will use
  message(name, ": ", fmt(nrow(d)), " transcripts written (", fmt(n_rna), " excluding control probes)")
  n_rna
}

runs <- list(all_cells = "all_cells", per_tissue_section = section_names)   # run names per approach
n_sent <- list()

## Per section first (each subset is a small copy)
n_sent[["per_tissue_section"]] <- sapply(section_names, function(sec)
  write_input(tx[section_id == sec, ..keep_cols],
              file.path(dir_ps, "per_tissue_section", "01_prepared_transcripts"), sec))

## All cells: edit tx in place instead of copying it (the full table is very large)
drop_cols <- setdiff(colnames(tx), keep_cols)          # e.g. section_id
if (length(drop_cols)) tx[, (drop_cols) := NULL]
setcolorder(tx, keep_cols)
n_sent[["all_cells"]] <- c(all_cells = write_input(tx,
                                                   file.path(dir_ps, "all_cells", "01_prepared_transcripts"), "all_cells"))
rm(tx); gc()


## =============================================================================
## CHUNK 5: run ProSeg (all cells once, then each section)
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
  a <- if (!is.na(z_col)) c(a, "--z-column", z_col) else c(a, "--z-column", "x_global_px", "--ignore-z-coord")
  if (has_compartment) a <- c(a, "--compartment-column", "CellComp", "--compartment-nuclear", "Nuclear")
  if (!is.null(max_transcript_nucleus_distance)) a <- c(a, "--max-transcript-nucleus-distance", max_transcript_nucleus_distance)
  if (!is.null(n_threads)) a <- c(a, "--nthreads", n_threads)
  a
}
quote_for_shell <- function(x) ifelse(x == "" | grepl("[^A-Za-z0-9_./=-]", x), shQuote(x), x)

for (a in approaches) for (nm in runs[[a]]) {
  out_path <- file.path(dir_ps, a, "02_proseg_raw_output", nm); dir.create(out_path, showWarnings = FALSE)
  tm_path  <- file.path(out_path, "transcript-metadata.csv")
  if (file.exists(tm_path) && !rerun_proseg) { message(a, "/", nm, ": output exists - reused"); next }
  message(a, "/", nm, ": running ProSeg (can take a long time)...")
  in_path  <- file.path(dir_ps, a, "01_prepared_transcripts", paste0(nm, "_tx_prepared.csv.gz"))
  log_lines <- system2("proseg", quote_for_shell(build_proseg_args(in_path, out_path)), stdout = TRUE, stderr = TRUE)
  writeLines(log_lines, file.path(out_path, "proseg_run_log.txt"))
  if (!file.exists(tm_path)) warning(a, "/", nm, ": ProSeg wrote no output - see proseg_run_log.txt")
}


## =============================================================================
## CHUNK 6: link ProSeg output back to cells and compute purity
## =============================================================================
score_run <- function(a, nm) {
  tm_path <- file.path(dir_ps, a, "02_proseg_raw_output", nm, "transcript-metadata.csv")
  if (!file.exists(tm_path)) { warning(a, "/", nm, ": no ProSeg output - skipped"); return(NULL) }
  ps <- fread(tm_path, select = c("transcript_id", "observed_x", "observed_y", "assignment", "background"))
  lookup <- readRDS(file.path(dir_ps, a, "01_prepared_transcripts", paste0(nm, "_row_id_lookup.RDS")))
  
  ## Join on ID, then confirm by coordinates (ProSeg position = pixel position x 0.12028)
  linked  <- merge(ps, lookup, by.x = "transcript_id", by.y = "row_id")
  max_dev <- max(pmax(abs(linked$observed_x - linked$x_global_px * pixel_size_um),
                      abs(linked$observed_y - linked$y_global_px * pixel_size_um)))
  if (nrow(linked) != nrow(ps) || max_dev > coord_tolerance_um) {
    warning(a, "/", nm, ": linking check failed (max deviation ", signif(max_dev, 3), " µm) - skipped")
    return(NULL)
  }
  
  ## Evaluable transcript: in an AtoMx cell, not background, assigned to a ProSeg cell
  bg <- linked$background
  if (!is.logical(bg)) bg <- tolower(as.character(bg)) %in% c("true", "1")
  ev <- linked[seurat_cell != "" & !bg & !is.na(assignment)]
  
  p <- ev[, .N, by = .(seurat_cell, assignment)][
    , .(n_transcripts_evaluated = sum(N), n_in_top_proseg_cell = max(N),
        n_proseg_cells_spanned = .N), by = seurat_cell]
  p[, purity := n_in_top_proseg_cell / n_transcripts_evaluated]
  p[n_transcripts_evaluated < min_transcripts_eval, purity := NA_real_]
  
  sent <- n_sent[[a]][[nm]]
  report <- data.table(approach = a, run = nm, transcripts_sent_excl_controls = sent,
                       transcripts_kept_by_proseg = nrow(ps),
                       pct_dropped_far_from_cells = round(100 * (1 - nrow(ps) / sent), 2),
                       transcripts_evaluable = nrow(ev),
                       max_coordinate_deviation_um = signif(max_dev, 3))
  message(a, "/", nm, ": linked and verified; ", fmt(nrow(p)), " cells scored")
  list(p = p, report = report)
}

## One per-cell table per approach (every cell in metadata, same order)
ps_results <- list(); link_reports <- list()
for (a in approaches) {
  sc <- lapply(runs[[a]], function(nm) score_run(a, nm))
  p  <- rbindlist(lapply(sc, `[[`, "p"))
  link_reports[[a]] <- rbindlist(lapply(sc, `[[`, "report"))
  i  <- match(metadata$cell_id, p$seurat_cell)
  d  <- data.table(cell_id = metadata$cell_id, section_id = metadata$section_id, fov = metadata$fov,
                   x_slide_mm = metadata$x_slide_mm, y_slide_mm = metadata$y_slide_mm,
                   nCount_RNA = metadata$nCount_RNA,
                   n_transcripts_evaluated = p$n_transcripts_evaluated[i],
                   n_in_top_proseg_cell    = p$n_in_top_proseg_cell[i],
                   n_proseg_cells_spanned  = p$n_proseg_cells_spanned[i],
                   purity                  = p$purity[i])
  d[, evaluated := !is.na(purity)]
  d[, flagged := evaluated & purity < purity_threshold]          # not evaluated -> not flagged
  d[, status := fifelse(!evaluated, "not evaluated",
                        fifelse(flagged, "flagged (purity < 0.5)", "pass"))]
  ps_results[[a]] <- d
}


## =============================================================================
## CHUNK 7: save results - per-cell CSV, per-section summary, linking report
## =============================================================================
summarise_ps <- function(d) {
  s <- d[, .(n_cells = .N, n_evaluated = sum(evaluated), n_not_evaluated = sum(!evaluated),
             median_purity = round(median(purity, na.rm = TRUE), 3),
             n_flagged = sum(flagged),
             pct_flagged_of_all_cells = round(100 * mean(flagged), 2),
             pct_flagged_of_evaluated = round(100 * sum(flagged) / sum(evaluated), 2)),
         by = section_id][order(section_id)]
  s <- rbind(s, d[, .(section_id = "All sections", n_cells = .N, n_evaluated = sum(evaluated),
                      n_not_evaluated = sum(!evaluated),
                      median_purity = round(median(purity, na.rm = TRUE), 3), n_flagged = sum(flagged),
                      pct_flagged_of_all_cells = round(100 * mean(flagged), 2),
                      pct_flagged_of_evaluated = round(100 * sum(flagged) / sum(evaluated), 2))])
  s[, flagged_label := paste0(fmt(n_flagged), " of ", fmt(n_cells), " cells flagged")]
  s
}

ps_summary <- list()
for (a in approaches) {
  res_dir <- file.path(dir_ps, a, "03_results")
  ps_summary[[a]] <- summarise_ps(ps_results[[a]])
  fwrite(ps_results[[a]],   file.path(res_dir, paste0("ProSeg_purity_per_cell_", a, ".csv")))
  fwrite(ps_results[[a]][flagged == TRUE], file.path(res_dir, paste0("ProSeg_flagged_cells_", a, ".csv")))
  fwrite(ps_summary[[a]],   file.path(res_dir, paste0("ProSeg_summary_per_section_", a, ".csv")))
  fwrite(link_reports[[a]], file.path(res_dir, paste0("ProSeg_linking_report_", a, ".csv")))
  saveRDS(setNames(ps_results[[a]]$flagged, ps_results[[a]]$cell_id),
          file.path(res_dir, paste0("flag_ProSeg_low_purity_", a, ".RDS")))
  message("\n", a, ":"); print(ps_summary[[a]]); print(link_reports[[a]])
}


## =============================================================================
## CHUNK 8: plots - slide (all 6 sections), grid, histograms, per section
## =============================================================================
st_cols <- c("pass" = "grey80", "not evaluated" = "grey40", "flagged (purity < 0.5)" = "red")
sec_pos <- as.data.table(metadata)[, .(x = mean(range(x_slide_mm)), y = max(y_slide_mm)), by = section_id]

for (a in approaches) {
  d   <- ps_results[[a]]; s <- ps_summary[[a]]; res_dir <- file.path(dir_ps, a, "03_results")
  d[, status := factor(status, levels = names(st_cols))]
  dp  <- d[order(status)]                                          # flagged drawn on top
  lab <- merge(sec_pos, s[, .(section_id, n_flagged, n_cells, pct_flagged_of_evaluated)], by = "section_id")
  lab[, label := paste0(section_id, ": ", fmt(n_flagged), " of ", fmt(n_cells),
                        " (", pct_flagged_of_evaluated, "% of evaluated)")]
  tot <- s[section_id == "All sections"]
  ttl <- paste0("ProSeg low purity - ", a, ": ", tot$flagged_label)
  
  ## Whole slide, sections in position, counts above each section
  p_slide <- ggplot(dp, aes(x_slide_mm, y_slide_mm, colour = status)) +
    geom_point(size = 0.05) + scale_colour_manual(values = st_cols, drop = FALSE) +
    geom_label(data = lab, aes(x, y, label = label), inherit.aes = FALSE, vjust = -0.3,
               size = 3, fontface = "bold", fill = "white", label.size = 0.2) +
    scale_y_continuous(expand = expansion(mult = c(0.03, 0.08))) +
    coord_fixed() + theme_bw() + guides(colour = guide_legend(override.aes = list(size = 3))) +
    labs(title = ttl, x = "x (mm)", y = "y (mm)", colour = NULL)
  ggsave(file.path(res_dir, paste0("01_ProSeg_flagged_cells_on_slide_", a, ".png")), p_slide, width = 12, height = 12, dpi = 300)
  print(p_slide)
  
  ## Grid: one panel per section
  pl <- setNames(paste0(s$section_id, ": ", s$flagged_label), s$section_id)
  dp[, panel := pl[section_id]]
  p_grid <- ggplot(dp, aes(x_slide_mm, y_slide_mm, colour = status)) +
    geom_point(size = 0.05) + scale_colour_manual(values = st_cols, drop = FALSE) +
    facet_wrap(~panel, nrow = 2, scales = "free") + theme_bw() +
    guides(colour = guide_legend(override.aes = list(size = 3))) +
    labs(title = ttl, x = "x (mm)", y = "y (mm)", colour = NULL)
  ggsave(file.path(res_dir, paste0("02_ProSeg_flagged_cells_per_section_grid_", a, ".png")), p_grid, width = 16, height = 10, dpi = 300)
  
  ## Purity histograms per section
  p_hist <- ggplot(d[evaluated == TRUE], aes(purity)) +
    geom_histogram(bins = 50, fill = "steelblue") +
    geom_vline(xintercept = purity_threshold, colour = "red", linetype = "dashed") +
    facet_wrap(~section_id, nrow = 2, scales = "free_y") + theme_classic() +
    labs(title = paste0("ProSeg purity (evaluated cells) - ", a), x = "ProSeg purity", y = "N cells")
  ggsave(file.path(res_dir, paste0("03_ProSeg_purity_histograms_", a, ".png")), p_hist, width = 14, height = 8, dpi = 300)
  
  ## One spatial plot per section
  dir_ind <- file.path(res_dir, "per_section_plots"); dir.create(dir_ind, showWarnings = FALSE)
  for (sec in section_names) {
    r <- s[section_id == sec]
    ggsave(file.path(dir_ind, paste0(sec, "_ProSeg_flagged_cells_", a, ".png")),
           ggplot(dp[section_id == sec], aes(x_slide_mm, y_slide_mm, colour = status)) +
             geom_point(size = 0.1) + scale_colour_manual(values = st_cols, drop = FALSE) +
             coord_fixed() + theme_bw() + guides(colour = guide_legend(override.aes = list(size = 3))) +
             labs(title = paste0(sec, " (", a, "): ", r$flagged_label),
                  subtitle = paste0(fmt(r$n_evaluated), " evaluated, median purity ", r$median_purity),
                  x = "x (mm)", y = "y (mm)", colour = NULL),
           width = 8, height = 8, dpi = 300)
  }
}


## =============================================================================
## CHUNK 9: results on screen
## =============================================================================
for (a in approaches) {
  message("\n", a, ":")
  for (sec in c(section_names, "All sections")) {
    r <- ps_summary[[a]][section_id == sec]
    message("  ", sec, ": ", r$flagged_label, " (", r$pct_flagged_of_evaluated, "% of ",
            fmt(r$n_evaluated), " evaluated)")
  }
}
message("ProSeg QC saved to ", dir_ps)



















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


## =============================================================================
## CHUNK 52: packages, inputs, folders
## =============================================================================
## Needed once:
# install.packages(c("scales", "spatstat.geom", "devtools"))
# devtools::install_github("https://github.com/Nanostring-Biostats/InSituCor")
library(scales)       # viridis_pal() (vignette)
library(data.table)

output_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"
## Fresh session only:
# metadata <- readRDS(file.path(output_dir, "metadata_unfiltered.RDS"))
# xy       <- readRDS(file.path(output_dir, "xy_unfiltered.RDS"))
stopifnot(identical(rownames(xy), metadata$cell_id))
section_names <- sort(unique(metadata$section_id))

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
## 11_add_condition_and_filter_seurat.R filter 03-09
## =============================================================================
## PURPOSE
##   1. Original Seurat object + 2 new metadata columns, saved under a NEW name
##      (seurat_with_sections.RDS itself is never overwritten or updated):
##        condition = "AAA" for every cell (TAA added later)
##        sample_id = AAA section of each cell (AAA2, AAA38, AAA40, AAA45,
##                    AAA46, AAA48), copied from section_id used in scripts 02-10
##   2. A new, QC-filtered Seurat object: the same counts and the same metadata
##      rows, for the cells KEPT by per-tissue-section QCs 03-09 (script 10, Part A).
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
flag_file  <- file.path(qc_dir_A, "final_QC_flag_per_tissue_section_03_to_09.RDS")  # script 10
flag_csv   <- file.path(qc_dir_A, "cells", "01_removed_cells_summary_per_section.csv")  # script 10

seurat_out_all      <- file.path(output_dir, "seurat_with_sections_condition.RDS")          # all cells + 2 columns
seurat_out_filtered <- file.path(output_dir, "seurat_QC_filtered_per_section_03_to_09.RDS")  # kept cells only
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

saveRDS(obj, seurat_out_all)                             # original + 2 columns, nothing else changed
message("Saved (all cells + condition, sample_id): ", seurat_out_all)


## =============================================================================
## CHUNK 3: which cells to keep (per-tissue-section QCs 03-09)
## =============================================================================
removed <- readRDS(flag_file)                            # TRUE = removed, named by cell_id
stopifnot(setequal(names(removed), colnames(obj)))
removed <- removed[colnames(obj)]                        # align by cell name
keep_cells <- colnames(obj)[!removed]                    # kept cells, original order

## Removed per section must equal script 10's summary (the plot numbers)
check <- data.table(section_id = unname(obj$sample_id), removed = unname(removed))[
  , .(n_cells = .N, n_removed = sum(removed), n_kept = sum(!removed)), by = section_id][order(section_id)]
s10 <- fread(flag_csv)[section_id != "All sections", .(section_id, n_removed_script10 = n_removed_independent)]
check <- merge(check, s10, by = "section_id")
print(check)
stopifnot("removed counts differ from script 10" = all(check$n_removed == check$n_removed_script10))


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
stopifnot(all(after$n_cells_after_QC == after$n_cells_before_QC - after$n_removed))
print(after)
obj_filtered


## =============================================================================
## CHUNK 6: save
## =============================================================================
saveRDS(obj_filtered, seurat_out_filtered)
fwrite(after, file.path(qc_dir_A, "cells_before_and_after_QC_per_section.csv"))
message("Saved (QC-filtered, per-section QCs 03-09): ", seurat_out_filtered)
message("Cells: ", fmt(ncol(obj)), " before, ", fmt(sum(removed)), " removed, ", fmt(ncol(obj_filtered)), " kept")

## Downstream examples (not run):
# obj_filtered <- readRDS(seurat_out_filtered)
# sub <- subset(obj_filtered, subset = sample_id %in% c("AAA48", "AAA38"))   # chosen AAA sections
# Idents(obj_filtered) <- "sample_id"                                         # compare sections





















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
## 12_Normalisation_PCA_UMAP_Leiden.R
## =============================================================================
## PURPOSE
## -------
## Continues from 11_add_condition_and_filter_seurat.R. Input = the QC-filtered
## Seurat object (per-tissue-section QCs removed), e.g.
##   seurat_QC_filtered_per_section_03_to_09.RDS
## The input file is only read, never changed.
##
## TWO ANALYSES, each in its own folder inside 12_Normalisation_PCA_UMAP_Leiden/
##   PRIMARY     (01_...): counts -> HVGs -> batch-aware Pearson-residual PCA
##                         -> neighbour graph -> UMAP   (layout of the graph)
##                                            -> Leiden (clusters the graph)
##   SENSITIVITY (02_...): counts -> HVGs -> batch-aware Pearson-residual PCA
##                         -> Harmony -> neighbour graph -> UMAP / Leiden
##   Shared steps (HVGs + batch-aware PCA) are computed ONCE (00_...) and reused,
##   so the two analyses differ ONLY by Harmony. Comparison in 03_...
##   Batch = sample_id (the 6 AAA tissue sections).
##
## CHANGING n PCs / n NEIGHBOURS / RESOLUTION (CHUNK 1)
##   The PCA is saved once and reused. Changing n_pcs or n_neighbors writes NEW
##   folders named ..._PC<n>_NN<k>, so earlier runs are never overwritten.
##   Resolutions: leiden_resolution (main) + resolution_sweep (all on one graph).
##
## RESOLUTION SWEEP VISUALS (in each branch folder)
##   Leiden_resolution_sweep_summary.tiff / _UMAP_panels.tiff / _cluster_tree.tiff
##
## GRAPH ARCHITECTURE - the same structure as scanpy
## ---------------------------------------------------------------------------
##   scanpy:  sc.pp.neighbors()  -> ONE graph (adata.obsp['connectivities'])
##            sc.tl.umap()        lays that graph out in 2D
##            sc.tl.leiden()      clusters that SAME graph
##   here:    STEP A  kNN search on the PCs (cosine, 15 neighbours incl. self)
##            STEP B  uwot::similarity_graph(kNN) -> ONE graph G
##                    (UMAP fuzzy connectivity weights = what sc.pp.neighbors
##                     stores as 'connectivities')
##            STEP C  uwot::optimize_graph_layout(G) -> UMAP coordinates
##            STEP D  igraph::cluster_leiden(G)       -> Leiden clusters
##   UMAP is NOT run to create the graph, and Leiden NEVER sees UMAP
##   coordinates. UMAP is only a 2D picture of G.
##
## PARAMETERS
## ---------------------------------------------------------------------------
##   HVGs: 2000, Seurat FindVariableFeatures (vignette: "2000 for 6k panel").
##   Batch-aware Pearson PCA: scPearsonPCA, exactly the vignette call
##     (scale.max = 10, do.scale = TRUE, do.center = TRUE; 50 PCs = package default),
##     batch = sample_id.
##   Neighbours / UMAP / Leiden: the vignette's Pearson-residual workflow
##     (sc.pp.neighbors n_neighbors = 15, n_pcs = 40, metric = cosine;
##      sc.tl.umap min_dist = 0.02, spread = 2;
##      sc.tl.leiden resolution = 0.5, flavor igraph, n_iterations = 2).
##   All are set in CHUNK 1 and saved to parameters_used.csv.
##
## Nothing is removed from the data in this script. Results are also added to a
## COPY of the Seurat object (reductions + cluster columns), saved under a new name.
## =============================================================================

## =============================================================================
## CHUNK 0: packages - checked, NOT silently installed
## =============================================================================
required_pkgs <- c("data.table", "Matrix", "ggplot2", "scales", "Seurat",
                   "scPearsonPCA", "RSpectra", "uwot", "RcppAnnoy", "igraph", "harmony")   # RSpectra: used inside scPearsonPCA
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall with:\n",
       "  install.packages(c('data.table','Matrix','ggplot2','scales','Seurat','RSpectra','uwot','RcppAnnoy','igraph','harmony'))\n",
       "  remotes::install_github('Nanostring-Biostats/CosMx-Analysis-Scratch-Space',\n",
       "                          subdir = '_code/scPearsonPCA', ref = 'Main')")
}
## functions this script depends on - stop now rather than fail halfway
if (!all(c("similarity_graph", "optimize_graph_layout") %in% getNamespaceExports("uwot"))) {
  stop("Your uwot version (", utils::packageVersion("uwot"), ") lacks similarity_graph() / ",
       "optimize_graph_layout(). Update it: install.packages('uwot')")
}
if (utils::packageVersion("harmony") < "1.0") {
  stop("harmony >= 1.0 needed (RunHarmony on a PC matrix, as in the vignette). install.packages('harmony')")
}
## the sweep cluster-tree uses the 'linewidth' aesthetic (ggplot2 >= 3.4)
if (utils::packageVersion("ggplot2") < "3.4.0") {
  stop("ggplot2 >= 3.4.0 needed for the resolution-sweep cluster tree. install.packages('ggplot2')")
}

library(data.table)
library(Matrix)
library(ggplot2)
library(Seurat)
library(scPearsonPCA)
library(RcppAnnoy)

## =============================================================================
## CHUNK 1: paths + parameters
## =============================================================================
output_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"
seurat_in  <- file.path(output_dir, "seurat_QC_filtered_per_section_03_to_09.RDS")   # script 11 (read only)
input_tag  <- sub("\\.RDS$", "", basename(seurat_in))    # folder name follows the input object

batch_col           <- "sample_id"    # the 6 AAA tissue sections (AAA2 ... AAA48)
n_hvgs              <- 2000           # vignette: 2000 for a 6k panel
n_pcs               <- 20             # vignette: sc.pp.neighbors(n_pcs = 40)            <- change here
n_neighbors         <- 250          # vignette: sc.pp.neighbors(n_neighbors = 15), incl. the cell itself  <- change here
rerun_pca           <- FALSE          # FALSE = reuse the saved HVGs + PCA; TRUE = recompute

dir_main    <- file.path(output_dir, "12_Normalisation_PCA_UMAP_Leiden", input_tag)
dir_shared  <- file.path(dir_main, "00_shared_HVGs_and_batch_aware_Pearson_PCA")
run_tag     <- paste0("PC", n_pcs, "_NN", n_neighbors)                                   # e.g. PC40_NN15
dir_primary <- file.path(dir_main, paste0("01_PRIMARY_batchawarePCA_", run_tag))
dir_sens    <- file.path(dir_main, paste0("02_SENSITIVITY_batchawarePCA_Harmony_", run_tag))
dir_compare <- file.path(dir_main, paste0("03_PRIMARY_vs_SENSITIVITY_", run_tag))
for (d in c(dir_main, dir_shared, dir_primary, dir_sens, dir_compare)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

annoy_n_trees       <- 50             # approximate kNN index size (uwot default)
umap_min_dist       <- 0.02           # vignette: sc.tl.umap(min_dist = 0.02)
umap_spread         <- 2              # vignette: sc.tl.umap(spread = 2)
leiden_resolution   <- 1            # vignette: sc.tl.leiden(resolution = 0.5)          <- change here
leiden_n_iterations <- 2              # vignette: sc.tl.leiden(n_iterations = 2)
resolution_sweep    <- c(0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0)  # resolutions compared on the same graph <- change here
dominant_section_flag_pct <- 90       # flag (for investigation only) clusters with > 90% of cells from one section
sweep_umap_max_cells   <- 200000      # PLOTTING ONLY: cells drawn per panel in the sweep UMAP grid (random subsample; all cells used for clustering and labels)
sweep_tree_min_in_prop <- 0.05        # PLOTTING ONLY: cluster-tree edges drawn if >= 5% of the child cluster came from that parent
seed    <- 42
dpi_out <- 300

params <- data.table(
  parameter = c("input_metadata", "input_counts", "batch_variable", "n_HVGs", "HVG_method",
                "PCA", "PCA_scale.max", "PCA_do.scale", "PCA_do.center", "totalcounts",
                "n_PCs_for_neighbours", "n_neighbors_incl_self", "neighbour_metric", "kNN_method",
                "graph", "UMAP", "UMAP_min_dist", "UMAP_spread", "UMAP_init",
                "clustering", "Leiden_objective", "Leiden_resolution", "Leiden_n_iterations",
                "Harmony (sensitivity only)", "seed"),
  value = c(basename(seurat_in), "RNA counts layer of the input Seurat object",
            batch_col, n_hvgs, "Seurat::FindVariableFeatures (vst)",
            "scPearsonPCA::sparse_quasipoisson_pca_seurat_batch (batch-aware Pearson residuals)",
            10, TRUE, TRUE, "per-cell sum over all genes (= nCount_RNA)",
            n_pcs, n_neighbors, "cosine", paste0("RcppAnnoy angular, ", annoy_n_trees, " trees"),
            "uwot::similarity_graph on the kNN (UMAP fuzzy connectivities, as sc.pp.neighbors)",
            "uwot::optimize_graph_layout on that SAME graph", umap_min_dist, umap_spread, "spectral",
            "igraph::cluster_leiden on that SAME graph (never on UMAP coordinates)", "modularity",
            leiden_resolution, leiden_n_iterations,
            paste0("harmony::RunHarmony on all batch-aware PCs, vars_use = ", batch_col), seed)
)
params <- rbind(params, data.table(
  parameter = c("Leiden_resolution_sweep", "sweep_UMAP_max_cells_per_panel (plot only)",
                "sweep_tree_min_in_prop (plot only)"),
  value     = c(paste(resolution_sweep, collapse = ", "), sweep_umap_max_cells, sweep_tree_min_in_prop)
))
fwrite(params, file.path(dir_main, paste0("parameters_used_", run_tag, ".csv")))

## =============================================================================
## CHUNK 2: load the QC-filtered Seurat object (script 11) + checks
## =============================================================================
stopifnot("input Seurat object not found - run script 11 first" = file.exists(seurat_in))
obj <- readRDS(seurat_in)
obj

metadata <- as.data.table(obj[[]])                       # one row per cell, same order as the object
stopifnot("cell_id column must match the object's cells" = identical(metadata$cell_id, colnames(obj)),
          all(c(batch_col, "x_slide_mm", "y_slide_mm", "fov") %in% names(metadata)))

## counts: genes x cells (as the vignette uses) and cells x genes
counts_t <- methods::as(LayerData(obj, assay = "RNA", layer = "counts"), "CsparseMatrix")
counts   <- Matrix::t(counts_t)
stopifnot(identical(rownames(counts), metadata$cell_id))

## per-cell total counts from the matrix (all genes) = nCount_RNA
totalcounts <- Matrix::rowSums(counts)
if (any(totalcounts == 0)) stop(sum(totalcounts == 0), " cells have zero counts.")
message("Cells where matrix total differs from nCount_RNA: ", sum(totalcounts != metadata$nCount_RNA))

## genes with zero counts in every kept cell carry no information
zero_genes <- colnames(counts)[Matrix::colSums(counts) == 0]
if (length(zero_genes) > 0) message(length(zero_genes), " genes have zero counts in all kept cells (not HVGs).")

sections <- sort(unique(metadata[[batch_col]]))
message("Loaded ", format(nrow(counts), big.mark = ","), " cells x ", ncol(counts), " genes; ",
        length(sections), " sections: ", paste(sections, collapse = ", "))
print(table(metadata[[batch_col]]))

xy <- data.table(x_mm = metadata$x_slide_mm, y_mm = metadata$y_slide_mm,    # same coordinates as all earlier plots
                 section = metadata[[batch_col]])
section_label_pos <- xy[, .(x_mm = mean(x_mm), y_mm = mean(y_mm)), by = section]  # label position only

## =============================================================================
## CHUNK 3: HVGs (shared) - Seurat FindVariableFeatures, as in the vignette
## =============================================================================
hvg_file <- file.path(dir_shared, "hvgs.RDS")
if (file.exists(hvg_file) && !rerun_pca) {
  hvgs <- readRDS(hvg_file)
  message("HVGs reused from ", hvg_file, ": ", length(hvgs))
} else {
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = n_hvgs)   # vst on the counts layer
  hvgs <- Seurat::VariableFeatures(obj)
  stopifnot(all(hvgs %in% rownames(counts_t)))
  message("\nHVGs selected: ", length(hvgs))
  
  ## batch-aware Pearson residuals use a gene frequency PER SECTION. An HVG with
  ## zero counts in a whole section has frequency 0 there (expected count 0, so
  ## the residual is undefined). scPearsonPCA's own input check requires
  ## "gene frequency > 0 for all genes in all batches, otherwise the gene should
  ## be removed" - but that check does not catch these genes (sparse zeros), so
  ## they would pass silently. Only such genes (if any) are removed, and listed.
  sec_factor <- factor(metadata[[batch_col]])
  det_by_sec <- sapply(levels(sec_factor), function(s)
    Matrix::colSums(counts[sec_factor == s, hvgs, drop = FALSE]) > 0)
  hvg_zero_in_a_section <- hvgs[rowSums(!det_by_sec) > 0]
  if (length(hvg_zero_in_a_section) > 0) {
    message("HVGs with ZERO counts in at least one whole section (removed, listed in CSV): ",
            length(hvg_zero_in_a_section))
    fwrite(data.table(gene = hvg_zero_in_a_section),
           file.path(dir_shared, "HVGs_removed_zero_counts_in_a_section.csv"))
    hvgs <- setdiff(hvgs, hvg_zero_in_a_section)
  } else {
    message("All HVGs have counts in every section - none removed.")
  }
  saveRDS(hvgs, hvg_file)
  fwrite(data.table(rank = seq_along(hvgs), gene = hvgs), file.path(dir_shared, "hvgs.csv"))
}

## =============================================================================
## CHUNK 4: batch-aware Pearson-residual PCA (shared) - the vignette's call
## =============================================================================
obs <- metadata[, .(cell_id, sample_id = get(batch_col))]
pca_file <- file.path(dir_shared, "pcaobj_batch_aware.RDS")
if (file.exists(pca_file) && !rerun_pca) {
  pcaobj <- readRDS(pca_file)
  message("Batch-aware PCA reused from ", pca_file)
} else {
  
  message("\nComputing per-section gene frequencies...")
  genefreq_batch <- scPearsonPCA::gene_frequency(
    x = counts_t,
    obs = obs,
    cellid_colname = "cell_id",
    batch_variable = "sample_id"
  )
  stopifnot("gene_frequency output missing HVGs" = all(hvgs %in% rownames(genefreq_batch)))
  gf_hvg <- genefreq_batch[hvgs, , drop = FALSE]
  if (any(!is.finite(as.matrix(gf_hvg))) || any(as.matrix(gf_hvg) <= 0)) {
    stop("Non-positive or non-finite gene frequencies for some HVGs - inspect genefreq_batch.")
  }
  saveRDS(genefreq_batch, file.path(dir_shared, "genefreq_batch.RDS"))
  
  message("Running batch-aware Pearson-residual PCA on ", length(hvgs), " HVGs...")
  set.seed(seed)
  pcaobj <- sparse_quasipoisson_pca_seurat_batch(
    counts_t[hvgs, ],
    totalcounts    = totalcounts,
    grate          = gf_hvg,
    obs            = obs,
    batch_variable = "sample_id",
    cellid_colname = "cell_id",
    scale.max      = 10,
    do.scale       = TRUE,
    do.center      = TRUE
  )
  saveRDS(pcaobj, pca_file)
}

if (is.null(pcaobj$reduction.data)) stop("pcaobj has no $reduction.data - fields: ", paste(names(pcaobj), collapse = ", "))
emb_pca <- pcaobj$reduction.data@cell.embeddings
if (is.null(rownames(emb_pca))) {
  stopifnot(nrow(emb_pca) == nrow(metadata))
  rownames(emb_pca) <- metadata$cell_id
}
emb_pca <- emb_pca[metadata$cell_id, , drop = FALSE]             # guarantee metadata order
stopifnot("PCA embeddings contain NaN/Inf" = all(is.finite(emb_pca)))
if (ncol(emb_pca) < n_pcs) stop("PCA returned ", ncol(emb_pca), " PCs, fewer than n_pcs = ", n_pcs)
message("PCA embeddings: ", nrow(emb_pca), " cells x ", ncol(emb_pca), " PCs")

## elbow plot (all PCs computed), with the n_pcs cut-off marked
pc_sd <- pcaobj$reduction.data@stdev
if (length(pc_sd) > 0) {
  elbow <- data.table(PC = seq_along(pc_sd), variance = pc_sd^2)
  elbow[, pct_variance_of_computed_PCs := 100 * variance / sum(variance)]
  fwrite(elbow, file.path(dir_shared, "PCA_variance_per_PC.csv"))
  p_elbow <- ggplot(elbow, aes(PC, pct_variance_of_computed_PCs)) +
    geom_point() + geom_line() +
    geom_vline(xintercept = n_pcs, colour = "red", linetype = "dashed") +
    labs(title = "Batch-aware Pearson PCA - variance per PC",
         subtitle = paste0("Red line = ", n_pcs, " PCs used for the neighbour graph"),
         x = "PC", y = "% of variance (of computed PCs)") +
    theme_minimal(base_size = 13) + theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(dir_shared, paste0("PCA_elbow_", run_tag, ".tiff")), p_elbow, width = 8, height = 5.5,
         dpi = dpi_out, device = "tiff", compression = "lzw", bg = "white")
}

## =============================================================================
## CHUNK 5: helper functions - neighbours, graph, UMAP, Leiden, outputs
## =============================================================================

## STEP A - kNN on PC space, cosine distance, the cell itself in column 1
## (sc.pp.neighbors' n_neighbors also counts the cell itself).
knn_cosine_annoy <- function(X, k, n_trees, seed) {
  n <- nrow(X)
  search_k <- 2L * k * n_trees                           # uwot's default search depth
  ann <- new(AnnoyAngular, ncol(X))
  ann$setSeed(seed)
  for (i in seq_len(n)) ann$addItem(i - 1L, X[i, ])
  ann$build(n_trees)
  idx  <- matrix(NA_integer_, n, k)
  dist <- matrix(NA_real_, n, k)
  for (i in seq_len(n)) {
    res <- ann$getNNsByItemList(i - 1L, k, search_k, TRUE)
    if (length(res$item) < k) stop("kNN search returned fewer than ", k, " neighbours for cell ", i)
    ids <- res$item + 1L
    d   <- res$distance
    self_pos <- match(i, ids)
    if (is.na(self_pos)) {                               # force self into column 1
      ids <- c(i, ids[-k]); d <- c(0, d[-k])
    } else if (self_pos != 1L) {
      ids <- c(i, ids[-self_pos]); d <- c(0, d[-self_pos])
    }
    idx[i, ] <- ids; dist[i, ] <- d
  }
  dist <- dist^2 / 2                                     # Annoy angular -> cosine distance (1 - cos), as uwot does
  dist[, 1] <- 0
  list(idx = idx, dist = dist)
}

## STEP D - Leiden on the graph; labels 0,1,2,... by decreasing size (scanpy convention)
run_leiden <- function(g, resolution) {
  args <- list(graph = g, objective_function = "modularity",
               weights = igraph::E(g)$weight, n_iterations = leiden_n_iterations)
  if ("resolution" %in% names(formals(igraph::cluster_leiden))) {
    args$resolution <- resolution                        # igraph >= 2.0
  } else {
    args$resolution_parameter <- resolution              # igraph 1.x
  }
  set.seed(seed)
  cl <- do.call(igraph::cluster_leiden, args)
  mem <- igraph::membership(cl)
  sizes <- sort(table(mem), decreasing = TRUE)
  relabel <- setNames(seq_along(sizes) - 1L, names(sizes))
  factor(unname(relabel[as.character(mem)]), levels = 0:(length(sizes) - 1L))
}

## VERY DISTINCT CLUSTER COLOURS
## Hand-picked, saturated, mid-lightness colours. No near-white (vanishes on
## the white background at alpha 0.3) and no black (clashes with the black
## cluster labels). Ordered so the first ~15 are the most distinct from each
## other - these go to the largest clusters (cluster 0 = largest).
distinct_palette <- c(
  "#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#42D4F4", "#F032E6",  # red, green, blue, orange, purple, cyan, magenta
  "#9A6324", "#469990", "#BFEF45", "#800000", "#FFE119", "#000075", "#808000",  # brown, teal, lime, maroon, yellow, navy, olive
  "#FF8FB1", "#A9A9A9", "#005300", "#009FFF", "#C77CFF", "#FF6F61", "#00FFBE",  # pink, grey, dark green, azure, violet, coral, aquamarine
  "#D4A017", "#7F7FFF", "#8B008B", "#2F4F4F", "#E9967A", "#556B2F", "#FF1493",  # mustard, periwinkle, dark magenta, slate, salmon, dark olive, deep pink
  "#20B2AA", "#B8860B", "#6A5ACD", "#CD853F", "#4B0082", "#00CED1", "#8FBC8F",  # sea green, goldenrod, slate blue, peru, indigo, dark turquoise, sage
  "#A0522D"                                                                      # sienna
)
stopifnot("duplicate colour in distinct_palette" = !anyDuplicated(toupper(distinct_palette)))

cluster_colours <- function(lv) {
  n <- length(lv)
  if (n <= length(distinct_palette)) {
    cols <- distinct_palette[seq_len(n)]
  } else if (requireNamespace("Polychrome", quietly = TRUE)) {
    ## more clusters than hand-picked colours: keep the hand-picked ones first,
    ## then let Polychrome add colours maximally distant from all of them
    set.seed(seed)
    cols <- unname(Polychrome::createPalette(n, seedcolors = distinct_palette, range = c(25, 80)))
  } else {
    warning(n, " clusters but only ", length(distinct_palette), " hand-picked colours; ",
            "extra colours are less distinct. install.packages('Polychrome') for better ones.")
    k <- n - length(distinct_palette)
    extra <- grDevices::hcl(h = seq(0, 360, length.out = k + 1)[-1], c = 90,
                            l = rep(c(40, 60, 75), length.out = k))   # vary lightness too
    cols <- c(distinct_palette, extra)
  }
  setNames(cols, lv)
}
section_colours <- setNames(scales::hue_pal()(length(sections)), sections)

save_tiff <- function(p, path, w, h) {
  ggsave(path, p, width = w, height = h, dpi = dpi_out, device = "tiff", compression = "lzw", bg = "white")
}

## RESOLUTION SWEEP VISUALS - plotting only, nothing is re-clustered here.
## sweep_cols: named list (leiden_res_<r>) of cluster factors from the sweep
## sweep:      the summary table written to Leiden_resolution_sweep_same_graph.csv
plot_resolution_sweep <- function(sweep, sweep_cols, cells, branch_dir, branch_name) {
  res_ord <- sort(resolution_sweep)
  get_cl  <- function(r) sweep_cols[[paste0("leiden_res_", r)]]
  panel_lab <- function(r) paste0("res ", r, "  (", nlevels(get_cl(r)), " clusters)",
                                  ifelse(r == leiden_resolution, "  [main]", ""))
  panel_levels <- vapply(res_ord, panel_lab, character(1))
  
  ## ---- 1. summary curves ---------------------------------------------------
  sw_long <- rbind(
    sweep[, .(resolution, metric = "Number of clusters", value = n_clusters)],
    sweep[, .(resolution, metric = "Clusters with < 50 cells", value = n_clusters_lt_50_cells)],
    sweep[, .(resolution, metric = "Smallest cluster (cells)", value = smallest_cluster)],
    sweep[, .(resolution, metric = paste0("ARI vs main resolution (", leiden_resolution, ")"),
              value = ARI_vs_resolution_main)]
  )
  sw_long[, metric := factor(metric, levels = unique(metric))]
  p_sum <- ggplot(sw_long, aes(resolution, value)) +
    geom_vline(xintercept = leiden_resolution, colour = "red", linetype = "dashed") +
    geom_line() + geom_point(size = 2.5) +
    scale_x_log10(breaks = res_ord, labels = res_ord) +
    facet_wrap(~ metric, ncol = 2, scales = "free_y") +
    theme_bw(base_size = 13) + theme(plot.title = element_text(face = "bold")) +
    labs(title = paste0(branch_name, ": Leiden resolution sweep (same graph)"),
         subtitle = "Red dashed line = main resolution; x axis on log scale",
         x = "Leiden resolution", y = NULL)
  save_tiff(p_sum, file.path(branch_dir, "Leiden_resolution_sweep_summary.tiff"), 11, 8)
  
  ## ---- 2. UMAP panels, one per resolution -----------------------------------
  ## (UMAP coordinates are identical in every panel; only the clustering changes.
  ##  Cluster numbers are by size WITHIN each resolution, so "cluster 3" in one
  ##  panel is not the same cells as "cluster 3" in another - use the tree for that.)
  max_k   <- max(sweep$n_clusters)
  all_lv  <- as.character(0:(max_k - 1L))
  cols_sw <- cluster_colours(all_lv)
  set.seed(seed)
  idx <- if (nrow(cells) > sweep_umap_max_cells) sample(nrow(cells), sweep_umap_max_cells) else sample(nrow(cells))
  um <- rbindlist(lapply(res_ord, function(r) {
    data.table(UMAP_1 = cells$UMAP_1[idx], UMAP_2 = cells$UMAP_2[idx],
               cluster = factor(as.character(get_cl(r)[idx]), levels = all_lv),
               panel = panel_lab(r))
  }))
  um[, panel := factor(panel, levels = panel_levels)]
  lab <- rbindlist(lapply(res_ord, function(r) {                  # label positions from ALL cells
    data.table(UMAP_1 = cells$UMAP_1, UMAP_2 = cells$UMAP_2, cluster = as.character(get_cl(r)))[
      , .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = cluster][, panel := panel_lab(r)]
  }))
  lab[, panel := factor(panel, levels = panel_levels)]
  p_um <- ggplot(um, aes(UMAP_1, UMAP_2, colour = cluster)) +
    geom_point(size = 0.03, alpha = 0.3, shape = 16) +
    geom_text(data = lab, aes(label = cluster), colour = "black", size = 2.3, fontface = "bold") +
    scale_colour_manual(values = cols_sw, guide = "none") +
    facet_wrap(~ panel, ncol = 4) + coord_fixed() + theme_bw(base_size = 12) +
    labs(title = paste0(branch_name, ": UMAP coloured by Leiden cluster at each sweep resolution"),
         subtitle = paste0("Same UMAP in every panel; ",
                           if (length(idx) < nrow(cells)) paste0(format(length(idx), big.mark = ","),
                                                                 " randomly sampled cells drawn per panel; ") else "",
                           "cluster numbers are by size within each resolution"))
  save_tiff(p_um, file.path(branch_dir, "Leiden_resolution_sweep_UMAP_panels.tiff"), 20, 11)
  
  ## ---- 3. cluster tree: how clusters split between consecutive resolutions --
  ## nodes = clusters at each resolution (size = n cells); edges = cells shared
  ## between a cluster and a cluster at the next resolution. in_prop = fraction
  ## of the lower-row (child) cluster's cells that came from the upper (parent).
  node_list <- vector("list", length(res_ord)); edge_list <- list(); x_prev <- NULL
  for (i in seq_along(res_ord)) {
    r  <- res_ord[i]; cl <- get_cl(r)
    nd <- as.data.table(table(cluster = as.character(cl)))
    setnames(nd, "N", "n_cells")
    if (i == 1L) {
      nd[, order_key := as.numeric(cluster)]
    } else {
      pr <- get_cl(res_ord[i - 1L])
      ct_i <- as.data.table(table(parent = as.character(pr), cluster = as.character(cl)))[N > 0]
      ct_i[, parent_x := x_prev[parent]]
      nd <- merge(nd, ct_i[, .(order_key = sum(N * parent_x) / sum(N)), by = cluster], by = "cluster")
      ct_i[, in_prop := N / sum(N), by = cluster]
      edge_list[[length(edge_list) + 1L]] <-
        ct_i[, .(res_from = res_ord[i - 1L], parent, res_to = r, cluster, n_cells = N, in_prop = round(in_prop, 4))]
    }
    setorder(nd, order_key)                              # children placed under their parents
    nd[, x := (seq_len(.N) - 0.5) / .N]
    nd[, `:=`(resolution = r, level = i)]
    x_prev <- setNames(nd$x, nd$cluster)
    node_list[[i]] <- nd[, .(resolution, level, cluster, n_cells, x)]
  }
  nodes <- rbindlist(node_list)
  edges <- rbindlist(edge_list)
  fwrite(nodes, file.path(branch_dir, "Leiden_resolution_sweep_cluster_tree_nodes.csv"))
  fwrite(edges, file.path(branch_dir, "Leiden_resolution_sweep_cluster_tree_edges.csv"))
  
  e_plot <- merge(edges[in_prop >= sweep_tree_min_in_prop],
                  nodes[, .(res_from = resolution, parent = cluster, x_from = x, y_from = -level)],
                  by = c("res_from", "parent"))
  e_plot <- merge(e_plot, nodes[, .(res_to = resolution, cluster, x_to = x, y_to = -level)],
                  by = c("res_to", "cluster"))
  main_level <- match(leiden_resolution, res_ord)
  
  p_tree <- ggplot() +
    { if (!is.na(main_level)) geom_hline(yintercept = -main_level, colour = "red",
                                         linetype = "dashed", alpha = 0.5) } +
    geom_segment(data = e_plot, aes(x = x_from, y = y_from, xend = x_to, yend = y_to,
                                    alpha = in_prop, linewidth = in_prop), colour = "grey25") +
    geom_point(data = nodes, aes(x, -level, size = n_cells, fill = factor(resolution)),
               shape = 21, colour = "black") +
    geom_text(data = nodes, aes(x, -level, label = cluster), size = 2.2) +
    scale_y_continuous(breaks = -seq_along(res_ord), labels = res_ord,
                       expand = expansion(add = 0.5)) +
    scale_x_continuous(expand = expansion(add = 0.02)) +
    scale_size_area(max_size = 12, labels = scales::comma) +
    scale_linewidth(range = c(0.2, 2.5), limits = c(0, 1)) +
    scale_alpha(range = c(0.15, 1), limits = c(0, 1)) +
    guides(fill = "none") +
    theme_minimal(base_size = 13) +
    theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(), plot.title = element_text(face = "bold")) +
    labs(title = paste0(branch_name, ": Leiden cluster tree across resolutions"),
         subtitle = paste0("Each row = one resolution (red dashed = main). Edges: fraction of the lower cluster's cells ",
                           "from the upper cluster (drawn if >= ", sweep_tree_min_in_prop, ").\n",
                           "One clean line in = cluster just splits off; several thick lines in = cells re-shuffled (unstable)."),
         x = NULL, y = "Leiden resolution", size = "Cells", linewidth = "in_prop", alpha = "in_prop")
  save_tiff(p_tree, file.path(branch_dir, "Leiden_resolution_sweep_cluster_tree.tiff"), 18, 11)
  
  message("  resolution sweep plots saved (summary, UMAP panels, cluster tree)")
  invisible(list(nodes = nodes, edges = edges))
}

run_branch <- function(emb_full, branch_dir, branch_name, emb_label) {
  message("\n", strrep("=", 78), "\n", branch_name, "\n", strrep("=", 78))
  stopifnot(identical(rownames(emb_full), metadata$cell_id))
  X <- emb_full[, seq_len(n_pcs), drop = FALSE]
  sec <- metadata[[batch_col]]
  
  ## PC1 vs PC2 of the input embedding, by section
  df_pc <- data.table(PC1 = emb_full[, 1], PC2 = emb_full[, 2], section = sec)
  set.seed(seed); df_pc <- df_pc[sample(.N)]
  p_pc <- ggplot(df_pc, aes(PC1, PC2, colour = section)) +
    geom_point(size = 0.05, alpha = 0.3, shape = 16) +
    scale_colour_manual(values = section_colours) + coord_fixed() + theme_bw(base_size = 13) +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    labs(title = paste0(emb_label, ": dimension 1 vs 2, coloured by section"), colour = "Section")
  save_tiff(p_pc, file.path(branch_dir, "PC1_vs_PC2_by_section.tiff"), 9, 8)
  
  ## ---- STEP A: kNN on PC space ------------------------------------------
  message("STEP A: kNN on ", n_pcs, " dims of ", emb_label, " (cosine, k = ", n_neighbors, " incl. self)...")
  nn <- knn_cosine_annoy(X, n_neighbors, annoy_n_trees, seed)
  stopifnot(all(nn$idx[, 1] == seq_len(nrow(X))))
  saveRDS(nn, file.path(branch_dir, "STEP_A_kNN_idx_dist.RDS"))
  
  ## ---- STEP B: ONE neighbour graph from the kNN -------------------------
  message("STEP B: building the neighbour graph (UMAP fuzzy connectivities) from the kNN...")
  G <- uwot::similarity_graph(X = NULL, nn_method = nn, set_op_mix_ratio = 1,
                              local_connectivity = 1, verbose = TRUE)
  G <- methods::as(G, "CsparseMatrix")
  stopifnot(nrow(G) == nrow(X), ncol(G) == nrow(X))
  stopifnot("neighbour graph is not symmetric" = Matrix::isSymmetric(G, tol = 1e-6))
  dimnames(G) <- list(metadata$cell_id, metadata$cell_id)
  saveRDS(G, file.path(branch_dir, "STEP_B_neighbour_graph_connectivities.RDS"))
  message("  graph: ", format(nrow(G), big.mark = ","), " cells, ",
          format(Matrix::nnzero(G) / 2, big.mark = ","), " edges")
  
  ## ---- STEP C: UMAP = layout of graph G ---------------------------------
  message("STEP C: UMAP layout of the neighbour graph...")
  graph_for_umap <- G
  set.seed(seed)
  umap_xy <- uwot::optimize_graph_layout(graph = graph_for_umap, X = NULL, n_components = 2,
                                         init = "spectral", spread = umap_spread,
                                         min_dist = umap_min_dist, verbose = TRUE)
  colnames(umap_xy) <- c("UMAP_1", "UMAP_2")
  rownames(umap_xy) <- metadata$cell_id
  
  ## ---- STEP D: Leiden on graph G (NOT on UMAP) --------------------------
  message("STEP D: Leiden clustering of the neighbour graph (resolution ", leiden_resolution, ")...")
  graph_for_leiden <- G
  stopifnot("UMAP and Leiden must use the same graph" = identical(graph_for_umap, graph_for_leiden))
  g <- igraph::graph_from_adjacency_matrix(graph_for_leiden, mode = "undirected", weighted = TRUE, diag = FALSE)
  message("  connected components in graph: ", igraph::components(g)$no)
  clusters <- run_leiden(g, leiden_resolution)
  message("  clusters at resolution ", leiden_resolution, ": ", nlevels(clusters))
  
  ## resolution sweep on the SAME graph
  sweep_cols <- list()
  sweep <- rbindlist(lapply(resolution_sweep, function(r) {
    cl_r <- if (r == leiden_resolution) clusters else run_leiden(g, r)
    sweep_cols[[paste0("leiden_res_", r)]] <<- cl_r
    tab <- table(cl_r)
    data.table(resolution = r, n_clusters = length(tab), smallest_cluster = min(tab),
               largest_cluster = max(tab), n_clusters_lt_50_cells = sum(tab < 50),
               ARI_vs_resolution_main = round(igraph::compare(as.integer(cl_r), as.integer(clusters),
                                                              method = "adjusted.rand"), 4))
  }))
  print(sweep)
  fwrite(sweep, file.path(branch_dir, "Leiden_resolution_sweep_same_graph.csv"))
  
  ## ---- tables --------------------------------------------------------------
  cells <- data.table(cell_id = metadata$cell_id, section_id = sec,
                      fov = metadata$fov, x_mm = xy$x_mm, y_mm = xy$y_mm,
                      totalcounts = totalcounts, umap_xy,
                      leiden = clusters)
  for (nm in names(sweep_cols)) cells[, (nm) := sweep_cols[[nm]]]
  fwrite(cells, file.path(branch_dir, "cells_UMAP_and_Leiden_clusters.csv"))
  
  sizes <- cells[, .(n_cells = .N), by = leiden][order(leiden)]
  sizes[, pct_of_cells := round(100 * n_cells / sum(n_cells), 2)]
  fwrite(sizes, file.path(branch_dir, "Leiden_cluster_sizes.csv"))
  
  comp <- cells[, .(n_cells = .N), by = .(leiden, section_id)]
  comp[, pct_of_cluster := round(100 * n_cells / sum(n_cells), 2), by = leiden]
  comp_wide <- dcast(comp, leiden ~ section_id, value.var = "n_cells", fill = 0)
  dom <- comp[order(-n_cells), .SD[1], by = leiden][, .(leiden, dominant_section = section_id,
                                                        pct_in_dominant_section = pct_of_cluster)]
  comp_wide <- merge(comp_wide, dom, by = "leiden")
  comp_wide[, flag_for_investigation := pct_in_dominant_section > dominant_section_flag_pct]
  fwrite(comp_wide[order(leiden)], file.path(branch_dir, "Leiden_cluster_composition_by_section.csv"))
  
  ## section mixing on the kNN graph: fraction of each cell's neighbours
  ## (self excluded) from its own section vs. random expectation (N_s - 1)/(N - 1)
  sec_int <- as.integer(factor(sec))
  nbr <- nn$idx[, -1, drop = FALSE]
  same_frac <- rowMeans(matrix(sec_int[nbr], nrow = nrow(nbr)) == sec_int)
  N  <- length(sec); Ns <- table(sec)
  mix <- data.table(section_id = sec, same_frac = same_frac,
                    expected = (as.numeric(Ns[sec]) - 1) / (N - 1))
  mix <- mix[, .(n_cells = .N,
                 mean_same_section_neighbour_fraction = round(mean(same_frac), 4),
                 expected_if_randomly_mixed = round(expected[1], 4),
                 ratio_observed_over_expected = round(mean(same_frac) / expected[1], 2)),
             by = section_id][order(section_id)]
  print(mix)
  fwrite(mix, file.path(branch_dir, "section_mixing_on_kNN_graph.csv"))
  
  ## ---- plots ---------------------------------------------------------------
  cl_cols <- cluster_colours(levels(clusters))
  set.seed(seed); plot_df <- cells[sample(.N)]
  lab_pos <- cells[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = leiden]
  
  p_umap_cl <- ggplot(plot_df, aes(UMAP_1, UMAP_2, colour = leiden)) +
    geom_point(size = 0.05, alpha = 0.3, shape = 16) +
    geom_label(data = lab_pos, aes(label = leiden), colour = "black", size = 3.5, alpha = 0.8) +
    scale_colour_manual(values = cl_cols) + coord_fixed() + theme_bw(base_size = 13) +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
    labs(title = paste0(branch_name, ": UMAP coloured by Leiden cluster"),
         subtitle = paste0("Leiden on the neighbour graph (res ", leiden_resolution, ", ",
                           nlevels(clusters), " clusters); UMAP = layout of the same graph"),
         colour = "Leiden")
  save_tiff(p_umap_cl, file.path(branch_dir, "UMAP_by_Leiden_cluster.tiff"), 11, 9)
  
  p_umap_sec <- ggplot(plot_df, aes(UMAP_1, UMAP_2, colour = section_id)) +
    geom_point(size = 0.05, alpha = 0.3, shape = 16) +
    scale_colour_manual(values = section_colours) + coord_fixed() + theme_bw(base_size = 13) +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    labs(title = paste0(branch_name, ": UMAP coloured by tissue section"), colour = "Section")
  save_tiff(p_umap_sec, file.path(branch_dir, "UMAP_by_section.tiff"), 11, 9)
  
  p_umap_sec_facet <- ggplot(plot_df, aes(UMAP_1, UMAP_2, colour = leiden)) +
    geom_point(size = 0.05, alpha = 0.3, shape = 16) +
    scale_colour_manual(values = cl_cols, guide = "none") + coord_fixed() + theme_bw(base_size = 12) +
    facet_wrap(~ section_id, ncol = 3) +
    labs(title = paste0(branch_name, ": UMAP split by section, coloured by Leiden cluster"))
  save_tiff(p_umap_sec_facet, file.path(branch_dir, "UMAP_by_Leiden_split_by_section.tiff"), 15, 10)
  
  p_umap_counts <- ggplot(plot_df, aes(UMAP_1, UMAP_2, colour = log10(totalcounts))) +
    geom_point(size = 0.05, alpha = 0.3, shape = 16) +
    scale_colour_viridis_c(option = "B") + coord_fixed() + theme_bw(base_size = 13) +
    labs(title = paste0(branch_name, ": UMAP coloured by log10 total counts"),
         subtitle = "Check that clusters are not driven by total counts", colour = "log10\ncounts")
  save_tiff(p_umap_counts, file.path(branch_dir, "UMAP_by_log10_total_counts.tiff"), 11, 9)
  
  p_xy <- ggplot(plot_df, aes(x_mm, y_mm, colour = leiden)) +
    geom_point(size = 0.05, alpha = 0.6, shape = 16) +
    geom_text(data = section_label_pos, aes(x_mm, y_mm, label = section), inherit.aes = FALSE,
              colour = "black", fontface = "bold", size = 5) +
    scale_colour_manual(values = cl_cols) + coord_fixed() + theme_bw(base_size = 13) +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
    labs(title = paste0(branch_name, ": Leiden clusters in tissue space"),
         x = "x (mm)", y = "y (mm)", colour = "Leiden")
  save_tiff(p_xy, file.path(branch_dir, "Spatial_Leiden_clusters_all_sections.tiff"), 13, 10)
  
  comp_plot <- copy(comp)
  p_comp <- ggplot(comp_plot, aes(x = leiden, y = pct_of_cluster, fill = section_id)) +
    geom_col(width = 0.8) +
    geom_hline(yintercept = dominant_section_flag_pct, linetype = "dashed") +
    scale_fill_manual(values = section_colours) + theme_minimal(base_size = 13) +
    labs(title = paste0(branch_name, ": section composition of each Leiden cluster"),
         subtitle = paste0("Dashed line = ", dominant_section_flag_pct,
                           "%; above it is a flag to investigate (section-specific biology OR residual batch)"),
         x = "Leiden cluster", y = "% of cluster", fill = "Section")
  save_tiff(p_comp, file.path(branch_dir, "Leiden_cluster_composition_by_section.tiff"), 12, 6)
  
  ## ---- resolution sweep visuals (plots only; clusters already computed above) --
  sweep_tree <- plot_resolution_sweep(sweep, sweep_cols, cells, branch_dir, branch_name)
  
  results <- list(branch = branch_name, embedding_used = X, knn = nn, graph = G,
                  umap = umap_xy, leiden = clusters, resolution_sweep = sweep,
                  resolution_sweep_tree = sweep_tree,
                  section_mixing = mix, cells = cells)
  saveRDS(results, file.path(branch_dir, "RESULTS_all_objects.RDS"))
  message(branch_name, " done -> ", branch_dir)
  results
}

## =============================================================================
## CHUNK 6: PRIMARY - batch-aware PCA -> neighbour graph -> UMAP / Leiden
## =============================================================================
res_primary <- run_branch(emb_pca, dir_primary,
                          "PRIMARY (batch-aware Pearson PCA)", "Batch-aware Pearson PCA")

## =============================================================================
## CHUNK 7: SENSITIVITY - batch-aware PCA -> Harmony -> neighbour graph -> UMAP / Leiden
## =============================================================================
message("\nRunning Harmony on all ", ncol(emb_pca), " batch-aware PCs (vars_use = ", batch_col, ")...")
set.seed(seed)
emb_harmony <- harmony::RunHarmony(
  emb_pca,                     # positional, as in the vignette (S3 generic dispatches on the 1st argument)
  meta_data = data.frame(sample_id = metadata[[batch_col]]),
  vars_use  = "sample_id"
)
emb_harmony <- as.matrix(emb_harmony)
stopifnot(nrow(emb_harmony) == nrow(emb_pca), ncol(emb_harmony) == ncol(emb_pca))
if (is.null(rownames(emb_harmony))) rownames(emb_harmony) <- rownames(emb_pca)
stopifnot(identical(rownames(emb_harmony), metadata$cell_id),
          "Harmony embeddings contain NaN/Inf" = all(is.finite(emb_harmony)))
saveRDS(emb_harmony, file.path(dir_sens, "harmony_embeddings.RDS"))

res_sens <- run_branch(emb_harmony, dir_sens,
                       "SENSITIVITY (batch-aware Pearson PCA + Harmony)", "Harmony (on batch-aware PCs)")

## =============================================================================
## CHUNK 8: PRIMARY vs SENSITIVITY - how do the clusterings correspond?
## =============================================================================
ari <- igraph::compare(as.integer(res_primary$leiden), as.integer(res_sens$leiden), method = "adjusted.rand")
message("\nAdjusted Rand Index, PRIMARY vs SENSITIVITY Leiden (res ", leiden_resolution, "): ", round(ari, 4))

ct <- as.data.table(table(primary = res_primary$leiden, sensitivity = res_sens$leiden))
setnames(ct, "N", "n_cells")
ct[, primary     := factor(primary,     levels = levels(res_primary$leiden))]   # numeric order, not 0,1,10,2...
ct[, sensitivity := factor(sensitivity, levels = levels(res_sens$leiden))]
ct[, pct_of_primary_cluster := round(100 * n_cells / sum(n_cells), 2), by = primary]
fwrite(dcast(ct, primary ~ sensitivity, value.var = "n_cells", fill = 0),
       file.path(dir_compare, "crosstab_cells_PRIMARY_rows_x_SENSITIVITY_cols.csv"))
fwrite(dcast(ct, primary ~ sensitivity, value.var = "pct_of_primary_cluster", fill = 0),
       file.path(dir_compare, "crosstab_rowpct_PRIMARY_rows_x_SENSITIVITY_cols.csv"))
fwrite(data.table(comparison = "PRIMARY vs SENSITIVITY Leiden", resolution = leiden_resolution,
                  n_clusters_primary = nlevels(res_primary$leiden),
                  n_clusters_sensitivity = nlevels(res_sens$leiden),
                  adjusted_rand_index = round(ari, 4)),
       file.path(dir_compare, "ARI_PRIMARY_vs_SENSITIVITY.csv"))

p_ct <- ggplot(ct, aes(x = sensitivity, y = primary, fill = pct_of_primary_cluster)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(pct_of_primary_cluster >= 5, sprintf("%.0f", pct_of_primary_cluster), "")), size = 2.8) +
  scale_fill_viridis_c() + theme_minimal(base_size = 12) +
  labs(title = "Leiden clusters: PRIMARY vs SENSITIVITY (Harmony)",
       subtitle = paste0("Each row = 100% of a PRIMARY cluster; ARI = ", round(ari, 3)),
       x = "SENSITIVITY cluster", y = "PRIMARY cluster", fill = "% of\nPRIMARY\ncluster")
save_tiff(p_ct, file.path(dir_compare, "crosstab_heatmap_PRIMARY_vs_SENSITIVITY.tiff"), 12, 10)

mix_both <- merge(res_primary$section_mixing[, .(section_id, primary_same_section_frac = mean_same_section_neighbour_fraction,
                                                 expected_if_randomly_mixed)],
                  res_sens$section_mixing[, .(section_id, sensitivity_same_section_frac = mean_same_section_neighbour_fraction)],
                  by = "section_id")
fwrite(mix_both, file.path(dir_compare, "section_mixing_PRIMARY_vs_SENSITIVITY.csv"))
print(mix_both)

## =============================================================================
## CHUNK 9: add results to a COPY of the Seurat object and save under a new name
## =============================================================================
## Reductions: Pearson PCA, Harmony, and the two UMAPs; metadata: Leiden clusters
## (main resolution + every sweep resolution) for both analyses.
obj_out <- obj
with_key <- function(m, key) { m <- as.matrix(m); colnames(m) <- paste0(key, seq_len(ncol(m))); m }  # column names = key_1, key_2, ...
obj_out[["pearson_pca"]]  <- CreateDimReducObject(embeddings = with_key(emb_pca, "PPC_"),       key = "PPC_",   assay = "RNA")
obj_out[["harmony"]]      <- CreateDimReducObject(embeddings = with_key(emb_harmony, "HARM_"),  key = "HARM_",  assay = "RNA")
obj_out[["umap_primary"]] <- CreateDimReducObject(embeddings = with_key(res_primary$umap, "UMAPP_"), key = "UMAPP_", assay = "RNA")
obj_out[["umap_harmony"]] <- CreateDimReducObject(embeddings = with_key(res_sens$umap, "UMAPH_"),    key = "UMAPH_", assay = "RNA")
for (nm in grep("^leiden", names(res_primary$cells), value = TRUE)) {
  obj_out[[paste0("primary_", run_tag, "_", nm)]] <- res_primary$cells[[nm]]
  obj_out[[paste0("harmony_", run_tag, "_", nm)]] <- res_sens$cells[[nm]]
}
stopifnot(identical(colnames(obj_out), metadata$cell_id))
seurat_out <- file.path(dir_main, paste0("seurat_", input_tag, "_PCA_UMAP_Leiden_", run_tag, ".RDS"))
saveRDS(obj_out, seurat_out)
message("Saved Seurat object with results: ", seurat_out)

## =============================================================================
## CHUNK 10: summary
## =============================================================================
writeLines(capture.output(sessionInfo()), file.path(dir_main, "sessionInfo.txt"))

message("\n", strrep("=", 78))
message("DONE. Everything is in: ", dir_main)
message("  00_shared  : hvgs.csv, genefreq_batch.RDS, pcaobj_batch_aware.RDS, PCA_elbow.tiff")
message("  01_PRIMARY : kNN (STEP A) -> graph (STEP B) -> UMAP of graph (STEP C) + Leiden on graph (STEP D)")
message("               ", nlevels(res_primary$leiden), " clusters at resolution ", leiden_resolution)
message("  02_SENSITIVITY (Harmony): same steps; ", nlevels(res_sens$leiden), " clusters")
message("  01_ and 02_ each also hold the resolution sweep: Leiden_resolution_sweep_same_graph.csv +")
message("               _summary.tiff, _UMAP_panels.tiff, _cluster_tree.tiff (+ tree nodes/edges CSVs)")
message("  03_comparison: ARI = ", round(ari, 4), ", crosstab + section mixing")
message("  parameters_used_", run_tag, ".csv + sessionInfo.txt for the methods section")
message(strrep("=", 78))




















## =============================================================================
## _UMAP_Leiden_nCount_RNA_and_tissue.R UMAP side by side UMAP tissue, and UMAP nCount_RNA 
## =============================================================================
## =============================================================================
## 13_UMAP_Leiden_nCount_RNA_and_tissue.R
##   UMAP side by side: Leiden vs nCount_RNA, and Leiden vs tissue section
## =============================================================================
## PURPOSE
##   Re-plot the UMAPs from script 12 at 20 PCs, 250 neighbours, Leiden
##   resolution 1, for BOTH analyses:
##     01_PRIMARY_batchawarePCA_PC20_NN250              (batch-aware Pearson PCA)
##     02_SENSITIVITY_batchawarePCA_Harmony_PC20_NN250  (batch-aware Pearson PCA + Harmony)
##   coloured three ways:
##     1. Leiden cluster (resolution 1)
##     2. total RNA counts per cell (nCount_RNA)
##     3. tissue section (the 6 AAA sections)
##   Inside EACH analysis folder:
##     UMAP_leiden_nCount_RNA/             : 1, 2, and 1 + 2 side by side
##     UMAP_leiden_color_cells_by_tissue/  : 1, 3, and 1 + 3 side by side
##
##   Nothing is recomputed. UMAP coordinates and clusters are read from
##   script 12's output. The UMAP depends only on PCs + neighbours; the
##   resolution changes only which Leiden cluster each cell is coloured as.
##   Leiden clusters of PRIMARY and SENSITIVITY are separate clusterings:
##   cluster "3" in one is not the same cells as cluster "3" in the other.
##
## REQUIRES: script 12 run with  n_pcs <- 20  and  n_neighbors <- 250.
##   Resolution 1 is in script 12's resolution_sweep, so "leiden_res_1"
##   exists in both cells tables.
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 1: packages, settings, paths
## =============================================================================
library(data.table)   # read the cells tables
library(ggplot2)      # plots
library(patchwork)    # put two plots side by side (also used in script 07)
library(scales)       # comma labels on the counts legend

output_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"
input_tag  <- "seurat_QC_filtered_per_section_03_to_09"   # same input as script 12

n_pcs       <- 20     # PCs used for the neighbour graph
n_neighbors <- 250    # neighbours (incl. the cell itself)
leiden_res  <- 1      # Leiden resolution to colour by

run_tag  <- paste0("PC", n_pcs, "_NN", n_neighbors)       # "PC20_NN250" (script 12 naming)
plot_tag <- paste0(run_tag, "_res", leiden_res)          # "PC20_NN250_res1" (file names)
dir_main <- file.path(output_dir, "12_Normalisation_PCA_UMAP_Leiden", input_tag)

## The two script 12 analysis folders; name = label used in plot subtitles
branches <- c(
  "Batch-aware Pearson PCA (PRIMARY)"               = file.path(dir_main, paste0("01_PRIMARY_batchawarePCA_", run_tag)),
  "Batch-aware Pearson PCA + Harmony (SENSITIVITY)" = file.path(dir_main, paste0("02_SENSITIVITY_batchawarePCA_Harmony_", run_tag)))
cells_files <- file.path(branches, "cells_UMAP_and_Leiden_clusters.csv")   # UMAP + clusters per cell
stopifnot("Run script 12 with n_pcs = 20 and n_neighbors = 250 first" = all(file.exists(cells_files)))

seed    <- 42    # same seed as script 12 (random plotting order)
dpi_out <- 300   # same resolution as script 12


## =============================================================================
## CHUNK 2: nCount_RNA per cell (shared by both analyses)
## =============================================================================
metadata <- readRDS(file.path(output_dir, "metadata_unfiltered.RDS"))
ncount   <- setNames(metadata$nCount_RNA, metadata$cell_id)   # named by cell, for matching
rm(metadata); invisible(gc())                                 # rest of metadata not needed


## =============================================================================
## CHUNK 3: colours (shared by both analyses)
## =============================================================================
## Leiden: same hand-picked palette and order as script 12, so cluster colours match
distinct_palette <- c(
  "#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#42D4F4", "#F032E6",
  "#9A6324", "#469990", "#BFEF45", "#800000", "#FFE119", "#000075", "#808000",
  "#FF8FB1", "#A9A9A9", "#005300", "#009FFF", "#C77CFF", "#FF6F61", "#00FFBE",
  "#D4A017", "#7F7FFF", "#8B008B", "#2F4F4F", "#E9967A", "#556B2F", "#FF1493",
  "#20B2AA", "#B8860B", "#6A5ACD", "#CD853F", "#4B0082", "#00CED1", "#8FBC8F",
  "#A0522D")
cluster_colours <- function(cl_lv) {
  if (length(cl_lv) <= length(distinct_palette)) {
    return(setNames(distinct_palette[seq_along(cl_lv)], cl_lv))
  }
  ## more clusters than colours: as script 12
  stopifnot("install.packages('Polychrome') for > 36 clusters" = requireNamespace("Polychrome", quietly = TRUE))
  set.seed(seed)
  setNames(unname(Polychrome::createPalette(length(cl_lv), seedcolors = distinct_palette,
                                            range = c(25, 80))), cl_lv)
}

## Tissue sections: 6 colours easy to tell apart by eye (same in both analyses)
tissue_cols <- c(AAA2  = "#2CA02C",   # green
                 AAA38 = "#1F4EDB",   # blue
                 AAA40 = "#E31A1C",   # red
                 AAA45 = "#FF8C00",   # orange
                 AAA46 = "#8E44AD",   # purple
                 AAA48 = "#E6B800")   # gold

## Save as TIFF, as in script 12
save_tiff <- function(p, path, w, h) {
  ggsave(path, p, width = w, height = h, dpi = dpi_out, device = "tiff", compression = "lzw", bg = "white")
  message("Saved: ", path)
}


## =============================================================================
## CHUNK 4: function - all plots for ONE analysis folder
## =============================================================================
## branch_dir   = script 12 analysis folder (PRIMARY or SENSITIVITY)
## branch_label = text shown in the plot subtitles
make_umap_plots <- function(branch_dir, branch_label) {
  message("\n==== ", branch_label, " ====")
  
  ## ---- load UMAP, Leiden (res 1) and nCount_RNA per cell ----
  cells <- fread(file.path(branch_dir, "cells_UMAP_and_Leiden_clusters.csv"))
  leiden_col <- paste0("leiden_res_", leiden_res)                      # "leiden_res_1"
  stopifnot("resolution not in script 12's sweep" = leiden_col %in% names(cells))
  cells[, cluster := factor(get(leiden_col), levels = sort(unique(get(leiden_col))))]  # 0,1,2,... largest first
  message("Leiden clusters at resolution ", leiden_res, ": ", nlevels(cells$cluster))
  
  cells[, nCount_RNA := unname(ncount[cell_id])]                       # match by cell name
  stopifnot("cells missing from metadata" = !anyNA(cells$nCount_RNA),
            "nCount_RNA differs from script 12 totals" = all(cells$nCount_RNA == cells$totalcounts))
  stopifnot("section names differ from tissue_cols" = setequal(names(tissue_cols), unique(cells$section_id)))
  
  ## ---- output folders inside this analysis folder ----
  dir_counts <- file.path(branch_dir, "UMAP_leiden_nCount_RNA")
  dir_tissue <- file.path(branch_dir, "UMAP_leiden_color_cells_by_tissue")
  dir.create(dir_counts, showWarnings = FALSE, recursive = TRUE)
  dir.create(dir_tissue, showWarnings = FALSE, recursive = TRUE)
  
  ## ---- colours and legend labels ----
  cl_cols <- cluster_colours(levels(cells$cluster))
  n_sec <- table(cells$section_id)[names(tissue_cols)]                 # cells per section
  tissue_labels <- setNames(paste0(names(tissue_cols), " (", format(as.integer(n_sec), big.mark = ","), " cells)"),
                            names(tissue_cols))
  
  ## ---- plot data ----
  set.seed(seed)
  plot_df <- cells[sample(.N)]                                         # random order: nothing hidden underneath
  lab_pos <- cells[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = cluster]   # label at cluster centre
  sub_txt <- paste0(branch_label, ", ", n_pcs, " PCs, ", n_neighbors, " neighbours")
  
  ## 1. Leiden clusters (same style as script 12)
  p_leiden <- ggplot(plot_df, aes(UMAP_1, UMAP_2, colour = cluster)) +
    geom_point(size = 0.05, alpha = 0.3, shape = 16) +
    geom_label(data = lab_pos, aes(label = cluster), colour = "black", size = 3.5, alpha = 0.8) +
    scale_colour_manual(values = cl_cols) +
    coord_fixed() + theme_bw(base_size = 13) +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
    labs(title = paste0("UMAP - Leiden clusters (resolution ", leiden_res, ", ", nlevels(cells$cluster), " clusters)"),
         subtitle = sub_txt, colour = "Leiden\ncluster")
  
  ## 2. Total RNA counts per cell; colour on log10 scale, legend shows real counts
  p_counts <- ggplot(plot_df, aes(UMAP_1, UMAP_2, colour = nCount_RNA)) +
    geom_point(size = 0.05, alpha = 0.3, shape = 16) +
    scale_colour_viridis_c(option = "B", trans = "log10", labels = label_comma()) +   # black = few, yellow = many
    coord_fixed() + theme_bw(base_size = 13) +
    guides(colour = guide_colourbar(barheight = unit(8, "cm"))) +                     # tall, readable legend
    labs(title = "UMAP - total RNA counts per cell (nCount_RNA)",
         subtitle = paste0(sub_txt, "; median ", format(median(cells$nCount_RNA), big.mark = ","), " counts per cell"),
         colour = "Total RNA counts\nper cell\n(nCount_RNA,\nlog10 colour scale)")
  
  ## 3. Tissue section of each cell
  p_tissue <- ggplot(plot_df, aes(UMAP_1, UMAP_2, colour = section_id)) +
    geom_point(size = 0.05, alpha = 0.3, shape = 16) +
    scale_colour_manual(values = tissue_cols, labels = tissue_labels) +
    coord_fixed() + theme_bw(base_size = 13) +
    guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    labs(title = "UMAP - cells coloured by tissue section", subtitle = sub_txt, colour = "Tissue section")
  
  ## ---- folder UMAP_leiden_nCount_RNA ----
  save_tiff(p_leiden, file.path(dir_counts, paste0("UMAP_Leiden_", plot_tag, ".tiff")), 11, 9)
  save_tiff(p_counts, file.path(dir_counts, paste0("UMAP_nCount_RNA_", plot_tag, ".tiff")), 11, 9)
  p_leiden_counts <- (p_leiden | p_counts) +                           # side by side
    plot_annotation(title = paste0("Leiden clusters vs total RNA counts - ", branch_label, " (", plot_tag, ")"),
                    theme = theme(plot.title = element_text(face = "bold", size = 16)))
  save_tiff(p_leiden_counts, file.path(dir_counts, paste0("UMAP_Leiden_and_nCount_RNA_side_by_side_", plot_tag, ".tiff")), 22, 10)
  
  ## ---- folder UMAP_leiden_color_cells_by_tissue ----
  save_tiff(p_leiden, file.path(dir_tissue, paste0("UMAP_Leiden_", plot_tag, ".tiff")), 11, 9)
  save_tiff(p_tissue, file.path(dir_tissue, paste0("UMAP_by_tissue_section_", plot_tag, ".tiff")), 11, 9)
  p_leiden_tissue <- (p_leiden | p_tissue) +                           # side by side
    plot_annotation(title = paste0("Leiden clusters vs tissue section - ", branch_label, " (", plot_tag, ")"),
                    theme = theme(plot.title = element_text(face = "bold", size = 16)))
  save_tiff(p_leiden_tissue, file.path(dir_tissue, paste0("UMAP_Leiden_and_tissue_section_side_by_side_", plot_tag, ".tiff")), 22, 10)
  
  message("Done: ", dir_counts, "\n      ", dir_tissue)
}


## =============================================================================
## CHUNK 5: run for PRIMARY and SENSITIVITY (Harmony)
## =============================================================================
for (b in names(branches)) make_umap_plots(branches[[b]], b)







.rs.restartR()
library(Seurat)
packageVersion("Seurat")        # expect 5.6.0
packageVersion("SeuratObject")



##install.packages('leidenbase')
##install.packages('glmGamPoi')
##BiocManager::install('glmGamPoi')
library(leidenbase)
library(glmGamPoi)


## =============================================================================
## 12b_Seurat_functions_PCA_UMAP_Leiden.R  (1) NO FINDALL MARKERS
## =============================================================================
## PURPOSE
##   Same analysis as 12_Normalisation_PCA_UMAP_Leiden.R, but the neighbour
##   search, UMAP, Leiden clustering and plots use SEURAT functions:
##     FindNeighbors()  -> kNN on the PCs (cosine, annoy)
##     RunUMAP()        -> UMAP from that kNN (nn.name = ...)
##     FindClusters()   -> Leiden (algorithm = 4) on the neighbour graph from
##                         the SAME kNN (never on UMAP coordinates)
##     ElbowPlot(), DimPlot(), FeaturePlot(), ImageDimPlot() -> plots
##
##   PCA is NOT recomputed. Script 12's batch-aware Pearson-residual PCA and
##   Harmony embeddings are loaded and added to the Seurat object as reductions.
##   (Seurat has no batch-aware Pearson-residual PCA. SCTransform also makes
##    Pearson residuals, but from a different model - regularised negative
##    binomial, not scPearsonPCA's per-section quasi-Poisson - so it is not used.)
##
##   PRIMARY     : pearson_pca          -> FindNeighbors -> RunUMAP / FindClusters
##   SENSITIVITY : pearson_pca -> harmony -> FindNeighbors -> RunUMAP / FindClusters
##   Settings: 20 PCs, 250 neighbours (incl. the cell itself), main resolution 1,
##   plus the resolution sweep on the same graph.
##
## ONE DIFFERENCE FROM SCRIPT 12 (Seurat's design, not a choice made here)
##   Script 12 (scanpy-style) clustered the UMAP fuzzy-connectivity graph.
##   Seurat's FindClusters clusters the shared-nearest-neighbour (SNN) graph that
##   FindNeighbors builds from the kNN. UMAP and Leiden still come from the SAME
##   kNN. Clusters can therefore differ slightly from script 12.
##   NOTE: with k = 250 the SNN graph is large (needs several GB of RAM).
##
## OUTPUT: <script 12 folder>/seurat_functions_used/PC20_NN250/
##   00_shared/ , 01_PRIMARY_.../ , 02_SENSITIVITY_.../ , 03_PRIMARY_vs_SENSITIVITY/
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 0: packages
## =============================================================================
library(data.table)   # tables
library(Matrix)       # sparse matrices (kNN check)
library(ggplot2)      # extra plot layers
library(patchwork)    # combines Seurat plots (the & operator)
library(Seurat)       # FindNeighbors, RunUMAP, FindClusters, DimPlot, ...
## FindClusters(algorithm = 4) needs the 'leidenbase' package (Seurat >= 5.1):
## install.packages("leidenbase") if FindClusters asks for it.


## =============================================================================
## CHUNK 1: paths + parameters (same values as script 12)
## =============================================================================
output_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"
seurat_in  <- file.path(output_dir, "seurat_QC_filtered_per_section_03_to_09.RDS")  # script 11 (read only)
input_tag  <- sub("\\.RDS$", "", basename(seurat_in))
dir_main   <- file.path(output_dir, "12_Normalisation_PCA_UMAP_Leiden", input_tag)  # script 12 folder

batch_col           <- "sample_id"    # the 6 AAA sections
n_pcs               <- 20             # PCs used by FindNeighbors
n_neighbors         <- 250            # k.param (Seurat counts the cell itself)
annoy_n_trees       <- 50             # annoy index size (as script 12)
umap_min_dist       <- 0.02           # as script 12
umap_spread         <- 2              # as script 12
leiden_resolution   <- 1              # main resolution
leiden_n_iterations <- 2              # as script 12
resolution_sweep    <- c(0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0)  # all on the same graph
dominant_section_flag_pct <- 90       # flag clusters with > 90% cells from one section
seed    <- 42
dpi_out <- 300
stopifnot(leiden_resolution %in% resolution_sweep)

run_tag     <- paste0("PC", n_pcs, "_NN", n_neighbors)                     # "PC20_NN250"
dir_seu     <- file.path(dir_main, "seurat_functions_used", run_tag)       # all outputs here
dir_shared  <- file.path(dir_seu, "00_shared")
dir_primary <- file.path(dir_seu, "01_PRIMARY_batchawarePCA")
dir_sens    <- file.path(dir_seu, "02_SENSITIVITY_batchawarePCA_Harmony")
dir_compare <- file.path(dir_seu, "03_PRIMARY_vs_SENSITIVITY")
for (d in c(dir_shared, dir_primary, dir_sens, dir_compare)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## script 12 inputs reused (not recomputed)
pca_file  <- file.path(dir_main, "00_shared_HVGs_and_batch_aware_Pearson_PCA", "pcaobj_batch_aware.RDS")
harm_file <- file.path(dir_main, paste0("02_SENSITIVITY_batchawarePCA_Harmony_", run_tag), "harmony_embeddings.RDS")
stopifnot("Run script 12 first (batch-aware PCA not found)" = file.exists(pca_file))

fwrite(data.table(
  parameter = c("input", "PCA", "Harmony", "n_PCs", "k.param (incl. self)", "nn.method", "annoy.metric",
                "n.trees", "UMAP", "min.dist", "spread", "clustering", "graph clustered", "prune.SNN",
                "Leiden n.iter", "main resolution", "resolution sweep", "seed"),
  value = c(basename(seurat_in), "script 12 scPearsonPCA (batch-aware Pearson residuals), loaded",
            "script 12 RunHarmony on all PCs, loaded", n_pcs, n_neighbors, "annoy", "cosine",
            annoy_n_trees, "Seurat::RunUMAP(nn.name = kNN, umap.method = 'uwot')", umap_min_dist,
            umap_spread, "Seurat::FindClusters(algorithm = 4, Leiden)", "SNN from FindNeighbors (same kNN)",
            "1/15 (Seurat default)", leiden_n_iterations, leiden_resolution,
            paste(resolution_sweep, collapse = ", "), seed)),
  file.path(dir_seu, "parameters_used.csv"))


## =============================================================================
## CHUNK 2: load object, add script 12's PCA and Harmony as Seurat reductions
## =============================================================================
obj <- readRDS(seurat_in)
cells_order <- colnames(obj)                                   # cell order of the object

## batch-aware Pearson PCA (all 50 PCs + their standard deviations)
pcaobj  <- readRDS(pca_file)
emb_pca <- pcaobj$reduction.data@cell.embeddings
if (is.null(rownames(emb_pca))) rownames(emb_pca) <- cells_order  # as in script 12
emb_pca <- emb_pca[cells_order, , drop = FALSE]                  # same order as the object
pc_sd   <- pcaobj$reduction.data@stdev                           # for ElbowPlot
rm(pcaobj); invisible(gc())

## Harmony on all PCs: load script 12's result; if absent, the identical call
if (file.exists(harm_file)) {
  emb_harmony <- readRDS(harm_file)
} else {
  set.seed(seed)
  emb_harmony <- as.matrix(harmony::RunHarmony(emb_pca, meta_data = data.frame(sample_id = obj[[batch_col]][, 1]),
                                               vars_use = "sample_id"))
  rownames(emb_harmony) <- rownames(emb_pca)
}
emb_harmony <- emb_harmony[cells_order, , drop = FALSE]
stopifnot(all(is.finite(emb_pca)), all(is.finite(emb_harmony)))

## Seurat needs column names = key + number (PPC_1, PPC_2, ...)
with_key <- function(m, key) { m <- as.matrix(m); colnames(m) <- paste0(key, seq_len(ncol(m))); m }
obj[["pearson_pca"]] <- CreateDimReducObject(embeddings = with_key(emb_pca, "PPC_"), stdev = pc_sd,
                                             key = "PPC_", assay = "RNA")
obj[["harmony"]]     <- CreateDimReducObject(embeddings = with_key(emb_harmony, "HARM_"),
                                             key = "HARM_", assay = "RNA")
rm(emb_pca, emb_harmony); invisible(gc())
obj   # now lists the two reductions

## Elbow plot (Seurat), n_pcs marked
p_elbow <- ElbowPlot(obj, reduction = "pearson_pca", ndims = length(pc_sd)) +
  geom_vline(xintercept = n_pcs, colour = "red", linetype = "dashed") +
  labs(title = "Batch-aware Pearson PCA", subtitle = paste0("Red line = ", n_pcs, " PCs used"))
ggsave(file.path(dir_shared, "ElbowPlot_pearson_pca.tiff"), p_elbow, width = 8, height = 5.5,
       dpi = dpi_out, compression = "lzw", bg = "white")


## =============================================================================
## CHUNK 3: helpers - colours and saving (same palette as script 12)
## =============================================================================
distinct_palette <- c(
  "#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#42D4F4", "#F032E6",
  "#9A6324", "#469990", "#BFEF45", "#800000", "#FFE119", "#000075", "#808000",
  "#FF8FB1", "#A9A9A9", "#005300", "#009FFF", "#C77CFF", "#FF6F61", "#00FFBE",
  "#D4A017", "#7F7FFF", "#8B008B", "#2F4F4F", "#E9967A", "#556B2F", "#FF1493",
  "#20B2AA", "#B8860B", "#6A5ACD", "#CD853F", "#4B0082", "#00CED1", "#8FBC8F",
  "#A0522D")
cluster_colours <- function(lv) {                 # one colour per cluster level
  if (length(lv) <= length(distinct_palette)) return(setNames(distinct_palette[seq_along(lv)], lv))
  stopifnot("install.packages('Polychrome') for > 36 clusters" = requireNamespace("Polychrome", quietly = TRUE))
  set.seed(seed)
  setNames(unname(Polychrome::createPalette(length(lv), seedcolors = distinct_palette, range = c(25, 80))), lv)
}
save_tiff <- function(p, path, w, h) {
  ggsave(path, p, width = w, height = h, dpi = dpi_out, device = "tiff", compression = "lzw", bg = "white")
  message("Saved: ", path)
}
img_name <- Images(obj)[1]                        # spatial FOV (cell centroids), for ImageDimPlot


## =============================================================================
## CHUNK 4: one analysis branch with Seurat functions
## =============================================================================
## reduction = "pearson_pca" (PRIMARY) or "harmony" (SENSITIVITY)
## prefix    = short name used for graphs, UMAP and cluster columns
run_branch_seurat <- function(obj, reduction, prefix, branch_dir, branch_name) {
  message("\n", strrep("=", 70), "\n", branch_name, "\n", strrep("=", 70))
  nn_name   <- paste0(prefix, "_neighbors")    # Seurat Neighbor object: kNN indices + distances
  nn_graph  <- paste0(prefix, "_nn")           # kNN as a graph
  snn_graph <- paste0(prefix, "_snn")          # shared-nearest-neighbour graph (clustered)
  umap_name <- paste0("umap_", prefix)
  cl_cols   <- paste0(prefix, "_leiden_res_", resolution_sweep)        # one column per resolution
  main_col  <- paste0(prefix, "_leiden_res_", leiden_resolution)
  sec       <- obj[[batch_col]][, 1]
  
  ## ---- PC1 vs PC2 of the input reduction, by section --------------------
  save_tiff(DimPlot(obj, reduction = reduction, dims = 1:2, group.by = batch_col, shuffle = TRUE, seed = seed) +
              labs(title = paste0(branch_name, ": dimension 1 vs 2")),
            file.path(branch_dir, "PC1_vs_PC2_by_section.tiff"), 9, 8)
  
  ## ---- STEP A: kNN (FindNeighbors) ---------------------------------------
  ## (a) kNN kept as a Neighbor object (indices + distances) -> used by RunUMAP
  obj <- FindNeighbors(obj, reduction = reduction, dims = 1:n_pcs, k.param = n_neighbors,
                       nn.method = "annoy", annoy.metric = "cosine", n.trees = annoy_n_trees,
                       return.neighbor = TRUE, graph.name = nn_name)
  ## (b) same search, returned as graphs (kNN + SNN) -> SNN used by FindClusters
  obj <- FindNeighbors(obj, reduction = reduction, dims = 1:n_pcs, k.param = n_neighbors,
                       nn.method = "annoy", annoy.metric = "cosine", n.trees = annoy_n_trees,
                       compute.SNN = TRUE, prune.SNN = 1/15, graph.name = c(nn_graph, snn_graph))
  
  ## check: (a) and (b) are the SAME kNN, and the cell itself is neighbour 1
  idx <- Indices(obj[[nn_name]])                                       # cells x k
  knn_adj <- sparseMatrix(i = rep(seq_len(nrow(idx)), ncol(idx)), j = as.vector(idx), x = 1,
                          dims = c(nrow(idx), nrow(idx)))
  pct_same <- 100 * mean(Matrix::rowSums(knn_adj * as(obj[[nn_graph]], "dgCMatrix")) == ncol(idx))
  message("Cells with identical kNN in UMAP input and clustering graph: ", round(pct_same, 2), "%")
  if (pct_same < 100) warning("kNN differs between the two FindNeighbors calls for some cells")
  message("Cells with themselves as neighbour 1: ", round(100 * mean(idx[, 1] == seq_len(nrow(idx))), 2), "%")
  saveRDS(obj[[nn_name]], file.path(branch_dir, "STEP_A_kNN_Seurat_Neighbor.RDS"))
  obj[[nn_graph]] <- NULL                                             # not needed further (memory)
  rm(knn_adj); invisible(gc())
  
  ## ---- STEP C: UMAP from the kNN (RunUMAP) -------------------------------
  obj <- RunUMAP(obj, nn.name = nn_name, umap.method = "uwot", n.neighbors = n_neighbors,
                 min.dist = umap_min_dist, spread = umap_spread, seed.use = seed,
                 reduction.name = umap_name, reduction.key = paste0(toupper(prefix), "UMAP_"))
  
  ## ---- STEP D: Leiden on the SNN graph (FindClusters), every resolution ---
  for (i in seq_along(resolution_sweep)) {
    obj <- FindClusters(obj, graph.name = snn_graph, algorithm = 4, resolution = resolution_sweep[i],
                        n.iter = leiden_n_iterations, random.seed = seed, cluster.name = cl_cols[i])
  }
  clusters <- obj[[main_col]][, 1]                                     # 0 = largest cluster
  message("Clusters at resolution ", leiden_resolution, ": ", nlevels(clusters))
  
  ## ---- resolution sweep table ---------------------------------------------
  sweep <- rbindlist(lapply(seq_along(resolution_sweep), function(i) {
    cl <- obj[[cl_cols[i]]][, 1]; tab <- table(cl)
    data.table(resolution = resolution_sweep[i], n_clusters = length(tab), smallest_cluster = min(tab),
               largest_cluster = max(tab), n_clusters_lt_50_cells = sum(tab < 50),
               ARI_vs_resolution_main = round(igraph::compare(as.integer(cl), as.integer(clusters),
                                                              method = "adjusted.rand"), 4))
  }))
  print(sweep)
  fwrite(sweep, file.path(branch_dir, "Leiden_resolution_sweep_same_graph.csv"))
  
  ## ---- tables ---------------------------------------------------------------
  cells <- data.table(cell_id = colnames(obj), section_id = sec, nCount_RNA = obj$nCount_RNA,
                      Embeddings(obj, umap_name), obj@meta.data[, cl_cols])
  fwrite(cells, file.path(branch_dir, "cells_UMAP_and_Leiden_clusters.csv"))
  
  sizes <- cells[, .(n_cells = .N), by = .(leiden = get(main_col))][order(leiden)]
  sizes[, pct_of_cells := round(100 * n_cells / sum(n_cells), 2)]
  fwrite(sizes, file.path(branch_dir, "Leiden_cluster_sizes.csv"))
  
  comp <- cells[, .(n_cells = .N), by = .(leiden = get(main_col), section_id)]
  comp[, pct_of_cluster := round(100 * n_cells / sum(n_cells), 2), by = leiden]
  comp_wide <- dcast(comp, leiden ~ section_id, value.var = "n_cells", fill = 0)
  dom <- comp[order(-n_cells), .SD[1], by = leiden][, .(leiden, dominant_section = section_id,
                                                        pct_in_dominant_section = pct_of_cluster)]
  comp_wide <- merge(comp_wide, dom, by = "leiden")
  comp_wide[, flag_for_investigation := pct_in_dominant_section > dominant_section_flag_pct]
  fwrite(comp_wide[order(leiden)], file.path(branch_dir, "Leiden_cluster_composition_by_section.csv"))
  
  ## section mixing on the kNN: own-section neighbours (self excluded) vs random expectation
  sec_int   <- as.integer(factor(sec))
  same_frac <- rowMeans(matrix(sec_int[idx[, -1]], nrow = nrow(idx)) == sec_int)
  Ns  <- table(sec); N <- length(sec)
  mix <- data.table(section_id = sec, same_frac, expected = (as.numeric(Ns[sec]) - 1) / (N - 1))[
    , .(n_cells = .N, mean_same_section_neighbour_fraction = round(mean(same_frac), 4),
        expected_if_randomly_mixed = round(expected[1], 4),
        ratio_observed_over_expected = round(mean(same_frac) / expected[1], 2)), by = section_id][order(section_id)]
  print(mix)
  fwrite(mix, file.path(branch_dir, "section_mixing_on_kNN_graph.csv"))
  
  ## ---- plots (Seurat) -----------------------------------------------------
  cols_main <- cluster_colours(levels(clusters))
  save_tiff(DimPlot(obj, reduction = umap_name, group.by = main_col, cols = cols_main, label = TRUE,
                    label.box = TRUE, shuffle = TRUE, seed = seed) +
              labs(title = paste0(branch_name, ": UMAP by Leiden cluster"),
                   subtitle = paste0("res ", leiden_resolution, ", ", nlevels(clusters), " clusters; ",
                                     n_pcs, " PCs, ", n_neighbors, " neighbours")),
            file.path(branch_dir, "UMAP_by_Leiden_cluster.tiff"), 11, 9)
  
  save_tiff(DimPlot(obj, reduction = umap_name, group.by = batch_col, shuffle = TRUE, seed = seed) +
              labs(title = paste0(branch_name, ": UMAP by tissue section")),
            file.path(branch_dir, "UMAP_by_section.tiff"), 11, 9)
  
  save_tiff(DimPlot(obj, reduction = umap_name, group.by = main_col, split.by = batch_col, ncol = 3,
                    cols = cols_main) + NoLegend() +
              labs(title = paste0(branch_name, ": UMAP split by section")),
            file.path(branch_dir, "UMAP_by_Leiden_split_by_section.tiff"), 15, 10)
  
  save_tiff(FeaturePlot(obj, features = "nCount_RNA", reduction = umap_name, order = FALSE) +
              scale_colour_viridis_c(option = "B", trans = "log10") +
              labs(title = paste0(branch_name, ": UMAP by total counts (nCount_RNA)")),
            file.path(branch_dir, "UMAP_by_nCount_RNA.tiff"), 11, 9)
  
  if (!is.na(img_name)) {                                                # clusters in tissue space
    save_tiff(ImageDimPlot(obj, fov = img_name, group.by = main_col, cols = cols_main, size = 0.3,
                           dark.background = FALSE, axes = TRUE) +
                labs(title = paste0(branch_name, ": Leiden clusters in tissue space")),
              file.path(branch_dir, "Spatial_Leiden_clusters_all_sections.tiff"), 13, 10)
  }
  
  p_comp <- ggplot(comp, aes(leiden, pct_of_cluster, fill = section_id)) +   # no Seurat function for this
    geom_col(width = 0.8) + geom_hline(yintercept = dominant_section_flag_pct, linetype = "dashed") +
    theme_minimal(base_size = 13) +
    labs(title = paste0(branch_name, ": section composition of each cluster"),
         x = "Leiden cluster", y = "% of cluster", fill = "Section")
  save_tiff(p_comp, file.path(branch_dir, "Leiden_cluster_composition_by_section.tiff"), 12, 6)
  
  ## ---- resolution sweep plots ---------------------------------------------
  p_sum <- ggplot(melt(sweep[, .(resolution, n_clusters, n_clusters_lt_50_cells, smallest_cluster,
                                 ARI_vs_resolution_main)], id.vars = "resolution"), aes(resolution, value)) +
    geom_vline(xintercept = leiden_resolution, colour = "red", linetype = "dashed") +
    geom_line() + geom_point(size = 2.5) + scale_x_log10(breaks = resolution_sweep) +
    facet_wrap(~ variable, ncol = 2, scales = "free_y") + theme_bw(base_size = 13) +
    labs(title = paste0(branch_name, ": Leiden resolution sweep (same graph)"), x = "Leiden resolution", y = NULL)
  save_tiff(p_sum, file.path(branch_dir, "Leiden_resolution_sweep_summary.tiff"), 11, 8)
  
  max_k <- max(sweep$n_clusters)                                        # one UMAP panel per resolution
  p_panels <- DimPlot(obj, reduction = umap_name, group.by = cl_cols, ncol = 4, label = TRUE, label.size = 2.5,
                      cols = cluster_colours(as.character(0:(max_k - 1)))) & NoLegend()
  save_tiff(p_panels, file.path(branch_dir, "Leiden_resolution_sweep_UMAP_panels.tiff"), 20, 11)
  
  if (requireNamespace("clustree", quietly = TRUE)) {                   # cluster tree across resolutions
    save_tiff(clustree::clustree(obj, prefix = paste0(prefix, "_leiden_res_")) +
                labs(title = paste0(branch_name, ": cluster tree across resolutions")),
              file.path(branch_dir, "Leiden_resolution_sweep_cluster_tree.tiff"), 18, 11)
  } else message("install.packages('clustree') for the cluster-tree plot")
  
  message(branch_name, " done -> ", branch_dir)
  obj
}


## =============================================================================
## CHUNK 5: PRIMARY - batch-aware Pearson PCA
## =============================================================================
obj <- run_branch_seurat(obj, "pearson_pca", "primary", dir_primary, "PRIMARY (batch-aware Pearson PCA)")


## =============================================================================
## CHUNK 6: SENSITIVITY - batch-aware Pearson PCA + Harmony
## =============================================================================
obj <- run_branch_seurat(obj, "harmony", "harmony", dir_sens, "SENSITIVITY (batch-aware Pearson PCA + Harmony)")


## =============================================================================
## CHUNK 7: PRIMARY vs SENSITIVITY at the main resolution
## =============================================================================
cl_p <- obj[[paste0("primary_leiden_res_", leiden_resolution)]][, 1]
cl_s <- obj[[paste0("harmony_leiden_res_", leiden_resolution)]][, 1]
ari  <- igraph::compare(as.integer(cl_p), as.integer(cl_s), method = "adjusted.rand")
message("ARI PRIMARY vs SENSITIVITY (res ", leiden_resolution, "): ", round(ari, 4))

ct <- as.data.table(table(primary = cl_p, sensitivity = cl_s))
setnames(ct, "N", "n_cells")
ct[, `:=`(primary = factor(primary, levels = levels(cl_p)), sensitivity = factor(sensitivity, levels = levels(cl_s)))]
ct[, pct_of_primary_cluster := round(100 * n_cells / sum(n_cells), 2), by = primary]
fwrite(dcast(ct, primary ~ sensitivity, value.var = "n_cells", fill = 0),
       file.path(dir_compare, "crosstab_cells_PRIMARY_rows_x_SENSITIVITY_cols.csv"))
fwrite(data.table(resolution = leiden_resolution, n_clusters_primary = nlevels(cl_p),
                  n_clusters_sensitivity = nlevels(cl_s), adjusted_rand_index = round(ari, 4)),
       file.path(dir_compare, "ARI_PRIMARY_vs_SENSITIVITY.csv"))
p_ct <- ggplot(ct, aes(sensitivity, primary, fill = pct_of_primary_cluster)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(pct_of_primary_cluster >= 5, sprintf("%.0f", pct_of_primary_cluster), "")), size = 2.8) +
  scale_fill_viridis_c() + theme_minimal(base_size = 12) +
  labs(title = "Leiden clusters: PRIMARY vs SENSITIVITY (Seurat functions)",
       subtitle = paste0("Each row = 100% of a PRIMARY cluster; ARI = ", round(ari, 3)),
       x = "SENSITIVITY cluster", y = "PRIMARY cluster", fill = "% of\nPRIMARY\ncluster")
save_tiff(p_ct, file.path(dir_compare, "crosstab_heatmap_PRIMARY_vs_SENSITIVITY.tiff"), 12, 10)

## side by side UMAPs, same colours per analysis
save_tiff(DimPlot(obj, reduction = "umap_primary", group.by = paste0("primary_leiden_res_", leiden_resolution),
                  label = TRUE, cols = cluster_colours(levels(cl_p))) + NoLegend() + labs(title = "PRIMARY") |
            DimPlot(obj, reduction = "umap_harmony", group.by = paste0("harmony_leiden_res_", leiden_resolution),
                    label = TRUE, cols = cluster_colours(levels(cl_s))) + NoLegend() + labs(title = "SENSITIVITY (Harmony)"),
          file.path(dir_compare, "UMAP_PRIMARY_and_SENSITIVITY_side_by_side.tiff"), 20, 9)


## =============================================================================
## CHUNK 8: save the Seurat object (reductions + all cluster columns)
## =============================================================================
## The kNN objects were saved in each branch folder; graphs and neighbours are
## dropped here because at k = 250 they would make the file several GB.
obj@graphs    <- list()
obj@neighbors <- list()
seurat_out <- file.path(dir_seu, paste0("seurat_", input_tag, "_SeuratFunctions_", run_tag, ".RDS"))
saveRDS(obj, seurat_out)
writeLines(capture.output(sessionInfo()), file.path(dir_seu, "sessionInfo.txt"))
message("Saved: ", seurat_out, "\nAll outputs in: ", dir_seu)



## =============================================================================
## CHUNK 9: side-by-side UMAPs (Leiden | counts) and (Leiden | tissue section)
## =============================================================================
## Uses the UMAPs and clusters already in obj (nothing recomputed).
## Saved inside each branch folder, in the same sub-folders as script 13:
##   UMAP_leiden_nCount_RNA/             : Leiden UMAP beside nCount_RNA UMAP
##   UMAP_leiden_color_cells_by_tissue/  : Leiden UMAP beside tissue-section UMAP
tissue_cols <- c(AAA2 = "#2CA02C", AAA38 = "#1F4EDB", AAA40 = "#E31A1C",   # same section colours as script 13
                 AAA45 = "#FF8C00", AAA46 = "#8E44AD", AAA48 = "#E6B800")
plot_tag <- paste0(run_tag, "_res", leiden_resolution)                       # "PC20_NN250_res1" (file names)

## one entry per analysis: UMAP reduction, cluster column, folder, label
branches_sbs <- list(
  list(umap = "umap_primary", col = paste0("primary_leiden_res_", leiden_resolution),
       dir = dir_primary, label = "PRIMARY (batch-aware Pearson PCA)"),
  list(umap = "umap_harmony", col = paste0("harmony_leiden_res_", leiden_resolution),
       dir = dir_sens,    label = "SENSITIVITY (batch-aware Pearson PCA + Harmony)"))

for (b in branches_sbs) {
  dir_counts <- file.path(b$dir, "UMAP_leiden_nCount_RNA")                  # output sub-folders
  dir_tissue <- file.path(b$dir, "UMAP_leiden_color_cells_by_tissue")
  dir.create(dir_counts, showWarnings = FALSE); dir.create(dir_tissue, showWarnings = FALSE)
  sub_txt <- paste0(n_pcs, " PCs, ", n_neighbors, " neighbours, resolution ", leiden_resolution)
  
  ## 1. Leiden clusters (labels on clusters, colours as in the other plots)
  cl <- obj[[b$col]][, 1]                                                   # cluster of each cell
  p_leiden <- DimPlot(obj, reduction = b$umap, group.by = b$col, cols = cluster_colours(levels(cl)),
                      label = TRUE, label.box = TRUE, shuffle = TRUE, seed = seed) +
    labs(title = paste0("Leiden clusters (", nlevels(cl), ")"), subtitle = sub_txt)
  
  ## 2. total RNA counts per cell, log10 colour scale
  p_counts <- FeaturePlot(obj, features = "nCount_RNA", reduction = b$umap) +
    scale_colour_viridis_c(option = "B", trans = "log10", labels = scales::label_comma()) +
    labs(title = "Total RNA counts per cell (nCount_RNA)", subtitle = sub_txt, colour = "nCount_RNA")
  
  ## 3. tissue section of each cell
  p_tissue <- DimPlot(obj, reduction = b$umap, group.by = batch_col, cols = tissue_cols,
                      shuffle = TRUE, seed = seed) +
    labs(title = "Tissue section", subtitle = sub_txt)
  
  ## side by side (patchwork |) and save
  save_tiff((p_leiden | p_counts) + plot_annotation(title = paste0(b$label, ": Leiden vs total counts")),
            file.path(dir_counts, paste0("UMAP_Leiden_and_nCount_RNA_side_by_side_", plot_tag, ".tiff")), 22, 10)
  save_tiff((p_leiden | p_tissue) + plot_annotation(title = paste0(b$label, ": Leiden vs tissue section")),
            file.path(dir_tissue, paste0("UMAP_Leiden_and_tissue_section_side_by_side_", plot_tag, ".tiff")), 22, 10)
}



































install.packages("openxlsx", dependencies = TRUE)


library(openxlsx) 



install.packages("reticulate")
reticulate::virtualenv_create("r-celltypist", packages = c("celltypist", "scanpy"))
# if no Python is found: reticulate::install_python(), then repeat the line above


library(reticulate)



## =============================================================================
## 12b_Seurat_functions_PCA_UMAP_Leiden.R (2) FINDALLMARKERS
## =============================================================================
## PURPOSE
##   Same analysis as 12_Normalisation_PCA_UMAP_Leiden.R, but the neighbour
##   search, UMAP, Leiden clustering and plots use SEURAT functions:
##     FindNeighbors()  -> kNN on the PCs (cosine, annoy)
##     RunUMAP()        -> UMAP from that kNN (nn.name = ...)
##     FindClusters()   -> Leiden (algorithm = 4) on the neighbour graph from
##                         the SAME kNN (never on UMAP coordinates)
##     ElbowPlot(), DimPlot(), FeaturePlot(), ImageDimPlot() -> plots
##
##   PCA is NOT recomputed. Script 12's batch-aware Pearson-residual PCA and
##   Harmony embeddings are loaded and added to the Seurat object as reductions.
##   (Seurat has no batch-aware Pearson-residual PCA. SCTransform also makes
##    Pearson residuals, but from a different model - regularised negative
##    binomial, not scPearsonPCA's per-section quasi-Poisson - so it is not used.)
##
##   PRIMARY     : pearson_pca          -> FindNeighbors -> RunUMAP / FindClusters
##   SENSITIVITY : pearson_pca -> harmony -> FindNeighbors -> RunUMAP / FindClusters
##   Settings: 20 PCs, 250 neighbours (incl. the cell itself), main resolution 1,
##   plus the resolution sweep on the same graph.
##
## ONE DIFFERENCE FROM SCRIPT 12 (Seurat's design, not a choice made here)
##   Script 12 (scanpy-style) clustered the UMAP fuzzy-connectivity graph.
##   Seurat's FindClusters clusters the shared-nearest-neighbour (SNN) graph that
##   FindNeighbors builds from the kNN. UMAP and Leiden still come from the SAME
##   kNN. Clusters can therefore differ slightly from script 12.
##   NOTE: with k = 250 the SNN graph is large (needs several GB of RAM).
##
## OUTPUT: <script 12 folder>/seurat_functions_used/PC20_NN250/
##   00_shared/ , 01_PRIMARY_.../ , 02_SENSITIVITY_.../ , 03_PRIMARY_vs_SENSITIVITY/
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 0: packages
## =============================================================================
library(data.table)   # tables
library(Matrix)       # sparse matrices (kNN check)
library(ggplot2)      # extra plot layers
library(patchwork)    # combines Seurat plots (the & operator)
library(Seurat)       # FindNeighbors, RunUMAP, FindClusters, DimPlot, ...
## FindClusters(algorithm = 4) needs the 'leidenbase' package (Seurat >= 5.1):
## install.packages("leidenbase") if FindClusters asks for it.


## =============================================================================
## CHUNK 1: paths + parameters (same values as script 12)
## =============================================================================
output_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"
seurat_in  <- file.path(output_dir, "seurat_QC_filtered_per_section_03_to_09.RDS")  # script 11 (read only)
input_tag  <- sub("\\.RDS$", "", basename(seurat_in))
dir_main   <- file.path(output_dir, "12_Normalisation_PCA_UMAP_Leiden", input_tag)  # script 12 folder

batch_col           <- "sample_id"    # the 6 AAA sections
n_pcs               <- 20             # PCs used by FindNeighbors
n_neighbors         <- 250            # k.param (Seurat counts the cell itself)
annoy_n_trees       <- 50             # annoy index size (as script 12)
umap_min_dist       <- 0.02           # as script 12
umap_spread         <- 2              # as script 12
leiden_resolution   <- 1              # main resolution
leiden_n_iterations <- 2              # as script 12
resolution_sweep    <- c(0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0)  # all on the same graph
dominant_section_flag_pct <- 90       # flag clusters with > 90% cells from one section
seed    <- 42
dpi_out <- 300
stopifnot(leiden_resolution %in% resolution_sweep)

run_tag     <- paste0("PC", n_pcs, "_NN", n_neighbors)                     # "PC20_NN250"
dir_seu     <- file.path(dir_main, "seurat_functions_used", run_tag)       # all outputs here
dir_shared  <- file.path(dir_seu, "00_shared")
dir_primary <- file.path(dir_seu, "01_PRIMARY_batchawarePCA")
dir_sens    <- file.path(dir_seu, "02_SENSITIVITY_batchawarePCA_Harmony")
dir_compare <- file.path(dir_seu, "03_PRIMARY_vs_SENSITIVITY")
for (d in c(dir_shared, dir_primary, dir_sens, dir_compare)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## script 12 inputs reused (not recomputed)
pca_file  <- file.path(dir_main, "00_shared_HVGs_and_batch_aware_Pearson_PCA", "pcaobj_batch_aware.RDS")
harm_file <- file.path(dir_main, paste0("02_SENSITIVITY_batchawarePCA_Harmony_", run_tag), "harmony_embeddings.RDS")
stopifnot("Run script 12 first (batch-aware PCA not found)" = file.exists(pca_file))

fwrite(data.table(
  parameter = c("input", "PCA", "Harmony", "n_PCs", "k.param (incl. self)", "nn.method", "annoy.metric",
                "n.trees", "UMAP", "min.dist", "spread", "clustering", "graph clustered", "prune.SNN",
                "Leiden n.iter", "main resolution", "resolution sweep", "seed"),
  value = c(basename(seurat_in), "script 12 scPearsonPCA (batch-aware Pearson residuals), loaded",
            "script 12 RunHarmony on all PCs, loaded", n_pcs, n_neighbors, "annoy", "cosine",
            annoy_n_trees, "Seurat::RunUMAP(nn.name = kNN, umap.method = 'uwot')", umap_min_dist,
            umap_spread, "Seurat::FindClusters(algorithm = 4, Leiden)", "SNN from FindNeighbors (same kNN)",
            "1/15 (Seurat default)", leiden_n_iterations, leiden_resolution,
            paste(resolution_sweep, collapse = ", "), seed)),
  file.path(dir_seu, "parameters_used.csv"))


## =============================================================================
## CHUNK 2: load object, add script 12's PCA and Harmony as Seurat reductions
## =============================================================================
obj <- readRDS(seurat_in)
cells_order <- colnames(obj)                                   # cell order of the object

## batch-aware Pearson PCA (all 50 PCs + their standard deviations)
pcaobj  <- readRDS(pca_file)
emb_pca <- pcaobj$reduction.data@cell.embeddings
if (is.null(rownames(emb_pca))) rownames(emb_pca) <- cells_order  # as in script 12
emb_pca <- emb_pca[cells_order, , drop = FALSE]                  # same order as the object
pc_sd   <- pcaobj$reduction.data@stdev                           # for ElbowPlot
rm(pcaobj); invisible(gc())

## Harmony on all PCs: load script 12's result; if absent, the identical call
if (file.exists(harm_file)) {
  emb_harmony <- readRDS(harm_file)
} else {
  set.seed(seed)
  emb_harmony <- as.matrix(harmony::RunHarmony(emb_pca, meta_data = data.frame(sample_id = obj[[batch_col]][, 1]),
                                               vars_use = "sample_id"))
  rownames(emb_harmony) <- rownames(emb_pca)
}
emb_harmony <- emb_harmony[cells_order, , drop = FALSE]
stopifnot(all(is.finite(emb_pca)), all(is.finite(emb_harmony)))

## Seurat needs column names = key + number (PPC_1, PPC_2, ...)
with_key <- function(m, key) { m <- as.matrix(m); colnames(m) <- paste0(key, seq_len(ncol(m))); m }
obj[["pearson_pca"]] <- CreateDimReducObject(embeddings = with_key(emb_pca, "PPC_"), stdev = pc_sd,
                                             key = "PPC_", assay = "RNA")
obj[["harmony"]]     <- CreateDimReducObject(embeddings = with_key(emb_harmony, "HARM_"),
                                             key = "HARM_", assay = "RNA")
rm(emb_pca, emb_harmony); invisible(gc())
obj   # now lists the two reductions

## Elbow plot (Seurat), n_pcs marked
p_elbow <- ElbowPlot(obj, reduction = "pearson_pca", ndims = length(pc_sd)) +
  geom_vline(xintercept = n_pcs, colour = "red", linetype = "dashed") +
  labs(title = "Batch-aware Pearson PCA", subtitle = paste0("Red line = ", n_pcs, " PCs used"))
ggsave(file.path(dir_shared, "ElbowPlot_pearson_pca.tiff"), p_elbow, width = 8, height = 5.5,
       dpi = dpi_out, compression = "lzw", bg = "white")


## =============================================================================
## CHUNK 3: helpers - colours and saving (same palette as script 12)
## =============================================================================
distinct_palette <- c(
  "#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#42D4F4", "#F032E6",
  "#9A6324", "#469990", "#BFEF45", "#800000", "#FFE119", "#000075", "#808000",
  "#FF8FB1", "#A9A9A9", "#005300", "#009FFF", "#C77CFF", "#FF6F61", "#00FFBE",
  "#D4A017", "#7F7FFF", "#8B008B", "#2F4F4F", "#E9967A", "#556B2F", "#FF1493",
  "#20B2AA", "#B8860B", "#6A5ACD", "#CD853F", "#4B0082", "#00CED1", "#8FBC8F",
  "#A0522D")
cluster_colours <- function(lv) {                 # one colour per cluster level
  if (length(lv) <= length(distinct_palette)) return(setNames(distinct_palette[seq_along(lv)], lv))
  stopifnot("install.packages('Polychrome') for > 36 clusters" = requireNamespace("Polychrome", quietly = TRUE))
  set.seed(seed)
  setNames(unname(Polychrome::createPalette(length(lv), seedcolors = distinct_palette, range = c(25, 80))), lv)
}
save_tiff <- function(p, path, w, h) {
  ggsave(path, p, width = w, height = h, dpi = dpi_out, device = "tiff", compression = "lzw", bg = "white")
  message("Saved: ", path)
}
img_name <- Images(obj)[1]                        # spatial FOV (cell centroids), for ImageDimPlot


## =============================================================================
## CHUNK 4: one analysis branch with Seurat functions
## =============================================================================
## reduction = "pearson_pca" (PRIMARY) or "harmony" (SENSITIVITY)
## prefix    = short name used for graphs, UMAP and cluster columns
run_branch_seurat <- function(obj, reduction, prefix, branch_dir, branch_name) {
  message("\n", strrep("=", 70), "\n", branch_name, "\n", strrep("=", 70))
  nn_name   <- paste0(prefix, "_neighbors")    # Seurat Neighbor object: kNN indices + distances
  nn_graph  <- paste0(prefix, "_nn")           # kNN as a graph
  snn_graph <- paste0(prefix, "_snn")          # shared-nearest-neighbour graph (clustered)
  umap_name <- paste0("umap_", prefix)
  cl_cols   <- paste0(prefix, "_leiden_res_", resolution_sweep)        # one column per resolution
  main_col  <- paste0(prefix, "_leiden_res_", leiden_resolution)
  sec       <- obj[[batch_col]][, 1]
  
  ## ---- PC1 vs PC2 of the input reduction, by section --------------------
  save_tiff(DimPlot(obj, reduction = reduction, dims = 1:2, group.by = batch_col, shuffle = TRUE, seed = seed) +
              labs(title = paste0(branch_name, ": dimension 1 vs 2")),
            file.path(branch_dir, "PC1_vs_PC2_by_section.tiff"), 9, 8)
  
  ## ---- STEP A: kNN (FindNeighbors) ---------------------------------------
  ## (a) kNN kept as a Neighbor object (indices + distances) -> used by RunUMAP
  obj <- FindNeighbors(obj, reduction = reduction, dims = 1:n_pcs, k.param = n_neighbors,
                       nn.method = "annoy", annoy.metric = "cosine", n.trees = annoy_n_trees,
                       return.neighbor = TRUE, graph.name = nn_name)
  ## (b) same search, returned as graphs (kNN + SNN) -> SNN used by FindClusters
  obj <- FindNeighbors(obj, reduction = reduction, dims = 1:n_pcs, k.param = n_neighbors,
                       nn.method = "annoy", annoy.metric = "cosine", n.trees = annoy_n_trees,
                       compute.SNN = TRUE, prune.SNN = 1/15, graph.name = c(nn_graph, snn_graph))
  
  ## check: (a) and (b) are the SAME kNN, and the cell itself is neighbour 1
  idx <- Indices(obj[[nn_name]])                                       # cells x k
  knn_adj <- sparseMatrix(i = rep(seq_len(nrow(idx)), ncol(idx)), j = as.vector(idx), x = 1,
                          dims = c(nrow(idx), nrow(idx)))
  pct_same <- 100 * mean(Matrix::rowSums(knn_adj * as(obj[[nn_graph]], "dgCMatrix")) == ncol(idx))
  message("Cells with identical kNN in UMAP input and clustering graph: ", round(pct_same, 2), "%")
  if (pct_same < 100) warning("kNN differs between the two FindNeighbors calls for some cells")
  message("Cells with themselves as neighbour 1: ", round(100 * mean(idx[, 1] == seq_len(nrow(idx))), 2), "%")
  saveRDS(obj[[nn_name]], file.path(branch_dir, "STEP_A_kNN_Seurat_Neighbor.RDS"))
  obj[[nn_graph]] <- NULL                                             # not needed further (memory)
  rm(knn_adj); invisible(gc())
  
  ## ---- STEP C: UMAP from the kNN (RunUMAP) -------------------------------
  obj <- RunUMAP(obj, nn.name = nn_name, umap.method = "uwot", n.neighbors = n_neighbors,
                 min.dist = umap_min_dist, spread = umap_spread, seed.use = seed,
                 reduction.name = umap_name, reduction.key = paste0(toupper(prefix), "UMAP_"))
  
  ## ---- STEP D: Leiden on the SNN graph (FindClusters), every resolution ---
  for (i in seq_along(resolution_sweep)) {
    obj <- FindClusters(obj, graph.name = snn_graph, algorithm = 4, resolution = resolution_sweep[i],
                        n.iter = leiden_n_iterations, random.seed = seed, cluster.name = cl_cols[i])
  }
  clusters <- obj[[main_col]][, 1]                                     # 0 = largest cluster
  message("Clusters at resolution ", leiden_resolution, ": ", nlevels(clusters))
  
  ## ---- resolution sweep table ---------------------------------------------
  sweep <- rbindlist(lapply(seq_along(resolution_sweep), function(i) {
    cl <- obj[[cl_cols[i]]][, 1]; tab <- table(cl)
    data.table(resolution = resolution_sweep[i], n_clusters = length(tab), smallest_cluster = min(tab),
               largest_cluster = max(tab), n_clusters_lt_50_cells = sum(tab < 50),
               ARI_vs_resolution_main = round(igraph::compare(as.integer(cl), as.integer(clusters),
                                                              method = "adjusted.rand"), 4))
  }))
  print(sweep)
  fwrite(sweep, file.path(branch_dir, "Leiden_resolution_sweep_same_graph.csv"))
  
  ## ---- tables ---------------------------------------------------------------
  cells <- data.table(cell_id = colnames(obj), section_id = sec, nCount_RNA = obj$nCount_RNA,
                      Embeddings(obj, umap_name), obj@meta.data[, cl_cols])
  fwrite(cells, file.path(branch_dir, "cells_UMAP_and_Leiden_clusters.csv"))
  
  sizes <- cells[, .(n_cells = .N), by = .(leiden = get(main_col))][order(leiden)]
  sizes[, pct_of_cells := round(100 * n_cells / sum(n_cells), 2)]
  fwrite(sizes, file.path(branch_dir, "Leiden_cluster_sizes.csv"))
  
  comp <- cells[, .(n_cells = .N), by = .(leiden = get(main_col), section_id)]
  comp[, pct_of_cluster := round(100 * n_cells / sum(n_cells), 2), by = leiden]
  comp_wide <- dcast(comp, leiden ~ section_id, value.var = "n_cells", fill = 0)
  dom <- comp[order(-n_cells), .SD[1], by = leiden][, .(leiden, dominant_section = section_id,
                                                        pct_in_dominant_section = pct_of_cluster)]
  comp_wide <- merge(comp_wide, dom, by = "leiden")
  comp_wide[, flag_for_investigation := pct_in_dominant_section > dominant_section_flag_pct]
  fwrite(comp_wide[order(leiden)], file.path(branch_dir, "Leiden_cluster_composition_by_section.csv"))
  
  ## section mixing on the kNN: own-section neighbours (self excluded) vs random expectation
  sec_int   <- as.integer(factor(sec))
  same_frac <- rowMeans(matrix(sec_int[idx[, -1]], nrow = nrow(idx)) == sec_int)
  Ns  <- table(sec); N <- length(sec)
  mix <- data.table(section_id = sec, same_frac, expected = (as.numeric(Ns[sec]) - 1) / (N - 1))[
    , .(n_cells = .N, mean_same_section_neighbour_fraction = round(mean(same_frac), 4),
        expected_if_randomly_mixed = round(expected[1], 4),
        ratio_observed_over_expected = round(mean(same_frac) / expected[1], 2)), by = section_id][order(section_id)]
  print(mix)
  fwrite(mix, file.path(branch_dir, "section_mixing_on_kNN_graph.csv"))
  
  ## ---- plots (Seurat) -----------------------------------------------------
  cols_main <- cluster_colours(levels(clusters))
  save_tiff(DimPlot(obj, reduction = umap_name, group.by = main_col, cols = cols_main, label = TRUE,
                    label.box = TRUE, shuffle = TRUE, seed = seed) +
              labs(title = paste0(branch_name, ": UMAP by Leiden cluster"),
                   subtitle = paste0("res ", leiden_resolution, ", ", nlevels(clusters), " clusters; ",
                                     n_pcs, " PCs, ", n_neighbors, " neighbours")),
            file.path(branch_dir, "UMAP_by_Leiden_cluster.tiff"), 11, 9)
  
  save_tiff(DimPlot(obj, reduction = umap_name, group.by = batch_col, shuffle = TRUE, seed = seed) +
              labs(title = paste0(branch_name, ": UMAP by tissue section")),
            file.path(branch_dir, "UMAP_by_section.tiff"), 11, 9)
  
  save_tiff(DimPlot(obj, reduction = umap_name, group.by = main_col, split.by = batch_col, ncol = 3,
                    cols = cols_main) + NoLegend() +
              labs(title = paste0(branch_name, ": UMAP split by section")),
            file.path(branch_dir, "UMAP_by_Leiden_split_by_section.tiff"), 15, 10)
  
  save_tiff(FeaturePlot(obj, features = "nCount_RNA", reduction = umap_name, order = FALSE) +
              scale_colour_viridis_c(option = "B", trans = "log10") +
              labs(title = paste0(branch_name, ": UMAP by total counts (nCount_RNA)")),
            file.path(branch_dir, "UMAP_by_nCount_RNA.tiff"), 11, 9)
  
  if (!is.na(img_name)) {                                                # clusters in tissue space
    save_tiff(ImageDimPlot(obj, fov = img_name, group.by = main_col, cols = cols_main, size = 0.3,
                           dark.background = FALSE, axes = TRUE) +
                labs(title = paste0(branch_name, ": Leiden clusters in tissue space")),
              file.path(branch_dir, "Spatial_Leiden_clusters_all_sections.tiff"), 13, 10)
  }
  
  p_comp <- ggplot(comp, aes(leiden, pct_of_cluster, fill = section_id)) +   # no Seurat function for this
    geom_col(width = 0.8) + geom_hline(yintercept = dominant_section_flag_pct, linetype = "dashed") +
    theme_minimal(base_size = 13) +
    labs(title = paste0(branch_name, ": section composition of each cluster"),
         x = "Leiden cluster", y = "% of cluster", fill = "Section")
  save_tiff(p_comp, file.path(branch_dir, "Leiden_cluster_composition_by_section.tiff"), 12, 6)
  
  ## ---- resolution sweep plots ---------------------------------------------
  p_sum <- ggplot(melt(sweep[, .(resolution, n_clusters, n_clusters_lt_50_cells, smallest_cluster,
                                 ARI_vs_resolution_main)], id.vars = "resolution"), aes(resolution, value)) +
    geom_vline(xintercept = leiden_resolution, colour = "red", linetype = "dashed") +
    geom_line() + geom_point(size = 2.5) + scale_x_log10(breaks = resolution_sweep) +
    facet_wrap(~ variable, ncol = 2, scales = "free_y") + theme_bw(base_size = 13) +
    labs(title = paste0(branch_name, ": Leiden resolution sweep (same graph)"), x = "Leiden resolution", y = NULL)
  save_tiff(p_sum, file.path(branch_dir, "Leiden_resolution_sweep_summary.tiff"), 11, 8)
  
  max_k <- max(sweep$n_clusters)                                        # one UMAP panel per resolution
  p_panels <- DimPlot(obj, reduction = umap_name, group.by = cl_cols, ncol = 4, label = TRUE, label.size = 2.5,
                      cols = cluster_colours(as.character(0:(max_k - 1)))) & NoLegend()
  save_tiff(p_panels, file.path(branch_dir, "Leiden_resolution_sweep_UMAP_panels.tiff"), 20, 11)
  
  if (requireNamespace("clustree", quietly = TRUE)) {                   # cluster tree across resolutions
    save_tiff(clustree::clustree(obj, prefix = paste0(prefix, "_leiden_res_")) +
                labs(title = paste0(branch_name, ": cluster tree across resolutions")),
              file.path(branch_dir, "Leiden_resolution_sweep_cluster_tree.tiff"), 18, 11)
  } else message("install.packages('clustree') for the cluster-tree plot")
  
  message(branch_name, " done -> ", branch_dir)
  obj
}


## =============================================================================
## CHUNK 5: PRIMARY - batch-aware Pearson PCA
## =============================================================================
obj <- run_branch_seurat(obj, "pearson_pca", "primary", dir_primary, "PRIMARY (batch-aware Pearson PCA)")


## =============================================================================
## CHUNK 6: SENSITIVITY - batch-aware Pearson PCA + Harmony
## =============================================================================
obj <- run_branch_seurat(obj, "harmony", "harmony", dir_sens, "SENSITIVITY (batch-aware Pearson PCA + Harmony)")


## =============================================================================
## CHUNK 7: PRIMARY vs SENSITIVITY at the main resolution
## =============================================================================
cl_p <- obj[[paste0("primary_leiden_res_", leiden_resolution)]][, 1]
cl_s <- obj[[paste0("harmony_leiden_res_", leiden_resolution)]][, 1]
ari  <- igraph::compare(as.integer(cl_p), as.integer(cl_s), method = "adjusted.rand")
message("ARI PRIMARY vs SENSITIVITY (res ", leiden_resolution, "): ", round(ari, 4))

ct <- as.data.table(table(primary = cl_p, sensitivity = cl_s))
setnames(ct, "N", "n_cells")
ct[, `:=`(primary = factor(primary, levels = levels(cl_p)), sensitivity = factor(sensitivity, levels = levels(cl_s)))]
ct[, pct_of_primary_cluster := round(100 * n_cells / sum(n_cells), 2), by = primary]
fwrite(dcast(ct, primary ~ sensitivity, value.var = "n_cells", fill = 0),
       file.path(dir_compare, "crosstab_cells_PRIMARY_rows_x_SENSITIVITY_cols.csv"))
fwrite(data.table(resolution = leiden_resolution, n_clusters_primary = nlevels(cl_p),
                  n_clusters_sensitivity = nlevels(cl_s), adjusted_rand_index = round(ari, 4)),
       file.path(dir_compare, "ARI_PRIMARY_vs_SENSITIVITY.csv"))
p_ct <- ggplot(ct, aes(sensitivity, primary, fill = pct_of_primary_cluster)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(pct_of_primary_cluster >= 5, sprintf("%.0f", pct_of_primary_cluster), "")), size = 2.8) +
  scale_fill_viridis_c() + theme_minimal(base_size = 12) +
  labs(title = "Leiden clusters: PRIMARY vs SENSITIVITY (Seurat functions)",
       subtitle = paste0("Each row = 100% of a PRIMARY cluster; ARI = ", round(ari, 3)),
       x = "SENSITIVITY cluster", y = "PRIMARY cluster", fill = "% of\nPRIMARY\ncluster")
save_tiff(p_ct, file.path(dir_compare, "crosstab_heatmap_PRIMARY_vs_SENSITIVITY.tiff"), 12, 10)

## side by side UMAPs, same colours per analysis
save_tiff(DimPlot(obj, reduction = "umap_primary", group.by = paste0("primary_leiden_res_", leiden_resolution),
                  label = TRUE, cols = cluster_colours(levels(cl_p))) + NoLegend() + labs(title = "PRIMARY") |
            DimPlot(obj, reduction = "umap_harmony", group.by = paste0("harmony_leiden_res_", leiden_resolution),
                    label = TRUE, cols = cluster_colours(levels(cl_s))) + NoLegend() + labs(title = "SENSITIVITY (Harmony)"),
          file.path(dir_compare, "UMAP_PRIMARY_and_SENSITIVITY_side_by_side.tiff"), 20, 9)


## =============================================================================
## CHUNK 8: save the Seurat object (reductions + all cluster columns)
## =============================================================================
## The kNN objects were saved in each branch folder; graphs and neighbours are
## dropped here because at k = 250 they would make the file several GB.
obj@graphs    <- list()
obj@neighbors <- list()
seurat_out <- file.path(dir_seu, paste0("seurat_", input_tag, "_SeuratFunctions_", run_tag, ".RDS"))
saveRDS(obj, seurat_out)
writeLines(capture.output(sessionInfo()), file.path(dir_seu, "sessionInfo.txt"))
message("Saved: ", seurat_out, "\nAll outputs in: ", dir_seu)


## =============================================================================
## CHUNK 9: side-by-side UMAPs (Leiden | counts) and (Leiden | tissue section)
## =============================================================================
## Uses the UMAPs and clusters already in obj (nothing recomputed).
## Saved inside each branch folder, in the same sub-folders as script 13:
##   UMAP_leiden_nCount_RNA/             : Leiden UMAP beside nCount_RNA UMAP
##   UMAP_leiden_color_cells_by_tissue/  : Leiden UMAP beside tissue-section UMAP
tissue_cols <- c(AAA2 = "#2CA02C", AAA38 = "#1F4EDB", AAA40 = "#E31A1C",   # same section colours as script 13
                 AAA45 = "#FF8C00", AAA46 = "#8E44AD", AAA48 = "#E6B800")
plot_tag <- paste0(run_tag, "_res", leiden_resolution)                       # "PC20_NN250_res1" (file names)

## one entry per analysis: UMAP reduction, cluster column, folder, label
branches_sbs <- list(
  list(umap = "umap_primary", col = paste0("primary_leiden_res_", leiden_resolution),
       dir = dir_primary, label = "PRIMARY (batch-aware Pearson PCA)"),
  list(umap = "umap_harmony", col = paste0("harmony_leiden_res_", leiden_resolution),
       dir = dir_sens,    label = "SENSITIVITY (batch-aware Pearson PCA + Harmony)"))

for (b in branches_sbs) {
  dir_counts <- file.path(b$dir, "UMAP_leiden_nCount_RNA")                  # output sub-folders
  dir_tissue <- file.path(b$dir, "UMAP_leiden_color_cells_by_tissue")
  dir.create(dir_counts, showWarnings = FALSE); dir.create(dir_tissue, showWarnings = FALSE)
  sub_txt <- paste0(n_pcs, " PCs, ", n_neighbors, " neighbours, resolution ", leiden_resolution)
  
  ## 1. Leiden clusters (labels on clusters, colours as in the other plots)
  cl <- obj[[b$col]][, 1]                                                   # cluster of each cell
  p_leiden <- DimPlot(obj, reduction = b$umap, group.by = b$col, cols = cluster_colours(levels(cl)),
                      label = TRUE, label.box = TRUE, shuffle = TRUE, seed = seed) +
    labs(title = paste0("Leiden clusters (", nlevels(cl), ")"), subtitle = sub_txt)
  
  ## 2. total RNA counts per cell, log10 colour scale
  p_counts <- FeaturePlot(obj, features = "nCount_RNA", reduction = b$umap) +
    scale_colour_viridis_c(option = "B", trans = "log10", labels = scales::label_comma()) +
    labs(title = "Total RNA counts per cell (nCount_RNA)", subtitle = sub_txt, colour = "nCount_RNA")
  
  ## 3. tissue section of each cell
  p_tissue <- DimPlot(obj, reduction = b$umap, group.by = batch_col, cols = tissue_cols,
                      shuffle = TRUE, seed = seed) +
    labs(title = "Tissue section", subtitle = sub_txt)
  
  ## side by side (patchwork |) and save
  save_tiff((p_leiden | p_counts) + plot_annotation(title = paste0(b$label, ": Leiden vs total counts")),
            file.path(dir_counts, paste0("UMAP_Leiden_and_nCount_RNA_side_by_side_", plot_tag, ".tiff")), 22, 10)
  save_tiff((p_leiden | p_tissue) + plot_annotation(title = paste0(b$label, ": Leiden vs tissue section")),
            file.path(dir_tissue, paste0("UMAP_Leiden_and_tissue_section_side_by_side_", plot_tag, ".tiff")), 22, 10)
}


## =============================================================================
## CHUNK 10: cluster markers - FindAllMarkers (PC20_NN250, resolution 1)
## =============================================================================
## For each Leiden cluster: genes HIGHER in that cluster than in all other cells
## (one cluster vs the rest), Wilcoxon rank-sum test, as in the Seurat vignette:
##   FindAllMarkers(obj, only.pos = TRUE, test.use = "wilcox")
## Clusters are NOT recomputed: they are the FindClusters result at resolution 1
## from chunks 5-6 (columns primary_leiden_res_1 / harmony_leiden_res_1).
##
## NORMALISATION: FindAllMarkers reads the RNA "data" layer. That layer still
## holds raw counts (the Pearson residuals were used only inside the PCA), so
## Seurat's standard log-normalisation is run first, as in the Seurat vignette:
##   NormalizeData() = log1p(counts / total counts per cell x 10,000)
## This changes neither the clusters nor the UMAP.
##
## SIGNIFICANT  = p_val_adj < 0.05 (Bonferroni, Seurat's default correction)
## RANKING      = within each cluster: pct.1 high -> low (fraction of the
##                cluster's cells expressing the gene); ties by avg_log2FC
## OUTPUT (inside each branch folder, sub-folder FindAllMarkers_res1/):
##   ..._top20_per_cluster.xlsx      one sheet per cluster, top 20 ranked genes
##   ..._all_significant_ranked.csv  every significant gene, all clusters
##   ..._all_results.RDS             FindAllMarkers output before filtering
##
## Install once: install.packages("openxlsx")
## Speed: with ~100k cells this is slow. Seurat v5 uses the 'presto' package
## automatically if installed (same Wilcoxon test, much faster):
##   remotes::install_github("immunogenomics/presto")
library(openxlsx)     # write.xlsx(): one Excel sheet per list element

## fresh session: run CHUNKS 0, 1, 3 first, then this reloads the chunk 8 object
seurat_out <- file.path(dir_seu, paste0("seurat_", input_tag, "_SeuratFunctions_", run_tag, ".RDS"))
if (!exists("obj") || !"umap_primary" %in% Reductions(obj)) obj <- readRDS(seurat_out)

marker_padj_max <- 0.05   # significance cut-off on adjusted p-value
marker_top_n    <- 20     # genes per cluster sheet

## log-normalise RNA counts -> "data" layer (Seurat defaults)
obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)

## one entry per analysis: cluster-column prefix, folder, label
branches_mk <- list(
  list(prefix = "primary", dir = dir_primary, label = "PRIMARY (batch-aware Pearson PCA)"),
  list(prefix = "harmony", dir = dir_sens,    label = "SENSITIVITY (batch-aware Pearson PCA + Harmony)"))

for (b in branches_mk) {
  main_col <- paste0(b$prefix, "_leiden_res_", leiden_resolution)            # e.g. "primary_leiden_res_1"
  stopifnot("cluster column missing - run chunks 5-6 first" = main_col %in% colnames(obj[[]]))
  dir_mk <- file.path(b$dir, paste0("FindAllMarkers_res", leiden_resolution))
  dir.create(dir_mk, showWarnings = FALSE)
  file_tag <- paste0(b$prefix, "_", plot_tag)                                 # e.g. "primary_PC20_NN250_res1"
  
  ## clusters to compare, ordered 0, 1, 2, ... (not 0, 1, 10, 11, ...)
  Idents(obj) <- main_col
  Idents(obj) <- factor(Idents(obj), levels = as.character(sort(as.integer(levels(Idents(obj))))))
  message("\n", b$label, ": ", nlevels(Idents(obj)), " clusters - running FindAllMarkers...")
  
  ## each cluster vs all other cells; only genes higher in the cluster
  markers <- FindAllMarkers(obj, assay = "RNA", only.pos = TRUE, test.use = "wilcox")
  markers <- as.data.table(markers)                                            # columns: p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster, gene
  saveRDS(markers, file.path(dir_mk, paste0("FindAllMarkers_", file_tag, "_all_results.RDS")))
  
  ## significant genes, ranked within each cluster by pct.1 (then avg_log2FC)
  sig <- markers[p_val_adj < marker_padj_max & avg_log2FC > 0]
  setorder(sig, cluster, -pct.1, -avg_log2FC)                                  # cluster is a factor: 0, 1, 2, ...
  sig[, rank_in_cluster := seq_len(.N), by = cluster]
  setcolorder(sig, c("cluster", "rank_in_cluster", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"))
  fwrite(sig, file.path(dir_mk, paste0("FindAllMarkers_", file_tag, "_all_significant_ranked.csv")))
  
  ## Excel: one sheet per cluster with its top 20 genes
  cl_levels <- levels(Idents(obj))
  sheets <- lapply(cl_levels, function(k) {
    top <- sig[cluster == k & rank_in_cluster <= marker_top_n]
    if (nrow(top) == 0) data.frame(note = "no significant positive markers") else as.data.frame(top)
  })
  names(sheets) <- paste0("cluster_", cl_levels)                               # sheet names
  write.xlsx(sheets, file.path(dir_mk, paste0("FindAllMarkers_", file_tag, "_top", marker_top_n, "_per_cluster.xlsx")),
             firstRow = TRUE, colWidths = "auto", overwrite = TRUE)
  
  ## summary: cells and significant markers per cluster
  n_cells <- table(Idents(obj))
  n_sig   <- table(factor(sig$cluster, levels = cl_levels))
  smry <- data.table(cluster = cl_levels, n_cells = as.integer(n_cells[cl_levels]),
                     n_significant_markers = as.integer(n_sig[cl_levels]),
                     top5_by_pct1 = sapply(cl_levels, function(k) paste(head(sig[cluster == k, gene], 5), collapse = ", ")))
  print(smry)
  fwrite(smry, file.path(dir_mk, paste0("FindAllMarkers_", file_tag, "_summary_per_cluster.csv")))
  message(b$label, " markers saved -> ", dir_mk)
}












## 1. Restart R first: Session > Restart R

## 2. Create the Python environment with CellTypist (takes a few minutes)
reticulate::virtualenv_create("r-celltypist", packages = c("celltypist", "scanpy"))

## 3. Check it exists (should print TRUE)
reticulate::virtualenv_exists("r-celltypist")

library(reticulate)




## =============================================================================
## CHUNK 11: CellTypist cell-type call for each Leiden cluster (via reticulate)
##           -> paste at the end of 12b_Seurat_functions_PCA_UMAP_Leiden.R
## =============================================================================
## HOW CELLTYPIST WORKS (why the top-20 lists are not fed in directly)
##   CellTypist is a classifier trained on whole single-cell profiles. It cannot
##   take a list of gene names; it reads each CELL's expression, log1p-normalised
##   to 10,000 counts per cell (exactly what NormalizeData() made in chunk 10),
##   and returns a cell type per cell.
##   One call per cluster, one cluster after the other:
##     over_clustering = our Leiden clusters (resolution 1), majority_voting = TRUE
##     -> each cluster is called the label that most of its cells received.
##   The top-20 marker lists (chunk 10) are then checked against that call:
##     overlap between a cluster's top 20 markers and the model's top 20 driving
##     genes for the called cell type (model.extract_top_markers).
##     High overlap = the cluster's own markers support the call.
##   Only genes in BOTH the 6k panel and the model are used (number saved).
##
## MODELS (CellTypist model list, checked): change ct_models to add or remove.
##   Adult_Human_Vascular.pkl : vascular populations from multiple adult human organs (42 types)
##   Healthy_Adult_Heart.pkl  : healthy adult human heart, 8 regions (75 types)
##   Immune_All_Low.pkl       : immune subtypes from 20 tissues (98 types)
##   All are scRNA-seq references, not aorta or CosMx: treat calls as a guide.
##
## ONE-TIME SETUP (Python environment for CellTypist, built from R):
##   install.packages("reticulate")
##   reticulate::virtualenv_create("r-celltypist", packages = c("celltypist", "scanpy"))
##   (if R says no Python was found: reticulate::install_python(), then repeat)
## R can choose its Python only once per session, so: restart R, run
## CHUNKS 0, 1 and 3, then this chunk.
##
## OUTPUT (inside each branch folder, sub-folder FindAllMarkers_res1/):
##   CellTypist_<branch>_..._per_cluster_calls.csv    one row per cluster, call per model
##   CellTypist_<branch>_..._per_cluster_details.csv  + % cells agreeing, confidence,
##                                                    2nd label, top-20 overlap
##   CellTypist_<branch>_..._per_cell.csv.gz          label and confidence of every cell
## =============================================================================
library(reticulate)
use_virtualenv("r-celltypist", required = TRUE)          # the Python environment made above
celltypist <- import("celltypist")                       # CellTypist (Python)
anndata    <- import("anndata")                          # AnnData = CellTypist's input format

ct_models <- c("Adult_Human_Vascular.pkl", "Healthy_Adult_Heart.pkl", "Immune_All_Low.pkl")
ct_top_n  <- 20                                          # same 20 as the marker sheets
celltypist$models$download_models(model = as.list(ct_models))   # downloads the models (stored locally)

## object with clusters (chunk 8 save) and log-normalised RNA "data" layer
plot_tag   <- paste0(run_tag, "_res", leiden_resolution)                     # "PC20_NN250_res1"
seurat_out <- file.path(dir_seu, paste0("seurat_", input_tag, "_SeuratFunctions_", run_tag, ".RDS"))
if (!exists("obj") || !"umap_primary" %in% Reductions(obj)) obj <- readRDS(seurat_out)
obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)

## expression for CellTypist: cells x genes, log1p(counts / total x 10,000)
X <- tryCatch(LayerData(obj, assay = "RNA", layer = "data"),
              error = function(e) GetAssayData(obj, assay = "RNA", slot = "data"))   # genes x cells
stopifnot("data layer is not log1p to 10,000" = abs(sum(expm1(X[, 1])) - 1e4) < 1)  # CellTypist's own check
Xt <- Matrix::t(X)                                       # cells x genes (CellTypist orientation)
rm(X); invisible(gc())

branches_ct <- list(
  list(prefix = "primary", dir = dir_primary, label = "PRIMARY (batch-aware Pearson PCA)"),
  list(prefix = "harmony", dir = dir_sens,    label = "SENSITIVITY (batch-aware Pearson PCA + Harmony)"))

for (b in branches_ct) {
  main_col <- paste0(b$prefix, "_leiden_res_", leiden_resolution)          # e.g. "primary_leiden_res_1"
  dir_mk   <- file.path(b$dir, paste0("FindAllMarkers_res", leiden_resolution))
  file_tag <- paste0(b$prefix, "_", plot_tag)
  sig_file <- file.path(dir_mk, paste0("FindAllMarkers_", file_tag, "_all_significant_ranked.csv"))
  stopifnot("run chunk 10 (FindAllMarkers) first" = file.exists(sig_file))
  message("\n", b$label, ": CellTypist ...")
  
  ## top-20 markers per cluster from chunk 10 (same ranking: pct.1, then avg_log2FC)
  sig <- fread(sig_file)
  sig[, cluster := as.character(cluster)]
  
  ## AnnData: expression + this branch's Leiden cluster of each cell
  cl    <- as.character(obj[[main_col]][, 1])
  adata <- anndata$AnnData(X   = Xt,
                           obs = data.frame(leiden = cl, row.names = colnames(obj)),
                           var = data.frame(row.names = rownames(obj)))
  
  per_cell    <- data.table(cell_id = colnames(obj), section_id = obj[[batch_col]][, 1], leiden = cl)
  per_cluster <- list()
  
  for (m in ct_models) {
    mname <- sub("\\.pkl$", "", m)
    model <- celltypist$models$Model$load(model = m)
    n_feat <- length(intersect(rownames(obj), as.character(model$features)))   # genes CellTypist can use
    message("  ", mname, ": ", n_feat, " panel genes found in the model")
    
    ## cell-by-cell prediction, then majority vote within each Leiden cluster
    pred <- celltypist$annotate(adata, model = m, majority_voting = TRUE, over_clustering = "leiden")
    lab  <- pred$predicted_labels                        # per cell: predicted_labels, over_clustering, majority_voting
    prob <- pred$probability_matrix                      # per cell: probability of every cell type
    stopifnot(identical(rownames(lab), colnames(obj)))
    d <- data.table(leiden = cl,
                    pred   = as.character(lab$predicted_labels),   # label of the single cell
                    call   = as.character(lab$majority_voting),    # label of its cluster
                    conf   = do.call(pmax, prob))                  # highest probability of the cell
    per_cell[, paste0(mname, c("_cell_label", "_cell_confidence", "_cluster_call")) :=
               list(d$pred, round(d$conf, 4), d$call)]
    
    ## one row per cluster
    s <- d[, {
      tab <- sort(table(pred), decreasing = TRUE)        # labels of the cluster's cells, most common first
      .(n_cells             = .N,
        celltypist_call     = call[1],
        pct_cells_with_call = round(100 * mean(pred == call[1]), 1),   # how unanimous the cluster is
        median_confidence   = round(median(conf), 3),
        second_label        = if (length(tab) > 1) names(tab)[2] else NA_character_,
        pct_second_label    = if (length(tab) > 1) round(100 * tab[[2]] / .N, 1) else NA_real_)
    }, by = leiden]
    
    ## top-20 markers of the cluster vs the model's top-20 genes for the call
    mk_top  <- lapply(s$leiden, function(k) sig[cluster == k & rank_in_cluster <= ct_top_n, gene])
    mod_top <- lapply(s$celltypist_call, function(x) as.character(model$extract_top_markers(x, as.integer(ct_top_n))))
    ov      <- mapply(intersect, mk_top, mod_top, SIMPLIFY = FALSE)
    s[, `:=`(model                          = mname,
             n_panel_genes_in_model         = n_feat,
             n_overlap_top20                = lengths(ov),
             overlap_genes                  = sapply(ov, paste, collapse = ", "),
             cluster_top20_markers          = sapply(mk_top, paste, collapse = ", "),
             celltypist_top20_genes_of_call = sapply(mod_top, paste, collapse = ", "))]
    per_cluster[[mname]] <- s
  }
  
  ## save
  details <- rbindlist(per_cluster)[order(as.integer(leiden), model)]
  calls   <- dcast(details, leiden ~ model, value.var = "celltypist_call")[order(as.integer(leiden))]
  calls   <- merge(calls, details[, .(n_cells = n_cells[1]), by = leiden], by = "leiden")[order(as.integer(leiden))]
  fwrite(calls,    file.path(dir_mk, paste0("CellTypist_", file_tag, "_per_cluster_calls.csv")))
  fwrite(details,  file.path(dir_mk, paste0("CellTypist_", file_tag, "_per_cluster_details.csv")))
  fwrite(per_cell, file.path(dir_mk, paste0("CellTypist_", file_tag, "_per_cell.csv.gz")))
  print(calls)
  message(b$label, ": CellTypist results saved -> ", dir_mk)
}
























## =============================================================================
## _PearsonResiduals_vs_SCTransform.R
## =============================================================================
## PURPOSE
##   How similar is batch-aware Pearson-residual normalisation (the model behind
##   scPearsonPCA in script 12) to batch-aware SCTransform?
##   "Batch-aware" is done the SAME way for both: each model is fitted
##   separately within each tissue section (sample_id).
##     Pearson residuals : expected mu = totalcounts(cell) x gene frequency(gene, section)
##                         residual    = (count - mu) / sqrt(mu)
##     SCTransform v2    : regularised negative binomial (glmGamPoi), one model
##                         per section; returns Pearson residuals of that model
##   Same genes (script 12 HVGs), same cells, same post-processing for both:
##   no centring, residuals clipped at +/- sqrt(number of cells in the section).
##
## PER GENE, compared between the two methods (x = Pearson, y = SCTransform):
##   01 mean residual across cells           (the plot you asked for)
##   02 residual variance across cells       (what both methods use to rank genes)
##   03 correlation of the gene's residuals across cells (cell-by-cell agreement)
##   Groups: each tissue section, and All cells (6 sections pooled).
##
## NOTE ON PLOT 01: a Pearson residual is (observed - expected) / sd, so each
##   gene's MEAN residual is close to 0 by design in both methods. Plot 01 shows
##   how well each model's expected value fits each gene; plots 02 and 03 show
##   how similar the normalised values are.
##
## Input is only read. Outputs:
##   <script 12 folder>/14_PearsonResiduals_vs_SCTransform/
## Run chunk by chunk (Cmd+Enter). SCTransform runs once per section (minutes each).
## =============================================================================


## =============================================================================
## CHUNK 0: packages
## =============================================================================
## Needed once: install.packages("glmGamPoi") or BiocManager::install("glmGamPoi")
stopifnot("install glmGamPoi (used by SCTransform v2)" = requireNamespace("glmGamPoi", quietly = TRUE))
library(data.table)   # tables
library(Matrix)       # sparse counts
library(ggplot2)      # plots
library(patchwork)    # panels side by side
library(Seurat)       # SCTransform, LayerData


## =============================================================================
## CHUNK 1: paths and settings
## =============================================================================
output_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"
seurat_in  <- file.path(output_dir, "seurat_QC_filtered_per_section_03_to_09.RDS")   # same input as script 12
input_tag  <- sub("\\.RDS$", "", basename(seurat_in))
dir_main   <- file.path(output_dir, "12_Normalisation_PCA_UMAP_Leiden", input_tag)  # script 12 folder
dir_shared <- file.path(dir_main, "00_shared_HVGs_and_batch_aware_Pearson_PCA")
hvg_file   <- file.path(dir_shared, "hvgs.RDS")                   # script 12 HVGs
gf_file    <- file.path(dir_shared, "genefreq_batch.RDS")         # script 12 gene frequencies (check only)
pca_file   <- file.path(dir_shared, "pcaobj_batch_aware.RDS")     # script 12 PCA (check only)
stopifnot("Run script 12 first" = file.exists(hvg_file))

dir_out <- file.path(dir_main, "14_PearsonResiduals_vs_SCTransform")
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

batch_col  <- "sample_id"   # the 6 tissue sections = batch
min_cells  <- 5             # SCTransform v2 needs a gene in >= 5 cells of a section
chunk_size <- 500           # genes processed at once (memory)
n_val      <- 1000          # cells per section kept for the check in CHUNK 5
seed       <- 42
dpi_out    <- 300
fmt <- function(n) format(n, big.mark = ",")


## =============================================================================
## CHUNK 2: load the object; is there a mean-expression column already?
## =============================================================================
obj <- readRDS(seurat_in)

## Cell metadata = one row per CELL, so it cannot hold a per-GENE mean.
grep("mean", colnames(obj[[]]), ignore.case = TRUE, value = TRUE)   # cell-level columns with "mean"

## Gene metadata = one row per GENE (RNA assay). A mean is here only if
## FindVariableFeatures results were saved (vf_vst_counts_mean = mean RAW count,
## not normalised). Either way the normalised means are computed below.
colnames(obj[["RNA"]][[]])

## Counts (genes x cells), total counts, sections
counts_t    <- methods::as(LayerData(obj, assay = "RNA", layer = "counts"), "CsparseMatrix")
totalcounts <- Matrix::colSums(counts_t)                              # = nCount_RNA
message("Cells where total counts differ from nCount_RNA: ", sum(totalcounts != obj$nCount_RNA))
sec          <- obj[[batch_col]][, 1]
sections     <- sort(unique(sec))
cells_by_sec <- split(colnames(obj), sec)[sections]                   # cell names per section
all_cells    <- colnames(obj)                                         # original order (for CHUNK 5)
rm(obj); invisible(gc())                                              # object no longer needed

## Genes: script 12 HVGs detected in >= min_cells cells in EVERY section
hvgs <- readRDS(hvg_file)
det  <- sapply(sections, function(s) Matrix::rowSums(counts_t[hvgs, cells_by_sec[[s]]] > 0))
genes <- hvgs[apply(det >= min_cells, 1, all)]
message("Genes compared: ", length(genes), " of ", length(hvgs), " HVGs")


## =============================================================================
## CHUNK 3: batch-aware Pearson residuals (scPearsonPCA model, per section)
## =============================================================================
## Gene frequency per section = gene's counts / all counts in that section
gene_freq <- sapply(sections, function(s) {
  i <- cells_by_sec[[s]]
  Matrix::rowSums(counts_t[genes, i]) / sum(totalcounts[i])
})                                                                    # genes x sections

## Check against script 12's gene_frequency() output (should be ~0)
if (file.exists(gf_file)) {
  gf12 <- as.matrix(readRDS(gf_file))
  if (all(sections %in% colnames(gf12)))
    message("Max difference from script 12 gene frequencies: ",
            signif(max(abs(gf12[genes, sections] - gene_freq)), 3))
}

## Pearson residual of a block of genes x cells, clipped at +/- clip
pearson_resid <- function(X, tc, gf, clip) {
  mu <- outer(gf, tc)                       # expected count = gene frequency x cell total
  r  <- (X - mu) / sqrt(mu)                 # (observed - expected) / sd (Poisson)
  r[r >  clip] <-  clip
  r[r < -clip] <- -clip
  r
}


## =============================================================================
## CHUNK 4: per section - SCTransform, Pearson residuals, per-gene sums
## =============================================================================
## Only running sums are kept (not the full residual matrices), so memory stays low.
stats_list <- list(); val_P <- list(); val_cells <- character(0)
set.seed(seed)
for (s in sections) {
  i    <- cells_by_sec[[s]]
  n_s  <- length(i)
  clip <- sqrt(n_s)                                                   # same clip for both methods
  message("\n", s, ": SCTransform v2 on ", fmt(n_s), " cells ...")
  
  ## SCTransform on this section only (fresh object = per-section model)
  so <- CreateSeuratObject(counts = counts_t[, i])
  so <- SCTransform(so, vst.flavor = "v2", residual.features = genes,
                    do.center = FALSE, clip.range = c(-clip, clip),
                    seed.use = seed, verbose = FALSE)
  S_all <- GetAssayData(so, assay = "SCT", layer = "scale.data")      # SCT Pearson residuals
  stopifnot("SCTransform did not return all genes" = all(genes %in% rownames(S_all)),
            identical(colnames(S_all), i))
  rm(so); invisible(gc())
  
  v  <- sample(i, min(n_val, n_s)); val_cells <- c(val_cells, v)      # cells kept for CHUNK 5
  vp <- list()
  for (g in split(genes, ceiling(seq_along(genes) / chunk_size))) {
    X <- as.matrix(counts_t[g, i])                                    # raw counts
    P <- pearson_resid(X, totalcounts[i], gene_freq[g, s], clip)      # Pearson residuals
    S <- as.matrix(S_all[g, i])                                       # SCT residuals, same genes/cells
    stats_list[[length(stats_list) + 1]] <- data.table(
      section_id = s, gene = g, n = n_s,
      sum_count = rowSums(X),                                         # for mean raw count
      sum_p = rowSums(P), sum_s = rowSums(S),                         # for means
      sum_p2 = rowSums(P^2), sum_s2 = rowSums(S^2),                   # for variances
      sum_ps = rowSums(P * S))                                        # for correlation
    vp[[length(vp) + 1]] <- P[, v, drop = FALSE]
  }
  val_P[[s]] <- do.call(rbind, vp)
  rm(S_all, X, P, S); invisible(gc())
  message(s, ": done")
}


## =============================================================================
## CHUNK 5: check - do these Pearson residuals reproduce script 12's PC1?
## =============================================================================
## Scale residuals as scPearsonPCA does (centre, scale, cap at 10), project onto
## script 12's PC1 loadings, correlate with script 12's PC1. |r| near 1 = same residuals.
stats_sec <- rbindlist(stats_list)
pool      <- stats_sec[, lapply(.SD, sum), by = gene, .SDcols = c("n", "sum_p", "sum_p2")]
mean_p    <- setNames(pool$sum_p / pool$n, pool$gene)
sd_p      <- setNames(sqrt((pool$sum_p2 - pool$n * mean_p^2) / (pool$n - 1)), pool$gene)

check_pc1 <- NA_real_
if (file.exists(pca_file)) {
  rd   <- readRDS(pca_file)$reduction.data
  load <- rd@feature.loadings
  emb  <- rd@cell.embeddings
  if (is.null(rownames(emb))) rownames(emb) <- all_cells              # as in script 12
  if (nrow(load) > 0) {
    g_v <- intersect(genes, rownames(load))
    Pv  <- do.call(cbind, val_P)[g_v, ]                               # genes x checked cells
    Z   <- (Pv - mean_p[g_v]) / sd_p[g_v]                             # centre + scale per gene
    Z[Z > 10] <- 10                                                   # scale.max = 10
    check_pc1 <- cor(colSums(Z * load[g_v, 1]), emb[colnames(Pv), 1])
  }
  rm(rd); invisible(gc())
}
message("Correlation with script 12 PC1 (", length(val_cells), " cells): ", round(check_pc1, 4))


## =============================================================================
## CHUNK 6: per-gene mean, variance, correlation - each section and all cells
## =============================================================================
pooled <- stats_sec[, lapply(.SD, sum), by = gene,
                    .SDcols = c("n", "sum_count", "sum_p", "sum_s", "sum_p2", "sum_s2", "sum_ps")]
pooled[, section_id := "All cells"]
gs <- rbind(stats_sec, pooled, use.names = TRUE)
gs[, `:=`(mean_count = sum_count / n,
          mean_p = sum_p / n, mean_s = sum_s / n)]
gs[, `:=`(var_p  = (sum_p2 - n * mean_p^2) / (n - 1),
          var_s  = (sum_s2 - n * mean_s^2) / (n - 1),
          cor_ps = (sum_ps - n * mean_p * mean_s) /
            sqrt((sum_p2 - n * mean_p^2) * (sum_s2 - n * mean_s^2)))]
groups <- c(sections, "All cells")
gs[, section_id := factor(section_id, levels = groups)]

## One summary row per group
summ <- gs[, .(n_cells = n[1], n_genes = .N,
               r_mean        = round(cor(mean_p, mean_s), 3),
               rho_mean      = round(cor(mean_p, mean_s, method = "spearman"), 3),
               r_log10_var   = round(cor(log10(var_p), log10(var_s)), 3),
               rho_var       = round(cor(var_p, var_s, method = "spearman"), 3),
               median_gene_cor = round(median(cor_ps), 3),
               pct_genes_cor_above_0.9 = round(100 * mean(cor_ps > 0.9), 1)), by = section_id]
summ[, check_cor_with_script12_PC1 := round(check_pc1, 4)]
print(summ)

fwrite(gs[, .(section_id, gene, n_cells = n, mean_raw_count = mean_count, mean_pearson = mean_p,
              mean_sct = mean_s, var_pearson = var_p, var_sct = var_s, cor_pearson_vs_sct = cor_ps)],
       file.path(dir_out, "per_gene_Pearson_vs_SCTransform.csv"))
fwrite(summ, file.path(dir_out, "summary_Pearson_vs_SCTransform.csv"))


## =============================================================================
## CHUNK 7: plots - one folder per section + All cells, and overview grids
## =============================================================================
save_tiff <- function(p, path, w, h) {
  ggsave(path, p, width = w, height = h, dpi = dpi_out, device = "tiff", compression = "lzw", bg = "white")
  message("Saved: ", path)
}
x_lab <- "Batch-aware Pearson residuals (scPearsonPCA model)"
y_lab <- "Batch-aware SCTransform v2"
col_scale <- scale_colour_viridis_c(name = "log10 mean\nraw count")   # colour = how highly expressed

## 01 mean residual per gene
plot_mean <- function(d, ttl, sub) {
  ggplot(d, aes(mean_p, mean_s, colour = log10(mean_count))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +   # y = x
    geom_point(size = 0.8, alpha = 0.7) + col_scale + theme_bw(base_size = 12) +
    labs(title = ttl, subtitle = sub, x = paste0("Mean residual per gene\n", x_lab),
         y = paste0("Mean residual per gene\n", y_lab))
}
## 02 residual variance per gene (log scales)
plot_var <- function(d, ttl, sub) {
  ggplot(d, aes(var_p, var_s, colour = log10(mean_count))) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey40") +
    geom_point(size = 0.8, alpha = 0.7) + scale_x_log10() + scale_y_log10() + col_scale +
    theme_bw(base_size = 12) +
    labs(title = ttl, subtitle = sub, x = paste0("Residual variance per gene\n", x_lab),
         y = paste0("Residual variance per gene\n", y_lab))
}
## 03 per-gene correlation of residuals across cells
plot_cor <- function(d, ttl, sub) {
  ggplot(d, aes(cor_ps)) +
    geom_histogram(bins = 50, fill = "steelblue") +
    geom_vline(xintercept = median(d$cor_ps), colour = "red", linetype = "dashed") +  # median
    theme_bw(base_size = 12) +
    labs(title = ttl, subtitle = sub, y = "Number of genes",
         x = "Correlation across cells of each gene's residuals\n(Pearson residuals vs SCTransform)")
}

## One folder per group, 3 plots + the 3 side by side
for (grp in groups) {
  d  <- gs[section_id == grp]; r <- summ[section_id == grp]
  od <- file.path(dir_out, gsub(" ", "_", grp)); dir.create(od, showWarnings = FALSE)
  tag <- paste0(grp, " (", fmt(r$n_cells), " cells, ", r$n_genes, " genes)")
  p1 <- plot_mean(d, paste0(tag, ": mean residual per gene"),
                  paste0("Pearson r = ", r$r_mean, ", Spearman rho = ", r$rho_mean))
  p2 <- plot_var(d, paste0(tag, ": residual variance per gene"),
                 paste0("Pearson r (log10) = ", r$r_log10_var, ", Spearman rho = ", r$rho_var))
  p3 <- plot_cor(d, paste0(tag, ": agreement per gene across cells"),
                 paste0("Median r = ", r$median_gene_cor, "; ", r$pct_genes_cor_above_0.9, "% of genes r > 0.9"))
  save_tiff(p1, file.path(od, paste0("01_mean_residual_Pearson_vs_SCT_", gsub(" ", "_", grp), ".tiff")), 8, 7)
  save_tiff(p2, file.path(od, paste0("02_residual_variance_Pearson_vs_SCT_", gsub(" ", "_", grp), ".tiff")), 8, 7)
  save_tiff(p3, file.path(od, paste0("03_per_gene_correlation_Pearson_vs_SCT_", gsub(" ", "_", grp), ".tiff")), 8, 6)
  save_tiff(p1 | p2 | p3, file.path(od, paste0("04_all_three_side_by_side_", gsub(" ", "_", grp), ".tiff")), 24, 7)
}

## Overview: all 6 sections + All cells in one figure per plot type
lab_mean <- summ[, .(section_id, label = paste0("r = ", r_mean, ", rho = ", rho_mean))]
lab_var  <- summ[, .(section_id, label = paste0("r(log10) = ", r_log10_var))]
lab_cor  <- summ[, .(section_id, label = paste0("median r = ", median_gene_cor))]
add_lab  <- function(lab) geom_text(data = lab, aes(x = -Inf, y = Inf, label = label), inherit.aes = FALSE,
                                    hjust = -0.05, vjust = 1.5, size = 3.2)
facet <- facet_wrap(~ section_id, ncol = 4, scales = "free")

save_tiff(plot_mean(gs, "Mean residual per gene: Pearson residuals vs SCTransform (both batch-aware)", NULL) +
            facet + add_lab(lab_mean),
          file.path(dir_out, "01_overview_mean_residual_all_groups.tiff"), 18, 9)
save_tiff(plot_var(gs, "Residual variance per gene: Pearson residuals vs SCTransform (both batch-aware)", NULL) +
            facet + add_lab(lab_var),
          file.path(dir_out, "02_overview_residual_variance_all_groups.tiff"), 18, 9)
save_tiff(ggplot(gs, aes(cor_ps)) + geom_histogram(bins = 50, fill = "steelblue") + facet +
            add_lab(lab_cor) + theme_bw(base_size = 12) +
            labs(title = "Per-gene correlation across cells: Pearson residuals vs SCTransform",
                 x = "Correlation of each gene's residuals across cells", y = "Number of genes"),
          file.path(dir_out, "03_overview_per_gene_correlation_all_groups.tiff"), 18, 9)

message("\nAll outputs in: ", dir_out)