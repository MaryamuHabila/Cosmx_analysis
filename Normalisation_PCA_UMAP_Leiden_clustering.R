## =============================================================================
## 12b_Seurat_functions_PCA_UMAP_Leiden.R
## =============================================================================
## PURPOSE
##   The SAME analysis as 12_Normalisation_PCA_UMAP_Leiden.R - same steps, same
##   parameters, same tables and plots - but each step uses the SEURAT function
##   wherever Seurat has one. Nothing is added that script 12 does not do.
##
##   step (script 12 chunk)        script 12                          12b (this script)
##   ----------------------------  ---------------------------------  ------------------------------------------
##   HVGs (3)                      Seurat::FindVariableFeatures       same (already Seurat)
##   batch-aware Pearson PCA (4)   scPearsonPCA                       same - Seurat has no batch-aware Pearson PCA;
##                                                                    stored as Seurat reduction "pearson_pca"
##   elbow plot (4)                ggplot (% variance)                Seurat::ElbowPlot (y axis = SD of each PC);
##                                                                    % variance still saved in the CSV
##   STEP A kNN, cosine, annoy     RcppAnnoy loop                     Seurat::FindNeighbors(nn.method = "annoy",
##                                                                    annoy.metric = "cosine", return.neighbor = TRUE)
##   STEP B neighbour graph        uwot::similarity_graph             same - Seurat has no function for this graph
##                                                                    (its own graph is SNN, a different graph);
##                                                                    stored as a Seurat Graph
##   STEP C UMAP of that graph     uwot::optimize_graph_layout        Seurat::RunUMAP(nn.name = STEP A kNN,
##                                                                    umap.method = "uwot") - uwot rebuilds the
##                                                                    same graph from the same kNN, then lays it out
##   STEP D Leiden on that graph   igraph::cluster_leiden             Seurat::FindClusters(graph.name = STEP B graph,
##                                                                    algorithm = 4) - never on UMAP coordinates
##   Harmony (7)                   harmony::RunHarmony(PC matrix)     harmony::RunHarmony(Seurat object) - same call
##   UMAP / PCA / spatial plots    ggplot                             DimPlot, FeaturePlot, ImageDimPlot
##   composition bars, sweep       ggplot                             same ggplot - Seurat has no function for these
##   summary, cluster tree,
##   crosstab heatmap
##
## EXPECTED DIFFERENCE FROM SCRIPT 12 (implementation, not method)
##   - kNN: annoy is approximate; Seurat uses annoy's default search depth,
##     script 12 used 2 x k x n_trees (FindNeighbors does not expose it).
##   - Leiden: Seurat runs it through leidenbase, script 12 through igraph.
##   So clusters/UMAP are near-identical to script 12, not bit-identical.
##
## OUTPUT: <script 12 folder>/seurat_functions_used/PC<n>_NN<k>/
##   00_shared/ , 01_PRIMARY_batchawarePCA/ , 02_SENSITIVITY_batchawarePCA_Harmony/ ,
##   03_PRIMARY_vs_SENSITIVITY/  - same file names as script 12 inside each.
##   Script 12's own results are never overwritten.
##
## Script 13 (side-by-side UMAPs) and FindAllMarkers are NOT part of script 12;
## they are in 13b_Seurat_UMAP_side_by_side_and_FindAllMarkers.R, run after this.
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 0: packages - checked, NOT silently installed
## =============================================================================
required_pkgs <- c("data.table", "Matrix", "ggplot2", "patchwork", "scales", "Seurat", "SeuratObject",
                   "scPearsonPCA", "RSpectra", "uwot", "igraph", "harmony", "leidenbase")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall with:\n",
       "  install.packages(c('data.table','Matrix','ggplot2','patchwork','scales','Seurat','RSpectra',\n",
       "                     'uwot','igraph','harmony','leidenbase'))\n",
       "  remotes::install_github('Nanostring-Biostats/CosMx-Analysis-Scratch-Space',\n",
       "                          subdir = '_code/scPearsonPCA', ref = 'Main')")
}
## functions this script depends on - stop now rather than fail halfway
if (utils::packageVersion("Seurat") < "5.1.0") {                 # FindClusters(cluster.name), Leiden via leidenbase
  stop("Seurat >= 5.1.0 needed. install.packages('Seurat')")
}
if (!"similarity_graph" %in% getNamespaceExports("uwot")) {       # STEP B
  stop("Your uwot version lacks similarity_graph(). install.packages('uwot')")
}
if (utils::packageVersion("harmony") < "1.0") {                   # RunHarmony Seurat method
  stop("harmony >= 1.0 needed. install.packages('harmony')")
}
if (utils::packageVersion("ggplot2") < "3.4.0") {                 # cluster tree uses 'linewidth' (as script 12)
  stop("ggplot2 >= 3.4.0 needed for the resolution-sweep cluster tree. install.packages('ggplot2')")
}

library(data.table)   # tables
library(Matrix)       # sparse matrices
library(ggplot2)      # plots Seurat has no function for (+ extra layers on Seurat plots)
library(patchwork)    # wrap_plots(): sweep UMAP panels
library(Seurat)       # FindVariableFeatures, FindNeighbors, RunUMAP, FindClusters, DimPlot, ...
library(scPearsonPCA) # batch-aware Pearson-residual PCA (vignette)


## =============================================================================
## CHUNK 1: paths + parameters (same values as script 12)
## =============================================================================
output_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"
seurat_in  <- file.path(output_dir, "seurat_QC_filtered_per_section_03_to_09.RDS")   # script 11 (read only)
input_tag  <- sub("\\.RDS$", "", basename(seurat_in))    # folder name follows the input object

batch_col           <- "sample_id"    # the 6 AAA tissue sections (AAA2 ... AAA48)
n_hvgs              <- 2000           # vignette: 2000 for a 6k panel
n_pcs               <- 20             # vignette: sc.pp.neighbors(n_pcs = 40)                              <- change here
n_neighbors         <- 250            # vignette: sc.pp.neighbors(n_neighbors = 15), incl. the cell itself  <- change here
rerun_pca           <- FALSE          # FALSE = reuse saved HVGs + PCA; TRUE = recompute

dir_main    <- file.path(output_dir, "12_Normalisation_PCA_UMAP_Leiden", input_tag)   # script 12 folder
dir_shared12 <- file.path(dir_main, "00_shared_HVGs_and_batch_aware_Pearson_PCA")    # script 12 HVGs + PCA (read only)
run_tag     <- paste0("PC", n_pcs, "_NN", n_neighbors)                                # e.g. PC20_NN250
dir_seu     <- file.path(dir_main, "seurat_functions_used", run_tag)                  # all 12b outputs here
dir_shared  <- file.path(dir_seu, "00_shared")
dir_primary <- file.path(dir_seu, "01_PRIMARY_batchawarePCA")
dir_sens    <- file.path(dir_seu, "02_SENSITIVITY_batchawarePCA_Harmony")
dir_compare <- file.path(dir_seu, "03_PRIMARY_vs_SENSITIVITY")
for (d in c(dir_shared, dir_primary, dir_sens, dir_compare)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

annoy_n_trees       <- 50             # annoy index size (as script 12)
umap_min_dist       <- 0.02           # vignette: sc.tl.umap(min_dist = 0.02)
umap_spread         <- 2              # vignette: sc.tl.umap(spread = 2)
leiden_resolution   <- 1              # vignette: sc.tl.leiden(resolution = 0.5)                  <- change here
leiden_n_iterations <- 2              # vignette: sc.tl.leiden(n_iterations = 2)
resolution_sweep    <- c(0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0)  # compared on the same graph <- change here
dominant_section_flag_pct <- 90       # flag (for investigation only) clusters with > 90% of cells from one section
sweep_umap_max_cells   <- 200000      # PLOTTING ONLY: cells drawn per panel in the sweep UMAP grid
sweep_tree_min_in_prop <- 0.05        # PLOTTING ONLY: cluster-tree edges drawn if >= 5% of the child cluster
seed    <- 42
dpi_out <- 300
stopifnot(leiden_resolution %in% resolution_sweep)   # main clusters are taken from the sweep

## parameters table - same rows as script 12, values name the Seurat functions
params <- data.table(
  parameter = c("input_metadata", "input_counts", "batch_variable", "n_HVGs", "HVG_method",
                "PCA", "PCA_scale.max", "PCA_do.scale", "PCA_do.center", "totalcounts",
                "n_PCs_for_neighbours", "n_neighbors_incl_self", "neighbour_metric", "kNN_method",
                "graph", "UMAP", "UMAP_min_dist", "UMAP_spread", "UMAP_init",
                "clustering", "Leiden_objective", "Leiden_resolution", "Leiden_n_iterations",
                "Harmony (sensitivity only)", "seed",
                "Leiden_resolution_sweep", "sweep_UMAP_max_cells_per_panel (plot only)",
                "sweep_tree_min_in_prop (plot only)"),
  value = c(basename(seurat_in), "RNA counts layer of the input Seurat object",
            batch_col, n_hvgs, "Seurat::FindVariableFeatures (vst)",
            "scPearsonPCA::sparse_quasipoisson_pca_seurat_batch (batch-aware Pearson residuals)",
            10, TRUE, TRUE, "per-cell sum over all genes (= nCount_RNA)",
            n_pcs, n_neighbors, "cosine",
            paste0("Seurat::FindNeighbors(nn.method = 'annoy', annoy.metric = 'cosine', n.trees = ", annoy_n_trees, ")"),
            "uwot::similarity_graph on the kNN (UMAP fuzzy connectivities, as sc.pp.neighbors), as Seurat Graph",
            "Seurat::RunUMAP(nn.name = kNN, umap.method = 'uwot') - layout of that SAME graph",
            umap_min_dist, umap_spread, "spectral",
            "Seurat::FindClusters(algorithm = 4, Leiden) on that SAME graph (never on UMAP coordinates)",
            "modularity (RBConfiguration)", leiden_resolution, leiden_n_iterations,
            paste0("harmony::RunHarmony (Seurat object) on all batch-aware PCs, group.by.vars = ", batch_col), seed,
            paste(resolution_sweep, collapse = ", "), sweep_umap_max_cells, sweep_tree_min_in_prop))
fwrite(params, file.path(dir_seu, paste0("parameters_used_", run_tag, ".csv")))


## =============================================================================
## CHUNK 2: load the QC-filtered Seurat object (script 11) + checks
## =============================================================================
stopifnot("input Seurat object not found - run script 11 first" = file.exists(seurat_in))
obj <- readRDS(seurat_in)
obj

stopifnot("cell_id column must match the object's cells" = identical(obj$cell_id, colnames(obj)),
          all(c(batch_col, "x_slide_mm", "y_slide_mm", "fov") %in% colnames(obj[[]])))

counts_t <- LayerData(obj, assay = "RNA", layer = "counts")   # genes x cells (as the vignette uses)

## per-cell total counts from the matrix (all genes) = nCount_RNA
totalcounts <- Matrix::colSums(counts_t)
if (any(totalcounts == 0)) stop(sum(totalcounts == 0), " cells have zero counts.")
message("Cells where matrix total differs from nCount_RNA: ", sum(totalcounts != obj$nCount_RNA))

## genes with zero counts in every kept cell carry no information
zero_genes <- rownames(counts_t)[Matrix::rowSums(counts_t) == 0]
if (length(zero_genes) > 0) message(length(zero_genes), " genes have zero counts in all kept cells (not HVGs).")

sec      <- obj[[batch_col]][, 1]                                  # section of each cell
sections <- sort(unique(sec))
message("Loaded ", format(ncol(obj), big.mark = ","), " cells x ", nrow(obj), " genes; ",
        length(sections), " sections: ", paste(sections, collapse = ", "))
print(table(sec))

## slide coordinates (mm) as a Seurat FOV, so ImageDimPlot draws the SAME
## coordinates as all earlier plots (x_slide_mm / y_slide_mm)
obj[["slide_mm"]] <- CreateFOV(coords = data.frame(x = obj$x_slide_mm, y = obj$y_slide_mm, cell = colnames(obj)),
                               type = "centroids", assay = "RNA")
section_label_pos <- data.table(x = obj$x_slide_mm, y = obj$y_slide_mm, section = sec)[
  , .(x = mean(x), y = mean(y)), by = section]                     # label position only


## =============================================================================
## CHUNK 3: HVGs (shared) - Seurat FindVariableFeatures, as in the vignette
## =============================================================================
## Reused if already saved: this script's 00_shared first, else script 12's
## (identical computation), unless rerun_pca = TRUE.
hvg_file <- file.path(dir_shared, "hvgs.RDS")
hvg_src  <- c(hvg_file, file.path(dir_shared12, "hvgs.RDS"))
hvg_src  <- hvg_src[file.exists(hvg_src)]
if (length(hvg_src) > 0 && !rerun_pca) {
  hvgs <- readRDS(hvg_src[1])
  message("HVGs reused from ", hvg_src[1], ": ", length(hvgs))
} else {
  obj  <- FindVariableFeatures(obj, nfeatures = n_hvgs)            # vst on the counts layer
  hvgs <- VariableFeatures(obj)
  message("\nHVGs selected: ", length(hvgs))
  
  ## remove HVGs with ZERO counts in a whole section (per-section frequency 0 ->
  ## Pearson residual undefined; scPearsonPCA's own check misses sparse zeros)
  det_by_sec <- sapply(sections, function(s) Matrix::rowSums(counts_t[hvgs, sec == s, drop = FALSE]) > 0)
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
VariableFeatures(obj) <- hvgs                                      # HVGs recorded in the object


## =============================================================================
## CHUNK 4: batch-aware Pearson-residual PCA (shared) - the vignette's call
## =============================================================================
## No Seurat equivalent (SCTransform uses a different model), so scPearsonPCA
## is kept; its result is stored as the Seurat reduction "pearson_pca".
obs <- data.table(cell_id = colnames(obj), sample_id = sec)
pca_file <- file.path(dir_shared, "pcaobj_batch_aware.RDS")
pca_src  <- c(pca_file, file.path(dir_shared12, "pcaobj_batch_aware.RDS"))
pca_src  <- pca_src[file.exists(pca_src)]
if (length(pca_src) > 0 && !rerun_pca) {
  pcaobj <- readRDS(pca_src[1])
  message("Batch-aware PCA reused from ", pca_src[1])
} else {
  message("\nComputing per-section gene frequencies...")
  genefreq_batch <- scPearsonPCA::gene_frequency(x = counts_t, obs = obs,
                                                 cellid_colname = "cell_id", batch_variable = "sample_id")
  stopifnot("gene_frequency output missing HVGs" = all(hvgs %in% rownames(genefreq_batch)))
  gf_hvg <- genefreq_batch[hvgs, , drop = FALSE]                   # per-section frequency of each HVG
  if (any(!is.finite(as.matrix(gf_hvg))) || any(as.matrix(gf_hvg) <= 0)) {
    stop("Non-positive or non-finite gene frequencies for some HVGs - inspect genefreq_batch.")
  }
  saveRDS(genefreq_batch, file.path(dir_shared, "genefreq_batch.RDS"))
  
  message("Running batch-aware Pearson-residual PCA on ", length(hvgs), " HVGs...")
  set.seed(seed)
  pcaobj <- sparse_quasipoisson_pca_seurat_batch(
    counts_t[hvgs, ], totalcounts = totalcounts, grate = gf_hvg, obs = obs,
    batch_variable = "sample_id", cellid_colname = "cell_id",
    scale.max = 10, do.scale = TRUE, do.center = TRUE)
  saveRDS(pcaobj, pca_file)
}
rm(counts_t); invisible(gc())                                      # raw matrix no longer needed

if (is.null(pcaobj$reduction.data)) stop("pcaobj has no $reduction.data - fields: ", paste(names(pcaobj), collapse = ", "))
emb_pca <- pcaobj$reduction.data@cell.embeddings
if (is.null(rownames(emb_pca))) { stopifnot(nrow(emb_pca) == ncol(obj)); rownames(emb_pca) <- colnames(obj) }
emb_pca <- emb_pca[colnames(obj), , drop = FALSE]                   # same cell order as the object
stopifnot("PCA embeddings contain NaN/Inf" = all(is.finite(emb_pca)))
if (ncol(emb_pca) < n_pcs) stop("PCA returned ", ncol(emb_pca), " PCs, fewer than n_pcs = ", n_pcs)
message("PCA embeddings: ", nrow(emb_pca), " cells x ", ncol(emb_pca), " PCs")

## PCA -> Seurat reduction (column names must be key + number: PPC_1, PPC_2, ...)
with_key <- function(m, key) { m <- as.matrix(m); colnames(m) <- paste0(key, seq_len(ncol(m))); m }
pc_sd <- pcaobj$reduction.data@stdev                               # SD per PC (ElbowPlot needs it)
obj[["pearson_pca"]] <- CreateDimReducObject(embeddings = with_key(emb_pca, "PPC_"), stdev = pc_sd,
                                             key = "PPC_", assay = "RNA")
rm(pcaobj, emb_pca); invisible(gc())

## elbow plot (all PCs computed) with the n_pcs cut-off - Seurat ElbowPlot
elbow <- data.table(PC = seq_along(pc_sd), variance = pc_sd^2)
elbow[, pct_variance_of_computed_PCs := 100 * variance / sum(variance)]
fwrite(elbow, file.path(dir_shared, "PCA_variance_per_PC.csv"))   # % variance table, as script 12
p_elbow <- ElbowPlot(obj, reduction = "pearson_pca", ndims = length(pc_sd)) +
  geom_vline(xintercept = n_pcs, colour = "red", linetype = "dashed") +
  labs(title = "Batch-aware Pearson PCA - elbow plot",
       subtitle = paste0("Red line = ", n_pcs, " PCs used for the neighbour graph"))
ggsave(file.path(dir_shared, paste0("PCA_elbow_", run_tag, ".tiff")), p_elbow, width = 8, height = 5.5,
       dpi = dpi_out, device = "tiff", compression = "lzw", bg = "white")


## =============================================================================
## CHUNK 5: helper functions - colours, saving, sweep plots, one analysis branch
## =============================================================================

## VERY DISTINCT CLUSTER COLOURS - same palette and order as script 12
distinct_palette <- c(
  "#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#42D4F4", "#F032E6",
  "#9A6324", "#469990", "#BFEF45", "#800000", "#FFE119", "#000075", "#808000",
  "#FF8FB1", "#A9A9A9", "#005300", "#009FFF", "#C77CFF", "#FF6F61", "#00FFBE",
  "#D4A017", "#7F7FFF", "#8B008B", "#2F4F4F", "#E9967A", "#556B2F", "#FF1493",
  "#20B2AA", "#B8860B", "#6A5ACD", "#CD853F", "#4B0082", "#00CED1", "#8FBC8F",
  "#A0522D")
stopifnot("duplicate colour in distinct_palette" = !anyDuplicated(toupper(distinct_palette)))

cluster_colours <- function(lv) {                                  # as script 12
  n <- length(lv)
  if (n <= length(distinct_palette)) {
    cols <- distinct_palette[seq_len(n)]
  } else if (requireNamespace("Polychrome", quietly = TRUE)) {
    set.seed(seed)
    cols <- unname(Polychrome::createPalette(n, seedcolors = distinct_palette, range = c(25, 80)))
  } else {
    warning(n, " clusters but only ", length(distinct_palette), " hand-picked colours; ",
            "extra colours are less distinct. install.packages('Polychrome') for better ones.")
    k <- n - length(distinct_palette)
    extra <- grDevices::hcl(h = seq(0, 360, length.out = k + 1)[-1], c = 90,
                            l = rep(c(40, 60, 75), length.out = k))
    cols <- c(distinct_palette, extra)
  }
  setNames(cols, lv)
}
section_colours <- setNames(scales::hue_pal()(length(sections)), sections)

save_tiff <- function(p, path, w, h) {
  ggsave(path, p, width = w, height = h, dpi = dpi_out, device = "tiff", compression = "lzw", bg = "white")
}

## Leiden labels 0,1,2,... by decreasing size (scanpy / script 12 convention)
relabel_by_size <- function(f) {
  sizes <- sort(table(f), decreasing = TRUE)
  new   <- setNames(seq_along(sizes) - 1L, names(sizes))
  factor(unname(new[as.character(f)]), levels = 0:(length(sizes) - 1L))
}

## RESOLUTION SWEEP VISUALS - plotting only, nothing is re-clustered here.
## 1. summary curves (ggplot, as script 12 - no Seurat function)
## 2. UMAP panels, one per resolution (Seurat DimPlot)
## 3. cluster tree + nodes/edges tables (ggplot, as script 12 - no Seurat function)
plot_resolution_sweep <- function(obj, sweep, cl_cols, umap_name, branch_dir, branch_name) {
  res_ord <- sort(resolution_sweep)
  col_of  <- setNames(cl_cols, resolution_sweep)                    # cluster column of each resolution
  get_cl  <- function(r) obj[[col_of[[as.character(r)]]]][, 1]
  panel_lab <- function(r) paste0("res ", r, "  (", nlevels(get_cl(r)), " clusters)",
                                  ifelse(r == leiden_resolution, "  [main]", ""))
  
  ## ---- 1. summary curves ---------------------------------------------------
  sw_long <- rbind(
    sweep[, .(resolution, metric = "Number of clusters", value = n_clusters)],
    sweep[, .(resolution, metric = "Clusters with < 50 cells", value = n_clusters_lt_50_cells)],
    sweep[, .(resolution, metric = "Smallest cluster (cells)", value = smallest_cluster)],
    sweep[, .(resolution, metric = paste0("ARI vs main resolution (", leiden_resolution, ")"),
              value = ARI_vs_resolution_main)])
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
  
  ## ---- 2. UMAP panels (DimPlot), one per resolution --------------------------
  ## Same UMAP in every panel; cluster numbers are by size WITHIN each resolution.
  ## Random subsample drawn only if > sweep_umap_max_cells (labels then sit at the
  ## median of the drawn cells; script 12 used all cells for label positions).
  max_k <- max(sweep$n_clusters)
  cols_sw <- cluster_colours(as.character(0:(max_k - 1L)))
  set.seed(seed)
  cells_draw <- if (ncol(obj) > sweep_umap_max_cells) sample(colnames(obj), sweep_umap_max_cells) else colnames(obj)
  panels <- lapply(res_ord, function(r) {
    DimPlot(obj, reduction = umap_name, group.by = col_of[[as.character(r)]], cells = cells_draw,
            cols = cols_sw, label = TRUE, label.size = 2.3, shuffle = TRUE, seed = seed,
            pt.size = 0.03, alpha = 0.3, raster = FALSE) +
      NoLegend() + coord_fixed() + labs(title = panel_lab(r))
  })
  p_um <- wrap_plots(panels, ncol = 4) +
    plot_annotation(title = paste0(branch_name, ": UMAP coloured by Leiden cluster at each sweep resolution"),
                    subtitle = paste0("Same UMAP in every panel; ",
                                      if (length(cells_draw) < ncol(obj)) paste0(format(length(cells_draw), big.mark = ","),
                                                                                 " randomly sampled cells drawn per panel; ") else "",
                                      "cluster numbers are by size within each resolution"))
  save_tiff(p_um, file.path(branch_dir, "Leiden_resolution_sweep_UMAP_panels.tiff"), 20, 11)
  
  ## ---- 3. cluster tree: how clusters split between consecutive resolutions --
  ## nodes = clusters at each resolution; edges = cells shared between a cluster
  ## and a cluster at the next resolution. in_prop = fraction of the child's cells
  ## that came from that parent. (exactly script 12's code)
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
    setorder(nd, order_key)                                         # children placed under their parents
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
    scale_y_continuous(breaks = -seq_along(res_ord), labels = res_ord, expand = expansion(add = 0.5)) +
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

## ONE ANALYSIS BRANCH (script 12's run_branch, with Seurat functions)
## reduction = "pearson_pca" (PRIMARY) or "harmony" (SENSITIVITY)
## prefix    = short name for the kNN, graph, UMAP and cluster columns in obj
## returns list(obj = updated object, results = same list as script 12)
run_branch_seurat <- function(obj, reduction, prefix, branch_dir, branch_name, emb_label) {
  message("\n", strrep("=", 78), "\n", branch_name, "\n", strrep("=", 78))
  nn_name    <- paste0(prefix, "_neighbors")                        # Seurat Neighbor object (STEP A kNN)
  graph_name <- paste0(prefix, "_connectivities")                   # Seurat Graph (STEP B)
  umap_name  <- paste0("umap_", prefix)                             # UMAP reduction (STEP C)
  cl_cols    <- paste0(prefix, "_leiden_res_", resolution_sweep)    # one cluster column per resolution (STEP D)
  main_col   <- paste0(prefix, "_leiden_res_", leiden_resolution)
  
  ## PC1 vs PC2 of the input embedding, by section (DimPlot)
  save_tiff(DimPlot(obj, reduction = reduction, dims = 1:2, group.by = batch_col, cols = section_colours,
                    shuffle = TRUE, seed = seed, pt.size = 0.05, alpha = 0.3, raster = FALSE) +
              coord_fixed() + theme_bw(base_size = 13) +
              guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
              labs(title = paste0(emb_label, ": dimension 1 vs 2, coloured by section"), colour = "Section"),
            file.path(branch_dir, "PC1_vs_PC2_by_section.tiff"), 9, 8)
  
  ## ---- STEP A: kNN on PC space (FindNeighbors) -----------------------------
  message("STEP A: FindNeighbors on ", n_pcs, " dims of ", emb_label, " (cosine, k = ", n_neighbors, " incl. self)...")
  obj <- FindNeighbors(obj, reduction = reduction, dims = 1:n_pcs, k.param = n_neighbors,
                       nn.method = "annoy", annoy.metric = "cosine", n.trees = annoy_n_trees,
                       return.neighbor = TRUE, graph.name = nn_name)
  idx <- Indices(obj[[nn_name]])                                    # cells x k: row numbers of neighbours
  n   <- nrow(idx); k <- ncol(idx)
  
  ## the cell itself in column 1 (as script 12 / sc.pp.neighbors); rare rows fixed
  bad <- which(idx[, 1] != seq_len(n))
  message("  cells without themselves as neighbour 1 (fixed, as script 12): ", length(bad))
  for (i in bad) idx[i, ] <- c(i, setdiff(idx[i, ], i))[seq_len(k)]
  stopifnot(all(idx[, 1] == seq_len(n)))
  
  ## cosine distance (1 - cos) to each neighbour = script 12's annoy angular^2 / 2;
  ## computed from the same PCs so it does not depend on how Seurat reports annoy distances
  X   <- Embeddings(obj, reduction)[, 1:n_pcs]
  Xn  <- X / sqrt(rowSums(X^2))                                     # unit-length rows
  dst <- vapply(seq_len(k), function(j) 1 - rowSums(Xn * Xn[idx[, j], , drop = FALSE]), numeric(n))
  dst[dst < 0] <- 0; dst[, 1] <- 0                                  # rounding; self = 0
  rm(Xn); invisible(gc())
  
  nb <- obj[[nn_name]]; nb@nn.idx <- idx; nb@nn.dist <- dst; obj[[nn_name]] <- nb   # back into the Neighbor object
  nn <- list(idx = idx, dist = dst)
  saveRDS(nn, file.path(branch_dir, "STEP_A_kNN_idx_dist.RDS"))     # same file as script 12
  
  ## ---- STEP B: ONE neighbour graph from the kNN -----------------------------
  message("STEP B: building the neighbour graph (UMAP fuzzy connectivities) from the kNN...")
  G <- uwot::similarity_graph(X = NULL, nn_method = nn, set_op_mix_ratio = 1,
                              local_connectivity = 1, verbose = TRUE)
  G <- methods::as(G, "CsparseMatrix")
  stopifnot(nrow(G) == n, ncol(G) == n)
  stopifnot("neighbour graph is not symmetric" = Matrix::isSymmetric(G, tol = 1e-6))
  dimnames(G) <- list(colnames(obj), colnames(obj))
  saveRDS(G, file.path(branch_dir, "STEP_B_neighbour_graph_connectivities.RDS"))
  message("  graph: ", format(n, big.mark = ","), " cells, ", format(Matrix::nnzero(G) / 2, big.mark = ","), " edges")
  g_obj <- as.Graph(G); DefaultAssay(g_obj) <- "RNA"                # Seurat Graph object
  obj[[graph_name]] <- g_obj
  rm(g_obj); invisible(gc())
  
  ## ---- STEP C: UMAP = layout of graph G (RunUMAP) ---------------------------
  ## RunUMAP(nn.name = ...) gives uwot the STEP A kNN; uwot builds the same
  ## fuzzy graph as STEP B (set_op_mix_ratio = 1, local_connectivity = 1) and
  ## lays it out (spectral init, min_dist, spread as script 12)
  message("STEP C: UMAP layout of the neighbour graph (RunUMAP)...")
  obj <- RunUMAP(obj, nn.name = nn_name, umap.method = "uwot", min.dist = umap_min_dist,
                 spread = umap_spread, seed.use = seed, reduction.name = umap_name,
                 reduction.key = paste0("UMAP", toupper(substr(prefix, 1, 1)), "_"))   # UMAPP_ / UMAPH_
  umap_xy <- Embeddings(obj, umap_name)
  colnames(umap_xy) <- c("UMAP_1", "UMAP_2")
  
  ## ---- STEP D: Leiden on graph G (FindClusters), NOT on UMAP ----------------
  message("STEP D: Leiden clustering of the neighbour graph (FindClusters)...")
  g <- igraph::graph_from_adjacency_matrix(G, mode = "undirected", weighted = TRUE, diag = FALSE)
  message("  connected components in graph: ", igraph::components(g)$no)
  rm(g); invisible(gc())
  
  ## every resolution on the SAME graph (main resolution is one of them)
  for (i in seq_along(resolution_sweep)) {
    obj <- FindClusters(obj, graph.name = graph_name, algorithm = 4, resolution = resolution_sweep[i],
                        n.iter = leiden_n_iterations, random.seed = seed, group.singletons = FALSE,
                        cluster.name = cl_cols[i], verbose = FALSE)
    obj[[cl_cols[i]]] <- relabel_by_size(obj[[cl_cols[i]]][, 1])    # 0 = largest (script 12 convention)
  }
  obj$seurat_clusters <- NULL                                       # FindClusters' extra copy
  Idents(obj) <- main_col
  clusters <- obj[[main_col]][, 1]
  message("  clusters at resolution ", leiden_resolution, ": ", nlevels(clusters))
  obj[[graph_name]] <- NULL; obj[[nn_name]] <- NULL; invisible(gc()) # saved above (STEP A/B files)
  
  ## resolution sweep table
  sweep <- rbindlist(lapply(seq_along(resolution_sweep), function(i) {
    cl_r <- obj[[cl_cols[i]]][, 1]; tab <- table(cl_r)
    data.table(resolution = resolution_sweep[i], n_clusters = length(tab), smallest_cluster = min(tab),
               largest_cluster = max(tab), n_clusters_lt_50_cells = sum(tab < 50),
               ARI_vs_resolution_main = round(igraph::compare(as.integer(cl_r), as.integer(clusters),
                                                              method = "adjusted.rand"), 4))
  }))
  print(sweep)
  fwrite(sweep, file.path(branch_dir, "Leiden_resolution_sweep_same_graph.csv"))
  
  ## ---- tables (same columns as script 12) -----------------------------------
  cells <- data.table(cell_id = colnames(obj), section_id = sec,
                      fov = obj$fov, x_mm = obj$x_slide_mm, y_mm = obj$y_slide_mm,
                      totalcounts = totalcounts, umap_xy, leiden = clusters)
  for (i in seq_along(cl_cols)) cells[, (paste0("leiden_res_", resolution_sweep[i])) := obj[[cl_cols[i]]][, 1]]
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
  
  ## section mixing on the kNN: own-section neighbours (self excluded) vs random (N_s - 1)/(N - 1)
  sec_int   <- as.integer(factor(sec))
  same_frac <- rowMeans(matrix(sec_int[idx[, -1]], nrow = n) == sec_int)
  N  <- length(sec); Ns <- table(sec)
  mix <- data.table(section_id = sec, same_frac = same_frac, expected = (as.numeric(Ns[sec]) - 1) / (N - 1))
  mix <- mix[, .(n_cells = .N,
                 mean_same_section_neighbour_fraction = round(mean(same_frac), 4),
                 expected_if_randomly_mixed = round(expected[1], 4),
                 ratio_observed_over_expected = round(mean(same_frac) / expected[1], 2)),
             by = section_id][order(section_id)]
  print(mix)
  fwrite(mix, file.path(branch_dir, "section_mixing_on_kNN_graph.csv"))
  
  ## ---- plots (Seurat where a function exists) --------------------------------
  cl_colours <- cluster_colours(levels(clusters))
  set.seed(seed); shuffled <- sample(colnames(obj))                 # random plotting order, as script 12
  
  ## UMAP by Leiden cluster (labels at cluster medians, as script 12)
  save_tiff(DimPlot(obj, reduction = umap_name, group.by = main_col, cols = cl_colours, label = TRUE,
                    label.box = TRUE, label.size = 3.5, shuffle = TRUE, seed = seed,
                    pt.size = 0.05, alpha = 0.3, raster = FALSE) +
              coord_fixed() + theme_bw(base_size = 13) +
              guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
              labs(title = paste0(branch_name, ": UMAP coloured by Leiden cluster"),
                   subtitle = paste0("Leiden on the neighbour graph (res ", leiden_resolution, ", ",
                                     nlevels(clusters), " clusters); UMAP = layout of the same graph"),
                   colour = "Leiden", x = "UMAP_1", y = "UMAP_2"),
            file.path(branch_dir, "UMAP_by_Leiden_cluster.tiff"), 11, 9)
  
  ## UMAP by tissue section
  save_tiff(DimPlot(obj, reduction = umap_name, group.by = batch_col, cols = section_colours,
                    shuffle = TRUE, seed = seed, pt.size = 0.05, alpha = 0.3, raster = FALSE) +
              coord_fixed() + theme_bw(base_size = 13) +
              guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
              labs(title = paste0(branch_name, ": UMAP coloured by tissue section"), colour = "Section",
                   x = "UMAP_1", y = "UMAP_2"),
            file.path(branch_dir, "UMAP_by_section.tiff"), 11, 9)
  
  ## UMAP split by section, coloured by Leiden cluster
  save_tiff(DimPlot(obj, reduction = umap_name, group.by = main_col, split.by = batch_col, ncol = 3,
                    cols = cl_colours, shuffle = TRUE, seed = seed, pt.size = 0.05, alpha = 0.3, raster = FALSE) +
              NoLegend() + coord_fixed() + theme_bw(base_size = 12) + theme(legend.position = "none") +
              labs(title = paste0(branch_name, ": UMAP split by section, coloured by Leiden cluster"),
                   x = "UMAP_1", y = "UMAP_2"),
            file.path(branch_dir, "UMAP_by_Leiden_split_by_section.tiff"), 15, 10)
  
  ## UMAP by total counts (FeaturePlot; log10 colour scale = script 12's log10(totalcounts) colours)
  save_tiff(FeaturePlot(obj, features = "nCount_RNA", reduction = umap_name, cells = shuffled,
                        pt.size = 0.05, alpha = 0.3, raster = FALSE) +
              scale_colour_viridis_c(option = "B", trans = "log10", labels = scales::label_comma()) +
              coord_fixed() + theme_bw(base_size = 13) +
              labs(title = paste0(branch_name, ": UMAP coloured by log10 total counts"),
                   subtitle = "Check that clusters are not driven by total counts",
                   colour = "total counts\n(log10 scale)", x = "UMAP_1", y = "UMAP_2"),
            file.path(branch_dir, "UMAP_by_log10_total_counts.tiff"), 11, 9)
  
  ## Leiden clusters in tissue space (ImageDimPlot on the slide-mm FOV)
  save_tiff(ImageDimPlot(obj, fov = "slide_mm", group.by = main_col, cols = cl_colours, size = 0.05,
                         alpha = 0.6, border.color = NA, dark.background = FALSE, axes = TRUE, flip_xy = FALSE) +
              geom_text(data = section_label_pos, aes(x, y, label = section), inherit.aes = FALSE,
                        colour = "black", fontface = "bold", size = 5) +
              coord_fixed() + theme_bw(base_size = 13) +
              guides(fill = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
              labs(title = paste0(branch_name, ": Leiden clusters in tissue space"),
                   x = "x (mm)", y = "y (mm)", fill = "Leiden"),
            file.path(branch_dir, "Spatial_Leiden_clusters_all_sections.tiff"), 13, 10)
  
  ## section composition of each cluster (ggplot, as script 12 - no Seurat function)
  p_comp <- ggplot(comp, aes(x = leiden, y = pct_of_cluster, fill = section_id)) +
    geom_col(width = 0.8) +
    geom_hline(yintercept = dominant_section_flag_pct, linetype = "dashed") +
    scale_fill_manual(values = section_colours) + theme_minimal(base_size = 13) +
    labs(title = paste0(branch_name, ": section composition of each Leiden cluster"),
         subtitle = paste0("Dashed line = ", dominant_section_flag_pct,
                           "%; above it is a flag to investigate (section-specific biology OR residual batch)"),
         x = "Leiden cluster", y = "% of cluster", fill = "Section")
  save_tiff(p_comp, file.path(branch_dir, "Leiden_cluster_composition_by_section.tiff"), 12, 6)
  
  ## resolution sweep visuals (plots only; clusters already computed above)
  sweep_tree <- plot_resolution_sweep(obj, sweep, cl_cols, umap_name, branch_dir, branch_name)
  
  ## same results list as script 12
  results <- list(branch = branch_name, embedding_used = X, knn = nn, graph = G,
                  umap = umap_xy, leiden = clusters, resolution_sweep = sweep,
                  resolution_sweep_tree = sweep_tree,
                  section_mixing = mix, cells = cells)
  saveRDS(results, file.path(branch_dir, "RESULTS_all_objects.RDS"))
  message(branch_name, " done -> ", branch_dir)
  list(obj = obj, results = results)
}


## =============================================================================
## CHUNK 6: PRIMARY - batch-aware PCA -> neighbour graph -> UMAP / Leiden
## =============================================================================
out <- run_branch_seurat(obj, "pearson_pca", "primary", dir_primary,
                         "PRIMARY (batch-aware Pearson PCA)", "Batch-aware Pearson PCA")
obj <- out$obj; res_primary <- out$results; rm(out)


## =============================================================================
## CHUNK 7: SENSITIVITY - batch-aware PCA -> Harmony -> neighbour graph -> UMAP / Leiden
## =============================================================================
message("\nRunning Harmony on all ", ncol(Embeddings(obj, "pearson_pca")), " batch-aware PCs (", batch_col, ")...")
set.seed(seed)
obj <- harmony::RunHarmony(obj, group.by.vars = batch_col, reduction.use = "pearson_pca",
                           dims.use = seq_len(ncol(Embeddings(obj, "pearson_pca"))),   # all PCs, as script 12
                           reduction.save = "harmony", project.dim = FALSE)           # no gene loadings to project
emb_harmony <- Embeddings(obj, "harmony")
stopifnot(identical(rownames(emb_harmony), colnames(obj)),
          ncol(emb_harmony) == ncol(Embeddings(obj, "pearson_pca")),
          "Harmony embeddings contain NaN/Inf" = all(is.finite(emb_harmony)))
saveRDS(emb_harmony, file.path(dir_sens, "harmony_embeddings.RDS"))
rm(emb_harmony)

out <- run_branch_seurat(obj, "harmony", "harmony", dir_sens,
                         "SENSITIVITY (batch-aware Pearson PCA + Harmony)", "Harmony (on batch-aware PCs)")
obj <- out$obj; res_sens <- out$results; rm(out)


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

## crosstab heatmap (ggplot, as script 12 - no Seurat function)
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
## CHUNK 9: save the Seurat object with the results under a new name
## =============================================================================
## Already in obj (added by the Seurat functions above):
##   reductions: pearson_pca, harmony, umap_primary, umap_harmony
##   metadata:   primary_leiden_res_<r>, harmony_leiden_res_<r> for every sweep resolution
## kNN and graphs are saved as RDS in each branch folder (as script 12), so they
## and the plotting-only FOV are dropped here (several GB at k = 250).
obj[["slide_mm"]] <- NULL
obj@graphs    <- list()
obj@neighbors <- list()
seurat_out <- file.path(dir_seu, paste0("seurat_", input_tag, "_SeuratFunctions_", run_tag, ".RDS"))
saveRDS(obj, seurat_out)
message("Saved Seurat object with results: ", seurat_out)


## =============================================================================
## CHUNK 10: summary
## =============================================================================
writeLines(capture.output(sessionInfo()), file.path(dir_seu, "sessionInfo.txt"))

message("\n", strrep("=", 78))
message("DONE. Everything is in: ", dir_seu)
message("  00_shared  : hvgs, PCA_variance_per_PC.csv, PCA_elbow_", run_tag, ".tiff (HVGs + PCA reused if saved)")
message("  01_PRIMARY : FindNeighbors kNN (STEP A) -> graph (STEP B) -> RunUMAP (STEP C) + FindClusters on graph (STEP D)")
message("               ", nlevels(res_primary$leiden), " clusters at resolution ", leiden_resolution)
message("  02_SENSITIVITY (Harmony): same steps; ", nlevels(res_sens$leiden), " clusters")
message("  01_ and 02_ each also hold the resolution sweep: Leiden_resolution_sweep_same_graph.csv +")
message("               _summary.tiff, _UMAP_panels.tiff, _cluster_tree.tiff (+ tree nodes/edges CSVs)")
message("  03_comparison: ARI = ", round(ari, 4), ", crosstab + section mixing")
message("  parameters_used_", run_tag, ".csv + sessionInfo.txt for the methods section")
message(strrep("=", 78))