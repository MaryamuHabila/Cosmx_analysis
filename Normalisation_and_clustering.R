## =============================================================================
## 12b_Seurat_functions_PCA_UMAP_Leiden.R
## =============================================================================
## PURPOSE
##   Same analysis, same outputs and same plots as 12_Normalisation_PCA_UMAP_Leiden.R
##   (+ the side-by-side UMAPs of script 13 + FindAllMarkers), but written with
##   SEURAT functions wherever Seurat has one:
##
##   step                         script 12                        here (Seurat)
##   ---------------------------  -------------------------------  --------------------------------------
##   HVGs                         FindVariableFeatures             FindVariableFeatures -> VariableFeatures(obj)
##   batch-aware Pearson PCA      scPearsonPCA                     scPearsonPCA (no Seurat equivalent),
##                                                                 stored as Seurat reduction "pearson_pca"
##   Harmony                      RunHarmony(matrix)               RunHarmony(Seurat object) -> reduction "harmony"
##   STEP A kNN (cosine, annoy)   RcppAnnoy loop                   FindNeighbors(return.neighbor = TRUE)
##   STEP B neighbour graph       uwot::similarity_graph           uwot::similarity_graph on the Seurat kNN,
##                                                                 stored as a Seurat Graph (as.Graph)
##   STEP C UMAP of the graph     uwot::optimize_graph_layout      RunUMAP(nn.name = kNN, umap.method = "uwot")
##   STEP D Leiden on the graph   igraph::cluster_leiden           FindClusters(graph.name = STEP B graph, algorithm = 4)
##   plots                        ggplot2                          ElbowPlot, DimPlot, FeaturePlot, ImageDimPlot
##   cluster tree                 hand-made ggplot                 clustree (built for Seurat cluster columns)
##
## WHY STEP B IS NOT FindNeighbors' SNN GRAPH
##   Script 12 (scanpy-style) clusters the UMAP fuzzy-connectivity graph - the
##   same graph UMAP lays out. Seurat's own graph (SNN) is a different graph.
##   To do EXACTLY what script 12 does, the fuzzy graph is built from Seurat's
##   kNN and handed to FindClusters. RunUMAP(nn.name = ...) builds this same
##   fuzzy graph internally (uwot, set_op_mix_ratio = 1, local_connectivity = 1)
##   from the same kNN, so UMAP and Leiden use one graph, as in script 12.
##   (Bonus: no SNN at k = 250, so far less RAM than the earlier 12b.)
##
## EXPECTED DIFFERENCE FROM SCRIPT 12
##   Annoy is an approximate search and Seurat's Leiden (leidenbase) is a
##   different implementation from igraph's. Method and parameters are the same,
##   so results are near-identical, not bit-identical. CHUNK 8 reports the ARI
##   against script 12's clusters as a check.
##
## OUTPUT: <script 12 folder>/seurat_functions_used/PC<n>_NN<k>/
##   00_shared/ , 01_PRIMARY_batchawarePCA/ , 02_SENSITIVITY_batchawarePCA_Harmony/ ,
##   03_PRIMARY_vs_SENSITIVITY/
##   Each branch folder holds the same file names as script 12, so script 13
##   also runs on these folders.
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 0: packages - checked, NOT silently installed
## =============================================================================
required_pkgs <- c("data.table", "Matrix", "ggplot2", "patchwork", "scales", "Seurat", "SeuratObject",
                   "scPearsonPCA", "RSpectra", "uwot", "igraph", "harmony", "leidenbase",
                   "clustree", "openxlsx")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing package(s): ", paste(missing_pkgs, collapse = ", "), "\nInstall with:\n",
       "  install.packages(c('data.table','Matrix','ggplot2','patchwork','scales','Seurat','RSpectra',\n",
       "                     'uwot','igraph','harmony','leidenbase','clustree','openxlsx'))\n",
       "  remotes::install_github('Nanostring-Biostats/CosMx-Analysis-Scratch-Space',\n",
       "                          subdir = '_code/scPearsonPCA', ref = 'Main')")
}
if (utils::packageVersion("Seurat") < "5.1.0") {                  # cluster.name + leidenbase Leiden
  stop("Seurat >= 5.1.0 needed (FindClusters(algorithm = 4) via leidenbase). install.packages('Seurat')")
}
if (!"similarity_graph" %in% getNamespaceExports("uwot")) {        # STEP B
  stop("uwot too old (no similarity_graph). install.packages('uwot')")
}
if (utils::packageVersion("harmony") < "1.0") {                    # RunHarmony Seurat method
  stop("harmony >= 1.0 needed. install.packages('harmony')")
}

library(data.table)   # tables
library(Matrix)       # sparse matrices
library(ggplot2)      # extra plot layers / plots Seurat has no function for
library(patchwork)    # combine plots ( | and wrap_plots )
library(Seurat)       # FindVariableFeatures, FindNeighbors, RunUMAP, FindClusters, DimPlot, ...
library(scPearsonPCA) # batch-aware Pearson-residual PCA (vignette)


## =============================================================================
## CHUNK 1: paths + parameters (same values as script 12)
## =============================================================================
output_dir <- "/Users/maryamuhabilausman/Desktop/Analysis_AAA_29th_Sept_Seurat_Results"
seurat_in  <- file.path(output_dir, "seurat_QC_filtered_per_section_03_to_09.RDS")  # script 11 (read only)
input_tag  <- sub("\\.RDS$", "", basename(seurat_in))                              # folder name follows the input
dir_main   <- file.path(output_dir, "12_Normalisation_PCA_UMAP_Leiden", input_tag)  # script 12 folder

batch_col           <- "sample_id"    # the 6 AAA sections
n_hvgs              <- 2000           # vignette: 2000 for a 6k panel
n_pcs               <- 20             # vignette: 40                        <- change here
n_neighbors         <- 250            # vignette: 15 (incl. the cell itself) <- change here
rerun_pca           <- FALSE          # FALSE = reuse saved HVGs + PCA (script 12's or this script's)
annoy_n_trees       <- 50             # annoy index size (as script 12)
umap_min_dist       <- 0.02           # vignette
umap_spread         <- 2              # vignette
leiden_resolution   <- 1              # main resolution                      <- change here
leiden_n_iterations <- 2              # vignette
resolution_sweep    <- c(0.1, 0.2, 0.3, 0.5, 0.8, 1.0, 1.5, 2.0)  # all on the same graph <- change here
dominant_section_flag_pct <- 90       # flag clusters with > 90% cells from one section
sweep_umap_max_cells   <- 200000      # PLOTTING ONLY: cells drawn per sweep UMAP panel
sweep_tree_min_in_prop <- 0.05        # PLOTTING ONLY: cluster-tree edges drawn if >= 5% of child
seed    <- 42
dpi_out <- 300
stopifnot(leiden_resolution %in% resolution_sweep)

run_tag     <- paste0("PC", n_pcs, "_NN", n_neighbors)                     # e.g. "PC20_NN250"
dir_seu     <- file.path(dir_main, "seurat_functions_used", run_tag)       # all outputs here
dir_shared  <- file.path(dir_seu, "00_shared")
dir_primary <- file.path(dir_seu, "01_PRIMARY_batchawarePCA")
dir_sens    <- file.path(dir_seu, "02_SENSITIVITY_batchawarePCA_Harmony")
dir_compare <- file.path(dir_seu, "03_PRIMARY_vs_SENSITIVITY")
for (d in c(dir_shared, dir_primary, dir_sens, dir_compare)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

## parameters table for the methods section
fwrite(data.table(
  parameter = c("input", "batch_variable", "n_HVGs", "HVG_method", "PCA", "PCA_scale.max", "PCA_do.scale",
                "PCA_do.center", "Harmony", "n_PCs_for_neighbours", "k.param (incl. self)", "kNN",
                "kNN distances", "graph", "UMAP", "UMAP_min_dist", "UMAP_spread", "UMAP_init", "clustering",
                "Leiden_n_iterations", "group.singletons", "Leiden_resolution", "Leiden_resolution_sweep",
                "sweep_UMAP_max_cells_per_panel (plot only)", "sweep_tree_min_in_prop (plot only)", "seed"),
  value = c(basename(seurat_in), batch_col, n_hvgs, "Seurat::FindVariableFeatures (vst)",
            "scPearsonPCA::sparse_quasipoisson_pca_seurat_batch (batch-aware Pearson residuals)",
            10, TRUE, TRUE, paste0("harmony::RunHarmony (Seurat method) on all batch-aware PCs, group.by.vars = ", batch_col),
            n_pcs, n_neighbors,
            paste0("Seurat::FindNeighbors(nn.method = 'annoy', annoy.metric = 'cosine', n.trees = ", annoy_n_trees, ")"),
            "exact cosine distance (1 - cos) on the PCs; self forced to neighbour 1",
            "uwot::similarity_graph on the kNN (UMAP fuzzy connectivities, as sc.pp.neighbors), stored as Seurat Graph",
            "Seurat::RunUMAP(nn.name = kNN, umap.method = 'uwot') - same fuzzy graph, laid out",
            umap_min_dist, umap_spread, "spectral (uwot default)",
            "Seurat::FindClusters(algorithm = 4, Leiden via leidenbase) on that SAME graph (never on UMAP)",
            leiden_n_iterations, FALSE, leiden_resolution, paste(resolution_sweep, collapse = ", "),
            sweep_umap_max_cells, sweep_tree_min_in_prop, seed)),
  file.path(dir_seu, paste0("parameters_used_", run_tag, ".csv")))


## =============================================================================
## CHUNK 2: load the QC-filtered object (script 11) + checks
## =============================================================================
stopifnot("input Seurat object not found - run script 11 first" = file.exists(seurat_in))
obj <- readRDS(seurat_in)
obj
stopifnot(all(c(batch_col, "x_slide_mm", "y_slide_mm", "fov", "nCount_RNA") %in% colnames(obj[[]])))
if ("cell_id" %in% colnames(obj[[]]))                                     # as script 12
  stopifnot("cell_id column must match the object's cells" = identical(obj$cell_id, colnames(obj)))

counts_t    <- LayerData(obj, assay = "RNA", layer = "counts")            # genes x cells (raw counts)
totalcounts <- Matrix::colSums(counts_t)                                  # per-cell total, all genes
if (any(totalcounts == 0)) stop(sum(totalcounts == 0), " cells have zero counts.")
message("Cells where matrix total differs from nCount_RNA: ", sum(totalcounts != obj$nCount_RNA))

## genes with zero counts in every kept cell carry no information (as script 12)
zero_genes <- rownames(counts_t)[Matrix::rowSums(counts_t) == 0]
if (length(zero_genes) > 0) message(length(zero_genes), " genes have zero counts in all kept cells (not HVGs).")

sec      <- obj[[batch_col]][, 1]                                         # section of each cell
sections <- sort(unique(sec))
message("Loaded ", format(ncol(obj), big.mark = ","), " cells x ", nrow(obj), " genes; ",
        length(sections), " sections: ", paste(sections, collapse = ", "))
print(table(sec))

## slide coordinates (mm) as a Seurat FOV -> ImageDimPlot uses the same
## coordinates as all earlier plots (x_slide_mm / y_slide_mm)
obj[["slide_mm"]] <- CreateFOV(coords = data.frame(x = obj$x_slide_mm, y = obj$y_slide_mm, cell = colnames(obj)),
                               type = "centroids", assay = "RNA")
section_label_pos <- data.table(x = obj$x_slide_mm, y = obj$y_slide_mm, section = sec)[
  , .(x = mean(x), y = mean(y)), by = section]                            # section label positions


## =============================================================================
## CHUNK 3: HVGs (Seurat FindVariableFeatures) + batch-aware Pearson PCA
## =============================================================================
## Reuse saved HVGs + PCA (this script's 00_shared, else script 12's) unless rerun_pca = TRUE.
reuse_dir <- NULL
if (!rerun_pca) for (d in c(dir_shared, file.path(dir_main, "00_shared_HVGs_and_batch_aware_Pearson_PCA"))) {
  if (is.null(reuse_dir) && all(file.exists(file.path(d, c("hvgs.RDS", "pcaobj_batch_aware.RDS"))))) reuse_dir <- d
}

if (!is.null(reuse_dir)) {
  hvgs   <- readRDS(file.path(reuse_dir, "hvgs.RDS"))                     # HVGs used for the PCA
  pcaobj <- readRDS(file.path(reuse_dir, "pcaobj_batch_aware.RDS"))       # PCA result
  message("HVGs (", length(hvgs), ") and batch-aware PCA reused from ", reuse_dir)
} else {
  ## HVGs - Seurat vst on the counts layer
  obj  <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = n_hvgs)
  hvgs <- VariableFeatures(obj)
  
  ## remove HVGs with zero counts in a whole section (per-section frequency 0 ->
  ## Pearson residual undefined; scPearsonPCA's input check misses sparse zeros)
  det_by_sec <- sapply(sections, function(s) Matrix::rowSums(counts_t[hvgs, sec == s, drop = FALSE]) > 0)
  hvg_zero <- hvgs[rowSums(!det_by_sec) > 0]
  if (length(hvg_zero) > 0) {
    message("HVGs with ZERO counts in a whole section (removed): ", length(hvg_zero))
    fwrite(data.table(gene = hvg_zero), file.path(dir_shared, "HVGs_removed_zero_counts_in_a_section.csv"))
    hvgs <- setdiff(hvgs, hvg_zero)
  } else message("All HVGs have counts in every section - none removed.")
  saveRDS(hvgs, file.path(dir_shared, "hvgs.RDS"))
  fwrite(data.table(rank = seq_along(hvgs), gene = hvgs), file.path(dir_shared, "hvgs.csv"))
  
  ## batch-aware Pearson-residual PCA - exactly the vignette call
  obs <- data.table(cell_id = colnames(obj), sample_id = sec)
  genefreq_batch <- scPearsonPCA::gene_frequency(x = counts_t, obs = obs,
                                                 cellid_colname = "cell_id", batch_variable = "sample_id")
  gf_hvg <- genefreq_batch[hvgs, , drop = FALSE]                          # per-section frequency of each HVG
  if (any(!is.finite(as.matrix(gf_hvg))) || any(as.matrix(gf_hvg) <= 0)) stop("Non-positive gene frequencies for some HVGs.")
  saveRDS(genefreq_batch, file.path(dir_shared, "genefreq_batch.RDS"))
  set.seed(seed)
  pcaobj <- sparse_quasipoisson_pca_seurat_batch(
    counts_t[hvgs, ], totalcounts = totalcounts, grate = gf_hvg, obs = obs,
    batch_variable = "sample_id", cellid_colname = "cell_id",
    scale.max = 10, do.scale = TRUE, do.center = TRUE)
  saveRDS(pcaobj, file.path(dir_shared, "pcaobj_batch_aware.RDS"))
}
VariableFeatures(obj) <- hvgs                                             # HVGs recorded in the object
rm(counts_t); invisible(gc())                                             # raw matrix no longer needed

## PCA -> Seurat reduction "pearson_pca" (columns must be key + number: PPC_1, PPC_2, ...)
with_key <- function(m, key) { m <- as.matrix(m); colnames(m) <- paste0(key, seq_len(ncol(m))); m }
emb_pca <- pcaobj$reduction.data@cell.embeddings
if (is.null(rownames(emb_pca))) { stopifnot(nrow(emb_pca) == ncol(obj)); rownames(emb_pca) <- colnames(obj) }
emb_pca <- emb_pca[colnames(obj), , drop = FALSE]                          # same cell order as the object
stopifnot("PCA embeddings contain NaN/Inf" = all(is.finite(emb_pca)))
if (ncol(emb_pca) < n_pcs) stop("PCA has ", ncol(emb_pca), " PCs, fewer than n_pcs = ", n_pcs)
pc_sd <- pcaobj$reduction.data@stdev                                      # needed by ElbowPlot
obj[["pearson_pca"]] <- CreateDimReducObject(embeddings = with_key(emb_pca, "PPC_"), stdev = pc_sd,
                                             key = "PPC_", assay = "RNA")
rm(pcaobj, emb_pca); invisible(gc())

## variance per PC (table) + Seurat ElbowPlot with the n_pcs cut-off
elbow <- data.table(PC = seq_along(pc_sd), variance = pc_sd^2)
elbow[, pct_variance_of_computed_PCs := 100 * variance / sum(variance)]
fwrite(elbow, file.path(dir_shared, "PCA_variance_per_PC.csv"))
p_elbow <- ElbowPlot(obj, reduction = "pearson_pca", ndims = length(pc_sd)) +
  geom_vline(xintercept = n_pcs, colour = "red", linetype = "dashed") +
  labs(title = "Batch-aware Pearson PCA - elbow plot",
       subtitle = paste0("Red line = ", n_pcs, " PCs used for the neighbour graph"))
ggsave(file.path(dir_shared, paste0("PCA_elbow_", run_tag, ".tiff")), p_elbow, width = 8, height = 5.5,
       dpi = dpi_out, device = "tiff", compression = "lzw", bg = "white")


## =============================================================================
## CHUNK 4: helpers - colours, saving, cluster relabelling
## =============================================================================
distinct_palette <- c(                                                    # same palette/order as script 12
  "#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#42D4F4", "#F032E6",
  "#9A6324", "#469990", "#BFEF45", "#800000", "#FFE119", "#000075", "#808000",
  "#FF8FB1", "#A9A9A9", "#005300", "#009FFF", "#C77CFF", "#FF6F61", "#00FFBE",
  "#D4A017", "#7F7FFF", "#8B008B", "#2F4F4F", "#E9967A", "#556B2F", "#FF1493",
  "#20B2AA", "#B8860B", "#6A5ACD", "#CD853F", "#4B0082", "#00CED1", "#8FBC8F",
  "#A0522D")
cluster_colours <- function(lv) {                                         # one colour per cluster level
  if (length(lv) <= length(distinct_palette)) return(setNames(distinct_palette[seq_along(lv)], lv))
  stopifnot("install.packages('Polychrome') for > 36 clusters" = requireNamespace("Polychrome", quietly = TRUE))
  set.seed(seed)
  setNames(unname(Polychrome::createPalette(length(lv), seedcolors = distinct_palette, range = c(25, 80))), lv)
}
section_colours <- setNames(scales::hue_pal()(length(sections)), sections) # as script 12

save_tiff <- function(p, path, w, h) {
  ggsave(path, p, width = w, height = h, dpi = dpi_out, device = "tiff", compression = "lzw", bg = "white")
  message("Saved: ", path)
}

## labels 0,1,2,... by decreasing cluster size (scanpy / script 12 convention)
relabel_by_size <- function(f) {
  sizes <- sort(table(f), decreasing = TRUE)
  new   <- setNames(seq_along(sizes) - 1L, names(sizes))
  factor(unname(new[as.character(f)]), levels = 0:(length(sizes) - 1L))
}


## =============================================================================
## CHUNK 5: one analysis branch with Seurat functions
## =============================================================================
## reduction = "pearson_pca" (PRIMARY) or "harmony" (SENSITIVITY)
## prefix    = short name used for kNN, graph, UMAP and cluster columns
run_branch_seurat <- function(obj, reduction, prefix, branch_dir, branch_name) {
  message("\n", strrep("=", 78), "\n", branch_name, "\n", strrep("=", 78))
  nn_name    <- paste0(prefix, "_neighbors")                              # Seurat Neighbor object (kNN)
  graph_name <- paste0(prefix, "_connectivities")                         # STEP B graph (as scanpy's name)
  umap_name  <- paste0("umap_", prefix)                                   # UMAP reduction
  cl_cols    <- paste0(prefix, "_leiden_res_", resolution_sweep)          # one cluster column per resolution
  main_col   <- paste0(prefix, "_leiden_res_", leiden_resolution)
  
  ## ---- PC1 vs PC2 of the input reduction, by section ---------------------
  save_tiff(DimPlot(obj, reduction = reduction, dims = 1:2, group.by = batch_col, cols = section_colours,
                    shuffle = TRUE, seed = seed, pt.size = 0.05, alpha = 0.3, raster = FALSE) +
              coord_fixed() + guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
              labs(title = paste0(branch_name, ": dimension 1 vs 2, coloured by section"), colour = "Section"),
            file.path(branch_dir, "PC1_vs_PC2_by_section.tiff"), 9, 8)
  
  ## ---- STEP A: kNN on the PCs (FindNeighbors) ------------------------------
  message("STEP A: FindNeighbors on ", n_pcs, " dims of ", reduction, " (cosine, k = ", n_neighbors, " incl. self)")
  obj <- FindNeighbors(obj, reduction = reduction, dims = 1:n_pcs, k.param = n_neighbors,
                       nn.method = "annoy", annoy.metric = "cosine", n.trees = annoy_n_trees,
                       return.neighbor = TRUE, graph.name = nn_name)
  idx <- Indices(obj[[nn_name]])                                          # cells x k, row numbers of neighbours
  n   <- nrow(idx); k <- ncol(idx)
  
  ## the cell itself must be neighbour 1 (sc.pp.neighbors / script 12); fix the rare rows where it is not
  bad <- which(idx[, 1] != seq_len(n))
  message("  cells without themselves as neighbour 1 (fixed): ", length(bad))
  for (i in bad) idx[i, ] <- c(i, setdiff(idx[i, ], i))[seq_len(k)]
  
  ## exact cosine distance (1 - cos) to each neighbour, from the same PCs
  ## (= what script 12 / uwot use; avoids depending on how Seurat reports annoy distances)
  X  <- Embeddings(obj, reduction)[, 1:n_pcs]
  Xn <- X / sqrt(rowSums(X^2))                                            # unit-length rows
  dst <- vapply(seq_len(k), function(j) 1 - rowSums(Xn * Xn[idx[, j], , drop = FALSE]), numeric(n))
  dst[dst < 0] <- 0; dst[, 1] <- 0                                        # rounding; self = 0
  rm(X, Xn); invisible(gc())
  
  nb <- obj[[nn_name]]; nb@nn.idx <- idx; nb@nn.dist <- dst; obj[[nn_name]] <- nb   # write back to the Neighbor object
  saveRDS(obj[[nn_name]], file.path(branch_dir, "STEP_A_kNN_Seurat_Neighbor.RDS"))
  
  ## ---- STEP B: ONE neighbour graph (UMAP fuzzy connectivities) -------------
  message("STEP B: neighbour graph from the kNN (uwot::similarity_graph), stored as a Seurat Graph")
  G <- uwot::similarity_graph(X = NULL, nn_method = list(idx = idx, dist = dst),
                              set_op_mix_ratio = 1, local_connectivity = 1)  # = RunUMAP's internal graph
  G <- methods::as(G, "CsparseMatrix")
  stopifnot("neighbour graph is not symmetric" = Matrix::isSymmetric(G, tol = 1e-6))
  dimnames(G) <- list(colnames(obj), colnames(obj))
  saveRDS(G, file.path(branch_dir, "STEP_B_neighbour_graph_connectivities.RDS"))
  message("  graph: ", format(n, big.mark = ","), " cells, ", format(Matrix::nnzero(G) / 2, big.mark = ","), " edges")
  g_obj <- as.Graph(G); DefaultAssay(g_obj) <- "RNA"                      # Seurat Graph object
  obj[[graph_name]] <- g_obj
  message("  connected components in graph: ",
          igraph::components(igraph::graph_from_adjacency_matrix(G, mode = "undirected", weighted = TRUE))$no)
  rm(g_obj); invisible(gc())                                              # G kept for RESULTS_all_objects.RDS
  
  ## ---- STEP C: UMAP = layout of that graph (RunUMAP from the kNN) ---------
  message("STEP C: RunUMAP from the kNN (uwot builds the same graph, then lays it out)")
  obj <- RunUMAP(obj, nn.name = nn_name, umap.method = "uwot", min.dist = umap_min_dist,
                 spread = umap_spread, seed.use = seed, reduction.name = umap_name,
                 reduction.key = paste0("UMAP", toupper(substr(prefix, 1, 1)), "_"))   # UMAPP_ / UMAPH_
  
  ## ---- STEP D: Leiden on the STEP B graph, every resolution (FindClusters) -
  message("STEP D: FindClusters (Leiden) on ", graph_name, " - resolutions ", paste(resolution_sweep, collapse = ", "))
  for (i in seq_along(resolution_sweep)) {
    obj <- FindClusters(obj, graph.name = graph_name, algorithm = 4, resolution = resolution_sweep[i],
                        n.iter = leiden_n_iterations, random.seed = seed, group.singletons = FALSE,
                        cluster.name = cl_cols[i], verbose = FALSE)
    obj[[cl_cols[i]]] <- relabel_by_size(obj[[cl_cols[i]]][, 1])          # 0 = largest
  }
  obj$seurat_clusters <- NULL                                             # FindClusters' extra copy, not needed
  Idents(obj) <- main_col
  clusters <- obj[[main_col]][, 1]
  message("  clusters at resolution ", leiden_resolution, ": ", nlevels(clusters))
  obj[[graph_name]] <- NULL                                               # graph saved above; free RAM
  
  ## ---- resolution sweep table -----------------------------------------------
  sweep <- rbindlist(lapply(seq_along(resolution_sweep), function(i) {
    cl <- obj[[cl_cols[i]]][, 1]; tab <- table(cl)
    data.table(resolution = resolution_sweep[i], n_clusters = length(tab), smallest_cluster = min(tab),
               largest_cluster = max(tab), n_clusters_lt_50_cells = sum(tab < 50),
               ARI_vs_resolution_main = round(igraph::compare(as.integer(cl), as.integer(clusters),
                                                              method = "adjusted.rand"), 4))
  }))
  print(sweep)
  fwrite(sweep, file.path(branch_dir, "Leiden_resolution_sweep_same_graph.csv"))
  
  ## ---- tables (same columns as script 12 -> script 13 runs on these too) ---
  um <- Embeddings(obj, umap_name)
  cells <- data.table(cell_id = colnames(obj), section_id = sec, fov = obj$fov,
                      x_mm = obj$x_slide_mm, y_mm = obj$y_slide_mm, totalcounts = obj$nCount_RNA,
                      UMAP_1 = um[, 1], UMAP_2 = um[, 2], leiden = clusters)
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
  Ns  <- table(sec)
  mix <- data.table(section_id = sec, same_frac, expected = (as.numeric(Ns[sec]) - 1) / (n - 1))[
    , .(n_cells = .N, mean_same_section_neighbour_fraction = round(mean(same_frac), 4),
        expected_if_randomly_mixed = round(expected[1], 4),
        ratio_observed_over_expected = round(mean(same_frac) / expected[1], 2)), by = section_id][order(section_id)]
  print(mix)
  fwrite(mix, file.path(branch_dir, "section_mixing_on_kNN_graph.csv"))
  obj[[nn_name]] <- NULL; invisible(gc())                                 # kNN saved above; idx/dst kept for RESULTS
  
  ## ---- plots (Seurat) -------------------------------------------------------
  cols_main <- cluster_colours(levels(clusters))
  
  save_tiff(DimPlot(obj, reduction = umap_name, group.by = main_col, cols = cols_main, label = TRUE,
                    label.box = TRUE, shuffle = TRUE, seed = seed, pt.size = 0.05, alpha = 0.3, raster = FALSE) +
              coord_fixed() + guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
              labs(title = paste0(branch_name, ": UMAP coloured by Leiden cluster"),
                   subtitle = paste0("Leiden on the neighbour graph (res ", leiden_resolution, ", ", nlevels(clusters),
                                     " clusters); UMAP = layout of the same graph"), colour = "Leiden"),
            file.path(branch_dir, "UMAP_by_Leiden_cluster.tiff"), 11, 9)
  
  save_tiff(DimPlot(obj, reduction = umap_name, group.by = batch_col, cols = section_colours,
                    shuffle = TRUE, seed = seed, pt.size = 0.05, alpha = 0.3, raster = FALSE) +
              coord_fixed() + guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1))) +
              labs(title = paste0(branch_name, ": UMAP coloured by tissue section"), colour = "Section"),
            file.path(branch_dir, "UMAP_by_section.tiff"), 11, 9)
  
  save_tiff(DimPlot(obj, reduction = umap_name, group.by = main_col, split.by = batch_col, ncol = 3,
                    cols = cols_main, shuffle = TRUE, seed = seed, pt.size = 0.05, alpha = 0.3, raster = FALSE) +
              NoLegend() + labs(title = paste0(branch_name, ": UMAP split by section, coloured by Leiden cluster")),
            file.path(branch_dir, "UMAP_by_Leiden_split_by_section.tiff"), 15, 10)
  
  save_tiff(FeaturePlot(obj, features = "nCount_RNA", reduction = umap_name, pt.size = 0.05, raster = FALSE) +
              scale_colour_viridis_c(option = "B", trans = "log10", labels = scales::label_comma()) +
              coord_fixed() +
              labs(title = paste0(branch_name, ": UMAP coloured by total counts (log10 colour scale)"),
                   subtitle = "Check that clusters are not driven by total counts", colour = "nCount_RNA"),
            file.path(branch_dir, "UMAP_by_log10_total_counts.tiff"), 11, 9)
  
  save_tiff(ImageDimPlot(obj, fov = "slide_mm", group.by = main_col, cols = cols_main, size = 0.1,
                         border.color = NA, dark.background = FALSE, axes = TRUE, flip_xy = FALSE) +
              geom_text(data = section_label_pos, aes(x, y, label = section), inherit.aes = FALSE,
                        colour = "black", fontface = "bold", size = 5) +
              guides(fill = guide_legend(override.aes = list(size = 4), ncol = 2)) +
              labs(title = paste0(branch_name, ": Leiden clusters in tissue space"),
                   x = "x (mm)", y = "y (mm)", fill = "Leiden"),
            file.path(branch_dir, "Spatial_Leiden_clusters_all_sections.tiff"), 13, 10)
  
  p_comp <- ggplot(comp, aes(leiden, pct_of_cluster, fill = section_id)) +   # no Seurat function for this
    geom_col(width = 0.8) + geom_hline(yintercept = dominant_section_flag_pct, linetype = "dashed") +
    scale_fill_manual(values = section_colours) + theme_minimal(base_size = 13) +
    labs(title = paste0(branch_name, ": section composition of each Leiden cluster"),
         subtitle = paste0("Dashed line = ", dominant_section_flag_pct,
                           "%; above it is a flag to investigate (section-specific biology OR residual batch)"),
         x = "Leiden cluster", y = "% of cluster", fill = "Section")
  save_tiff(p_comp, file.path(branch_dir, "Leiden_cluster_composition_by_section.tiff"), 12, 6)
  
  ## ---- resolution sweep plots -----------------------------------------------
  ## 1. summary curves (ggplot; Seurat has no function for this)
  sw_long <- rbind(
    sweep[, .(resolution, metric = "Number of clusters", value = n_clusters)],
    sweep[, .(resolution, metric = "Clusters with < 50 cells", value = n_clusters_lt_50_cells)],
    sweep[, .(resolution, metric = "Smallest cluster (cells)", value = smallest_cluster)],
    sweep[, .(resolution, metric = paste0("ARI vs main resolution (", leiden_resolution, ")"), value = ARI_vs_resolution_main)])
  sw_long[, metric := factor(metric, levels = unique(metric))]
  p_sum <- ggplot(sw_long, aes(resolution, value)) +
    geom_vline(xintercept = leiden_resolution, colour = "red", linetype = "dashed") +
    geom_line() + geom_point(size = 2.5) +
    scale_x_log10(breaks = resolution_sweep, labels = resolution_sweep) +
    facet_wrap(~ metric, ncol = 2, scales = "free_y") + theme_bw(base_size = 13) +
    labs(title = paste0(branch_name, ": Leiden resolution sweep (same graph)"),
         subtitle = "Red dashed line = main resolution; x axis on log scale", x = "Leiden resolution", y = NULL)
  save_tiff(p_sum, file.path(branch_dir, "Leiden_resolution_sweep_summary.tiff"), 11, 8)
  
  ## 2. one UMAP panel per resolution (DimPlot); random subsample of cells drawn if very many
  set.seed(seed)
  cells_draw <- if (ncol(obj) > sweep_umap_max_cells) sample(colnames(obj), sweep_umap_max_cells) else colnames(obj)
  max_k <- max(sweep$n_clusters)
  panels <- DimPlot(obj, reduction = umap_name, group.by = cl_cols, cells = cells_draw,
                    cols = cluster_colours(as.character(0:(max_k - 1L))), label = TRUE, label.size = 2.3,
                    pt.size = 0.03, alpha = 0.3, raster = FALSE, combine = FALSE)
  panels <- Map(function(p, r, cc) p + NoLegend() + coord_fixed() +
                  labs(title = paste0("res ", r, "  (", nlevels(obj[[cc]][, 1]), " clusters)",
                                      ifelse(r == leiden_resolution, "  [main]", ""))),
                panels, resolution_sweep, cl_cols)
  save_tiff(wrap_plots(panels, ncol = 4) +
              plot_annotation(title = paste0(branch_name, ": UMAP coloured by Leiden cluster at each sweep resolution"),
                              subtitle = "Same UMAP in every panel; cluster numbers are by size within each resolution"),
            file.path(branch_dir, "Leiden_resolution_sweep_UMAP_panels.tiff"), 20, 11)
  
  ## 3. cluster tree across resolutions (clustree) + its nodes/edges tables
  tree_df <- obj[[]][, cl_cols]                                           # one column per resolution
  tree_prefix <- paste0(prefix, "_leiden_res_")
  p_tree <- clustree::clustree(tree_df, prefix = tree_prefix, prop_filter = sweep_tree_min_in_prop) +
    labs(title = paste0(branch_name, ": Leiden cluster tree across resolutions"),
         subtitle = paste0("Edges: fraction of the lower cluster's cells from the upper cluster (drawn if >= ",
                           sweep_tree_min_in_prop, ")"))
  save_tiff(p_tree, file.path(branch_dir, "Leiden_resolution_sweep_cluster_tree.tiff"), 18, 11)
  
  ## tree tables - same columns as script 12 (clustree's own tables have different columns)
  ## nodes: resolution, level, cluster, n_cells, x ; edges: res_from, parent, res_to, cluster, n_cells, in_prop
  ## in_prop = fraction of the child cluster's cells that came from that parent
  res_ord <- sort(resolution_sweep)
  get_cl  <- function(r) as.character(obj[[paste0(prefix, "_leiden_res_", r)]][, 1])
  node_list <- vector("list", length(res_ord)); edge_list <- list(); x_prev <- NULL
  for (i in seq_along(res_ord)) {
    r  <- res_ord[i]; cl <- get_cl(r)
    nd <- as.data.table(table(cluster = cl)); setnames(nd, "N", "n_cells")
    if (i == 1L) {
      nd[, order_key := as.numeric(cluster)]
    } else {
      ct_i <- as.data.table(table(parent = get_cl(res_ord[i - 1L]), cluster = cl))[N > 0]
      ct_i[, parent_x := x_prev[parent]]
      nd <- merge(nd, ct_i[, .(order_key = sum(N * parent_x) / sum(N)), by = cluster], by = "cluster")
      ct_i[, in_prop := N / sum(N), by = cluster]
      edge_list[[length(edge_list) + 1L]] <-
        ct_i[, .(res_from = res_ord[i - 1L], parent, res_to = r, cluster, n_cells = N, in_prop = round(in_prop, 4))]
    }
    setorder(nd, order_key)                                               # children placed under their parents
    nd[, x := (seq_len(.N) - 0.5) / .N]
    nd[, `:=`(resolution = r, level = i)]
    x_prev <- setNames(nd$x, nd$cluster)
    node_list[[i]] <- nd[, .(resolution, level, cluster, n_cells, x)]
  }
  tree_nodes <- rbindlist(node_list); tree_edges <- rbindlist(edge_list)
  fwrite(tree_nodes, file.path(branch_dir, "Leiden_resolution_sweep_cluster_tree_nodes.csv"))
  fwrite(tree_edges, file.path(branch_dir, "Leiden_resolution_sweep_cluster_tree_edges.csv"))
  
  ## per-branch results list - same elements as script 12's RESULTS_all_objects.RDS
  saveRDS(list(branch = branch_name, embedding_used = Embeddings(obj, reduction)[, 1:n_pcs],
               knn = list(idx = idx, dist = dst), graph = G, umap = um, leiden = clusters,
               resolution_sweep = sweep, resolution_sweep_tree = list(nodes = tree_nodes, edges = tree_edges),
               section_mixing = mix, cells = cells),
          file.path(branch_dir, "RESULTS_all_objects.RDS"))
  rm(idx, dst, G); invisible(gc())
  message(branch_name, " done -> ", branch_dir)
  obj
}


## =============================================================================
## CHUNK 6: PRIMARY - batch-aware Pearson PCA -> graph -> UMAP / Leiden
## =============================================================================
obj <- run_branch_seurat(obj, "pearson_pca", "primary", dir_primary, "PRIMARY (batch-aware Pearson PCA)")


## =============================================================================
## CHUNK 7: SENSITIVITY - Harmony (Seurat method) -> graph -> UMAP / Leiden
## =============================================================================
message("\nRunning Harmony on all ", ncol(Embeddings(obj, "pearson_pca")), " batch-aware PCs (", batch_col, ")...")
set.seed(seed)
obj <- harmony::RunHarmony(obj, group.by.vars = batch_col, reduction.use = "pearson_pca",
                           dims.use = seq_len(ncol(Embeddings(obj, "pearson_pca"))),  # all PCs, as script 12
                           reduction.save = "harmony", project.dim = FALSE)          # no gene loadings to project
stopifnot("Harmony embeddings contain NaN/Inf" = all(is.finite(Embeddings(obj, "harmony"))))
saveRDS(Embeddings(obj, "harmony"), file.path(dir_sens, "harmony_embeddings.RDS"))

obj <- run_branch_seurat(obj, "harmony", "harmony", dir_sens, "SENSITIVITY (batch-aware Pearson PCA + Harmony)")


## =============================================================================
## CHUNK 8: PRIMARY vs SENSITIVITY at the main resolution (+ check vs script 12)
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
fwrite(dcast(ct, primary ~ sensitivity, value.var = "pct_of_primary_cluster", fill = 0),
       file.path(dir_compare, "crosstab_rowpct_PRIMARY_rows_x_SENSITIVITY_cols.csv"))
fwrite(data.table(comparison = "PRIMARY vs SENSITIVITY Leiden", resolution = leiden_resolution,
                  n_clusters_primary = nlevels(cl_p), n_clusters_sensitivity = nlevels(cl_s),
                  adjusted_rand_index = round(ari, 4)),
       file.path(dir_compare, "ARI_PRIMARY_vs_SENSITIVITY.csv"))

p_ct <- ggplot(ct, aes(sensitivity, primary, fill = pct_of_primary_cluster)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = ifelse(pct_of_primary_cluster >= 5, sprintf("%.0f", pct_of_primary_cluster), "")), size = 2.8) +
  scale_fill_viridis_c() + theme_minimal(base_size = 12) +
  labs(title = "Leiden clusters: PRIMARY vs SENSITIVITY (Harmony)",
       subtitle = paste0("Each row = 100% of a PRIMARY cluster; ARI = ", round(ari, 3)),
       x = "SENSITIVITY cluster", y = "PRIMARY cluster", fill = "% of\nPRIMARY\ncluster")
save_tiff(p_ct, file.path(dir_compare, "crosstab_heatmap_PRIMARY_vs_SENSITIVITY.tiff"), 12, 10)

## section mixing, both analyses side by side
mix_p <- fread(file.path(dir_primary, "section_mixing_on_kNN_graph.csv"))
mix_s <- fread(file.path(dir_sens,    "section_mixing_on_kNN_graph.csv"))
mix_both <- merge(mix_p[, .(section_id, primary_same_section_frac = mean_same_section_neighbour_fraction,
                            expected_if_randomly_mixed)],
                  mix_s[, .(section_id, sensitivity_same_section_frac = mean_same_section_neighbour_fraction)],
                  by = "section_id")
fwrite(mix_both, file.path(dir_compare, "section_mixing_PRIMARY_vs_SENSITIVITY.csv"))
print(mix_both)

## UMAPs side by side
save_tiff(DimPlot(obj, reduction = "umap_primary", group.by = paste0("primary_leiden_res_", leiden_resolution),
                  label = TRUE, cols = cluster_colours(levels(cl_p)), pt.size = 0.05, alpha = 0.3, raster = FALSE) +
            NoLegend() + coord_fixed() + labs(title = "PRIMARY") |
            DimPlot(obj, reduction = "umap_harmony", group.by = paste0("harmony_leiden_res_", leiden_resolution),
                    label = TRUE, cols = cluster_colours(levels(cl_s)), pt.size = 0.05, alpha = 0.3, raster = FALSE) +
            NoLegend() + coord_fixed() + labs(title = "SENSITIVITY (Harmony)"),
          file.path(dir_compare, "UMAP_PRIMARY_and_SENSITIVITY_side_by_side.tiff"), 20, 9)

## check: agreement with script 12's clusters (same method, different implementation)
chk <- rbindlist(lapply(list(c("primary", paste0("01_PRIMARY_batchawarePCA_", run_tag)),
                             c("harmony", paste0("02_SENSITIVITY_batchawarePCA_Harmony_", run_tag))), function(b) {
                               f12 <- file.path(dir_main, b[2], "cells_UMAP_and_Leiden_clusters.csv")   # script 12 output
                               if (!file.exists(f12)) return(NULL)
                               c12 <- fread(f12, select = c("cell_id", "leiden"))
                               c12 <- c12[match(colnames(obj), cell_id)]
                               data.table(analysis = b[1], resolution = leiden_resolution,
                                          ARI_Seurat_vs_script12 = round(igraph::compare(as.integer(obj[[paste0(b[1], "_leiden_res_", leiden_resolution)]][, 1]),
                                                                                         c12$leiden + 1L, method = "adjusted.rand"), 4))
                             }))
if (nrow(chk) > 0) { print(chk); fwrite(chk, file.path(dir_compare, "ARI_Seurat_functions_vs_script12.csv")) }


## =============================================================================
## CHUNK 9: save the Seurat object (reductions + all cluster columns)
## =============================================================================
## Reductions: pearson_pca, harmony, umap_primary, umap_harmony.
## Metadata: primary_/harmony_leiden_res_<r> for every sweep resolution.
## kNN and graphs were saved as RDS in each branch folder and removed from obj.
obj@graphs    <- list()
obj@neighbors <- list()
seurat_out <- file.path(dir_seu, paste0("seurat_", input_tag, "_SeuratFunctions_", run_tag, ".RDS"))
saveRDS(obj, seurat_out)
writeLines(capture.output(sessionInfo()), file.path(dir_seu, "sessionInfo.txt"))
message("Saved: ", seurat_out, "\nAll outputs in: ", dir_seu)

## summary (as script 12 CHUNK 10)
message("\n", strrep("=", 78))
message("DONE. Everything is in: ", dir_seu)
message("  00_shared  : hvgs, PCA, PCA_variance_per_PC.csv, PCA_elbow_", run_tag, ".tiff")
message("  01_PRIMARY : FindNeighbors kNN (STEP A) -> graph (STEP B) -> RunUMAP (STEP C) + FindClusters (STEP D)")
message("               ", nlevels(cl_p), " clusters at resolution ", leiden_resolution)
message("  02_SENSITIVITY (Harmony): same steps; ", nlevels(cl_s), " clusters")
message("  01_ and 02_ each also hold the resolution sweep: Leiden_resolution_sweep_same_graph.csv +")
message("               _summary.tiff, _UMAP_panels.tiff, _cluster_tree.tiff (+ tree nodes/edges CSVs)")
message("  03_comparison: ARI = ", round(ari, 4), ", crosstab + section mixing")
message("  parameters_used_", run_tag, ".csv + sessionInfo.txt for the methods section")
message(strrep("=", 78))


## =============================================================================
## CHUNK 10: side-by-side UMAPs (Leiden | counts) and (Leiden | tissue section) - as script 13
## =============================================================================
tissue_cols <- c(AAA2 = "#2CA02C", AAA38 = "#1F4EDB", AAA40 = "#E31A1C",   # same section colours as script 13
                 AAA45 = "#FF8C00", AAA46 = "#8E44AD", AAA48 = "#E6B800")
stopifnot("section names differ from tissue_cols" = setequal(names(tissue_cols), sections))
n_sec <- table(sec)[names(tissue_cols)]                                    # cells per section, for the legend
tissue_labels <- setNames(paste0(names(tissue_cols), " (", format(as.integer(n_sec), big.mark = ","), " cells)"),
                          names(tissue_cols))
plot_tag <- paste0(run_tag, "_res", leiden_resolution)                     # e.g. "PC20_NN250_res1"

branches_sbs <- list(
  list(umap = "umap_primary", col = paste0("primary_leiden_res_", leiden_resolution),
       dir = dir_primary, label = "Batch-aware Pearson PCA (PRIMARY)"),
  list(umap = "umap_harmony", col = paste0("harmony_leiden_res_", leiden_resolution),
       dir = dir_sens,    label = "Batch-aware Pearson PCA + Harmony (SENSITIVITY)"))

for (b in branches_sbs) {
  dir_counts <- file.path(b$dir, "UMAP_leiden_nCount_RNA")                 # output sub-folders (as script 13)
  dir_tissue <- file.path(b$dir, "UMAP_leiden_color_cells_by_tissue")
  dir.create(dir_counts, showWarnings = FALSE); dir.create(dir_tissue, showWarnings = FALSE)
  sub_txt <- paste0(b$label, ", ", n_pcs, " PCs, ", n_neighbors, " neighbours")
  cl <- obj[[b$col]][, 1]                                                  # cluster of each cell
  
  ## 1. Leiden clusters
  p_leiden <- DimPlot(obj, reduction = b$umap, group.by = b$col, cols = cluster_colours(levels(cl)),
                      label = TRUE, label.box = TRUE, shuffle = TRUE, seed = seed,
                      pt.size = 0.05, alpha = 0.3, raster = FALSE) +
    coord_fixed() + guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
    labs(title = paste0("UMAP - Leiden clusters (resolution ", leiden_resolution, ", ", nlevels(cl), " clusters)"),
         subtitle = sub_txt, colour = "Leiden\ncluster")
  
  ## 2. total RNA counts per cell; log10 colour scale, legend in real counts
  p_counts <- FeaturePlot(obj, features = "nCount_RNA", reduction = b$umap, pt.size = 0.05, raster = FALSE) +
    scale_colour_viridis_c(option = "B", trans = "log10", labels = scales::label_comma()) +
    coord_fixed() + guides(colour = guide_colourbar(barheight = unit(8, "cm"))) +
    labs(title = "UMAP - total RNA counts per cell (nCount_RNA)",
         subtitle = paste0(sub_txt, "; median ", format(median(obj$nCount_RNA), big.mark = ","), " counts per cell"),
         colour = "Total RNA counts\nper cell\n(nCount_RNA,\nlog10 colour scale)")
  
  ## 3. tissue section of each cell (legend shows cell numbers)
  p_tissue <- DimPlot(obj, reduction = b$umap, group.by = batch_col, shuffle = TRUE, seed = seed,
                      pt.size = 0.05, alpha = 0.3, raster = FALSE) +
    scale_colour_manual(values = tissue_cols, labels = tissue_labels) +
    coord_fixed() + guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) +
    labs(title = "UMAP - cells coloured by tissue section", subtitle = sub_txt, colour = "Tissue section")
  
  ## single plots + side by side (patchwork |)
  save_tiff(p_leiden, file.path(dir_counts, paste0("UMAP_Leiden_", plot_tag, ".tiff")), 11, 9)
  save_tiff(p_counts, file.path(dir_counts, paste0("UMAP_nCount_RNA_", plot_tag, ".tiff")), 11, 9)
  save_tiff((p_leiden | p_counts) +
              plot_annotation(title = paste0("Leiden clusters vs total RNA counts - ", b$label, " (", plot_tag, ")"),
                              theme = theme(plot.title = element_text(face = "bold", size = 16))),
            file.path(dir_counts, paste0("UMAP_Leiden_and_nCount_RNA_side_by_side_", plot_tag, ".tiff")), 22, 10)
  
  save_tiff(p_leiden, file.path(dir_tissue, paste0("UMAP_Leiden_", plot_tag, ".tiff")), 11, 9)
  save_tiff(p_tissue, file.path(dir_tissue, paste0("UMAP_by_tissue_section_", plot_tag, ".tiff")), 11, 9)
  save_tiff((p_leiden | p_tissue) +
              plot_annotation(title = paste0("Leiden clusters vs tissue section - ", b$label, " (", plot_tag, ")"),
                              theme = theme(plot.title = element_text(face = "bold", size = 16))),
            file.path(dir_tissue, paste0("UMAP_Leiden_and_tissue_section_side_by_side_", plot_tag, ".tiff")), 22, 10)
}


## =============================================================================
## CHUNK 11: cluster markers - FindAllMarkers (main resolution)
## =============================================================================
## Each Leiden cluster vs all other cells, genes HIGHER in the cluster, Wilcoxon
## (Seurat vignette: FindAllMarkers(only.pos = TRUE)). Clusters NOT recomputed.
## FindAllMarkers reads the RNA "data" layer, which still holds raw counts
## (Pearson residuals were used only inside the PCA) -> NormalizeData first
## (log1p(counts / total x 10,000)); this changes neither clusters nor UMAP.
## SIGNIFICANT = p_val_adj < 0.05 (Bonferroni). RANKING = pct.1 high -> low, ties by avg_log2FC.
## Faster with presto: remotes::install_github("immunogenomics/presto")
library(openxlsx)     # write.xlsx(): one Excel sheet per cluster

## fresh session: run CHUNKS 0, 1, 4 first; this reloads the CHUNK 9 object
seurat_out <- file.path(dir_seu, paste0("seurat_", input_tag, "_SeuratFunctions_", run_tag, ".RDS"))
if (!exists("obj") || !"umap_primary" %in% Reductions(obj)) obj <- readRDS(seurat_out)
plot_tag <- paste0(run_tag, "_res", leiden_resolution)

marker_padj_max <- 0.05   # significance cut-off on adjusted p-value
marker_top_n    <- 20     # genes per cluster sheet

obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)

branches_mk <- list(
  list(prefix = "primary", dir = dir_primary, label = "PRIMARY (batch-aware Pearson PCA)"),
  list(prefix = "harmony", dir = dir_sens,    label = "SENSITIVITY (batch-aware Pearson PCA + Harmony)"))

for (b in branches_mk) {
  main_col <- paste0(b$prefix, "_leiden_res_", leiden_resolution)          # e.g. "primary_leiden_res_1"
  stopifnot("cluster column missing - run chunks 6-7 first" = main_col %in% colnames(obj[[]]))
  dir_mk <- file.path(b$dir, paste0("FindAllMarkers_res", leiden_resolution))
  dir.create(dir_mk, showWarnings = FALSE)
  file_tag <- paste0(b$prefix, "_", plot_tag)                               # e.g. "primary_PC20_NN250_res1"
  
  Idents(obj) <- main_col                                                   # levels already 0, 1, 2, ...
  message("\n", b$label, ": ", nlevels(Idents(obj)), " clusters - running FindAllMarkers...")
  
  markers <- as.data.table(FindAllMarkers(obj, assay = "RNA", only.pos = TRUE, test.use = "wilcox"))
  saveRDS(markers, file.path(dir_mk, paste0("FindAllMarkers_", file_tag, "_all_results.RDS")))
  
  ## significant genes, ranked within each cluster
  sig <- markers[p_val_adj < marker_padj_max & avg_log2FC > 0]
  setorder(sig, cluster, -pct.1, -avg_log2FC)
  sig[, rank_in_cluster := seq_len(.N), by = cluster]
  setcolorder(sig, c("cluster", "rank_in_cluster", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"))
  fwrite(sig, file.path(dir_mk, paste0("FindAllMarkers_", file_tag, "_all_significant_ranked.csv")))
  
  ## Excel: one sheet per cluster, top 20 genes
  cl_levels <- levels(Idents(obj))
  sheets <- lapply(cl_levels, function(k) {
    top <- sig[cluster == k & rank_in_cluster <= marker_top_n]
    if (nrow(top) == 0) data.frame(note = "no significant positive markers") else as.data.frame(top)
  })
  names(sheets) <- paste0("cluster_", cl_levels)
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