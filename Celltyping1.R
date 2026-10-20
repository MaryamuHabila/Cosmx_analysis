## install.packages("reticulate")
library(reticulate)
## =============================================================================
## 14_FindAllMarkers_from_clustering_path.R
## =============================================================================
## PURPOSE
##   FindAllMarkers on the Leiden clusters (main resolution) of ONE clustering
##   branch saved by 12b. Set ONE path in CHUNK 1 (the branch folder, or any file
##   inside it, e.g. its UMAP_by_Leiden_cluster.tiff); everything else (QC choice,
##   PC/NN run, PRIMARY vs SENSITIVITY, Seurat object, resolution) is read from it.
##   Nothing is recomputed (no PCA, neighbours, UMAP or clustering).
## OUTPUT: <output_dir>/Results/<QC choice>_final_All_markers/
##   e.g. .../AAA/Results/QC03_04_05_06_07_08_09_final_All_markers/
##   File names carry branch + run + resolution, so other runs do not overwrite.
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 1: the ONLY thing to change - path to the clustering results of choice
## =============================================================================
library(data.table)   # tables
library(Seurat)       # NormalizeData, FindAllMarkers
library(openxlsx)     # write.xlsx(): one Excel sheet per cluster


clustering_path <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Clustering_PCA_UMAP_Leiden_Markers/QC03_04_05_06_07_08_09/PC20_NN250_rep3/02_SENSITIVITY_batchawarePCA_Harmony/UMAP_by_Leiden_cluster.tiff"   # <- change here
##clustering_path <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Clustering_PCA_UMAP_Leiden_Markers/QC03_04_06_07_08_09/PC20_NN250_rep3/02_SENSITIVITY_batchawarePCA_Harmony/UMAP_by_Leiden_cluster.tiff"   # <- change here

leiden_resolution <- NULL   # NULL = main resolution used in 12b (from its parameters CSV); or e.g. 0.5 (must be in 12b's sweep)
marker_padj_max   <- 0.05   # significance cut-off on adjusted p-value
marker_top_n      <- 20     # genes per cluster sheet


## =============================================================================
## CHUNK 2: everything else is derived from the path (12b folder layout)
## =============================================================================
## 12b layout: <output_dir>/Clustering_PCA_UMAP_Leiden_Markers/<qc_tag>/<run_tag>/<branch folder>/
branch_dir <- if (dir.exists(clustering_path)) clustering_path else dirname(clustering_path)   # file given -> its folder
branch_dir <- normalizePath(branch_dir, mustWork = TRUE)                    # stops if the folder does not exist
dir_seu    <- dirname(branch_dir)                                           # .../PC20_NN250_rep3
run_tag    <- basename(dir_seu)                                             # "PC20_NN250_rep3"
qc_tag     <- basename(dirname(dir_seu))                                    # "QC03_04_05_06_07_08_09"
stopifnot("path is not inside Clustering_PCA_UMAP_Leiden_Markers/<QC>/<run>/<branch>/" =
            basename(dirname(dirname(dir_seu))) == "Clustering_PCA_UMAP_Leiden_Markers")
output_dir <- dirname(dirname(dirname(dir_seu)))                            # .../AAA

## PRIMARY or SENSITIVITY branch -> prefix of the Leiden columns / UMAP in the object (as 12b)
branch_name <- basename(branch_dir)
prefix <- if (grepl("^01_PRIMARY", branch_name)) "primary" else
  if (grepl("^02_SENSITIVITY", branch_name)) "harmony" else
    stop("branch folder must start with 01_PRIMARY or 02_SENSITIVITY: ", branch_name)
branch_label <- if (prefix == "primary") "Batch-aware Pearson PCA (PRIMARY)" else
  "Batch-aware Pearson PCA + Harmony (SENSITIVITY)"

## main resolution: from 12b's parameters table unless set above
if (is.null(leiden_resolution)) {
  params <- fread(file.path(dir_seu, paste0("parameters_used_", run_tag, ".csv")))
  leiden_resolution <- as.numeric(params[parameter == "Leiden_resolution", value])
}

## Seurat object saved by 12b CHUNK 9 (exactly one per run folder)
seurat_out <- list.files(dir_seu, pattern = paste0("^seurat_.*_SeuratFunctions_", run_tag, "\\.RDS$"),
                         full.names = TRUE)
stopifnot("expected exactly one 12b Seurat object in the run folder" = length(seurat_out) == 1)

## Results folder named by the QC combination
dir_mk <- file.path(output_dir, "Results", paste0(qc_tag, "_final_All_markers"))
dir.create(dir_mk, showWarnings = FALSE, recursive = TRUE)
file_tag <- paste0(prefix, "_", run_tag, "_res", leiden_resolution)        # e.g. "harmony_PC20_NN250_rep3_res1"

message("QC choice     : ", qc_tag)
message("Run           : ", run_tag)
message("Branch        : ", branch_label)
message("Resolution    : ", leiden_resolution)
message("Seurat object : ", seurat_out)
message("Results to    : ", dir_mk)

## copy the UMAP of these clusters next to the markers (renamed so runs do not overwrite)
umap_tiff <- file.path(branch_dir, "UMAP_by_Leiden_cluster.tiff")
if (file.exists(umap_tiff)) {
  file.copy(umap_tiff, file.path(dir_mk, paste0("UMAP_by_Leiden_cluster_", file_tag, ".tiff")), overwrite = TRUE)
}


## =============================================================================
## CHUNK 3: load the object + checks
## =============================================================================
obj <- readRDS(seurat_out)                                                  # reductions + Leiden columns from 12b
main_col <- paste0(prefix, "_leiden_res_", leiden_resolution)              # e.g. "harmony_leiden_res_1"
stopifnot("resolution not in 12b's sweep" = main_col %in% colnames(obj[[]]))


## =============================================================================
## CHUNK 4: cluster markers - FindAllMarkers (main resolution)
## =============================================================================
## Each Leiden cluster vs all other cells, genes HIGHER in the cluster, Wilcoxon
## (Seurat vignette: FindAllMarkers(only.pos = TRUE)). Clusters NOT recomputed.
## FindAllMarkers reads the RNA "data" layer, which still holds raw counts
## -> NormalizeData first (log1p(counts / total x 10,000)); this changes
## neither clusters nor UMAP.
## SIGNIFICANT = p_val_adj < 0.05 (Bonferroni). RANKING = pct.1 high -> low, ties by avg_log2FC.
## Faster with presto: remotes::install_github("immunogenomics/presto")
obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)

Idents(obj) <- main_col                                                     # levels already 0, 1, 2, ... (by size)
message("\n", branch_label, ": ", nlevels(Idents(obj)), " clusters - running FindAllMarkers...")

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

## record of what was run (for the methods section)
fwrite(data.table(parameter = c("clustering_path", "qc_tag", "run_tag", "branch", "leiden_resolution",
                                "seurat_object", "test", "only.pos", "padj_cutoff", "top_n_per_sheet"),
                  value = c(branch_dir, qc_tag, run_tag, branch_label, leiden_resolution,
                            seurat_out, "wilcox", TRUE, marker_padj_max, marker_top_n)),
       file.path(dir_mk, paste0("FindAllMarkers_", file_tag, "_parameters.csv")))
message(branch_label, " markers saved -> ", dir_mk)



































































## =============================================================================
## 15_CellTypist_clusters.R

library(reticulate)                                              # R <-> Python
use_condaenv("celltypist_env", required = TRUE)                  # same env as script 15
ct_models <- import("celltypist.models", convert = FALSE)        # CellTypist model functions

## 1. every model available, with a one-line description (as in the vignette)
ct_models$download_models(force_update = FALSE)                  # fetch any models not yet downloaded (one-off, all models)
mods <- py_to_r(ct_models$models_description(on_the_fly = TRUE)) # table: model name + description
print(mods, right = FALSE)                                       # read this to pick candidate models

## 2. cell types inside the models you are considering
show_types <- c("Adult_Human_Vascular.pkl", "Immune_All_Low.pkl", "Healthy_Adult_Heart.pkl")  # <- add/remove names
for (f in show_types) {
  m   <- ct_models$Model$load(model = f)                         # load one model
  cts <- unlist(py_to_r(m$cell_types$tolist()))                  # its cell-type names
  cat("\n====", f, "-", length(cts), "cell types ====\n")
  print(cts)                                                     # one cell type per element
}




## =============================================================================
## =============================================================================
## 15_CellTypist_top20_top50.R
## =============================================================================
## WHAT IT DOES
##   Names each Leiden cluster with CellTypist's human reference models:
##     Adult_Human_Vascular.pkl (Barnett 2024, Nat Med; 42 endothelial + mural types)
##     Immune_All_Low.pkl       (Dominguez Conde 2022, Science; 98 immune types)
##     Healthy_Adult_Heart.pkl  (Kanemaru 2023, Nature; 75 adult heart types, incl. fibroblasts, SMC, adipocytes)
##   CellTypist classifies CELLS from expression; it cannot read a gene list.
##   So for each cluster, its marker genes are compared with the genes each
##   CellTypist cell type is built on (model.extract_top_markers, as in the
##   CellTypist vignette). The best-matching cell type is the call, if the
##   overlap is larger than chance (hypergeometric test, BH-adjusted p < 0.05).
##     top 20 markers  (sheets of the top-20 xlsx)  -> celltypist_top20.xlsx
##     top 50 markers  (rank_in_cluster <= 50, CSV) -> celltypist_top50.xlsx
##   Both files are saved next to the input files.
##
## ONE-TIME SETUP (Terminal, not R):
##   conda create -n celltypist_env -c conda-forge python=3.11 -y
##   conda activate celltypist_env
##   pip install celltypist
##   python -c "import celltypist; print(celltypist.__version__)"
## Then in R: install.packages("reticulate") if needed, and RESTART R
## (reticulate binds to one Python per R session).
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 1: settings - the two files from script 14
## =============================================================================
library(data.table)   # tables
library(openxlsx)     # read + write xlsx
library(reticulate)   # R <-> Python (CellTypist)

## xlsx_top20 <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_top20_per_cluster.xlsx"   # <- change here
xlsx_top20 <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_05_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_top20_per_cluster.xlsx"   # <- change here
## csv_ranked <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_all_significant_ranked.csv"  # <- change here
csv_ranked <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_05_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_all_significant_ranked.csv"  # <- change here


models_use <- c(Vascular = "Adult_Human_Vascular.pkl",   # vasculature (endothelial + mural)
                Immune   = "Immune_All_Low.pkl",         # immune (latest CellTypist immune model, v2)
                Heart    = "Healthy_Adult_Heart.pkl")    # adult cardiovascular stroma incl. fibroblasts, SMC, adipocytes
ref_n      <- 100        # genes per CellTypist cell type used for matching (its top driving genes on your panel)
padj_max   <- 0.05       # a call is made only if BH-adjusted p < this
conda_env  <- "celltypist_env"

stopifnot("top-20 xlsx not found" = file.exists(xlsx_top20),
          "ranked CSV not found"  = file.exists(csv_ranked))
out_dir   <- dirname(xlsx_top20)                                            # save next to the inputs
rds_all   <- sub("_top20_per_cluster\\.xlsx$", "_all_results.RDS", xlsx_top20)   # script 14's full results
stopifnot("script 14's _all_results.RDS not found next to the xlsx" = file.exists(rds_all))


## =============================================================================
## CHUNK 2: CellTypist models (download latest, load)
## =============================================================================
use_condaenv(conda_env, required = TRUE)                                    # must come before any Python call
ct_models <- import("celltypist.models", convert = FALSE)                   # convert = FALSE: convert explicitly
ct_models$download_models(force_update = TRUE,                              # fetch latest versions of the models in models_use only
                          model = r_to_py(as.list(unname(models_use))))
models <- lapply(models_use, function(f) ct_models$Model$load(model = f))   # as in the vignette
model_types <- lapply(models, function(m) unlist(py_to_r(m$cell_types$tolist())))   # cell-type names
model_genes <- lapply(models, function(m) unlist(py_to_r(m$features$tolist())))     # genes each model uses
for (nm in names(models)) message(nm, ": ", length(model_types[[nm]]), " cell types, ",
                                  length(model_genes[[nm]]), " genes")


## =============================================================================
## CHUNK 3: your marker lists (top 20 from xlsx, top 50 from CSV)
## =============================================================================
## CosMx combined probes ("IGHG1/2", "IGLC1/2", "XCL1/2", "IFNA4/10/17") are split
## into separate genes so they can match model genes ("IGHG1", "IGHG2", ...).
split_probe <- function(g) {
  parts <- strsplit(g, "/", fixed = TRUE)[[1]]                              # "IGHG1/2" -> "IGHG1", "2"
  if (length(parts) == 1) return(g)                                         # ordinary gene: unchanged
  stem <- sub("[0-9]+$", "", parts[1])                                      # "IGHG1" -> "IGHG"
  c(parts[1], paste0(stem, parts[-1]))                                      # "IGHG1", "IGHG2"
}
split_genes <- function(gs) unique(unlist(lapply(gs, split_probe)))

## genes that could appear in a marker list = all genes FindAllMarkers tested (script 14)
tested_genes <- split_genes(unique(as.character(readRDS(rds_all)$gene)))

## top 20: one sheet per cluster ("cluster_0", "cluster_1", ...)
sheets <- getSheetNames(xlsx_top20)
top20 <- lapply(sheets, function(s) {
  d <- read.xlsx(xlsx_top20, sheet = s)
  if ("gene" %in% names(d)) d$gene else character(0)                        # sheet may say "no significant markers"
})
names(top20) <- sub("^cluster_", "", sheets)                                # "0", "1", ...

## top 50: ranks 1-50 of each cluster in the ranked CSV (same ranking as the xlsx)
ranked <- fread(csv_ranked)
ranked[, cluster := as.character(cluster)]
top50 <- lapply(names(top20), function(k) ranked[cluster == k & rank_in_cluster <= 50, gene])
names(top50) <- names(top20)

## check 1: xlsx top 20 = CSV ranks 1-20 (both from the same script-14 run)
for (k in names(top20)) {
  stopifnot("xlsx and CSV disagree - are they from the same run?" =
              setequal(top20[[k]], ranked[cluster == k & rank_in_cluster <= 20, gene]))
}
message("Check passed: ", length(top20), " clusters; xlsx top 20 = CSV ranks 1-20")
print(data.table(cluster = names(top20), n_top20 = lengths(top20), n_top50 = lengths(top50)))


## =============================================================================
## CHUNK 4: reference genes of every CellTypist cell type (on your panel)
## =============================================================================
## For each model: universe = genes in the model AND tested on your panel.
## Reference genes of a cell type = its top `ref_n` positive driving genes within
## that universe (model.extract_top_markers, ranked by the model's coefficients).
universe <- lapply(model_genes, function(g) intersect(g, tested_genes))
for (nm in names(universe)) message(nm, ": ", length(universe[[nm]]), " model genes are on your panel")

ref_genes <- lapply(names(models), function(nm) {
  m <- models[[nm]]
  all_n <- as.integer(length(model_genes[[nm]]))                            # ask for the full ranking
  setNames(lapply(model_types[[nm]], function(ct) {
    ranked_ct <- unlist(py_to_r(m$extract_top_markers(ct, top_n = all_n, only_positive = TRUE)$tolist()))
    head(ranked_ct[ranked_ct %in% universe[[nm]]], ref_n)                   # top ref_n that are on the panel
  }), model_types[[nm]])
})
names(ref_genes) <- names(models)


## =============================================================================
## CHUNK 5: match each cluster's markers to each cell type (hypergeometric)
## =============================================================================
## For one model and one cluster:
##   n = cluster markers that are in the universe
##   r = reference genes of the cell type, x = markers shared with them
##   p = chance of sharing >= x genes at random; BH across the model's cell types
match_clusters <- function(marker_lists, nm) {
  U <- universe[[nm]]; N <- length(U)
  rbindlist(lapply(names(marker_lists), function(k) {
    q <- intersect(split_genes(marker_lists[[k]]), U)                       # this cluster's usable markers
    res <- rbindlist(lapply(names(ref_genes[[nm]]), function(ct) {
      ref <- ref_genes[[nm]][[ct]]
      hit <- intersect(q, ref)
      data.table(cluster = k, cell_type = ct,
                 n_markers = length(marker_lists[[k]]),                     # markers given
                 n_markers_in_model = length(q),                            # of those, known to the model
                 n_shared = length(hit),
                 p_value = if (length(q) > 0 && length(ref) > 0)
                   phyper(length(hit) - 1, length(ref), N - length(ref), length(q), lower.tail = FALSE) else NA_real_,
                 shared_genes = paste(hit, collapse = ", "))
    }))
    res[, p_adj := p.adjust(p_value, method = "BH")]                        # per cluster, across cell types
    res[order(p_value, -n_shared, na.last = TRUE)]
  }))
}

## one row per cluster: best match, its genes, and the next two candidates
best_per_cluster <- function(all, nm) {
  out <- all[, .(n_markers = n_markers[1], n_markers_in_model = n_markers_in_model[1],
                 best_match = cell_type[1], n_shared = n_shared[1], shared_genes = shared_genes[1],
                 p_adj = signif(p_adj[1], 3),
                 call = if (!is.na(p_adj[1]) && p_adj[1] < padj_max) cell_type[1] else "no significant match",
                 second = paste0(cell_type[2], " (", n_shared[2], " shared, padj ", signif(p_adj[2], 2), ")"),
                 third  = paste0(cell_type[3], " (", n_shared[3], " shared, padj ", signif(p_adj[3], 2), ")")),
             by = cluster]
  out[order(as.numeric(cluster))]                                           # 0, 1, 2, ... (not 0, 1, 10)
}

run_and_save <- function(marker_lists, label, out_file) {
  per_model <- lapply(names(models), function(nm) best_per_cluster(match_clusters(marker_lists, nm), nm))
  names(per_model) <- names(models)
  
  ## Summary: one row per cluster, the call from each model side by side
  smry <- data.table(cluster = per_model[[1]]$cluster, n_markers = per_model[[1]]$n_markers)
  for (nm in names(per_model)) {
    p <- per_model[[nm]]
    smry[, (paste0(nm, "_call"))         := p$call]
    smry[, (paste0(nm, "_shared_genes")) := ifelse(p$call == "no significant match", "", p$shared_genes)]
    smry[, (paste0(nm, "_padj"))         := p$p_adj]
  }
  print(smry[, c("cluster", paste0(names(models), "_call")), with = FALSE])
  
  how_to_read <- data.frame(item = c(
    "Markers used", "Models", "Method", "call", "best_match", "n_markers_in_model",
    "n_shared / shared_genes", "p_adj", "no significant match"),
    meaning = c(
      label,
      paste(paste0(names(models_use), " = ", models_use), collapse = "; "),
      paste0("Each cluster's markers vs each CellTypist cell type's top ", ref_n,
             " driving genes (model.extract_top_markers) on the panel; hypergeometric test, BH across cell types"),
      paste0("best_match if p_adj < ", padj_max, ", otherwise 'no significant match'"),
      "cell type sharing the most markers beyond chance (reported even when not significant)",
      "how many of the cluster's markers the model knows (only these can match)",
      "markers shared with the best-matching cell type",
      "BH-adjusted p-value of the overlap",
      "markers do not resemble any cell type in that model (or too few markers are known to it)"))
  
  sheets_out <- c(list(Summary = as.data.frame(smry)),
                  lapply(per_model, as.data.frame),                         # one detailed sheet per model
                  list(How_to_read = how_to_read))
  write.xlsx(sheets_out, out_file, firstRow = TRUE, colWidths = "auto", overwrite = TRUE)
  message("Saved: ", out_file)
}


## =============================================================================
## CHUNK 6: run - top 20 and top 50
## =============================================================================
run_and_save(top20, "top 20 per cluster (xlsx sheets; ranked by pct.1, then avg_log2FC)",
             file.path(out_dir, "celltypist_top20.xlsx"))
run_and_save(top50, "top 50 per cluster (ranked CSV, rank_in_cluster <= 50)",
             file.path(out_dir, "celltypist_top50.xlsx"))










## 1. best match per cluster: genes shared vs genes expected by chance
for (nm in names(models)) {
  a <- match_clusters(top20, nm)                                    # every cluster x cell-type test (top 20)
  b <- a[, .SD[1], by = cluster]                                    # best match per cluster (already sorted by p)
  b[, ref_size := sapply(cell_type, function(ct) length(ref_genes[[nm]][[ct]]))]   # reference genes of that type
  b[, expected_by_chance := round(n_markers_in_model * ref_size / length(universe[[nm]]), 2)]  # n x r / N
  cat("\n====", nm, "====\n")
  print(b[order(as.numeric(cluster)),
          .(cluster, best = cell_type, n_markers_in_model, n_shared, expected_by_chance,
            p_adj = signif(p_adj, 2), shared_genes)])
}

## 2. cluster 7 (immunoglobulin genes): why no plasma-cell match?
g7 <- split_genes(top20[["7"]])                                     # cluster 7 markers, probes split
print(data.table(gene = g7,
                 in_Immune_model = g7 %in% model_genes$Immune,      # does the model use this gene at all?
                 on_panel_and_model = g7 %in% universe$Immune,      # usable for matching?
                 in_Plasma_ref = g7 %in% ref_genes$Immune[["Plasma cells"]],   # among plasma reference genes?
                 in_Heart_B_plasma_ref = g7 %in% ref_genes$Heart[["B_plasma"]]))
print(head(ref_genes$Immune[["Plasma cells"]], 30))                 # what the model's plasma genes actually are








## =============================================================================
## 16_CellTypist_per_cell_majority_voting.R
## =============================================================================
## WHAT IT DOES (the CellTypist vignette method, on your own Seurat object)
##   1. log1p-normalise every cell to 10,000 counts (what CellTypist requires)
##   2. celltypist.annotate() labels EVERY CELL using all panel genes the model knows
##   3. majority_voting = TRUE with over_clustering = your Leiden clusters
##      -> one label per Leiden cluster (the label most of its cells received)
##   Models: Vascular (Adult_Human_Vascular), Immune (Immune_All_Low), Heart (Healthy_Adult_Heart)
## OUTPUT (same Results folder as script 14):
##   celltypist_majority_voting_<tag>.xlsx  - Summary sheet + one sheet per model
##   celltypist_per_cell_<model>_<tag>.csv.gz - the label of every cell
## SETUP: same conda env as script 15 (celltypist_env). RESTART R before running.
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 1: settings - same path you give script 14
## =============================================================================
library(data.table)   # tables
library(Seurat)       # NormalizeData, LayerData
library(Matrix)       # sparse matrices
library(openxlsx)     # write xlsx
library(reticulate)   # R <-> Python (CellTypist)

clustering_path <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Clustering_PCA_UMAP_Leiden_Markers/QC03_04_05_06_07_08_09/PC20_NN250_rep3/02_SENSITIVITY_batchawarePCA_Harmony/UMAP_by_Leiden_cluster.tiff"   # <- change here
leiden_resolution <- NULL   # NULL = main resolution of 12b (as script 14)

models_use <- c(Vascular = "Adult_Human_Vascular.pkl",   # vasculature (endothelial + mural)
                Immune   = "Immune_All_Low.pkl",         # immune
                Heart    = "Healthy_Adult_Heart.pkl")    # adult cardiovascular stroma incl. fibroblasts, SMC, adipocytes
conda_env <- "celltypist_env"


## =============================================================================
## CHUNK 2: paths from clustering_path (identical to script 14) + load object
## =============================================================================
branch_dir <- if (dir.exists(clustering_path)) clustering_path else dirname(clustering_path)  # file -> its folder
branch_dir <- normalizePath(branch_dir, mustWork = TRUE)                    # stop if missing
dir_seu    <- dirname(branch_dir)                                           # .../PC20_NN250_rep3
run_tag    <- basename(dir_seu)                                             # "PC20_NN250_rep3"
qc_tag     <- basename(dirname(dir_seu))                                    # "QC03_04_06_07_08_09"
stopifnot("path is not inside Clustering_PCA_UMAP_Leiden_Markers/<QC>/<run>/<branch>/" =
            basename(dirname(dirname(dir_seu))) == "Clustering_PCA_UMAP_Leiden_Markers")
output_dir <- dirname(dirname(dirname(dir_seu)))                            # .../AAA
prefix <- if (grepl("^01_PRIMARY", basename(branch_dir))) "primary" else
  if (grepl("^02_SENSITIVITY", basename(branch_dir))) "harmony" else
    stop("branch folder must start with 01_PRIMARY or 02_SENSITIVITY")
if (is.null(leiden_resolution)) {                                           # main resolution from 12b
  params <- fread(file.path(dir_seu, paste0("parameters_used_", run_tag, ".csv")))
  leiden_resolution <- as.numeric(params[parameter == "Leiden_resolution", value])
}
file_tag <- paste0(prefix, "_", run_tag, "_res", leiden_resolution)        # "harmony_PC20_NN250_rep3_res1"
dir_mk   <- file.path(output_dir, "Results", paste0(qc_tag, "_final_All_markers"))   # script 14's folder
dir.create(dir_mk, showWarnings = FALSE, recursive = TRUE)

seurat_out <- list.files(dir_seu, pattern = paste0("^seurat_.*_SeuratFunctions_", run_tag, "\\.RDS$"),
                         full.names = TRUE)
stopifnot("expected exactly one 12b Seurat object" = length(seurat_out) == 1)
obj <- readRDS(seurat_out)                                                  # 12b object
main_col <- paste0(prefix, "_leiden_res_", leiden_resolution)              # "harmony_leiden_res_1"
stopifnot("resolution not in 12b's sweep" = main_col %in% colnames(obj[[]]))
clusters <- as.character(obj[[main_col]][, 1])                             # Leiden cluster of each cell
message("Object: ", ncol(obj), " cells x ", nrow(obj), " genes; ", length(unique(clusters)), " clusters")


## =============================================================================
## CHUNK 3: expression in the form CellTypist needs (vignette: log1p, 10,000 per cell)
## =============================================================================
obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize",
                     scale.factor = 10000)                                  # log1p(counts / total x 10,000)
X <- t(LayerData(obj, assay = "RNA", layer = "data"))                       # cells x genes (AnnData orientation)

## check (as the vignette's adata.X.expm1().sum(axis = 1)): every cell sums to 10,000
chk <- Matrix::rowSums(expm1(X[seq_len(min(1000, nrow(X))), ]))             # first 1,000 cells
stopifnot("expression is not log1p-normalised to 10,000" = all(abs(chk - 10000) < 1))
message("Check passed: expm1(expression) sums to 10,000 per cell")


## =============================================================================
## CHUNK 4: build the AnnData (cells x genes + your Leiden clusters)
## =============================================================================
use_condaenv(conda_env, required = TRUE)                                    # must come before any Python call
ct        <- import("celltypist", convert = FALSE)                          # celltypist.annotate
ct_models <- import("celltypist.models", convert = FALSE)                   # Model.load, download_models
anndata   <- import("anndata", convert = FALSE)
pd        <- import("pandas", convert = FALSE)
message("CellTypist version: ", py_to_r(ct$`__version__`))
ct_models$download_models(force_update = FALSE,                             # only downloads models not yet present
                          model = r_to_py(as.list(unname(models_use))))

obs <- pd$DataFrame(data = r_to_py(list(leiden = clusters)),                # Leiden cluster per cell
                    index = r_to_py(colnames(obj)))                         # cell ids
var <- pd$DataFrame(index = r_to_py(rownames(obj)))                         # gene names
ad  <- anndata$AnnData(X = r_to_py(X)$tocsr(), obs = obs, var = var)        # sparse matrix -> AnnData
rm(X); invisible(gc())
message("AnnData: ", py_to_r(ad$n_obs), " cells x ", py_to_r(ad$n_vars), " genes")


## =============================================================================
## CHUNK 5: annotate every cell, majority vote within each Leiden cluster
## =============================================================================
per_model <- list()                                                         # one cluster table per model
for (nm in names(models_use)) {
  f <- models_use[[nm]]
  n_genes <- sum(unlist(py_to_r(ct_models$Model$load(model = f)$features$tolist())) %in% rownames(obj))
  message("\n", nm, " (", f, "): ", n_genes, " of its genes are on your panel. Annotating...")
  
  pred <- ct$annotate(ad, model = f, majority_voting = TRUE,               # vignette call, but votes are taken
                      over_clustering = "leiden")                           # within YOUR Leiden clusters
  pl   <- py_to_r(pred$predicted_labels)                                    # predicted_labels, over_clustering, majority_voting
  conf <- as.numeric(py_to_r(pred$probability_matrix$max(axis = 1L)$to_numpy()))   # = vignette's conf_score
  stopifnot("cell order changed" = identical(rownames(pl), colnames(obj)))
  
  cells <- data.table(cell_id = colnames(obj), cluster = clusters,
                      predicted = as.character(pl$predicted_labels),        # label of the cell on its own
                      majority_vote = as.character(pl$majority_voting),     # label of its Leiden cluster
                      conf_score = round(conf, 4))
  fwrite(cells, file.path(dir_mk, paste0("celltypist_per_cell_", nm, "_", file_tag, ".csv.gz")))
  
  ## per cluster: the majority label, how many cells really got it, confidence, runners-up
  tab <- cells[, {
    fr <- sort(table(predicted), decreasing = TRUE) / .N                    # fraction of cells per label
    mv <- majority_vote[1]                                                  # same for all cells of the cluster
    .(n_cells = .N,
      majority_label = mv,
      pct_cells_agreeing = round(100 * mean(predicted == mv), 1),           # % cells individually given that label
      median_conf_agreeing = round(median(conf_score[predicted == mv]), 3), # their median confidence
      top3_labels = paste0(names(fr)[1:min(3, length(fr))], " (",
                           round(100 * fr[1:min(3, length(fr))], 1), "%)", collapse = "; "))
  }, by = cluster][order(as.numeric(cluster))]
  tab[, genes_used := n_genes]
  print(tab[, .(cluster, n_cells, majority_label, pct_cells_agreeing, median_conf_agreeing)])
  per_model[[nm]] <- tab
}


## =============================================================================
## CHUNK 6: save - Summary (all models side by side) + one sheet per model
## =============================================================================
smry <- per_model[[1]][, .(cluster, n_cells)]
for (nm in names(per_model)) {
  p <- per_model[[nm]]
  smry[, (paste0(nm, "_label"))     := p$majority_label]
  smry[, (paste0(nm, "_pct_agree")) := p$pct_cells_agreeing]
  smry[, (paste0(nm, "_conf"))      := p$median_conf_agreeing]
}
print(smry)

how_to_read <- data.frame(item = c("Method", "label", "pct_agree", "conf", "top3_labels", "genes_used", "Caution"),
                          meaning = c(
                            "celltypist.annotate on log1p-normalised expression of every cell; majority_voting = TRUE with over_clustering = Leiden clusters",
                            "majority label of the cluster (the label most of its cells received)",
                            "% of the cluster's cells that were individually given that label (low = mixed or uncertain cluster)",
                            "median CellTypist confidence (0-1) of those agreeing cells",
                            "the three most frequent per-cell labels in the cluster",
                            "model genes present on the CosMx panel (the model was trained on whole-transcriptome data)",
                            "CellTypist always returns a label, even for cells unlike any reference type. Check labels with canonical markers; ignore cardiomyocyte/neural labels from the Heart model."))

out_xlsx <- file.path(dir_mk, paste0("celltypist_majority_voting_", file_tag, ".xlsx"))
write.xlsx(c(list(Summary = as.data.frame(smry)), lapply(per_model, as.data.frame),
             list(How_to_read = how_to_read)),
           out_xlsx, firstRow = TRUE, colWidths = "auto", overwrite = TRUE)
message("Saved: ", out_xlsx)




























install.packages('clusterProfiler')


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
## =============================================================================
## 17_CellMarker_top20_top50.R
## =============================================================================
## WHAT IT DOES
##   Names each Leiden cluster with CellMarker 2.0 (Hu et al. 2023, NAR),
##   the method shown in the clusterProfiler book (enricher + CellMarker):
##     each cluster's marker genes vs each CellMarker cell type's marker genes,
##     hypergeometric test (clusterProfiler::enricher), BH-adjusted per cluster.
##   CellMarker is restricted to HUMAN, NORMAL cells from VASCULAR tissues
##   (aorta / artery / vessel) + BLOOD (immune). Neural cell types are excluded.
##     top 20 markers (sheets of the top-20 xlsx)  -> cellmarker_top20.xlsx
##     top 50 markers (rank_in_cluster <= 50, CSV) -> cellmarker_top50.xlsx
##   Both saved next to the input files. No Python / reticulate needed.
##
## ONE-TIME SETUP (in R):
##   install.packages(c("BiocManager", "readxl", "openxlsx", "data.table"))
##   BiocManager::install("clusterProfiler")
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 1: settings
## =============================================================================
library(data.table)       # tables
library(openxlsx)         # read the top-20 sheets + write results
library(readxl)           # read the CellMarker xlsx
library(clusterProfiler)  # enricher(): hypergeometric test against CellMarker sets

## xlsx_top20 <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_top20_per_cluster.xlsx"   # <- change here
xlsx_top20 <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_05_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_top20_per_cluster.xlsx"   # <- change here
## csv_ranked <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_all_significant_ranked.csv"  # <- change here
csv_ranked <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_05_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_all_significant_ranked.csv"  # <- change here

## CellMarker 2.0 human file (URL from the clusterProfiler book); downloaded once into the Results folder
cm_file <- file.path(dirname(xlsx_top20), "Cell_marker_Human.xlsx")
cm_url  <- "http://www.bio-bigdata.center/CellMarker_download_files/file/Cell_marker_Human.xlsx"

tissue_regex  <- "aort|arter|vessel|vascul|blood"             # vascular tissues + blood (immune)   <- change here
exclude_regex <- "neur|glia|schwann|astrocyt|oligodendro|microglia"   # no neural cell types       <- change here
min_set  <- 3      # a CellMarker cell type is tested only if >= 3 of its markers are on your panel
padj_max <- 0.05   # a call is made only if BH-adjusted p < this

stopifnot("top-20 xlsx not found" = file.exists(xlsx_top20),
          "ranked CSV not found"  = file.exists(csv_ranked))
out_dir <- dirname(xlsx_top20)                                              # save next to the inputs
rds_all <- sub("_top20_per_cluster\\.xlsx$", "_all_results.RDS", xlsx_top20)   # script 14's full results
stopifnot("script 14's _all_results.RDS not found next to the xlsx" = file.exists(rds_all))


## =============================================================================
## CHUNK 2: CellMarker 2.0 - download once, keep vascular + blood, normal cells
## =============================================================================
if (!file.exists(cm_file)) {                                                # download only if not yet there
  ok <- tryCatch({ download.file(cm_url, cm_file, mode = "wb"); TRUE }, error = function(e) FALSE)
  if (!ok || file.size(cm_file) < 1e5) {
    stop("Download failed. Download 'Cell_marker_Human.xlsx' by hand from the CellMarker 2.0 Download page\n",
         "(http://bio-bigdata.hrbmu.edu.cn/CellMarker or http://117.50.127.228/CellMarker/) and save it as:\n", cm_file)
  }
}
cm <- as.data.table(read_excel(cm_file))
print(names(cm))                                                            # columns of the file
gene_col <- if ("Symbol" %in% names(cm)) "Symbol" else "marker"             # official gene symbol if present
stopifnot("unexpected CellMarker columns" = all(c("tissue_class", "tissue_type", "cell_name", gene_col) %in% names(cm)))

if ("cell_type" %in% names(cm)) {                                           # "Normal cell" vs "Cancer cell"
  print(table(cm$cell_type))
  cm <- cm[cell_type == "Normal cell"]                                      # normal cells only
}
cm_vasc <- cm[grepl(tissue_regex, tissue_class, ignore.case = TRUE) |
                grepl(tissue_regex, tissue_type,  ignore.case = TRUE)]        # vascular + blood entries
cm_vasc <- cm_vasc[!grepl(exclude_regex, cell_name, ignore.case = TRUE)]    # drop neural cell types

## CHECK: which tissues were kept (edit tissue_regex in CHUNK 1 if something is wrong)
print(cm_vasc[, .(n_entries = .N), by = .(tissue_class, tissue_type)][order(-n_entries)], nrows = 200)

## TERM2GENE: one row per (cell type, marker gene)
t2g <- unique(cm_vasc[!is.na(get(gene_col)) & get(gene_col) != "",
                      .(term = cell_name, gene = get(gene_col))])
message(length(unique(t2g$term)), " CellMarker cell types, ", length(unique(t2g$gene)), " marker genes kept")


## =============================================================================
## CHUNK 3: your marker lists (top 20 from xlsx, top 50 from CSV)
## =============================================================================
## CosMx combined probes ("IGHG1/2", "IGLC1/2") -> separate genes ("IGHG1", "IGHG2")
split_probe <- function(g) {
  parts <- strsplit(g, "/", fixed = TRUE)[[1]]                              # "IGHG1/2" -> "IGHG1", "2"
  if (length(parts) == 1) return(g)                                         # ordinary gene: unchanged
  stem <- sub("[0-9]+$", "", parts[1])                                      # "IGHG1" -> "IGHG"
  c(parts[1], paste0(stem, parts[-1]))                                      # "IGHG1", "IGHG2"
}
split_genes <- function(gs) unique(unlist(lapply(gs, split_probe)))

## background = every gene FindAllMarkers tested on your panel (script 14)
universe <- split_genes(unique(as.character(readRDS(rds_all)$gene)))
message(length(universe), " panel genes in the background; ",
        length(intersect(unique(t2g$gene), universe)), " CellMarker genes are on your panel")

sheets <- getSheetNames(xlsx_top20)                                         # "cluster_0", "cluster_1", ...
top20 <- lapply(sheets, function(s) {
  d <- read.xlsx(xlsx_top20, sheet = s)
  if ("gene" %in% names(d)) d$gene else character(0)                        # sheet may say "no significant markers"
})
names(top20) <- sub("^cluster_", "", sheets)

ranked <- fread(csv_ranked)
ranked[, cluster := as.character(cluster)]
top50 <- lapply(names(top20), function(k) ranked[cluster == k & rank_in_cluster <= 50, gene])
names(top50) <- names(top20)

## check: xlsx top 20 = CSV ranks 1-20 (same script-14 run)
for (k in names(top20)) {
  stopifnot("xlsx and CSV disagree - are they from the same run?" =
              setequal(top20[[k]], ranked[cluster == k & rank_in_cluster <= 20, gene]))
}
message("Check passed: ", length(top20), " clusters; xlsx top 20 = CSV ranks 1-20")


## =============================================================================
## CHUNK 4: enricher per cluster (clusterProfiler + CellMarker, as in the book)
## =============================================================================
run_cluster <- function(genes) {
  g <- intersect(split_genes(genes), universe)                              # markers on the panel
  if (length(g) == 0) return(NULL)
  e <- enricher(g, TERM2GENE = t2g, universe = universe,                    # hypergeometric test, BH per cluster
                minGSSize = min_set, maxGSSize = 5000,
                pvalueCutoff = 1, qvalueCutoff = 1)                         # keep all; the call uses padj_max below
  if (is.null(e) || nrow(as.data.frame(e)) == 0) return(NULL)
  r <- as.data.table(as.data.frame(e))
  r[, n_markers_on_panel := length(g)]
  r
}

run_and_save <- function(marker_lists, label, out_file) {
  all <- rbindlist(lapply(names(marker_lists), function(k) {
    r <- run_cluster(marker_lists[[k]])
    if (is.null(r)) NULL else cbind(cluster = k, r)
  }), fill = TRUE)
  setorder(all, cluster, pvalue)
  
  ## one row per cluster: best CellMarker cell type + the next two
  smry <- rbindlist(lapply(names(marker_lists), function(k) {
    r <- all[cluster == k]
    if (nrow(r) == 0) return(data.table(cluster = k, call = "no CellMarker set tested", best_match = NA_character_))
    data.table(cluster = k,
               n_markers_on_panel = r$n_markers_on_panel[1],
               call = if (r$p.adjust[1] < padj_max) r$ID[1] else "no significant match",
               best_match = r$ID[1],
               n_shared = r$Count[1],
               shared_genes = gsub("/", ", ", r$geneID[1]),
               p_adj = signif(r$p.adjust[1], 3),
               second = if (nrow(r) > 1) paste0(r$ID[2], " (", r$Count[2], " shared, padj ", signif(r$p.adjust[2], 2), ")") else NA,
               third  = if (nrow(r) > 2) paste0(r$ID[3], " (", r$Count[3], " shared, padj ", signif(r$p.adjust[3], 2), ")") else NA)
  }), fill = TRUE)[order(as.numeric(cluster))]
  print(smry[, .(cluster, call, n_shared, p_adj, shared_genes)])
  
  keep <- intersect(c("cluster", "ID", "GeneRatio", "BgRatio", "Count", "pvalue", "p.adjust", "geneID",
                      "n_markers_on_panel"), names(all))                    # columns differ between versions
  sets <- t2g[gene %in% universe, .(n_markers_on_panel = .N,
                                    markers_on_panel = paste(sort(gene), collapse = ", ")), by = term][order(term)]
  how_to_read <- data.frame(item = c("Markers used", "Reference", "Method", "call", "best_match",
                                     "n_shared / shared_genes", "p_adj"),
                            meaning = c(label,
                                        paste0("CellMarker 2.0 human, normal cells, tissue matching '", tissue_regex,
                                               "', neural types excluded ('", exclude_regex, "')"),
                                        paste0("clusterProfiler::enricher (hypergeometric), background = panel genes tested by FindAllMarkers; ",
                                               "cell types with < ", min_set, " markers on the panel not tested"),
                                        paste0("best_match if p_adj < ", padj_max, ", otherwise 'no significant match'"),
                                        "CellMarker cell type sharing the most markers beyond chance (shown even if not significant)",
                                        "cluster markers that are CellMarker markers of that cell type",
                                        "BH-adjusted p-value across the CellMarker cell types tested for that cluster"))
  write.xlsx(list(Summary = as.data.frame(smry), All_tests = as.data.frame(all[, ..keep]),
                  CellMarker_sets_used = as.data.frame(sets), How_to_read = how_to_read),
             out_file, firstRow = TRUE, colWidths = "auto", overwrite = TRUE)
  message("Saved: ", out_file)
}


## =============================================================================
## CHUNK 5: run - top 20 and top 50
## =============================================================================
run_and_save(top20, "top 20 per cluster (xlsx sheets; ranked by pct.1, then avg_log2FC)",
             file.path(out_dir, "cellmarker_top20.xlsx"))
run_and_save(top50, "top 50 per cluster (ranked CSV, rank_in_cluster <= 50)",
             file.path(out_dir, "cellmarker_top50.xlsx"))





















## =============================================================================
## 17_CellMarker_top20_top50.R
## =============================================================================
## WHAT IT DOES
##   Names each Leiden cluster with CellMarker 2.0 (Hu et al. 2023, NAR),
##   the method shown in the clusterProfiler book (enricher + CellMarker):
##     each cluster's marker genes vs each CellMarker cell type's marker genes,
##     hypergeometric test (clusterProfiler::enricher), BH-adjusted per cluster.
##   CellMarker is restricted to HUMAN, NORMAL cells from VASCULAR tissues
##   (aorta / artery / vessel) + BLOOD (immune). Neural cell types are excluded.
##     top 20 markers (sheets of the top-20 xlsx)  -> cellmarker_top20.xlsx
##     top 50 markers (rank_in_cluster <= 50, CSV) -> cellmarker_top50.xlsx
##   CHUNK 6 (simple counting rule, no statistics): a cluster gets EVERY cell type
##   for which >= 5 of its top 10 genes are CellMarker markers -> cellmarker_top10_min5.xlsx
##   All saved next to the input files. No Python / reticulate needed.
##
## ONE-TIME SETUP (in R):
##   install.packages(c("BiocManager", "readxl", "openxlsx", "data.table"))
##   BiocManager::install("clusterProfiler")
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 1: settings
## =============================================================================
library(data.table)       # tables
library(openxlsx)         # read the top-20 sheets + write results
library(readxl)           # read the CellMarker xlsx
library(clusterProfiler)  # enricher(): hypergeometric test against CellMarker sets

## xlsx_top20 <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_top20_per_cluster.xlsx"   # <- change here
xlsx_top20 <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_05_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_top20_per_cluster.xlsx"   # <- change here
## csv_ranked <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_all_significant_ranked.csv"  # <- change here
csv_ranked <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_05_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_all_significant_ranked.csv"  # <- change here

## CellMarker 2.0 human file (URL from the clusterProfiler book); downloaded once into the Results folder
cm_file <- file.path(dirname(xlsx_top20), "Cell_marker_Human.xlsx")
cm_url  <- "http://www.bio-bigdata.center/CellMarker_download_files/file/Cell_marker_Human.xlsx"

tissue_regex  <- "aort|arter|vessel|vascul|blood"             # vascular tissues + blood (immune)   <- change here
exclude_regex <- "neur|glia|schwann|astrocyt|oligodendro|microglia"   # no neural cell types       <- change here
min_set  <- 3      # a CellMarker cell type is tested only if >= 3 of its markers are on your panel
padj_max <- 0.05   # a call is made only if BH-adjusted p < this

## CHUNK 6 rule: a cluster is named a cell type if >= min_hits of its top_k genes are that type's CellMarker markers
top_k    <- 10     # top genes per cluster (rank_in_cluster 1-10)                   <- change here
min_hits <- 5      # how many of them must be markers of the cell type              <- change here

stopifnot("top-20 xlsx not found" = file.exists(xlsx_top20),
          "ranked CSV not found"  = file.exists(csv_ranked))
out_dir <- dirname(xlsx_top20)                                              # save next to the inputs
rds_all <- sub("_top20_per_cluster\\.xlsx$", "_all_results.RDS", xlsx_top20)   # script 14's full results
stopifnot("script 14's _all_results.RDS not found next to the xlsx" = file.exists(rds_all))


## =============================================================================
## CHUNK 2: CellMarker 2.0 - download once, keep vascular + blood, normal cells
## =============================================================================
if (!file.exists(cm_file)) {                                                # download only if not yet there
  ok <- tryCatch({ download.file(cm_url, cm_file, mode = "wb"); TRUE }, error = function(e) FALSE)
  if (!ok || file.size(cm_file) < 1e5) {
    stop("Download failed. Download 'Cell_marker_Human.xlsx' by hand from the CellMarker 2.0 Download page\n",
         "(http://bio-bigdata.hrbmu.edu.cn/CellMarker or http://117.50.127.228/CellMarker/) and save it as:\n", cm_file)
  }
}
cm <- as.data.table(read_excel(cm_file))
print(names(cm))                                                            # columns of the file
gene_col <- if ("Symbol" %in% names(cm)) "Symbol" else "marker"             # official gene symbol if present
stopifnot("unexpected CellMarker columns" = all(c("tissue_class", "tissue_type", "cell_name", gene_col) %in% names(cm)))

if ("cell_type" %in% names(cm)) {                                           # "Normal cell" vs "Cancer cell"
  print(table(cm$cell_type))
  cm <- cm[cell_type == "Normal cell"]                                      # normal cells only
}
cm_vasc <- cm[grepl(tissue_regex, tissue_class, ignore.case = TRUE) |
                grepl(tissue_regex, tissue_type,  ignore.case = TRUE)]        # vascular + blood entries
cm_vasc <- cm_vasc[!grepl(exclude_regex, cell_name, ignore.case = TRUE)]    # drop neural cell types

## CHECK: which tissues were kept (edit tissue_regex in CHUNK 1 if something is wrong)
print(cm_vasc[, .(n_entries = .N), by = .(tissue_class, tissue_type)][order(-n_entries)], nrows = 200)

## TERM2GENE: one row per (cell type, marker gene)
t2g <- unique(cm_vasc[!is.na(get(gene_col)) & get(gene_col) != "",
                      .(term = cell_name, gene = get(gene_col))])
message(length(unique(t2g$term)), " CellMarker cell types, ", length(unique(t2g$gene)), " marker genes kept")


## =============================================================================
## CHUNK 3: your marker lists (top 20 from xlsx, top 50 from CSV)
## =============================================================================
## CosMx combined probes ("IGHG1/2", "IGLC1/2") -> separate genes ("IGHG1", "IGHG2")
split_probe <- function(g) {
  parts <- strsplit(g, "/", fixed = TRUE)[[1]]                              # "IGHG1/2" -> "IGHG1", "2"
  if (length(parts) == 1) return(g)                                         # ordinary gene: unchanged
  stem <- sub("[0-9]+$", "", parts[1])                                      # "IGHG1" -> "IGHG"
  c(parts[1], paste0(stem, parts[-1]))                                      # "IGHG1", "IGHG2"
}
split_genes <- function(gs) unique(unlist(lapply(gs, split_probe)))

## background = every gene FindAllMarkers tested on your panel (script 14)
universe <- split_genes(unique(as.character(readRDS(rds_all)$gene)))
message(length(universe), " panel genes in the background; ",
        length(intersect(unique(t2g$gene), universe)), " CellMarker genes are on your panel")

sheets <- getSheetNames(xlsx_top20)                                         # "cluster_0", "cluster_1", ...
top20 <- lapply(sheets, function(s) {
  d <- read.xlsx(xlsx_top20, sheet = s)
  if ("gene" %in% names(d)) d$gene else character(0)                        # sheet may say "no significant markers"
})
names(top20) <- sub("^cluster_", "", sheets)

ranked <- fread(csv_ranked)
ranked[, cluster := as.character(cluster)]
top50 <- lapply(names(top20), function(k) ranked[cluster == k & rank_in_cluster <= 50, gene])
names(top50) <- names(top20)

## check: xlsx top 20 = CSV ranks 1-20 (same script-14 run)
for (k in names(top20)) {
  stopifnot("xlsx and CSV disagree - are they from the same run?" =
              setequal(top20[[k]], ranked[cluster == k & rank_in_cluster <= 20, gene]))
}
message("Check passed: ", length(top20), " clusters; xlsx top 20 = CSV ranks 1-20")


## =============================================================================
## CHUNK 4: enricher per cluster (clusterProfiler + CellMarker, as in the book)
## =============================================================================
run_cluster <- function(genes) {
  g <- intersect(split_genes(genes), universe)                              # markers on the panel
  if (length(g) == 0) return(NULL)
  e <- enricher(g, TERM2GENE = t2g, universe = universe,                    # hypergeometric test, BH per cluster
                minGSSize = min_set, maxGSSize = 5000,
                pvalueCutoff = 1, qvalueCutoff = 1)                         # keep all; the call uses padj_max below
  if (is.null(e) || nrow(as.data.frame(e)) == 0) return(NULL)
  r <- as.data.table(as.data.frame(e))
  r[, n_markers_on_panel := length(g)]
  r
}

run_and_save <- function(marker_lists, label, out_file) {
  all <- rbindlist(lapply(names(marker_lists), function(k) {
    r <- run_cluster(marker_lists[[k]])
    if (is.null(r)) NULL else cbind(cluster = k, r)
  }), fill = TRUE)
  setorder(all, cluster, pvalue)
  
  ## one row per cluster: best CellMarker cell type + the next two
  smry <- rbindlist(lapply(names(marker_lists), function(k) {
    r <- all[cluster == k]
    if (nrow(r) == 0) return(data.table(cluster = k, call = "no CellMarker set tested", best_match = NA_character_))
    data.table(cluster = k,
               n_markers_on_panel = r$n_markers_on_panel[1],
               call = if (r$p.adjust[1] < padj_max) r$ID[1] else "no significant match",
               best_match = r$ID[1],
               n_shared = r$Count[1],
               shared_genes = gsub("/", ", ", r$geneID[1]),
               p_adj = signif(r$p.adjust[1], 3),
               second = if (nrow(r) > 1) paste0(r$ID[2], " (", r$Count[2], " shared, padj ", signif(r$p.adjust[2], 2), ")") else NA,
               third  = if (nrow(r) > 2) paste0(r$ID[3], " (", r$Count[3], " shared, padj ", signif(r$p.adjust[3], 2), ")") else NA)
  }), fill = TRUE)[order(as.numeric(cluster))]
  print(smry[, .(cluster, call, n_shared, p_adj, shared_genes)])
  
  keep <- intersect(c("cluster", "ID", "GeneRatio", "BgRatio", "Count", "pvalue", "p.adjust", "geneID",
                      "n_markers_on_panel"), names(all))                    # columns differ between versions
  sets <- t2g[gene %in% universe, .(n_markers_on_panel = .N,
                                    markers_on_panel = paste(sort(gene), collapse = ", ")), by = term][order(term)]
  how_to_read <- data.frame(item = c("Markers used", "Reference", "Method", "call", "best_match",
                                     "n_shared / shared_genes", "p_adj"),
                            meaning = c(label,
                                        paste0("CellMarker 2.0 human, normal cells, tissue matching '", tissue_regex,
                                               "', neural types excluded ('", exclude_regex, "')"),
                                        paste0("clusterProfiler::enricher (hypergeometric), background = panel genes tested by FindAllMarkers; ",
                                               "cell types with < ", min_set, " markers on the panel not tested"),
                                        paste0("best_match if p_adj < ", padj_max, ", otherwise 'no significant match'"),
                                        "CellMarker cell type sharing the most markers beyond chance (shown even if not significant)",
                                        "cluster markers that are CellMarker markers of that cell type",
                                        "BH-adjusted p-value across the CellMarker cell types tested for that cluster"))
  write.xlsx(list(Summary = as.data.frame(smry), All_tests = as.data.frame(all[, ..keep]),
                  CellMarker_sets_used = as.data.frame(sets), How_to_read = how_to_read),
             out_file, firstRow = TRUE, colWidths = "auto", overwrite = TRUE)
  message("Saved: ", out_file)
}


## =============================================================================
## CHUNK 5: run - top 20 and top 50
## =============================================================================
run_and_save(top20, "top 20 per cluster (xlsx sheets; ranked by pct.1, then avg_log2FC)",
             file.path(out_dir, "cellmarker_top20.xlsx"))
run_and_save(top50, "top 50 per cluster (ranked CSV, rank_in_cluster <= 50)",
             file.path(out_dir, "cellmarker_top50.xlsx"))


## =============================================================================
## CHUNK 6: counting rule - >= min_hits of the top_k genes are markers of a cell type
## =============================================================================
## No test and no p-value: a cluster gets EVERY cell type that passes, so it can
## have several. expected_by_chance shows how many hits a random set of top_k
## panel genes would give for that cell type (large CellMarker lists catch many
## genes by chance - compare n_hits with it).
sets_on_panel <- split(t2g[gene %in% universe, gene], t2g[gene %in% universe, term])   # markers per cell type, panel only
N <- length(universe)                                                       # panel genes in the background

count_rule <- rbindlist(lapply(names(top20), function(k) {
  g <- ranked[cluster == k & rank_in_cluster <= top_k][order(rank_in_cluster), gene]   # top_k genes, in rank order
  rbindlist(lapply(names(sets_on_panel), function(ct) {
    s   <- sets_on_panel[[ct]]
    hit <- g[vapply(g, function(x) any(split_probe(x) %in% s), logical(1))] # a probe like IGHG1/2 counts once
    data.table(cluster = k, cell_type = ct, n_top_genes = length(g), n_hits = length(hit),
               hit_genes = paste(hit, collapse = ", "),
               markers_of_type_on_panel = length(s),
               expected_by_chance = round(length(g) * length(s) / N, 2))   # random top_k genes would hit this many
  }))
}))
setorder(count_rule, cluster, -n_hits, expected_by_chance)

passed <- count_rule[n_hits >= min_hits]                                   # every cell type that passes the rule

## one row per cluster: all cell types that pass (most hits first), plus the best near-miss
calls <- rbindlist(lapply(names(top20), function(k) {
  p <- passed[cluster == k]
  b <- count_rule[cluster == k][1]                                          # best cell type even if it failed
  data.table(cluster = k,
             top_genes = paste(ranked[cluster == k & rank_in_cluster <= top_k][order(rank_in_cluster), gene],
                               collapse = ", "),
             n_cell_types_called = nrow(p),
             cell_types_called = if (nrow(p) > 0) paste0(p$cell_type, " (", p$n_hits, "/", p$n_top_genes, ")",
                                                         collapse = "; ") else paste0("none (< ", min_hits, " hits)"),
             best_cell_type = b$cell_type, best_n_hits = b$n_hits, best_hit_genes = b$hit_genes)
}))[order(as.numeric(cluster))]
print(calls[, .(cluster, n_cell_types_called, cell_types_called)])

how_to_read_rule <- data.frame(item = c("Rule", "Reference", "cell_types_called", "best_cell_type",
                                        "n_hits / hit_genes", "expected_by_chance", "Caution"),
                               meaning = c(paste0("a cluster is named a cell type if >= ", min_hits, " of its top ", top_k,
                                                  " genes (rank_in_cluster 1-", top_k, ") are CellMarker markers of that type; several types allowed"),
                                           paste0("CellMarker 2.0 human, normal cells, tissue matching '", tissue_regex, "', neural types excluded"),
                                           "every cell type passing the rule, with hits/top genes in brackets",
                                           "cell type with the most hits, shown even when it does not pass (near-miss)",
                                           "how many / which of the top genes are markers of that cell type",
                                           "hits expected from a random set of the same number of panel genes",
                                           "No statistical test. A cell type with a long CellMarker list can reach the cut-off by chance: trust a call only if n_hits is well above expected_by_chance and the hit genes are canonical for that type."))
write.xlsx(list(Calls = as.data.frame(calls),
                Passing_cell_types = as.data.frame(passed[order(as.numeric(cluster), -n_hits)]),
                All_counts = as.data.frame(count_rule[n_hits > 0][order(as.numeric(cluster), -n_hits)]),
                How_to_read = how_to_read_rule),
           file.path(out_dir, paste0("cellmarker_top", top_k, "_min", min_hits, ".xlsx")),
           firstRow = TRUE, colWidths = "auto", overwrite = TRUE)
message("Saved: ", file.path(out_dir, paste0("cellmarker_top", top_k, "_min", min_hits, ".xlsx")))

































## =============================================================================
## CELLxGENE_enricher_top20_top50.R
## =============================================================================

##dir_res <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_06_07_08_09_final_All_markers"
dir_res <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_05_06_07_08_09_final_All_markers"

split_probe <- function(g) {                                                # "IGHG1/2" -> "IGHG1", "IGHG2"
  parts <- strsplit(g, "/", fixed = TRUE)[[1]]
  if (length(parts) == 1) return(g)
  c(parts[1], paste0(sub("[0-9]+$", "", parts[1]), parts[-1]))
}
genes <- unique(as.character(readRDS(file.path(dir_res,
                                               "FindAllMarkers_harmony_PC20_NN250_rep3_res1_all_results.RDS"))$gene))   # panel genes tested
writeLines(unique(unlist(lapply(genes, split_probe))),
           file.path(dir_res, "panel_genes_for_CELLxGENE.txt"))             # one gene per line







## =============================================================================
## 18A_CELLxGENE_enricher_top20_top50.R
## =============================================================================
## WHAT IT DOES (same approach as script 17 CHUNKS 4-5, reference = CELLxGENE)
##   Reference: marker genes per cell type from CZ CELLxGENE Census
##     (human vasculature, normal donors), built on YOUR panel genes by
##     18_CELLxGENE_reference_markers.py (run automatically below if missing).
##   Test: each cluster's markers vs each cell type's markers, hypergeometric
##     (clusterProfiler::enricher), BH-adjusted per cluster.
##     top 20 (xlsx sheets)          -> CELLxGENE_top20.xlsx
##     top 50 (CSV, rank <= 50)      -> CELLxGENE_top50.xlsx
##
## ONE-TIME SETUP
##   Terminal: conda activate celltypist_env ; pip install cellxgene-census
##   R:        BiocManager::install("clusterProfiler") ; install.packages("reticulate")
## RESTART R before running (reticulate binds one Python per session). Needs internet.
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 1: settings
## =============================================================================
library(data.table)       # tables
library(openxlsx)         # read the top-20 sheets + write results
library(clusterProfiler)  # enricher()
library(reticulate)       # runs the Python Census step

## xlsx_top20 <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_top20_per_cluster.xlsx"   # <- change here
xlsx_top20 <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_05_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_top20_per_cluster.xlsx"   # <- change here
## csv_ranked <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_all_significant_ranked.csv"  # <- change here
csv_ranked <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_05_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_all_significant_ranked.csv"  # <- change here
py_file    <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/18_CELLxGENE_reference_markers.py"   # <- where you saved it
conda_env  <- "celltypist_env"

## reference settings (must be the same in 18A and 18B)
tissue_general <- "vasculature"   # Census coarse tissue
normal_only    <- TRUE            # healthy donors only
max_cells      <- 3000            # cells sampled per reference cell type
min_cells      <- 50              # reference cell types with fewer cells dropped
ref_top_n      <- 100             # marker genes per reference cell type

min_set  <- 3      # a reference cell type is tested only if >= 3 of its markers are in the background
padj_max <- 0.05   # a call is made only if BH-adjusted p < this

stopifnot("top-20 xlsx not found" = file.exists(xlsx_top20), "ranked CSV not found" = file.exists(csv_ranked))
out_dir <- dirname(xlsx_top20)                                              # save next to the inputs
rds_all <- sub("_top20_per_cluster\\.xlsx$", "_all_results.RDS", xlsx_top20)   # script 14's full results
stopifnot("script 14's _all_results.RDS not found next to the xlsx" = file.exists(rds_all))
ref_csv <- file.path(out_dir, paste0("CELLxGENE_reference_markers_", tissue_general, ".csv"))


## =============================================================================
## CHUNK 2: CELLxGENE reference - built once (skipped if the CSV already exists)
## =============================================================================
## CosMx combined probes ("IGHG1/2") -> separate genes ("IGHG1", "IGHG2")
split_probe <- function(g) {
  parts <- strsplit(g, "/", fixed = TRUE)[[1]]                              # "IGHG1/2" -> "IGHG1", "2"
  if (length(parts) == 1) return(g)                                         # ordinary gene: unchanged
  stem <- sub("[0-9]+$", "", parts[1])                                      # "IGHG1" -> "IGHG"
  c(parts[1], paste0(stem, parts[-1]))                                      # "IGHG1", "IGHG2"
}
split_genes <- function(gs) unique(unlist(lapply(gs, split_probe)))
panel <- split_genes(unique(as.character(readRDS(rds_all)$gene)))          # genes FindAllMarkers tested

if (!file.exists(ref_csv)) {
  stopifnot("Python file not found - set py_file in CHUNK 1" = file.exists(py_file))
  panel_txt <- file.path(out_dir, "panel_genes_for_CELLxGENE.txt")
  writeLines(panel, panel_txt)                                              # gene list for Census
  use_condaenv(conda_env, required = TRUE)
  source_python(py_file)                                                    # makes build_reference() available
  build_reference(panel_txt, ref_csv, tissue_general = tissue_general, normal_only = normal_only,
                  max_cells_per_type = as.integer(max_cells), min_cells_per_type = as.integer(min_cells),
                  top_n = as.integer(ref_top_n))
}
ref <- fread(ref_csv)                                                       # cell_type, rank, gene, ...
genes_in_ref <- readLines(sub("\\.csv$", "_genes_in_reference.txt", ref_csv))
universe <- intersect(panel, genes_in_ref)                                  # background: panel genes Census also has
t2g <- unique(ref[gene %in% universe, .(term = cell_type, gene)])           # TERM2GENE for enricher
print(ref[, .(n_cells_used = n_cells_used[1], n_datasets = n_datasets[1], n_markers = .N),
          by = cell_type], nrows = 200)                                     # reference cell types + their size
message(length(unique(t2g$term)), " reference cell types; background ", length(universe), " genes")


## =============================================================================
## CHUNK 3: your marker lists (top 20 from xlsx, top 50 from CSV)
## =============================================================================
sheets <- getSheetNames(xlsx_top20)                                         # "cluster_0", "cluster_1", ...
top20 <- lapply(sheets, function(s) {
  d <- read.xlsx(xlsx_top20, sheet = s)
  if ("gene" %in% names(d)) d$gene else character(0)                        # sheet may say "no significant markers"
})
names(top20) <- sub("^cluster_", "", sheets)

ranked <- fread(csv_ranked)
ranked[, cluster := as.character(cluster)]
top50 <- lapply(names(top20), function(k) ranked[cluster == k & rank_in_cluster <= 50, gene])
names(top50) <- names(top20)

for (k in names(top20)) {                                                   # check: same script-14 run
  stopifnot("xlsx and CSV disagree - are they from the same run?" =
              setequal(top20[[k]], ranked[cluster == k & rank_in_cluster <= 20, gene]))
}
message("Check passed: ", length(top20), " clusters; xlsx top 20 = CSV ranks 1-20")


## =============================================================================
## CHUNK 4: enricher per cluster
## =============================================================================
run_cluster <- function(genes) {
  g <- intersect(split_genes(genes), universe)                              # markers in the background
  if (length(g) == 0) return(NULL)
  e <- enricher(g, TERM2GENE = t2g, universe = universe, minGSSize = min_set, maxGSSize = 5000,
                pvalueCutoff = 1, qvalueCutoff = 1)                         # keep all; call uses padj_max
  if (is.null(e) || nrow(as.data.frame(e)) == 0) return(NULL)
  r <- as.data.table(as.data.frame(e))
  r[, n_markers_in_background := length(g)]
  r
}

run_and_save <- function(marker_lists, label, out_file) {
  all <- rbindlist(lapply(names(marker_lists), function(k) {
    r <- run_cluster(marker_lists[[k]])
    if (is.null(r)) NULL else cbind(cluster = k, r)
  }), fill = TRUE)
  setorder(all, cluster, pvalue)
  
  smry <- rbindlist(lapply(names(marker_lists), function(k) {             # best match + next two per cluster
    r <- all[cluster == k]
    if (nrow(r) == 0) return(data.table(cluster = k, call = "nothing tested"))
    data.table(cluster = k, n_markers_in_background = r$n_markers_in_background[1],
               call = if (r$p.adjust[1] < padj_max) r$ID[1] else "no significant match",
               best_match = r$ID[1], n_shared = r$Count[1],
               shared_genes = gsub("/", ", ", r$geneID[1]), p_adj = signif(r$p.adjust[1], 3),
               second = if (nrow(r) > 1) paste0(r$ID[2], " (", r$Count[2], " shared, padj ", signif(r$p.adjust[2], 2), ")") else NA,
               third  = if (nrow(r) > 2) paste0(r$ID[3], " (", r$Count[3], " shared, padj ", signif(r$p.adjust[3], 2), ")") else NA)
  }), fill = TRUE)[order(as.numeric(cluster))]
  print(smry[, .(cluster, call, n_shared, p_adj, shared_genes)])
  
  keep <- intersect(c("cluster", "ID", "GeneRatio", "BgRatio", "Count", "pvalue", "p.adjust", "geneID",
                      "n_markers_in_background"), names(all))
  how_to_read <- data.frame(item = c("Markers used", "Reference", "Method", "call", "n_shared / shared_genes", "p_adj"),
                            meaning = c(label,
                                        paste0("CELLxGENE Census, human ", tissue_general, if (normal_only) ", normal donors" else "",
                                               "; top ", ref_top_n, " Wilcoxon markers per cell type computed on your panel genes"),
                                        "clusterProfiler::enricher (hypergeometric); background = panel genes present in Census",
                                        paste0("best_match if p_adj < ", padj_max, ", otherwise 'no significant match'"),
                                        "cluster markers that are reference markers of that cell type",
                                        "BH-adjusted p-value across the reference cell types tested for that cluster"))
  write.xlsx(list(Summary = as.data.frame(smry), All_tests = as.data.frame(all[, ..keep]),
                  Reference_markers = as.data.frame(ref), How_to_read = how_to_read),
             out_file, firstRow = TRUE, colWidths = "auto", overwrite = TRUE)
  message("Saved: ", out_file)
}


## =============================================================================
## CHUNK 5: run - top 20 and top 50
## =============================================================================
run_and_save(top20, "top 20 per cluster (xlsx sheets; ranked by pct.1, then avg_log2FC)",
             file.path(out_dir, "CELLxGENE_top20.xlsx"))
run_and_save(top50, "top 50 per cluster (ranked CSV, rank_in_cluster <= 50)",
             file.path(out_dir, "CELLxGENE_top50.xlsx"))





















## =============================================================================
## 18B_CELLxGENE_count_top10_min5.R
## =============================================================================
## WHAT IT DOES (same approach as script 17 CHUNK 6, reference = CELLxGENE)
##   Reference: marker genes per cell type from CZ CELLxGENE Census
##     (human vasculature, normal donors), built on YOUR panel genes by
##     18_CELLxGENE_reference_markers.py (shared with 18A; built here if missing).
##   Rule (no statistics): a cluster is named EVERY cell type for which
##     >= min_hits of its top_k genes are reference markers of that type.
##     -> CELLxGENE_top10_min5.xlsx
##
## ONE-TIME SETUP: as 18A. RESTART R before running. Run chunk by chunk.
## =============================================================================


## =============================================================================
## CHUNK 1: settings
## =============================================================================
library(data.table)   # tables
library(openxlsx)     # write results
library(reticulate)   # runs the Python Census step (only if the reference is missing)

## csv_ranked <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_all_significant_ranked.csv"  # <- change here
csv_ranked <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_05_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_all_significant_ranked.csv"  # <- change here
py_file    <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/18_CELLxGENE_reference_markers.py"   # <- where you saved it
conda_env  <- "celltypist_env"

## reference settings (must be the same in 18A and 18B)
tissue_general <- "vasculature"
normal_only    <- TRUE
max_cells      <- 3000
min_cells      <- 50
ref_top_n      <- 100

top_k    <- 10     # top genes per cluster (rank_in_cluster 1-10)        <- change here
min_hits <- 5      # how many must be reference markers of the cell type  <- change here

stopifnot("ranked CSV not found" = file.exists(csv_ranked))
out_dir <- dirname(csv_ranked)                                              # save next to the inputs
rds_all <- sub("_all_significant_ranked\\.csv$", "_all_results.RDS", csv_ranked)   # script 14's full results
stopifnot("script 14's _all_results.RDS not found next to the CSV" = file.exists(rds_all))
ref_csv <- file.path(out_dir, paste0("CELLxGENE_reference_markers_", tissue_general, ".csv"))


## =============================================================================
## CHUNK 2: CELLxGENE reference (same file as 18A; built only if missing)
## =============================================================================
split_probe <- function(g) {
  parts <- strsplit(g, "/", fixed = TRUE)[[1]]                              # "IGHG1/2" -> "IGHG1", "2"
  if (length(parts) == 1) return(g)
  stem <- sub("[0-9]+$", "", parts[1])
  c(parts[1], paste0(stem, parts[-1]))                                      # "IGHG1", "IGHG2"
}
split_genes <- function(gs) unique(unlist(lapply(gs, split_probe)))
panel <- split_genes(unique(as.character(readRDS(rds_all)$gene)))          # genes FindAllMarkers tested

if (!file.exists(ref_csv)) {
  stopifnot("Python file not found - set py_file in CHUNK 1" = file.exists(py_file))
  panel_txt <- file.path(out_dir, "panel_genes_for_CELLxGENE.txt")
  writeLines(panel, panel_txt)
  use_condaenv(conda_env, required = TRUE)
  source_python(py_file)
  build_reference(panel_txt, ref_csv, tissue_general = tissue_general, normal_only = normal_only,
                  max_cells_per_type = as.integer(max_cells), min_cells_per_type = as.integer(min_cells),
                  top_n = as.integer(ref_top_n))
}
ref <- fread(ref_csv)
universe <- intersect(panel, readLines(sub("\\.csv$", "_genes_in_reference.txt", ref_csv)))
sets <- split(ref[gene %in% universe, gene], ref[gene %in% universe, cell_type])   # markers per reference type
N <- length(universe)
message(length(sets), " reference cell types; background ", N, " genes")


## =============================================================================
## CHUNK 3: counting rule - >= min_hits of the top_k genes are markers of a type
## =============================================================================
ranked <- fread(csv_ranked)
ranked[, cluster := as.character(cluster)]
clusters <- as.character(sort(unique(as.numeric(ranked$cluster))))         # 0, 1, 2, ...

counts <- rbindlist(lapply(clusters, function(k) {
  g <- ranked[cluster == k & rank_in_cluster <= top_k][order(rank_in_cluster), gene]   # top_k genes in rank order
  rbindlist(lapply(names(sets), function(ct) {
    s   <- sets[[ct]]
    hit <- g[vapply(g, function(x) any(split_probe(x) %in% s), logical(1))]   # IGHG1/2 counts once
    data.table(cluster = k, cell_type = ct, n_top_genes = length(g), n_hits = length(hit),
               hit_genes = paste(hit, collapse = ", "), markers_of_type = length(s),
               expected_by_chance = round(length(g) * length(s) / N, 2))     # random genes would hit this many
  }))
}))
setorder(counts, cluster, -n_hits, expected_by_chance)
passed <- counts[n_hits >= min_hits]                                        # every type passing the rule

calls <- rbindlist(lapply(clusters, function(k) {
  p <- passed[cluster == k]
  b <- counts[cluster == k][1]                                              # best type even if it failed
  data.table(cluster = k,
             top_genes = paste(ranked[cluster == k & rank_in_cluster <= top_k][order(rank_in_cluster), gene], collapse = ", "),
             n_cell_types_called = nrow(p),
             cell_types_called = if (nrow(p) > 0) paste0(p$cell_type, " (", p$n_hits, "/", p$n_top_genes, ")", collapse = "; ")
             else paste0("none (< ", min_hits, " hits)"),
             best_cell_type = b$cell_type, best_n_hits = b$n_hits, best_hit_genes = b$hit_genes)
}))
print(calls[, .(cluster, n_cell_types_called, cell_types_called)])


## =============================================================================
## CHUNK 4: save
## =============================================================================
how_to_read <- data.frame(item = c("Rule", "Reference", "cell_types_called", "best_cell_type",
                                   "expected_by_chance", "Caution"),
                          meaning = c(paste0("a cluster is named a cell type if >= ", min_hits, " of its top ", top_k,
                                             " genes are reference markers of that type; several types allowed"),
                                      paste0("CELLxGENE Census, human ", tissue_general, if (normal_only) ", normal donors" else "",
                                             "; top ", ref_top_n, " Wilcoxon markers per cell type on your panel genes"),
                                      "every cell type passing the rule, hits/top genes in brackets",
                                      "cell type with the most hits, shown even when it does not pass (near-miss)",
                                      "hits expected from the same number of random panel genes",
                                      "No statistical test: trust a call only if n_hits is well above expected_by_chance and the hit genes are canonical for that type."))
out_file <- file.path(out_dir, paste0("CELLxGENE_top", top_k, "_min", min_hits, ".xlsx"))
write.xlsx(list(Calls = as.data.frame(calls),
                Passing_cell_types = as.data.frame(passed[order(as.numeric(cluster), -n_hits)]),
                All_counts = as.data.frame(counts[n_hits > 0][order(as.numeric(cluster), -n_hits)]),
                Reference_markers = as.data.frame(ref), How_to_read = how_to_read),
           out_file, firstRow = TRUE, colWidths = "auto", overwrite = TRUE)
message("Saved: ", out_file)



























































## =============================================================================
## 19_HCA_aorta_build_reference.R
## =============================================================================
## =============================================================================
## 19_HCA_aorta_build_reference.R
## =============================================================================
## WHAT IT DOES
##   Builds a human AORTA reference from the HCA Data Explorer project
##   "AneurysmalHumanAorticTissue" (Li et al. 2020 Circulation; GEO GSE155468;
##   ascending aorta, 3 control + 8 aneurysm donors, 10X scRNA-seq), using the
##   AUTHORS' OWN cell labels (AorticTissue_cell_types_annotations.csv):
##     celltype  = major cell types (e.g. SMC1)
##     celltype2 = cell types/states (e.g. Contractile SMC)
##   1. read the per-donor count tables (GSE155468_RAW.tar: GSM..._Con4.txt.gz ...)
##   2. keep the annotated cells, attach the authors' labels
##   3. markers of each labelled type, ON YOUR PANEL GENES ONLY
##      (FindAllMarkers, Wilcoxon), for both label levels:
##        HCA_reference_markers_aorta_celltype.csv
##        HCA_reference_markers_aorta_celltype2.csv
##   19A and 19B use those CSVs to name your CosMx clusters.
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 1: settings
## =============================================================================
library(data.table)   # tables
library(Matrix)       # sparse matrices
library(Seurat)       # normalisation + FindAllMarkers

hca_dir <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/HCA_aorta_reference"   # folder holding the curl download
##dir_res <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_06_07_08_09_final_All_markers"
dir_res <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_05_06_07_08_09_final_All_markers"
rds_all <- file.path(dir_res, "FindAllMarkers_harmony_PC20_NN250_rep3_res1_all_results.RDS")   # your panel genes

ref_top_n <- 100      # marker genes kept per reference cell type
seed      <- 42


## =============================================================================
## CHUNK 2: find the files, unpack the tar, read the authors' annotation
## =============================================================================
tar_file <- list.files(hca_dir, pattern = "GSE155468_RAW\\.tar$", recursive = TRUE, full.names = TRUE)
ann_file <- list.files(hca_dir, pattern = "cell_types_annotations\\.csv$", recursive = TRUE, full.names = TRUE)
stopifnot("GSE155468_RAW.tar not found under hca_dir" = length(tar_file) == 1,
          "annotation CSV not found under hca_dir"     = length(ann_file) == 1)
unpack <- file.path(hca_dir, "unpacked")
if (!dir.exists(unpack)) untar(tar_file, exdir = unpack)                    # 11 files: GSM..._<donor>.txt.gz
txts <- list.files(unpack, pattern = "\\.txt\\.gz$", full.names = TRUE)
print(basename(txts))

ann <- fread(ann_file)                                                      # one row per annotated cell
setnames(ann, 1, "cell_id")                                                 # first column has no header
bc_core <- function(x) sub("^([ACGT]+).*$", "\\1", x)                       # barcode letters only ("AAAC...GACAT")
ann[, key := paste0(orig.ident, "|", bc_core(cell_id))]                     # donor + barcode
stopifnot("donor+barcode not unique in annotation" = !anyDuplicated(ann$key))
print(table(ann$orig.ident))                                                # cells per donor
print(ann[, .N, by = celltype][order(-N)])                                  # major types
print(ann[, .N, by = .(celltype, celltype2)][order(celltype, -N)], nrows = 200)   # types/states


## =============================================================================
## CHUNK 3: read each donor's counts, keep the annotated cells
## =============================================================================
mats <- list()
for (f in txts) {
  s <- sub("^GSM[0-9]+_", "", sub("\\.txt\\.gz$", "", basename(f)))       # donor in file name, e.g. "Con4", "TAA1"
  s <- sub("^TAA", "ATAA", s)                                               # annotation calls them ATAA1..ATAA8
  d <- fread(f)                                                             # genes in rows, cells in columns
  genes <- make.unique(as.character(d[[1]]))                                # first column = gene symbol
  m <- Matrix(as.matrix(d[, -1]), sparse = TRUE); rownames(m) <- genes
  rm(d); invisible(gc())
  key  <- paste0(s, "|", bc_core(colnames(m)))
  keep <- key %in% ann$key                                                  # cells the authors annotated
  m <- m[, keep, drop = FALSE]
  colnames(m) <- ann$cell_id[match(key[keep], ann$key)]                     # use the annotation's cell id
  message(s, ": ", nrow(m), " genes; ", ncol(m), " of ", length(keep), " cells annotated")
  mats[[s]] <- m
}

## what kind of values are these? (decides whether to normalise in CHUNK 4)
x1 <- mats[[1]]
is_counts <- all(sapply(mats, function(m) all(m@x == round(m@x))))         # whole numbers = raw counts
tot_raw   <- Matrix::colSums(x1)                                            # per-cell sum as stored
tot_expm1 <- Matrix::colSums(expm1(x1))                                     # per-cell sum if values are log1p
print(summary(x1@x))                                                        # range of non-zero values
message("whole numbers: ", is_counts,
        "; median per-cell sum: ", round(median(tot_raw), 1),
        "; median per-cell sum of expm1(values): ", round(median(tot_expm1), 1),
        "; authors' median nCount_RNA: ", median(ann[orig.ident == names(mats)[1], nCount_RNA]))
value_type <- if (is_counts) "counts" else
  if (abs(median(tot_expm1) - 1e4) / 1e4 < 0.01) "lognorm_1e4" else "unknown"
message("Values treated as: ", value_type)
stopifnot("values are neither raw counts nor log1p(10,000) - paste the output above to Claude" =
            value_type != "unknown")
stopifnot("matrices look transposed (cells in rows)" = all(sapply(mats, nrow) > 5000))
n_matched <- sum(sapply(mats, ncol))
message("Matched ", n_matched, " of ", nrow(ann), " annotated cells")
stopifnot("fewer than 90% of annotated cells matched - check barcodes" = n_matched > 0.9 * nrow(ann))


## =============================================================================
## CHUNK 4: one Seurat object with the authors' labels, normalised
## =============================================================================
genes  <- Reduce(intersect, lapply(mats, rownames))                         # genes measured in every donor
counts <- do.call(cbind, lapply(mats, function(m) m[genes, , drop = FALSE]))
rm(mats); invisible(gc())
md <- as.data.frame(ann[match(colnames(counts), cell_id), .(orig.ident, stim, celltype, celltype2)])
rownames(md) <- colnames(counts)
ref <- CreateSeuratObject(counts = counts, meta.data = md)
rm(counts); invisible(gc())
if (value_type == "counts") {
  ref <- NormalizeData(ref)                                                 # log1p(counts / total x 10,000)
} else {
  LayerData(ref, assay = "RNA", layer = "data") <- LayerData(ref, assay = "RNA", layer = "counts")   # already log1p(10,000): use as is
}
message("Reference: ", ncol(ref), " cells x ", nrow(ref), " genes")


## =============================================================================
## CHUNK 5: markers of each labelled type, on YOUR panel genes only (both levels)
## =============================================================================
split_probe <- function(g) {                                                # "IGHG1/2" -> "IGHG1", "IGHG2"
  parts <- strsplit(g, "/", fixed = TRUE)[[1]]
  if (length(parts) == 1) return(g)
  c(parts[1], paste0(sub("[0-9]+$", "", parts[1]), parts[-1]))
}
panel <- unique(unlist(lapply(unique(as.character(readRDS(rds_all)$gene)), split_probe)))
feat  <- intersect(panel, rownames(ref))                                    # panel genes the reference measured
message(length(feat), " of ", length(panel), " panel genes are in the reference")

for (lvl in c("celltype", "celltype2")) {
  Idents(ref) <- lvl
  mk <- as.data.table(FindAllMarkers(ref, features = feat, only.pos = TRUE, test.use = "wilcox",
                                     logfc.threshold = 0.5, min.pct = 0.10,
                                     max.cells.per.ident = 3000, random.seed = seed))
  mk <- mk[p_val_adj < 0.05][order(cluster, -avg_log2FC)]
  mk <- mk[, head(.SD, ref_top_n), by = cluster]                            # top ref_top_n per type
  mk[, rank := seq_len(.N), by = cluster]
  info <- as.data.table(ref[[c(lvl, "orig.ident")]])
  setnames(info, lvl, "cell_type")
  info <- info[, .(n_cells_used = .N, n_samples = uniqueN(orig.ident)), by = cell_type]
  out <- merge(mk[, .(cell_type = as.character(cluster), rank, gene, log2FC = round(avg_log2FC, 3),
                      pct_in_type = pct.1, pct_other_types = pct.2, padj = signif(p_val_adj, 3))],
               info, by = "cell_type")[order(cell_type, rank)]
  out_csv <- file.path(dir_res, paste0("HCA_reference_markers_aorta_", lvl, ".csv"))
  fwrite(out, out_csv)
  writeLines(feat, sub("\\.csv$", "_genes_in_reference.txt", out_csv))     # background genes for 19A/19B
  print(out[, .(n_markers = .N, n_cells_used = n_cells_used[1]), by = cell_type])
  message("Saved: ", out_csv)
}
saveRDS(ref, file.path(hca_dir, "HCA_aorta_reference_seurat.RDS"))         # reference kept for later checks








## =============================================================================
## 19A_HCA_enricher_top20_top50.R
## =============================================================================
## WHAT IT DOES (same approach as script 17 CHUNKS 4-5, reference = HCA aorta)
##   Reference: marker genes per cell type from the HCA Data Explorer human aorta
##     dataset (Li et al. 2020), computed on your panel genes by
##     19_HCA_aorta_build_reference.R (run that first).
##   Test: each cluster's markers vs each reference cell type's markers,
##     hypergeometric (clusterProfiler::enricher), BH-adjusted per cluster.
##     top 20 (xlsx sheets)      -> HCA_DataExplorer_top20_<ref_level>.xlsx
##     top 50 (CSV, rank <= 50)  -> HCA_DataExplorer_top50_<ref_level>.xlsx
##   ref_level (CHUNK 1): "celltype" = authors' major types, "celltype2" = types/states.
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 1: settings
## =============================================================================
library(data.table)       # tables
library(openxlsx)         # read the top-20 sheets + write results
library(clusterProfiler)  # enricher()

## xlsx_top20 <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_top20_per_cluster.xlsx"   # <- change here
xlsx_top20 <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_05_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_top20_per_cluster.xlsx"   # <- change here
## csv_ranked <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_all_significant_ranked.csv"  # <- change here
csv_ranked <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_05_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_all_significant_ranked.csv"  # <- change here

ref_level <- "celltype"   # "celltype" = authors' major types; "celltype2" = authors' types/states   <- change here
min_set  <- 3      # a reference cell type is tested only if >= 3 of its markers are in the background
padj_max <- 0.05   # a call is made only if BH-adjusted p < this

stopifnot("top-20 xlsx not found" = file.exists(xlsx_top20), "ranked CSV not found" = file.exists(csv_ranked))
out_dir <- dirname(xlsx_top20)                                              # save next to the inputs
ref_csv <- file.path(out_dir, paste0("HCA_reference_markers_aorta_", ref_level, ".csv"))
stopifnot("run 19_HCA_aorta_build_reference.R first" = file.exists(ref_csv))


## =============================================================================
## CHUNK 2: HCA aorta reference
## =============================================================================
split_probe <- function(g) {                                                # "IGHG1/2" -> "IGHG1", "IGHG2"
  parts <- strsplit(g, "/", fixed = TRUE)[[1]]
  if (length(parts) == 1) return(g)
  c(parts[1], paste0(sub("[0-9]+$", "", parts[1]), parts[-1]))
}
split_genes <- function(gs) unique(unlist(lapply(gs, split_probe)))

ref <- fread(ref_csv)                                                       # cell_type, rank, gene, ...
universe <- readLines(sub("\\.csv$", "_genes_in_reference.txt", ref_csv))   # panel genes the reference measured
t2g <- unique(ref[gene %in% universe, .(term = cell_type, gene)])           # TERM2GENE for enricher
print(ref[, .(n_markers = .N, n_cells_used = n_cells_used[1], n_samples = n_samples[1]), by = cell_type])
message(length(unique(t2g$term)), " reference cell types; background ", length(universe), " genes")


## =============================================================================
## CHUNK 3: your marker lists (top 20 from xlsx, top 50 from CSV)
## =============================================================================
sheets <- getSheetNames(xlsx_top20)                                         # "cluster_0", "cluster_1", ...
top20 <- lapply(sheets, function(s) {
  d <- read.xlsx(xlsx_top20, sheet = s)
  if ("gene" %in% names(d)) d$gene else character(0)                        # sheet may say "no significant markers"
})
names(top20) <- sub("^cluster_", "", sheets)

ranked <- fread(csv_ranked)
ranked[, cluster := as.character(cluster)]
top50 <- lapply(names(top20), function(k) ranked[cluster == k & rank_in_cluster <= 50, gene])
names(top50) <- names(top20)

for (k in names(top20)) {                                                   # check: same script-14 run
  stopifnot("xlsx and CSV disagree - are they from the same run?" =
              setequal(top20[[k]], ranked[cluster == k & rank_in_cluster <= 20, gene]))
}
message("Check passed: ", length(top20), " clusters; xlsx top 20 = CSV ranks 1-20")


## =============================================================================
## CHUNK 4: enricher per cluster
## =============================================================================
run_cluster <- function(genes) {
  g <- intersect(split_genes(genes), universe)                              # markers in the background
  if (length(g) == 0) return(NULL)
  e <- enricher(g, TERM2GENE = t2g, universe = universe, minGSSize = min_set, maxGSSize = 5000,
                pvalueCutoff = 1, qvalueCutoff = 1)                         # keep all; call uses padj_max
  if (is.null(e) || nrow(as.data.frame(e)) == 0) return(NULL)
  r <- as.data.table(as.data.frame(e))
  r[, n_markers_in_background := length(g)]
  r
}

run_and_save <- function(marker_lists, label, out_file) {
  all <- rbindlist(lapply(names(marker_lists), function(k) {
    r <- run_cluster(marker_lists[[k]])
    if (is.null(r)) NULL else cbind(cluster = k, r)
  }), fill = TRUE)
  setorder(all, cluster, pvalue)
  
  smry <- rbindlist(lapply(names(marker_lists), function(k) {             # best match + next two per cluster
    r <- all[cluster == k]
    if (nrow(r) == 0) return(data.table(cluster = k, call = "nothing tested"))
    data.table(cluster = k, n_markers_in_background = r$n_markers_in_background[1],
               call = if (r$p.adjust[1] < padj_max) r$ID[1] else "no significant match",
               best_match = r$ID[1], n_shared = r$Count[1],
               shared_genes = gsub("/", ", ", r$geneID[1]), p_adj = signif(r$p.adjust[1], 3),
               second = if (nrow(r) > 1) paste0(r$ID[2], " (", r$Count[2], " shared, padj ", signif(r$p.adjust[2], 2), ")") else NA,
               third  = if (nrow(r) > 2) paste0(r$ID[3], " (", r$Count[3], " shared, padj ", signif(r$p.adjust[3], 2), ")") else NA)
  }), fill = TRUE)[order(as.numeric(cluster))]
  print(smry[, .(cluster, call, n_shared, p_adj, shared_genes)])
  
  keep <- intersect(c("cluster", "ID", "GeneRatio", "BgRatio", "Count", "pvalue", "p.adjust", "geneID",
                      "n_markers_in_background"), names(all))
  how_to_read <- data.frame(item = c("Markers used", "Reference", "Method", "call", "n_shared / shared_genes", "p_adj"),
                            meaning = c(label,
                                        paste0("HCA Data Explorer project AneurysmalHumanAorticTissue (Li et al. 2020 Circulation; GSE155468), ascending aorta, aneurysm + control; authors' labels (", ref_level, "); top Wilcoxon markers per cell type on your panel genes"),
                                        "clusterProfiler::enricher (hypergeometric); background = panel genes measured in the reference",
                                        paste0("best_match if p_adj < ", padj_max, ", otherwise 'no significant match'"),
                                        "cluster markers that are reference markers of that cell type",
                                        "BH-adjusted p-value across the reference cell types tested for that cluster"))
  write.xlsx(list(Summary = as.data.frame(smry), All_tests = as.data.frame(all[, ..keep]),
                  Reference_markers = as.data.frame(ref), How_to_read = how_to_read),
             out_file, firstRow = TRUE, colWidths = "auto", overwrite = TRUE)
  message("Saved: ", out_file)
}


## =============================================================================
## CHUNK 5: run - top 20 and top 50
## =============================================================================
run_and_save(top20, "top 20 per cluster (xlsx sheets; ranked by pct.1, then avg_log2FC)",
             file.path(out_dir, paste0("HCA_DataExplorer_top20_", ref_level, ".xlsx")))
run_and_save(top50, "top 50 per cluster (ranked CSV, rank_in_cluster <= 50)",
             file.path(out_dir, paste0("HCA_DataExplorer_top50_", ref_level, ".xlsx")))





















## =============================================================================
## 19A_HCA_enricher_top20_top50.R
## =============================================================================
## WHAT IT DOES (same approach as script 17 CHUNKS 4-5, reference = HCA aorta)
##   Reference: marker genes per cell type from the HCA Data Explorer human aorta
##     dataset (Li et al. 2020), computed on your panel genes by
##     19_HCA_aorta_build_reference.R (run that first).
##   Test: each cluster's markers vs each reference cell type's markers,
##     hypergeometric (clusterProfiler::enricher), BH-adjusted per cluster.
##     top 20 (xlsx sheets)      -> HCA_DataExplorer_top20_<ref_level>.xlsx
##     top 50 (CSV, rank <= 50)  -> HCA_DataExplorer_top50_<ref_level>.xlsx
##   ref_level (CHUNK 1): "celltype" = authors' major types, "celltype2" = types/states.
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 1: settings
## =============================================================================
library(data.table)       # tables
library(openxlsx)         # read the top-20 sheets + write results
library(clusterProfiler)  # enricher()

## xlsx_top20 <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_top20_per_cluster.xlsx"   # <- change here
xlsx_top20 <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_05_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_top20_per_cluster.xlsx"   # <- change here
## csv_ranked <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_all_significant_ranked.csv"  # <- change here
csv_ranked <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_05_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_all_significant_ranked.csv"  # <- change here

ref_level <- "celltype"   # "celltype" = authors' major types; "celltype2" = authors' types/states   <- change here
min_set  <- 3      # a reference cell type is tested only if >= 3 of its markers are in the background
padj_max <- 0.05   # a call is made only if BH-adjusted p < this

stopifnot("top-20 xlsx not found" = file.exists(xlsx_top20), "ranked CSV not found" = file.exists(csv_ranked))
out_dir <- dirname(xlsx_top20)                                              # save next to the inputs
ref_csv <- file.path(out_dir, paste0("HCA_reference_markers_aorta_", ref_level, ".csv"))
stopifnot("run 19_HCA_aorta_build_reference.R first" = file.exists(ref_csv))


## =============================================================================
## CHUNK 2: HCA aorta reference
## =============================================================================
split_probe <- function(g) {                                                # "IGHG1/2" -> "IGHG1", "IGHG2"
  parts <- strsplit(g, "/", fixed = TRUE)[[1]]
  if (length(parts) == 1) return(g)
  c(parts[1], paste0(sub("[0-9]+$", "", parts[1]), parts[-1]))
}
split_genes <- function(gs) unique(unlist(lapply(gs, split_probe)))

ref <- fread(ref_csv)                                                       # cell_type, rank, gene, ...
universe <- readLines(sub("\\.csv$", "_genes_in_reference.txt", ref_csv))   # panel genes the reference measured
t2g <- unique(ref[gene %in% universe, .(term = cell_type, gene)])           # TERM2GENE for enricher
print(ref[, .(n_markers = .N, n_cells_used = n_cells_used[1], n_samples = n_samples[1]), by = cell_type])
message(length(unique(t2g$term)), " reference cell types; background ", length(universe), " genes")


## =============================================================================
## CHUNK 3: your marker lists (top 20 from xlsx, top 50 from CSV)
## =============================================================================
sheets <- getSheetNames(xlsx_top20)                                         # "cluster_0", "cluster_1", ...
top20 <- lapply(sheets, function(s) {
  d <- read.xlsx(xlsx_top20, sheet = s)
  if ("gene" %in% names(d)) d$gene else character(0)                        # sheet may say "no significant markers"
})
names(top20) <- sub("^cluster_", "", sheets)

ranked <- fread(csv_ranked)
ranked[, cluster := as.character(cluster)]
top50 <- lapply(names(top20), function(k) ranked[cluster == k & rank_in_cluster <= 50, gene])
names(top50) <- names(top20)

for (k in names(top20)) {                                                   # check: same script-14 run
  stopifnot("xlsx and CSV disagree - are they from the same run?" =
              setequal(top20[[k]], ranked[cluster == k & rank_in_cluster <= 20, gene]))
}
message("Check passed: ", length(top20), " clusters; xlsx top 20 = CSV ranks 1-20")


## =============================================================================
## CHUNK 4: enricher per cluster
## =============================================================================
run_cluster <- function(genes) {
  g <- intersect(split_genes(genes), universe)                              # markers in the background
  if (length(g) == 0) return(NULL)
  e <- enricher(g, TERM2GENE = t2g, universe = universe, minGSSize = min_set, maxGSSize = 5000,
                pvalueCutoff = 1, qvalueCutoff = 1)                         # keep all; call uses padj_max
  if (is.null(e) || nrow(as.data.frame(e)) == 0) return(NULL)
  r <- as.data.table(as.data.frame(e))
  r[, n_markers_in_background := length(g)]
  r
}

run_and_save <- function(marker_lists, label, out_file) {
  all <- rbindlist(lapply(names(marker_lists), function(k) {
    r <- run_cluster(marker_lists[[k]])
    if (is.null(r)) return(NULL)
    r[, cluster := k]; setcolorder(r, "cluster"); r                         # cluster id as first column
  }), fill = TRUE)
  setorder(all, cluster, pvalue)
  
  smry <- rbindlist(lapply(names(marker_lists), function(k) {             # best match + next two per cluster
    r <- all[cluster == k]
    if (nrow(r) == 0) return(data.table(cluster = k, call = "nothing tested"))
    data.table(cluster = k, n_markers_in_background = r$n_markers_in_background[1],
               call = if (r$p.adjust[1] < padj_max) r$ID[1] else "no significant match",
               best_match = r$ID[1], n_shared = r$Count[1],
               shared_genes = gsub("/", ", ", r$geneID[1]), p_adj = signif(r$p.adjust[1], 3),
               second = if (nrow(r) > 1) paste0(r$ID[2], " (", r$Count[2], " shared, padj ", signif(r$p.adjust[2], 2), ")") else NA,
               third  = if (nrow(r) > 2) paste0(r$ID[3], " (", r$Count[3], " shared, padj ", signif(r$p.adjust[3], 2), ")") else NA)
  }), fill = TRUE)[order(as.numeric(cluster))]
  print(smry[, .(cluster, call, n_shared, p_adj, shared_genes)])
  
  keep <- intersect(c("cluster", "ID", "GeneRatio", "BgRatio", "Count", "pvalue", "p.adjust", "geneID",
                      "n_markers_in_background"), names(all))
  how_to_read <- data.frame(item = c("Markers used", "Reference", "Method", "call", "n_shared / shared_genes", "p_adj"),
                            meaning = c(label,
                                        paste0("HCA Data Explorer project AneurysmalHumanAorticTissue (Li et al. 2020 Circulation; GSE155468), ascending aorta, aneurysm + control; authors' labels (", ref_level, "); top Wilcoxon markers per cell type on your panel genes"),
                                        "clusterProfiler::enricher (hypergeometric); background = panel genes measured in the reference",
                                        paste0("best_match if p_adj < ", padj_max, ", otherwise 'no significant match'"),
                                        "cluster markers that are reference markers of that cell type",
                                        "BH-adjusted p-value across the reference cell types tested for that cluster"))
  write.xlsx(list(Summary = as.data.frame(smry), All_tests = as.data.frame(all[, ..keep]),
                  Reference_markers = as.data.frame(ref), How_to_read = how_to_read),
             out_file, firstRow = TRUE, colWidths = "auto", overwrite = TRUE)
  message("Saved: ", out_file)
}


## =============================================================================
## CHUNK 5: run - top 20 and top 50
## =============================================================================
run_and_save(top20, "top 20 per cluster (xlsx sheets; ranked by pct.1, then avg_log2FC)",
             file.path(out_dir, paste0("HCA_DataExplorer_top20_", ref_level, ".xlsx")))
run_and_save(top50, "top 50 per cluster (ranked CSV, rank_in_cluster <= 50)",
             file.path(out_dir, paste0("HCA_DataExplorer_top50_", ref_level, ".xlsx")))
























## =============================================================================
## 19B_HCA_count_top10_min5.R
## =============================================================================
## WHAT IT DOES (same approach as script 17 CHUNK 6, reference = HCA aorta)
##   Reference: marker genes per cell type from the HCA Data Explorer human aorta
##     dataset (Li et al. 2020), computed on your panel genes by
##     19_HCA_aorta_build_reference.R (run that first).
##   Rule (no statistics): a cluster is named EVERY reference cell type for which
##     >= min_hits of its top_k genes are reference markers of that type.
##     -> HCA_DataExplorer_top10_min5_<ref_level>.xlsx
##   ref_level (CHUNK 1): "celltype" = authors' major types, "celltype2" = types/states.
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 1: settings
## =============================================================================
library(data.table)   # tables
library(openxlsx)     # write results


## csv_ranked <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_all_significant_ranked.csv"  # <- change here
csv_ranked <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_05_06_07_08_09_final_All_markers/FindAllMarkers_harmony_PC20_NN250_rep3_res1_all_significant_ranked.csv"  # <- change here
ref_level <- "celltype"   # "celltype" = authors' major types; "celltype2" = authors' types/states   <- change here
top_k    <- 10     # top genes per cluster (rank_in_cluster 1-10)        <- change here
min_hits <- 5      # how many must be reference markers of the cell type  <- change here

stopifnot("ranked CSV not found" = file.exists(csv_ranked))
out_dir <- dirname(csv_ranked)                                              # save next to the inputs
ref_csv <- file.path(out_dir, paste0("HCA_reference_markers_aorta_", ref_level, ".csv"))
stopifnot("run 19_HCA_aorta_build_reference.R first" = file.exists(ref_csv))


## =============================================================================
## CHUNK 2: HCA aorta reference
## =============================================================================
split_probe <- function(g) {                                                # "IGHG1/2" -> "IGHG1", "IGHG2"
  parts <- strsplit(g, "/", fixed = TRUE)[[1]]
  if (length(parts) == 1) return(g)
  c(parts[1], paste0(sub("[0-9]+$", "", parts[1]), parts[-1]))
}
ref <- fread(ref_csv)
universe <- readLines(sub("\\.csv$", "_genes_in_reference.txt", ref_csv))   # panel genes the reference measured
sets <- split(ref[gene %in% universe, gene], ref[gene %in% universe, cell_type])   # markers per reference type
N <- length(universe)
message(length(sets), " reference cell types; background ", N, " genes")


## =============================================================================
## CHUNK 3: counting rule - >= min_hits of the top_k genes are markers of a type
## =============================================================================
ranked <- fread(csv_ranked)
ranked[, cluster := as.character(cluster)]
clusters <- as.character(sort(unique(as.numeric(ranked$cluster))))         # 0, 1, 2, ...

counts <- rbindlist(lapply(clusters, function(k) {
  g <- ranked[cluster == k & rank_in_cluster <= top_k][order(rank_in_cluster), gene]   # top_k genes in rank order
  rbindlist(lapply(names(sets), function(ct) {
    s   <- sets[[ct]]
    hit <- g[vapply(g, function(x) any(split_probe(x) %in% s), logical(1))]   # IGHG1/2 counts once
    data.table(cluster = k, cell_type = ct, n_top_genes = length(g), n_hits = length(hit),
               hit_genes = paste(hit, collapse = ", "), markers_of_type = length(s),
               expected_by_chance = round(length(g) * length(s) / N, 2))     # random genes would hit this many
  }))
}))
setorder(counts, cluster, -n_hits, expected_by_chance)
passed <- counts[n_hits >= min_hits]                                        # every type passing the rule

calls <- rbindlist(lapply(clusters, function(k) {
  p <- passed[cluster == k]
  b <- counts[cluster == k][1]                                              # best type even if it failed
  data.table(cluster = k,
             top_genes = paste(ranked[cluster == k & rank_in_cluster <= top_k][order(rank_in_cluster), gene], collapse = ", "),
             n_cell_types_called = nrow(p),
             cell_types_called = if (nrow(p) > 0) paste0(p$cell_type, " (", p$n_hits, "/", p$n_top_genes, ")", collapse = "; ")
             else paste0("none (< ", min_hits, " hits)"),
             best_cell_type = b$cell_type, best_n_hits = b$n_hits, best_hit_genes = b$hit_genes)
}))
print(calls[, .(cluster, n_cell_types_called, cell_types_called)])


## =============================================================================
## CHUNK 4: save
## =============================================================================
how_to_read <- data.frame(item = c("Rule", "Reference", "cell_types_called", "best_cell_type",
                                   "expected_by_chance", "Caution"),
                          meaning = c(paste0("a cluster is named a cell type if >= ", min_hits, " of its top ", top_k,
                                             " genes are reference markers of that type; several types allowed"),
                                      paste0("HCA Data Explorer project AneurysmalHumanAorticTissue (Li et al. 2020 Circulation; GSE155468), ascending aorta; authors' labels (", ref_level, "); top Wilcoxon markers per reference cell type on your panel genes"),
                                      "every cell type passing the rule, hits/top genes in brackets",
                                      "cell type with the most hits, shown even when it does not pass (near-miss)",
                                      "hits expected from the same number of random panel genes",
                                      "No statistical test: trust a call only if n_hits is well above expected_by_chance and the hit genes are canonical for that type."))
out_file <- file.path(out_dir, paste0("HCA_DataExplorer_top", top_k, "_min", min_hits, "_", ref_level, ".xlsx"))
write.xlsx(list(Calls = as.data.frame(calls),
                Passing_cell_types = as.data.frame(passed[order(as.numeric(cluster), -n_hits)]),
                All_counts = as.data.frame(counts[n_hits > 0][order(as.numeric(cluster), -n_hits)]),
                Reference_markers = as.data.frame(ref), How_to_read = how_to_read),
           out_file, firstRow = TRUE, colWidths = "auto", overwrite = TRUE)
message("Saved: ", out_file)




































## =============================================================================
## 15_InSituType_diagnostics.R
## =============================================================================
## =============================================================================
## 15_InSituType_diagnostics.R
## =============================================================================
## =============================================================================
## 15_InSituType_diagnostics.R
## =============================================================================

# Install the required installation tools
install.packages("remotes")

# Update key dependencies as recommended by the developers
install.packages("Matrix", type = "source")
install.packages("irlba", type = "source")

# Install InSituType from GitHub
remotes::install_github("Nanostring-Biostats/InSituType")


install.packages("Matrix", type = "source"); install.packages("irlba", type = "source")
devtools::install_github("https://github.com/Nanostring-Biostats/InSituType")

## =============================================================================
## 20_InSituType_unsupervised_clustering_markers.R
## =============================================================================
## WHAT IT DOES
##   1. Unsupervised InSituType (Danaher et al. 2022) on the QC-filtered cells:
##        insitutype(x = raw counts (cells x genes), neg = mean negprobe counts per
##                   cell, n_clusts = range -> best number chosen by InSituType)
##      InSituType has NO PCs / neighbours / resolution: it fits a negative-
##      binomial mixture to raw counts with negprobe background. n_clusts is the
##      only clustering choice. (It uses 20 PCs internally only to pick
##      representative cell subsets, not to define clusters.)
##   2. Plots: the 12b UMAP (20 PCs, 250 neighbours) coloured by InSituType
##      cluster, cells in tissue, InSituType flightpath, InSituType vs Leiden.
##   3. Genes enriched in each cluster vs all other cells: Seurat FindAllMarkers
##      (Wilcoxon, only.pos = TRUE); top 20 per cluster.
##   Clusters are renumbered 0, 1, 2, ... by size (as Leiden in 12b); the
##   InSituType letters are kept in a mapping table.
##   Output file names match script 14, so scripts 17 / 18A-B / 19A-B can be
##   pointed at them to name the clusters.
##
## INSTALL (once, in R):
##   install.packages("Matrix", type = "source"); install.packages("irlba", type = "source")
##   devtools::install_github("https://github.com/Nanostring-Biostats/InSituType")
## OUTPUT: <AAA>/Unsupervised_clustering_folder/<QC choice>/InSituType_<run>/
## Run chunk by chunk (Cmd+Enter). CHUNK 4 can take a long time (tens of minutes).
## =============================================================================

## 1. Bioconductor dependencies + irlba
install.packages("BiocManager")
BiocManager::install(c("SummarizedExperiment", "sparseMatrixStats", "SingleCellExperiment"),
                     update = FALSE, ask = FALSE)
install.packages("irlba", type = "source")

## 2. compile InSituType as C++14 (temporary setting, only for this install)
mk <- tempfile(fileext = ".mk")
writeLines("CXX11STD = -std=gnu++14", mk)       # makes the "C++11" build use C++14
Sys.setenv(R_MAKEVARS_USER = mk)
remotes::install_github("Nanostring-Biostats/InSituType", upgrade = "never")
Sys.unsetenv("R_MAKEVARS_USER")                 # back to normal for other packages

## 3. check
library(InSituType)
packageVersion("InSituType")                    # expect 2.0

BiocManager::install(c("SummarizedExperiment", "sparseMatrixStats", "SingleCellExperiment"), update = FALSE, ask = FALSE)


install.packages('SingleCellExperiment')

library(SingleCellExperiment)   # check: should load with no error
r
mk <- tempfile(fileext = ".mk"); writeLines("CXX11STD = -std=gnu++14", mk); Sys.setenv(R_MAKEVARS_USER = mk)
r
remotes::install_github("Nanostring-Biostats/InSituType", upgrade = "never")
r
Sys.unsetenv("R_MAKEVARS_USER")
library(InSituType); packageVersion("InSituType") 


renv::install(c("bioc::SummarizedExperiment", "bioc::sparseMatrixStats", "bioc::SingleCellExperiment"))

library(InSituType); packageVersion("InSituType")   # expect 2.0, no error
library(InSituType); packageVersion("InSituType")   # expect 2.0, no error

## =============================================================================
## CHUNK 1: settings
## =============================================================================
library(data.table)   # tables
library(Matrix)       # sparse counts
library(Seurat)       # object, DimPlot, FindAllMarkers
library(InSituType)   # insitutype(), flightpath_layout()
library(ggplot2)      # plots
library(openxlsx)     # one Excel sheet per cluster

## same path as script 14: the 12b clustering whose cells + UMAP are used
##clustering_path <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Clustering_PCA_UMAP_Leiden_Markers/QC03_04_05_06_07_08_09/PC20_NN250_rep3/02_SENSITIVITY_batchawarePCA_Harmony/UMAP_by_Leiden_cluster.tiff"   # <- change here
clustering_path <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Clustering_PCA_UMAP_Leiden_Markers/QC03_04_06_07_08_09/PC20_NN250_rep3/02_SENSITIVITY_batchawarePCA_Harmony/UMAP_by_Leiden_cluster.tiff" 


n_clusts  <- 14:20    # range: InSituType chooses the best number; or one number, e.g. 18   <- change here
n_starts  <- 10       # InSituType default (random starts in phase 1)
seed      <- 42
marker_padj_max <- 0.05        # significance cut-off (Bonferroni, as script 14)
marker_top_n    <- 20          # genes per cluster sheet
rank_by   <- "avg_log2FC"      # "avg_log2FC" (most specific first) or "pct.1" (as script 14)   <- change here


## =============================================================================
## CHUNK 2: paths from clustering_path (as script 14)
## =============================================================================
branch_dir <- if (dir.exists(clustering_path)) clustering_path else dirname(clustering_path)
branch_dir <- normalizePath(branch_dir, mustWork = TRUE)                    # stop if missing
dir_seu    <- dirname(branch_dir)                                           # .../PC20_NN250_rep3
run_tag    <- basename(dir_seu)                                             # "PC20_NN250_rep3"
qc_tag     <- basename(dirname(dir_seu))                                    # "QC03_04_05_06_07_08_09"
stopifnot("path is not inside Clustering_PCA_UMAP_Leiden_Markers/<QC>/<run>/<branch>/" =
            basename(dirname(dirname(dir_seu))) == "Clustering_PCA_UMAP_Leiden_Markers")
output_dir <- dirname(dirname(dirname(dir_seu)))                            # .../AAA
prefix <- if (grepl("^01_PRIMARY", basename(branch_dir))) "primary" else
  if (grepl("^02_SENSITIVITY", basename(branch_dir))) "harmony" else
    stop("branch folder must start with 01_PRIMARY or 02_SENSITIVITY")
params <- fread(file.path(dir_seu, paste0("parameters_used_", run_tag, ".csv")))
leiden_resolution <- as.numeric(params[parameter == "Leiden_resolution", value])
leiden_col <- paste0(prefix, "_leiden_res_", leiden_resolution)            # Leiden clusters, for comparison only
umap_name  <- paste0("umap_", prefix)                                       # 12b UMAP (20 PCs, 250 NN)

seurat_out <- list.files(dir_seu, pattern = paste0("^seurat_.*_SeuratFunctions_", run_tag, "\\.RDS$"),
                         full.names = TRUE)
stopifnot("expected exactly one 12b Seurat object" = length(seurat_out) == 1)

dir_out  <- file.path(output_dir, "Unsupervised_clustering_folder", qc_tag, paste0("InSituType_", run_tag))
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
file_tag <- paste0("insitutype_", run_tag)                                  # e.g. "insitutype_PC20_NN250_rep3"
message("QC choice: ", qc_tag, " | UMAP: ", umap_name, " | results -> ", dir_out)


## =============================================================================
## CHUNK 3: raw counts (cells x genes) and negprobe background
## =============================================================================
obj <- readRDS(seurat_out)                                                  # QC-filtered cells + 12b UMAPs
stopifnot("UMAP not in object" = umap_name %in% Reductions(obj),
          "Leiden column not in object" = leiden_col %in% colnames(obj[[]]))

x <- t(LayerData(obj, assay = "RNA", layer = "counts"))                     # cells x genes, raw counts
stopifnot("RNA counts must be raw whole numbers" = all(x@x == round(x@x)),
          "cells with 0 counts - InSituType stops on these" = all(Matrix::rowSums(x) > 0))

## neg = mean negative-probe count per cell (InSituType's background input)
if ("negmean" %in% colnames(obj[[]])) {
  neg <- obj$negmean                                                        # saved by the load script
  message("neg = metadata column negmean")
} else {
  neg_assay <- grep("neg", Assays(obj), ignore.case = TRUE, value = TRUE)
  stopifnot("no negmean column and no negprobe assay found" = length(neg_assay) == 1)
  neg <- Matrix::colMeans(LayerData(obj, assay = neg_assay, layer = "counts"))
  message("neg = mean of assay ", neg_assay)
}
neg <- neg[rownames(x)]                                                     # same cell order as x
stopifnot(!anyNA(neg))
message(nrow(x), " cells x ", ncol(x), " genes; median neg = ", signif(median(neg), 3))


## =============================================================================
## CHUNK 4: unsupervised InSituType
## =============================================================================
set.seed(seed)                                                              # InSituType samples cells
unsup <- insitutype(x = x, neg = neg, assay_type = "rna",
                    n_clusts = n_clusts,                                    # range -> best number chosen
                    n_starts = n_starts)                                    # other arguments: package defaults
saveRDS(unsup, file.path(dir_out, paste0("InSituType_result_", file_tag, ".RDS")))   # clust, prob, logliks, profiles

## renumber clusters 0, 1, 2, ... by size (0 = largest), keep the InSituType letters
cl_letter <- unsup$clust[rownames(x)]                                       # letter per cell, object order
stopifnot(length(cl_letter) == ncol(obj), !anyNA(cl_letter))
sizes   <- sort(table(cl_letter), decreasing = TRUE)
map     <- data.table(insitutype_letter = names(sizes), cluster = seq_along(sizes) - 1L,
                      n_cells = as.integer(sizes))
cl_num  <- factor(map$cluster[match(cl_letter, map$insitutype_letter)], levels = map$cluster)
obj$insitutype_cluster <- cl_num
obj$insitutype_prob    <- unsup$prob[rownames(x)]                           # confidence per cell
message("InSituType clusters: ", nlevels(cl_num))
print(map)
fwrite(map, file.path(dir_out, paste0("InSituType_cluster_letter_to_number_", file_tag, ".csv")))

## cluster mean profiles (expected counts per gene), columns renamed to numbers
prof <- unsup$profiles[, map$insitutype_letter, drop = FALSE]
colnames(prof) <- paste0("cluster_", map$cluster)
fwrite(as.data.table(prof, keep.rownames = "gene"), file.path(dir_out, paste0("InSituType_profiles_", file_tag, ".csv")))

## per-cell table
cells <- data.table(cell_id = colnames(obj), sample_id = obj$sample_id,
                    x_slide_mm = obj$x_slide_mm, y_slide_mm = obj$y_slide_mm,
                    insitutype_cluster = cl_num, insitutype_letter = cl_letter,
                    insitutype_prob = round(obj$insitutype_prob, 4),
                    leiden = obj[[leiden_col]][, 1], Embeddings(obj, umap_name)[, 1:2])
fwrite(cells, file.path(dir_out, paste0("InSituType_cells_", file_tag, ".csv")))


## =============================================================================
## CHUNK 5: plots
## =============================================================================
pal <- c("#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#42D4F4", "#F032E6",
         "#9A6324", "#469990", "#BFEF45", "#800000", "#FFE119", "#000075", "#808000",
         "#FF8FB1", "#A9A9A9", "#005300", "#009FFF", "#C77CFF", "#FF6F61", "#00FFBE",
         "#D4A017", "#7F7FFF", "#8B008B", "#2F4F4F", "#E9967A", "#556B2F", "#FF1493")  # as 12b
stopifnot("more clusters than colours - add colours to pal" = nlevels(cl_num) <= length(pal))
cl_cols <- setNames(pal[seq_len(nlevels(cl_num))], levels(cl_num))
save_tiff <- function(p, f, w, h) ggsave(file.path(dir_out, f), p, width = w, height = h, dpi = 300,
                                         device = "tiff", compression = "lzw", bg = "white")

## 1. 12b UMAP (20 PCs, 250 NN) coloured by InSituType cluster
p_umap <- DimPlot(obj, reduction = umap_name, group.by = "insitutype_cluster", cols = cl_cols,
                  label = TRUE, label.box = TRUE, shuffle = TRUE, seed = seed,
                  pt.size = 0.05, alpha = 0.3, raster = FALSE) +
  coord_fixed() + theme_bw(base_size = 13) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
  labs(title = paste0("InSituType unsupervised: ", nlevels(cl_num), " clusters"),
       subtitle = paste0("UMAP from 12b (", run_tag, ", ", prefix, "); clusters from InSituType, not from the UMAP"),
       colour = "InSituType\ncluster")
save_tiff(p_umap, paste0("01_UMAP_by_InSituType_cluster_", file_tag, ".tiff"), 11, 9)

## 2. clusters in tissue
p_sp <- ggplot(cells[sample(.N)], aes(x_slide_mm, y_slide_mm, colour = insitutype_cluster)) +
  geom_point(size = 0.05, alpha = 0.6) + scale_colour_manual(values = cl_cols) +
  coord_fixed() + theme_bw(base_size = 13) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
  labs(title = "InSituType clusters in tissue", x = "x (mm)", y = "y (mm)", colour = "InSituType\ncluster")
save_tiff(p_sp, paste0("02_Spatial_InSituType_clusters_", file_tag, ".tiff"), 13, 10)

## 3. InSituType flightpath (cells placed by their cluster log-likelihoods)
fp <- tryCatch(flightpath_layout(logliks = unsup$logliks, profiles = unsup$profiles),
               error = function(e) { message("flightpath skipped: ", conditionMessage(e)); NULL })
if (!is.null(fp)) {
  fpd <- data.table(x = fp$cellpos[, 1], y = fp$cellpos[, 2],
                    cluster = cl_num[match(rownames(fp$cellpos), rownames(x))])
  p_fp <- ggplot(fpd[sample(.N)], aes(x, y, colour = cluster)) + geom_point(size = 0.05, alpha = 0.4) +
    scale_colour_manual(values = cl_cols) + theme_bw(base_size = 13) +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
    labs(title = "InSituType flightpath", subtitle = "Cells between two clusters = uncertain assignment",
         x = NULL, y = NULL, colour = "InSituType\ncluster")
  save_tiff(p_fp, paste0("03_Flightpath_", file_tag, ".tiff"), 11, 9)
}

## 4. confidence per cluster
p_prob <- ggplot(cells, aes(insitutype_cluster, insitutype_prob, fill = insitutype_cluster)) +
  geom_violin(scale = "width") + scale_fill_manual(values = cl_cols) + theme_bw(base_size = 13) +
  theme(legend.position = "none") + labs(title = "InSituType confidence per cell", x = "cluster", y = "posterior probability")
save_tiff(p_prob, paste0("04_Confidence_per_cluster_", file_tag, ".tiff"), 12, 5)

## 5. InSituType vs Leiden (12b): % of each InSituType cluster in each Leiden cluster
ari <- igraph::compare(as.integer(cl_num), as.integer(factor(cells$leiden)), method = "adjusted.rand")
ct  <- cells[, .N, by = .(insitutype_cluster, leiden)][, pct := round(100 * N / sum(N), 1), by = insitutype_cluster]
p_ct <- ggplot(ct, aes(factor(leiden, levels = sort(unique(as.numeric(as.character(leiden))))), insitutype_cluster, fill = pct)) +
  geom_tile(colour = "white") + geom_text(aes(label = ifelse(pct >= 5, pct, "")), size = 2.8) +
  scale_fill_viridis_c() + theme_minimal(base_size = 12) +
  labs(title = "InSituType vs Leiden (12b)", subtitle = paste0("Each row = 100% of an InSituType cluster; ARI = ", round(ari, 3)),
       x = "Leiden cluster", y = "InSituType cluster", fill = "% of\nInSituType\ncluster")
save_tiff(p_ct, paste0("05_InSituType_vs_Leiden_", file_tag, ".tiff"), 12, 10)
fwrite(dcast(ct, insitutype_cluster ~ leiden, value.var = "N", fill = 0),
       file.path(dir_out, paste0("InSituType_vs_Leiden_cells_", file_tag, ".csv")))


## =============================================================================
## CHUNK 6: genes enriched in each InSituType cluster vs all other cells
## =============================================================================
## FindAllMarkers reads the "data" layer -> log-normalise first (as script 14);
## this does not change the InSituType clusters (they came from raw counts).
obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
Idents(obj) <- "insitutype_cluster"
markers <- as.data.table(FindAllMarkers(obj, assay = "RNA", only.pos = TRUE, test.use = "wilcox"))
saveRDS(markers, file.path(dir_out, paste0("FindAllMarkers_", file_tag, "_all_results.RDS")))

sig <- markers[p_val_adj < marker_padj_max & avg_log2FC > 0]
if (rank_by == "pct.1") setorder(sig, cluster, -pct.1, -avg_log2FC) else setorder(sig, cluster, -avg_log2FC, -pct.1)
sig[, rank_in_cluster := seq_len(.N), by = cluster]
setcolorder(sig, c("cluster", "rank_in_cluster", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"))
fwrite(sig, file.path(dir_out, paste0("FindAllMarkers_", file_tag, "_all_significant_ranked.csv")))

cl_levels <- levels(cl_num)
sheets <- lapply(cl_levels, function(k) {                                   # one sheet per cluster, top 20
  top <- sig[cluster == k & rank_in_cluster <= marker_top_n]
  if (nrow(top) == 0) data.frame(note = "no significant positive markers") else as.data.frame(top)
})
names(sheets) <- paste0("cluster_", cl_levels)
write.xlsx(sheets, file.path(dir_out, paste0("FindAllMarkers_", file_tag, "_top", marker_top_n, "_per_cluster.xlsx")),
           firstRow = TRUE, colWidths = "auto", overwrite = TRUE)

smry <- data.table(cluster = cl_levels,
                   insitutype_letter = map$insitutype_letter,
                   n_cells = map$n_cells,
                   median_prob = sapply(cl_levels, function(k) round(median(cells[insitutype_cluster == k, insitutype_prob]), 3)),
                   n_significant_markers = sapply(cl_levels, function(k) nrow(sig[cluster == k])),
                   top10 = sapply(cl_levels, function(k) paste(head(sig[cluster == k, gene], 10), collapse = ", ")))
print(smry)
fwrite(smry, file.path(dir_out, paste0("FindAllMarkers_", file_tag, "_summary_per_cluster.csv")))

## dot plot: top 3 markers per cluster
top3 <- unique(sig[rank_in_cluster <= 3, gene])
p_dot <- DotPlot(obj, features = top3, group.by = "insitutype_cluster") + RotatedAxis() +
  labs(title = "Top 3 enriched genes per InSituType cluster")
save_tiff(p_dot, paste0("06_DotPlot_top3_markers_", file_tag, ".tiff"), max(12, length(top3) * 0.3), 8)

fwrite(data.table(parameter = c("input_object", "qc_tag", "method", "n_clusts_range", "n_clusts_chosen",
                                "n_starts", "neg", "seed", "markers", "rank_by", "padj_cutoff"),
                  value = c(seurat_out, qc_tag, "InSituType::insitutype (unsupervised, raw counts, negprobe background)",
                            paste(range(n_clusts), collapse = "-"), nlevels(cl_num), n_starts,
                            "mean negprobe counts per cell", seed,
                            "Seurat::FindAllMarkers(only.pos = TRUE, wilcox) on log-normalised RNA",
                            rank_by, marker_padj_max)),
       file.path(dir_out, paste0("parameters_", file_tag, ".csv")))
writeLines(capture.output(sessionInfo()), file.path(dir_out, "sessionInfo.txt"))
message("Done -> ", dir_out)





## 1b. UMAP built from InSituType's own output (per-cell log-likelihoods for each cluster)
##     = the InSituType counterpart of the 12b UMAP for Leiden (separated by construction)
## 1b. UMAP built from InSituType's own output (per-cell log-likelihoods for each cluster)
##     = the InSituType counterpart of the 12b UMAP for Leiden (separated by construction)
##     Layout settings only - the InSituType clusters are NOT changed by any of these.
umap_ist_nn        <- 30     # n_neighbors: lower = crisper, more separate islands (250 = everything pulled together)  <- change here
umap_ist_min_dist  <- 0.3    # higher = looser points within each cluster                                             <- change here
umap_ist_spread    <- 2      # overall scale of the layout (as 12b)                                                   <- change here
umap_ist_repulsion <- 5      # higher = clusters pushed further apart (as 12b)                                        <- change here
umap_ist_tag <- paste0("NN", umap_ist_nn, "_md", umap_ist_min_dist, "_sp", umap_ist_spread, "_rep", umap_ist_repulsion)

ll <- unsup$logliks[rownames(x), map$insitutype_letter, drop = FALSE]  # cells x clusters, object order
ll <- ll - apply(ll, 1, max)                          # relative to each cell's best cluster (best = 0)
ll[ll < -50] <- -50                                   # cap very unlikely clusters so they do not dominate
set.seed(seed)
um_ist <- uwot::umap(ll, n_neighbors = umap_ist_nn,                      # input = 14 cluster scores per cell, not PCs
                     min_dist = umap_ist_min_dist, spread = umap_ist_spread,
                     repulsion_strength = umap_ist_repulsion,
                     metric = "euclidean")
rownames(um_ist) <- rownames(x); colnames(um_ist) <- c("UMAPIST_1", "UMAPIST_2")
obj[["umap_insitutype"]] <- CreateDimReducObject(embeddings = um_ist, key = "UMAPIST_", assay = "RNA")

p_umap_ist <- DimPlot(obj, reduction = "umap_insitutype", group.by = "insitutype_cluster", cols = cl_cols,
                      label = TRUE, label.box = TRUE, shuffle = TRUE, seed = seed,
                      pt.size = 0.05, alpha = 0.3, raster = FALSE) +
  coord_fixed() + theme_bw(base_size = 13) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
  labs(title = paste0("InSituType unsupervised: ", nlevels(cl_num), " clusters"),
       subtitle = paste0("UMAP of InSituType log-likelihoods (layout only): n_neighbors ", umap_ist_nn,
                         ", min_dist ", umap_ist_min_dist, ", spread ", umap_ist_spread,
                         ", repulsion ", umap_ist_repulsion),
       colour = "InSituType\ncluster")
save_tiff(p_umap_ist, paste0("01b_UMAP_of_InSituType_logliks_", umap_ist_tag, "_", file_tag, ".tiff"), 11, 9)

## record the layout settings for the methods section
fwrite(data.table(parameter = c("UMAP_insitutype_input", "n_neighbors", "min_dist", "spread", "repulsion_strength", "seed"),
                  value = c("uwot::umap on per-cell InSituType log-likelihoods relative to best cluster, capped at -50",
                            umap_ist_nn, umap_ist_min_dist, umap_ist_spread, umap_ist_repulsion, seed)),
       file.path(dir_out, paste0("parameters_UMAP_insitutype_", umap_ist_tag, ".csv")))
























## =============================================================================
## CHUNK 1: settings
## =============================================================================

install.packages('pheatmap')
library(data.table)   # tables
library(Matrix)       # sparse counts
library(Seurat)       # object, DimPlot, FindAllMarkers
library(InSituType)   # insitutype(), fastCohorting(), refineClusters(), flightpath_layout()
library(ggplot2)      # plots
library(openxlsx)     # one Excel sheet per cluster
library(pheatmap)

##clustering_path <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Clustering_PCA_UMAP_Leiden_Markers/QC03_04_05_06_07_08_09/PC20_NN250_rep3/02_SENSITIVITY_batchawarePCA_Harmony/UMAP_by_Leiden_cluster.tiff"   # <- change here
clustering_path <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Clustering_PCA_UMAP_Leiden_Markers/QC03_04_06_07_08_09/PC20_NN250_rep3/02_SENSITIVITY_batchawarePCA_Harmony/UMAP_by_Leiden_cluster.tiff" 

n_clusts  <- 14:20    # range: InSituType chooses the best number; or one number, e.g. 18   <- change here
n_starts  <- 10       # InSituType default (random starts in phase 1)
seed      <- 42
conf_min  <- 0.8      # FAQ: cells with posterior < 0.8 called "unclassified"   <- change here
marker_padj_max <- 0.05        # significance cut-off (Bonferroni, as script 14)
marker_top_n    <- 20          # genes per cluster sheet
rank_by   <- "avg_log2FC"      # "avg_log2FC" (most specific first) or "pct.1" (as script 14)   <- change here

## refineClusters: leave FALSE for the first run. Set TRUE only AFTER inspecting plots 01, 03, 07
## and the markers; then the saved InSituType result is re-loaded (not re-run) and refined.
## Use InSituType LETTERS (see InSituType_cluster_letter_to_number_*.csv).
run_refine        <- FALSE                        # <- change here
refine_merges     <- NULL                         # e.g. c("a" = "SMC", "b" = "SMC")   <- change here
refine_delete     <- NULL                         # e.g. c("d")                        <- change here
refine_subcluster <- NULL                         # e.g. list("e" = 2)                 <- change here


## =============================================================================
## CHUNK 2: paths from clustering_path (as script 14)
## =============================================================================
branch_dir <- if (dir.exists(clustering_path)) clustering_path else dirname(clustering_path)
branch_dir <- normalizePath(branch_dir, mustWork = TRUE)                    # stop if missing
dir_seu    <- dirname(branch_dir)                                           # .../PC20_NN250_rep3
run_tag    <- basename(dir_seu)                                             # "PC20_NN250_rep3"
qc_tag     <- basename(dirname(dir_seu))                                    # "QC03_04_05_06_07_08_09"
stopifnot("path is not inside Clustering_PCA_UMAP_Leiden_Markers/<QC>/<run>/<branch>/" =
            basename(dirname(dirname(dir_seu))) == "Clustering_PCA_UMAP_Leiden_Markers")
output_dir <- dirname(dirname(dirname(dir_seu)))                            # .../AAA
prefix <- if (grepl("^01_PRIMARY", basename(branch_dir))) "primary" else
  if (grepl("^02_SENSITIVITY", basename(branch_dir))) "harmony" else
    stop("branch folder must start with 01_PRIMARY or 02_SENSITIVITY")
params <- fread(file.path(dir_seu, paste0("parameters_used_", run_tag, ".csv")))
leiden_resolution <- as.numeric(params[parameter == "Leiden_resolution", value])
leiden_col <- paste0(prefix, "_leiden_res_", leiden_resolution)            # Leiden clusters, for comparison only
umap_name  <- paste0("umap_", prefix)                                       # 12b UMAP (20 PCs, 250 NN)

seurat_out <- list.files(dir_seu, pattern = paste0("^seurat_.*_SeuratFunctions_", run_tag, "\\.RDS$"),
                         full.names = TRUE)
stopifnot("expected exactly one 12b Seurat object" = length(seurat_out) == 1)

dir_out  <- file.path(output_dir, "Unsupervised_clustering_folder", qc_tag, paste0("InSituType_", run_tag))
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)
file_tag <- paste0("insitutype_", run_tag)                                  # e.g. "insitutype_PC20_NN250_rep3"
unsup_file <- file.path(dir_out, paste0("InSituType_result_", file_tag, ".RDS"))   # raw InSituType result
message("QC choice: ", qc_tag, " | UMAP: ", umap_name, " | results -> ", dir_out)


## =============================================================================
## CHUNK 3: raw counts (cells x genes) and negprobe background
## =============================================================================
obj <- readRDS(seurat_out)                                                  # QC-filtered cells + 12b UMAPs
stopifnot("UMAP not in object" = umap_name %in% Reductions(obj),
          "Leiden column not in object" = leiden_col %in% colnames(obj[[]]))

x <- t(LayerData(obj, assay = "RNA", layer = "counts"))                     # cells x genes, raw counts
stopifnot("RNA counts must be raw whole numbers" = all(x@x == round(x@x)),
          "cells with 0 counts - InSituType stops on these" = all(Matrix::rowSums(x) > 0))

## remove CosMx control probes from x (case-sensitive: real genes such as NEGR1 are kept)
## this object: "NegativeAdd" is removed, "NEGR1" is kept
ctrl <- grep("^Negative|^NegPrb|^SystemControl|^SystemCtrl|^FalseCode", colnames(x), value = TRUE)  # control-probe names only
message("control probes removed from x: ", length(ctrl), if (length(ctrl)) paste0(" (", paste(ctrl, collapse = ", "), ")"))
x <- x[, setdiff(colnames(x), ctrl), drop = FALSE]                          # genes only
stopifnot("NEGR1 must stay in x" = !("NEGR1" %in% rownames(obj)) || "NEGR1" %in% colnames(x),
          "cells with 0 counts after removing controls" = all(Matrix::rowSums(x) > 0))

## neg = mean negative-probe count per cell (InSituType's background input)
if ("negmean" %in% colnames(obj[[]])) {
  neg <- obj$negmean                                                        # saved by the load script
  message("neg = metadata column negmean")
} else {
  neg_assay <- grep("neg", Assays(obj), ignore.case = TRUE, value = TRUE)
  stopifnot("no negmean column and no negprobe assay found" = length(neg_assay) == 1)
  neg <- Matrix::colMeans(LayerData(obj, assay = neg_assay, layer = "counts"))
  message("neg = mean of assay ", neg_assay)
}
neg <- neg[rownames(x)]                                                     # same cell order as x
stopifnot(!anyNA(neg))
message(nrow(x), " cells x ", ncol(x), " genes; median neg = ", signif(median(neg), 3))


## =============================================================================
## CHUNK 4: unsupervised InSituType (as vignette: counts, negmean, cohort, n_clusts)
## =============================================================================
## cohorts from immunofluorescence + area (vignette: fastCohorting)
if_cols <- grep("^Mean\\.|^Area", colnames(obj[[]]), value = TRUE)         # AtoMx IF means + cell area   <- check names
message("cohort features: ", paste(if_cols, collapse = ", "))
stopifnot("no IF/Area columns in metadata" = length(if_cols) > 0)
immunofluordata <- as.matrix(obj[[]][rownames(x), if_cols])                 # cells x features, same order as x
stopifnot("NA in IF/Area columns" = !anyNA(immunofluordata))

if (run_refine) {
  ## refinement: re-use the saved first run (clusters must not change before refining)
  stopifnot("run the unrefined script first (run_refine = FALSE)" = file.exists(unsup_file))
  unsup <- readRDS(unsup_file)                                              # original InSituType result
  message("loaded saved InSituType result: ", unsup_file)
} else {
  set.seed(seed)                                                            # fastCohorting uses random starts
  cohort <- fastCohorting(immunofluordata, gaussian_transform = TRUE)       # one cohort label per cell
  names(cohort) <- rownames(x)
  print(table(cohort))                                                      # cohort sizes
  
  set.seed(seed)                                                            # InSituType subsamples cells
  unsup <- insitutype(x = x,                                                # raw counts, cells x genes
                      neg = neg,                                            # mean negprobe per cell
                      cohort = cohort,                                      # IF-based cohorts (vignette)
                      assay_type = "rna",                                   # InSituType >= 2.0 only
                      reference_profiles = NULL,                            # NULL = unsupervised
                      n_clusts = n_clusts,                                  # range -> best number chosen
                      n_starts = n_starts)                                  # phases/iters: package defaults
  saveRDS(unsup, unsup_file)                                                # clust, prob, logliks, profiles
}

## optional refineClusters (FAQ: generous n_clusts, then condense)
if (run_refine) {
  res <- refineClusters(logliks    = unsup$logliks,                         # from the original run
                        merges     = refine_merges,                         # letters -> new shared name
                        to_delete  = refine_delete,                         # letters to drop
                        subcluster = refine_subcluster,                     # letters to split
                        counts     = x,                                     # raw counts, for new profiles
                        neg        = neg)                                   # background, for new profiles
  file_tag <- paste0(file_tag, "_refined")                                  # refined outputs do not overwrite
  saveRDS(res, file.path(dir_out, paste0("InSituType_result_", file_tag, ".RDS")))
  message("refined clusters: ", paste(names(table(res$clust)), collapse = ", "))
} else {
  res <- unsup                                                              # everything below uses res
}

## renumber clusters 0, 1, 2, ... by size (0 = largest), keep the InSituType letters (or refined names)
cl_letter <- res$clust[rownames(x)]                                         # letter per cell, object order
stopifnot(length(cl_letter) == ncol(obj), !anyNA(cl_letter))
sizes   <- sort(table(cl_letter), decreasing = TRUE)
map     <- data.table(insitutype_letter = names(sizes), cluster = seq_along(sizes) - 1L,
                      n_cells = as.integer(sizes))
cl_num  <- factor(map$cluster[match(cl_letter, map$insitutype_letter)], levels = map$cluster)
obj$insitutype_cluster <- cl_num
obj$insitutype_prob    <- res$prob[rownames(x)]                             # confidence per cell
message("InSituType clusters: ", nlevels(cl_num))
print(map)
fwrite(map, file.path(dir_out, paste0("InSituType_cluster_letter_to_number_", file_tag, ".csv")))

## FAQ: posteriors are overconfident -> flag cells < conf_min (cluster kept, flag added)
obj$insitutype_confident <- obj$insitutype_prob >= conf_min                # TRUE = confident call
message("cells with prob < ", conf_min, ": ", sum(!obj$insitutype_confident))

## FAQ: cells below conf_min called "unclassified" (separate column; original clusters kept)
obj$insitutype_cluster_conf <- as.character(obj$insitutype_cluster)        # copy cluster numbers
obj$insitutype_cluster_conf[!obj$insitutype_confident] <- "unclassified"   # low-confidence cells
obj$insitutype_cluster_conf <- factor(obj$insitutype_cluster_conf,
                                      levels = c(levels(cl_num), "unclassified"))   # keep cluster order
print(table(obj$insitutype_cluster_conf))                                  # cells per label

## cluster mean profiles (expected counts per gene), columns renamed to numbers
prof <- res$profiles[, map$insitutype_letter, drop = FALSE]
colnames(prof) <- paste0("cluster_", map$cluster)
fwrite(as.data.table(prof, keep.rownames = "gene"), file.path(dir_out, paste0("InSituType_profiles_", file_tag, ".csv")))

## per-cell table
cells <- data.table(cell_id = colnames(obj), sample_id = obj$sample_id,
                    x_slide_mm = obj$x_slide_mm, y_slide_mm = obj$y_slide_mm,
                    insitutype_cluster = cl_num, insitutype_letter = cl_letter,
                    insitutype_prob = round(obj$insitutype_prob, 4),
                    insitutype_confident = obj$insitutype_confident,        # prob >= conf_min
                    insitutype_cluster_conf = obj$insitutype_cluster_conf,  # cluster or "unclassified"
                    leiden = obj[[leiden_col]][, 1], Embeddings(obj, umap_name)[, 1:2])
fwrite(cells, file.path(dir_out, paste0("InSituType_cells_", file_tag, ".csv")))


## =============================================================================
## CHUNK 5: plots
## =============================================================================
pal <- c("#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#42D4F4", "#F032E6",
         "#9A6324", "#469990", "#BFEF45", "#800000", "#FFE119", "#000075", "#808000",
         "#FF8FB1", "#A9A9A9", "#005300", "#009FFF", "#C77CFF", "#FF6F61", "#00FFBE",
         "#D4A017", "#7F7FFF", "#8B008B", "#2F4F4F", "#E9967A", "#556B2F", "#FF1493")  # as 12b
stopifnot("more clusters than colours - add colours to pal" = nlevels(cl_num) <= length(pal))
cl_cols      <- setNames(pal[seq_len(nlevels(cl_num))], levels(cl_num))
cl_cols_conf <- c(cl_cols, unclassified = "grey80")                         # + colour for low-confidence cells
save_tiff <- function(p, f, w, h) ggsave(file.path(dir_out, f), p, width = w, height = h, dpi = 300,
                                         device = "tiff", compression = "lzw", bg = "white")

## 1. 12b UMAP (20 PCs, 250 NN) coloured by InSituType cluster
p_umap <- DimPlot(obj, reduction = umap_name, group.by = "insitutype_cluster", cols = cl_cols,
                  label = TRUE, label.box = TRUE, shuffle = TRUE, seed = seed,
                  pt.size = 0.05, alpha = 0.3, raster = FALSE) +
  coord_fixed() + theme_bw(base_size = 13) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
  labs(title = paste0("InSituType unsupervised: ", nlevels(cl_num), " clusters"),
       subtitle = paste0("UMAP from 12b (", run_tag, ", ", prefix, "); clusters from InSituType, not from the UMAP"),
       colour = "InSituType\ncluster")
save_tiff(p_umap, paste0("01_UMAP_by_InSituType_cluster_", file_tag, ".tiff"), 11, 9)

## 1c. same UMAP, cells with prob < conf_min shown as "unclassified" (grey)
p_umap_conf <- DimPlot(obj, reduction = umap_name, group.by = "insitutype_cluster_conf", cols = cl_cols_conf,
                       label = TRUE, label.box = TRUE, shuffle = TRUE, seed = seed,
                       pt.size = 0.05, alpha = 0.3, raster = FALSE) +
  coord_fixed() + theme_bw(base_size = 13) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
  labs(title = paste0("InSituType clusters, prob < ", conf_min, " = unclassified"),
       subtitle = paste0(sum(!obj$insitutype_confident), " of ", ncol(obj), " cells unclassified"),
       colour = "InSituType\ncluster")
save_tiff(p_umap_conf, paste0("01c_UMAP_by_InSituType_cluster_conf", conf_min, "_", file_tag, ".tiff"), 11, 9)

## 2. clusters in tissue
p_sp <- ggplot(cells[sample(.N)], aes(x_slide_mm, y_slide_mm, colour = insitutype_cluster)) +
  geom_point(size = 0.05, alpha = 0.6) + scale_colour_manual(values = cl_cols) +
  coord_fixed() + theme_bw(base_size = 13) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
  labs(title = "InSituType clusters in tissue", x = "x (mm)", y = "y (mm)", colour = "InSituType\ncluster")
save_tiff(p_sp, paste0("02_Spatial_InSituType_clusters_", file_tag, ".tiff"), 13, 10)

## 3. InSituType flightpath (cells placed by their cluster log-likelihoods)
fp <- tryCatch(flightpath_layout(logliks = res$logliks, profiles = res$profiles),
               error = function(e) { message("flightpath skipped: ", conditionMessage(e)); NULL })
if (!is.null(fp)) {
  fpd <- data.table(x = fp$cellpos[, 1], y = fp$cellpos[, 2],
                    cluster = cl_num[match(rownames(fp$cellpos), rownames(x))])
  p_fp <- ggplot(fpd[sample(.N)], aes(x, y, colour = cluster)) + geom_point(size = 0.05, alpha = 0.4) +
    scale_colour_manual(values = cl_cols) + theme_bw(base_size = 13) +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
    labs(title = "InSituType flightpath", subtitle = "Cells between two clusters = uncertain assignment",
         x = NULL, y = NULL, colour = "InSituType\ncluster")
  save_tiff(p_fp, paste0("03_Flightpath_", file_tag, ".tiff"), 11, 9)
}

## 4. confidence per cluster (dashed line = conf_min)
p_prob <- ggplot(cells, aes(insitutype_cluster, insitutype_prob, fill = insitutype_cluster)) +
  geom_violin(scale = "width") + scale_fill_manual(values = cl_cols) +
  geom_hline(yintercept = conf_min, linetype = "dashed") + theme_bw(base_size = 13) +
  theme(legend.position = "none") + labs(title = "InSituType confidence per cell", x = "cluster", y = "posterior probability")
save_tiff(p_prob, paste0("04_Confidence_per_cluster_", file_tag, ".tiff"), 12, 5)

## 5. InSituType vs Leiden (12b): % of each InSituType cluster in each Leiden cluster
ari <- igraph::compare(as.integer(cl_num), as.integer(factor(cells$leiden)), method = "adjusted.rand")
ct  <- cells[, .N, by = .(insitutype_cluster, leiden)][, pct := round(100 * N / sum(N), 1), by = insitutype_cluster]
p_ct <- ggplot(ct, aes(factor(leiden, levels = sort(unique(as.numeric(as.character(leiden))))), insitutype_cluster, fill = pct)) +
  geom_tile(colour = "white") + geom_text(aes(label = ifelse(pct >= 5, pct, "")), size = 2.8) +
  scale_fill_viridis_c() + theme_minimal(base_size = 12) +
  labs(title = "InSituType vs Leiden (12b)", subtitle = paste0("Each row = 100% of an InSituType cluster; ARI = ", round(ari, 3)),
       x = "Leiden cluster", y = "InSituType cluster", fill = "% of\nInSituType\ncluster")
save_tiff(p_ct, paste0("05_InSituType_vs_Leiden_", file_tag, ".tiff"), 12, 10)
fwrite(dcast(ct, insitutype_cluster ~ leiden, value.var = "N", fill = 0),
       file.path(dir_out, paste0("InSituType_vs_Leiden_cells_", file_tag, ".csv")))

## 7. cluster profiles heatmap (FAQ code): each gene scaled to its max across clusters
mat <- prof                                                                 # genes x clusters (renumbered)
mat <- sweep(mat, 1, pmax(apply(mat, 1, max), 0.1), "/")                    # 0-1 per gene; floor 0.1 for low genes
pheatmap::pheatmap(mat, col = colorRampPalette(c("white", "darkblue"))(100),
                   fontsize_row = 5,                                        # as FAQ (names overlap with 6K genes)
                   filename = file.path(dir_out, paste0("07_Profiles_heatmap_", file_tag, ".tiff")),
                   width = 6, height = 20)                                  # inches, as FAQ


## =============================================================================
## CHUNK 6: genes enriched in each InSituType cluster vs all other cells
## =============================================================================
## FindAllMarkers reads the "data" layer -> log-normalise first (as script 14);
## this does not change the InSituType clusters (they came from raw counts).
## Markers use all cells (insitutype_cluster), not the "unclassified" column.
obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
Idents(obj) <- "insitutype_cluster"
markers <- as.data.table(FindAllMarkers(obj, assay = "RNA", only.pos = TRUE, test.use = "wilcox",
                                        features = colnames(x)))            # genes only (control probes excluded)
saveRDS(markers, file.path(dir_out, paste0("FindAllMarkers_", file_tag, "_all_results.RDS")))

sig <- markers[p_val_adj < marker_padj_max & avg_log2FC > 0]
if (rank_by == "pct.1") setorder(sig, cluster, -pct.1, -avg_log2FC) else setorder(sig, cluster, -avg_log2FC, -pct.1)
sig[, rank_in_cluster := seq_len(.N), by = cluster]
setcolorder(sig, c("cluster", "rank_in_cluster", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"))
fwrite(sig, file.path(dir_out, paste0("FindAllMarkers_", file_tag, "_all_significant_ranked.csv")))

cl_levels <- levels(cl_num)
sheets <- lapply(cl_levels, function(k) {                                   # one sheet per cluster, top 20
  top <- sig[cluster == k & rank_in_cluster <= marker_top_n]
  if (nrow(top) == 0) data.frame(note = "no significant positive markers") else as.data.frame(top)
})
names(sheets) <- paste0("cluster_", cl_levels)
write.xlsx(sheets, file.path(dir_out, paste0("FindAllMarkers_", file_tag, "_top", marker_top_n, "_per_cluster.xlsx")),
           firstRow = TRUE, colWidths = "auto", overwrite = TRUE)

smry <- data.table(cluster = cl_levels,
                   insitutype_letter = map$insitutype_letter,
                   n_cells = map$n_cells,
                   median_prob = sapply(cl_levels, function(k) round(median(cells[insitutype_cluster == k, insitutype_prob]), 3)),
                   pct_unclassified = sapply(cl_levels, function(k) round(100 * mean(!cells[insitutype_cluster == k, insitutype_confident]), 1)),  # % cells < conf_min
                   n_significant_markers = sapply(cl_levels, function(k) nrow(sig[cluster == k])),
                   top10 = sapply(cl_levels, function(k) paste(head(sig[cluster == k, gene], 10), collapse = ", ")))
print(smry)
fwrite(smry, file.path(dir_out, paste0("FindAllMarkers_", file_tag, "_summary_per_cluster.csv")))

## dot plot: top 3 markers per cluster
top3 <- unique(sig[rank_in_cluster <= 3, gene])
p_dot <- DotPlot(obj, features = top3, group.by = "insitutype_cluster") + RotatedAxis() +
  labs(title = "Top 3 enriched genes per InSituType cluster")
save_tiff(p_dot, paste0("06_DotPlot_top3_markers_", file_tag, ".tiff"), max(12, length(top3) * 0.3), 8)


## =============================================================================
## CHUNK 7: 1b. UMAP built from InSituType's own output (per-cell log-likelihoods)
## =============================================================================
## Layout only: clusters separate by construction, so NOT evidence that clusters are real
## (use plot 01 and the flightpath, plot 03, for that). Clusters are NOT changed here.
umap_ist_nn        <- 30     # n_neighbors: lower = crisper, more separate islands   <- change here
umap_ist_min_dist  <- 0.3    # higher = looser points within each cluster           <- change here
umap_ist_spread    <- 2      # overall scale of the layout (as 12b)                 <- change here
umap_ist_repulsion <- 5      # higher = clusters pushed further apart (as 12b)      <- change here
umap_ist_tag <- paste0("NN", umap_ist_nn, "_md", umap_ist_min_dist, "_sp", umap_ist_spread, "_rep", umap_ist_repulsion)

ll <- res$logliks[rownames(x), map$insitutype_letter, drop = FALSE]        # cells x clusters, object order
ll <- ll - apply(ll, 1, max)                                               # relative to each cell's best cluster (best = 0)
ll[ll < -50] <- -50                                                        # cap very unlikely clusters
set.seed(seed)
um_ist <- uwot::umap(ll, n_neighbors = umap_ist_nn,                        # input = cluster scores per cell, not PCs
                     min_dist = umap_ist_min_dist, spread = umap_ist_spread,
                     repulsion_strength = umap_ist_repulsion,
                     metric = "euclidean")
rownames(um_ist) <- rownames(x); colnames(um_ist) <- c("UMAPIST_1", "UMAPIST_2")
obj[["umap_insitutype"]] <- CreateDimReducObject(embeddings = um_ist, key = "UMAPIST_", assay = "RNA")

p_umap_ist <- DimPlot(obj, reduction = "umap_insitutype", group.by = "insitutype_cluster", cols = cl_cols,
                      label = TRUE, label.box = TRUE, shuffle = TRUE, seed = seed,
                      pt.size = 0.05, alpha = 0.3, raster = FALSE) +
  coord_fixed() + theme_bw(base_size = 13) +
  guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
  labs(title = paste0("InSituType unsupervised: ", nlevels(cl_num), " clusters"),
       subtitle = paste0("UMAP of InSituType log-likelihoods (layout only): n_neighbors ", umap_ist_nn,
                         ", min_dist ", umap_ist_min_dist, ", spread ", umap_ist_spread,
                         ", repulsion ", umap_ist_repulsion),
       colour = "InSituType\ncluster")
save_tiff(p_umap_ist, paste0("01b_UMAP_of_InSituType_logliks_", umap_ist_tag, "_", file_tag, ".tiff"), 11, 9)


## =============================================================================
## CHUNK 8: parameters + session info (methods section)
## =============================================================================
fwrite(data.table(parameter = c("input_object", "qc_tag", "method", "control_probes_removed", "n_genes_used",
                                "cohort_features", "n_clusts_range",
                                "n_clusts_chosen", "n_starts", "neg", "seed", "conf_min",
                                "n_unclassified", "refined", "refine_merges", "refine_delete",
                                "refine_subcluster", "markers", "rank_by", "padj_cutoff"),
                  value = c(seurat_out, qc_tag,
                            "InSituType::insitutype (unsupervised, raw counts, negprobe background, fastCohorting cohorts)",
                            if (length(ctrl)) paste(ctrl, collapse = "; ") else "none", ncol(x),
                            paste(if_cols, collapse = "; "), paste(range(n_clusts), collapse = "-"),
                            nlevels(cl_num), n_starts, "mean negprobe counts per cell", seed, conf_min,
                            sum(!obj$insitutype_confident), run_refine,
                            if (is.null(refine_merges)) "none" else paste(names(refine_merges), refine_merges, sep = "->", collapse = "; "),
                            if (is.null(refine_delete)) "none" else paste(refine_delete, collapse = "; "),
                            if (is.null(refine_subcluster)) "none" else paste(names(refine_subcluster), unlist(refine_subcluster), sep = "x", collapse = "; "),
                            "Seurat::FindAllMarkers(only.pos = TRUE, wilcox) on log-normalised RNA",
                            rank_by, marker_padj_max)),
       file.path(dir_out, paste0("parameters_", file_tag, ".csv")))
fwrite(data.table(parameter = c("UMAP_insitutype_input", "n_neighbors", "min_dist", "spread", "repulsion_strength", "seed"),
                  value = c("uwot::umap on per-cell InSituType log-likelihoods relative to best cluster, capped at -50",
                            umap_ist_nn, umap_ist_min_dist, umap_ist_spread, umap_ist_repulsion, seed)),
       file.path(dir_out, paste0("parameters_UMAP_insitutype_", umap_ist_tag, "_", file_tag, ".csv")))
writeLines(capture.output(sessionInfo()), file.path(dir_out, paste0("sessionInfo_", file_tag, ".txt")))
message("Done -> ", dir_out)













































## =============================================================================
## 20_InSituType_cluster_annotation.R
## =============================================================================
## WHAT IT DOES
##   Repeats scripts 15-19 (naming clusters) on the InSituType clusters instead of
##   the Leiden clusters. Input = the marker files written by the InSituType script
##   (its CHUNK 6). Methods, references and settings are the same as scripts 15-19:
##     CHUNK 4  (15)    CellTypist marker overlap, hypergeometric  -> celltypist_top20/top50_<tag>.xlsx
##     CHUNK 5  (16)    CellTypist per cell + majority vote        -> celltypist_majority_voting_<tag>.xlsx
##     CHUNK 6  (17)    CellMarker 2.0 enricher + counting rule    -> cellmarker_top20/top50/top10_min5_<tag>.xlsx
##     CHUNK 7  (18A/B) CELLxGENE Census enricher + counting rule  -> CELLxGENE_top20/top50/top10_min5_<tag>.xlsx
##     CHUNK 8  (19A/B) HCA aorta enricher + counting rule         -> HCA_DataExplorer_*_<level>_<tag>.xlsx
##   Reference files are NOT rebuilt: they are read from the Results folder where
##   scripts 17-19 saved them (same 6K panel genes). Everything is saved in the
##   InSituType folder; <tag> = file_tag, so the Leiden results are never overwritten.
## Packages already installed. use_condaenv() (CHUNK 4) must be the first Python
## call of the R session -> if Python was already used for something else, restart R.
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 1: settings
## =============================================================================
library(data.table)       # tables
library(Matrix)           # sparse matrices
library(Seurat)           # NormalizeData, LayerData (CHUNK 5)
library(openxlsx)         # read the top-20 sheets + write results
library(readxl)           # read the CellMarker xlsx
library(clusterProfiler)  # enricher()
library(reticulate)       # R <-> Python (CellTypist, CELLxGENE)

## InSituType folder (dir_out of the InSituType script) and its file_tag
ist_dir  <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Unsupervised_clustering_folder/QC03_04_06_07_08_09/InSituType_PC20_NN250_rep3"   # <- change here
##ist_dir  <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Unsupervised_clustering_folder/QC03_04_05_06_07_08_09/InSituType_PC20_NN250_rep3"

file_tag <- "insitutype_PC20_NN250_rep3"   # add "_refined" if refineClusters was used   <- change here

## folder where scripts 17-19 saved the reference files (Cell_marker_Human.xlsx, CELLxGENE + HCA CSVs)
ref_dir <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Results/QC03_04_06_07_08_09_final_All_markers"   # <- change here
py_file   <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/18_CELLxGENE_reference_markers.py"   # only used if the CELLxGENE CSV is missing
conda_env <- "celltypist_env"

## CellTypist (as scripts 15, 16)
models_use <- c(Vascular = "Adult_Human_Vascular.pkl",
                Immune   = "Immune_All_Low.pkl",
                Heart    = "Healthy_Adult_Heart.pkl")
ref_n <- 100          # genes per CellTypist cell type used for matching (script 15)

## CellMarker (as script 17)
tissue_regex  <- "aort|arter|vessel|vascul|blood"
exclude_regex <- "neur|glia|schwann|astrocyt|oligodendro|microglia"

## CELLxGENE reference settings (as 18A/18B; must match the saved reference)
tissue_general <- "vasculature"
normal_only    <- TRUE
max_cells      <- 3000
min_cells      <- 50
ref_top_n      <- 100

## HCA aorta label levels (as 19A/19B)
hca_levels <- c("celltype", "celltype2")   # authors' major types; types/states

## tests (as scripts 15-19)
min_set  <- 3      # a reference cell type is tested only if >= 3 of its markers are in the background
padj_max <- 0.05   # a call is made only if BH-adjusted p < this
top_k    <- 10     # counting rule: top genes per cluster
min_hits <- 5      # counting rule: how many must be markers of the cell type

## input files written by the InSituType script
xlsx_top20 <- file.path(ist_dir, paste0("FindAllMarkers_", file_tag, "_top20_per_cluster.xlsx"))
csv_ranked <- file.path(ist_dir, paste0("FindAllMarkers_", file_tag, "_all_significant_ranked.csv"))
rds_all    <- file.path(ist_dir, paste0("FindAllMarkers_", file_tag, "_all_results.RDS"))
cells_csv  <- file.path(ist_dir, paste0("InSituType_cells_", file_tag, ".csv"))
params_csv <- file.path(ist_dir, paste0("parameters_", file_tag, ".csv"))
stopifnot("InSituType marker / cell / parameter files not found - check ist_dir and file_tag" =
            all(file.exists(c(xlsx_top20, csv_ranked, rds_all, cells_csv, params_csv))))
out_dir <- ist_dir                                                          # results saved in the InSituType folder

## how the InSituType markers were ranked (rank_by in the InSituType script)
ist_par    <- fread(params_csv)
rank_by    <- ist_par[parameter == "rank_by", value]
seurat_out <- ist_par[parameter == "input_object", value]                   # 12b Seurat object (CHUNK 5)
rank_label <- if (rank_by == "pct.1") "ranked by pct.1, then avg_log2FC" else "ranked by avg_log2FC, then pct.1"
label20 <- paste0("top 20 per InSituType cluster (xlsx sheets; ", rank_label, ")")
label50 <- paste0("top 50 per InSituType cluster (ranked CSV, rank_in_cluster <= 50; ", rank_label, ")")
message("InSituType markers ", rank_label, " | results -> ", out_dir)


## =============================================================================
## CHUNK 2: InSituType marker lists (top 20 from xlsx, top 50 from CSV)
## =============================================================================
## CosMx combined probes ("IGHG1/2") -> separate genes ("IGHG1", "IGHG2")
split_probe <- function(g) {
  parts <- strsplit(g, "/", fixed = TRUE)[[1]]                              # "IGHG1/2" -> "IGHG1", "2"
  if (length(parts) == 1) return(g)                                         # ordinary gene: unchanged
  stem <- sub("[0-9]+$", "", parts[1])                                      # "IGHG1" -> "IGHG"
  c(parts[1], paste0(stem, parts[-1]))                                      # "IGHG1", "IGHG2"
}
split_genes <- function(gs) unique(unlist(lapply(gs, split_probe)))

## genes FindAllMarkers tested on the InSituType clusters (control probes already excluded)
panel <- split_genes(unique(as.character(readRDS(rds_all)$gene)))

## top 20: one sheet per cluster ("cluster_0", "cluster_1", ...)
sheets <- getSheetNames(xlsx_top20)
top20 <- lapply(sheets, function(s) {
  d <- read.xlsx(xlsx_top20, sheet = s)
  if ("gene" %in% names(d)) d$gene else character(0)                        # sheet may say "no significant markers"
})
names(top20) <- sub("^cluster_", "", sheets)                                # "0", "1", ...
clusters <- names(top20)[order(as.numeric(names(top20)))]                   # 0, 1, 2, ... (not 0, 1, 10)

## top 50: ranks 1-50 of each cluster in the ranked CSV
ranked <- fread(csv_ranked)
ranked[, cluster := as.character(cluster)]
top50 <- lapply(names(top20), function(k) ranked[cluster == k & rank_in_cluster <= 50, gene])
names(top50) <- names(top20)

## check: xlsx top 20 = CSV ranks 1-20 (same InSituType run)
for (k in names(top20)) {
  stopifnot("xlsx and CSV disagree - are they from the same run?" =
              setequal(top20[[k]], ranked[cluster == k & rank_in_cluster <= 20, gene]))
}
message("Check passed: ", length(top20), " InSituType clusters; xlsx top 20 = CSV ranks 1-20")
print(data.table(cluster = names(top20), n_top20 = lengths(top20), n_top50 = lengths(top50)))


## =============================================================================
## CHUNK 3: shared functions (identical to scripts 17-19)
## =============================================================================
## --- enricher (17 CHUNK 4, 18A, 19A): hypergeometric, BH per cluster ---
run_cluster <- function(genes, t2g, universe) {
  g <- intersect(split_genes(genes), universe)                              # markers in the background
  if (length(g) == 0) return(NULL)
  e <- enricher(g, TERM2GENE = t2g, universe = universe, minGSSize = min_set, maxGSSize = 5000,
                pvalueCutoff = 1, qvalueCutoff = 1)                         # keep all; call uses padj_max
  if (is.null(e) || nrow(as.data.frame(e)) == 0) return(NULL)
  r <- as.data.table(as.data.frame(e))
  r[, n_markers_in_background := length(g)]
  r
}

enricher_and_save <- function(marker_lists, label, out_file, t2g, universe,
                              ref_text, bg_text, extra_sheets) {
  all <- rbindlist(lapply(names(marker_lists), function(k) {
    r <- run_cluster(marker_lists[[k]], t2g, universe)
    if (is.null(r)) return(NULL)
    r[, cluster := k]; setcolorder(r, "cluster"); r                         # cluster id as first column
  }), fill = TRUE)
  setorder(all, cluster, pvalue)                                            # best p first within each cluster
  
  smry <- rbindlist(lapply(clusters, function(k) {                          # best match + next two per cluster
    r <- all[cluster == k]
    if (nrow(r) == 0) return(data.table(cluster = k, call = "nothing tested"))
    data.table(cluster = k, n_markers_in_background = r$n_markers_in_background[1],
               call = if (r$p.adjust[1] < padj_max) r$ID[1] else "no significant match",
               best_match = r$ID[1], n_shared = r$Count[1],
               shared_genes = gsub("/", ", ", r$geneID[1]), p_adj = signif(r$p.adjust[1], 3),
               second = if (nrow(r) > 1) paste0(r$ID[2], " (", r$Count[2], " shared, padj ", signif(r$p.adjust[2], 2), ")") else NA,
               third  = if (nrow(r) > 2) paste0(r$ID[3], " (", r$Count[3], " shared, padj ", signif(r$p.adjust[3], 2), ")") else NA)
  }), fill = TRUE)
  print(smry[, .(cluster, call, n_shared, p_adj, shared_genes)])
  
  keep <- intersect(c("cluster", "ID", "GeneRatio", "BgRatio", "Count", "pvalue", "p.adjust", "geneID",
                      "n_markers_in_background"), names(all))               # columns differ between versions
  how_to_read <- data.frame(item = c("Clusters", "Markers used", "Reference", "Method", "call",
                                     "n_shared / shared_genes", "p_adj"),
                            meaning = c(paste0("InSituType unsupervised clusters (", file_tag, ")"),
                                        label, ref_text,
                                        paste0("clusterProfiler::enricher (hypergeometric); background = ", bg_text,
                                               "; reference types with < ", min_set, " markers in the background not tested"),
                                        paste0("best_match if p_adj < ", padj_max, ", otherwise 'no significant match'"),
                                        "cluster markers that are reference markers of that cell type",
                                        "BH-adjusted p-value across the reference cell types tested for that cluster"))
  write.xlsx(c(list(Summary = as.data.frame(smry), All_tests = as.data.frame(all[, ..keep])),
               extra_sheets, list(How_to_read = how_to_read)),
             out_file, firstRow = TRUE, colWidths = "auto", overwrite = TRUE)
  message("Saved: ", out_file)
}

## --- counting rule (17 CHUNK 6, 18B, 19B): >= min_hits of the top_k genes are markers of a type ---
count_and_save <- function(sets, N, ref_text, out_file, extra_sheets) {
  counts <- rbindlist(lapply(clusters, function(k) {
    g <- ranked[cluster == k & rank_in_cluster <= top_k][order(rank_in_cluster), gene]   # top_k genes, rank order
    rbindlist(lapply(names(sets), function(ct) {
      s   <- sets[[ct]]
      hit <- g[vapply(g, function(x) any(split_probe(x) %in% s), logical(1))]   # IGHG1/2 counts once
      data.table(cluster = k, cell_type = ct, n_top_genes = length(g), n_hits = length(hit),
                 hit_genes = paste(hit, collapse = ", "), markers_of_type = length(s),
                 expected_by_chance = round(length(g) * length(s) / N, 2))   # random genes would hit this many
    }))
  }))
  setorder(counts, cluster, -n_hits, expected_by_chance)
  passed <- counts[n_hits >= min_hits]                                      # every type passing the rule
  
  calls <- rbindlist(lapply(clusters, function(k) {
    p <- passed[cluster == k]
    b <- counts[cluster == k][1]                                            # best type even if it failed
    data.table(cluster = k,
               top_genes = paste(ranked[cluster == k & rank_in_cluster <= top_k][order(rank_in_cluster), gene], collapse = ", "),
               n_cell_types_called = nrow(p),
               cell_types_called = if (nrow(p) > 0) paste0(p$cell_type, " (", p$n_hits, "/", p$n_top_genes, ")", collapse = "; ")
               else paste0("none (< ", min_hits, " hits)"),
               best_cell_type = b$cell_type, best_n_hits = b$n_hits, best_hit_genes = b$hit_genes)
  }))
  print(calls[, .(cluster, n_cell_types_called, cell_types_called)])
  
  how_to_read <- data.frame(item = c("Clusters", "Rule", "Reference", "cell_types_called", "best_cell_type",
                                     "expected_by_chance", "Caution"),
                            meaning = c(paste0("InSituType unsupervised clusters (", file_tag, "); top genes ", rank_label),
                                        paste0("a cluster is named a cell type if >= ", min_hits, " of its top ", top_k,
                                               " genes are reference markers of that type; several types allowed"),
                                        ref_text,
                                        "every cell type passing the rule, hits/top genes in brackets",
                                        "cell type with the most hits, shown even when it does not pass (near-miss)",
                                        "hits expected from the same number of random panel genes",
                                        "No statistical test: trust a call only if n_hits is well above expected_by_chance and the hit genes are canonical for that type."))
  write.xlsx(c(list(Calls = as.data.frame(calls),
                    Passing_cell_types = as.data.frame(passed[order(as.numeric(cluster), -n_hits)]),
                    All_counts = as.data.frame(counts[n_hits > 0][order(as.numeric(cluster), -n_hits)])),
               extra_sheets, list(How_to_read = how_to_read)),
             out_file, firstRow = TRUE, colWidths = "auto", overwrite = TRUE)
  message("Saved: ", out_file)
}


## =============================================================================
## CHUNK 4: (script 15) CellTypist marker overlap - top 20 and top 50
## =============================================================================
use_condaenv(conda_env, required = TRUE)                                    # first Python call of the session
ct_models <- import("celltypist.models", convert = FALSE)
ct_models$download_models(force_update = FALSE,                             # FALSE: same model versions as the Leiden run
                          model = r_to_py(as.list(unname(models_use))))
models <- lapply(models_use, function(f) ct_models$Model$load(model = f))
model_types <- lapply(models, function(m) unlist(py_to_r(m$cell_types$tolist())))   # cell-type names
model_genes <- lapply(models, function(m) unlist(py_to_r(m$features$tolist())))     # genes each model uses

## universe per model = model genes tested on your panel
ct_universe <- lapply(model_genes, function(g) intersect(g, panel))
for (nm in names(ct_universe)) message(nm, ": ", length(ct_universe[[nm]]), " model genes are on your panel")

## reference genes of each cell type = its top ref_n positive driving genes in the universe
ref_genes <- lapply(names(models), function(nm) {
  m <- models[[nm]]
  all_n <- as.integer(length(model_genes[[nm]]))                            # ask for the full ranking
  setNames(lapply(model_types[[nm]], function(ty) {
    ranked_ty <- unlist(py_to_r(m$extract_top_markers(ty, top_n = all_n, only_positive = TRUE)$tolist()))
    head(ranked_ty[ranked_ty %in% ct_universe[[nm]]], ref_n)                # top ref_n that are on the panel
  }), model_types[[nm]])
})
names(ref_genes) <- names(models)

## every cluster x cell-type test for one model (hypergeometric, BH across cell types)
match_clusters <- function(marker_lists, nm) {
  U <- ct_universe[[nm]]; N <- length(U)
  rbindlist(lapply(names(marker_lists), function(k) {
    q <- intersect(split_genes(marker_lists[[k]]), U)                       # this cluster's usable markers
    res <- rbindlist(lapply(names(ref_genes[[nm]]), function(ty) {
      ref <- ref_genes[[nm]][[ty]]
      hit <- intersect(q, ref)
      data.table(cluster = k, cell_type = ty,
                 n_markers = length(marker_lists[[k]]),
                 n_markers_in_model = length(q),
                 n_shared = length(hit),
                 p_value = if (length(q) > 0 && length(ref) > 0)
                   phyper(length(hit) - 1, length(ref), N - length(ref), length(q), lower.tail = FALSE) else NA_real_,
                 shared_genes = paste(hit, collapse = ", "))
    }))
    res[, p_adj := p.adjust(p_value, method = "BH")]
    res[order(p_value, -n_shared, na.last = TRUE)]
  }))
}

## one row per cluster: best match, its genes, and the next two candidates
best_per_cluster <- function(all) {
  out <- all[, .(n_markers = n_markers[1], n_markers_in_model = n_markers_in_model[1],
                 best_match = cell_type[1], n_shared = n_shared[1], shared_genes = shared_genes[1],
                 p_adj = signif(p_adj[1], 3),
                 call = if (!is.na(p_adj[1]) && p_adj[1] < padj_max) cell_type[1] else "no significant match",
                 second = paste0(cell_type[2], " (", n_shared[2], " shared, padj ", signif(p_adj[2], 2), ")"),
                 third  = paste0(cell_type[3], " (", n_shared[3], " shared, padj ", signif(p_adj[3], 2), ")")),
             by = cluster]
  out[order(as.numeric(cluster))]
}

celltypist_overlap_save <- function(marker_lists, label, out_file) {
  per_model <- lapply(names(models), function(nm) best_per_cluster(match_clusters(marker_lists, nm)))
  names(per_model) <- names(models)
  smry <- data.table(cluster = per_model[[1]]$cluster, n_markers = per_model[[1]]$n_markers)
  for (nm in names(per_model)) {                                            # calls of each model side by side
    p <- per_model[[nm]]
    smry[, (paste0(nm, "_call"))         := p$call]
    smry[, (paste0(nm, "_shared_genes")) := ifelse(p$call == "no significant match", "", p$shared_genes)]
    smry[, (paste0(nm, "_padj"))         := p$p_adj]
  }
  print(smry[, c("cluster", paste0(names(models), "_call")), with = FALSE])
  how_to_read <- data.frame(item = c("Clusters", "Markers used", "Models", "Method", "call", "best_match",
                                     "n_markers_in_model", "n_shared / shared_genes", "p_adj", "no significant match"),
                            meaning = c(paste0("InSituType unsupervised clusters (", file_tag, ")"), label,
                                        paste(paste0(names(models_use), " = ", models_use), collapse = "; "),
                                        paste0("Each cluster's markers vs each CellTypist cell type's top ", ref_n,
                                               " driving genes (model.extract_top_markers) on the panel; hypergeometric test, BH across cell types"),
                                        paste0("best_match if p_adj < ", padj_max, ", otherwise 'no significant match'"),
                                        "cell type sharing the most markers beyond chance (reported even when not significant)",
                                        "how many of the cluster's markers the model knows (only these can match)",
                                        "markers shared with the best-matching cell type",
                                        "BH-adjusted p-value of the overlap",
                                        "markers do not resemble any cell type in that model (or too few markers are known to it)"))
  write.xlsx(c(list(Summary = as.data.frame(smry)), lapply(per_model, as.data.frame),
               list(How_to_read = how_to_read)),
             out_file, firstRow = TRUE, colWidths = "auto", overwrite = TRUE)
  message("Saved: ", out_file)
}

celltypist_overlap_save(top20, label20, file.path(out_dir, paste0("celltypist_top20_", file_tag, ".xlsx")))
celltypist_overlap_save(top50, label50, file.path(out_dir, paste0("celltypist_top50_", file_tag, ".xlsx")))


## =============================================================================
## CHUNK 5: (script 16) CellTypist per cell, majority vote within InSituType clusters
## =============================================================================
obj <- readRDS(seurat_out)                                                  # 12b object (same cells as InSituType)
cells_ist <- fread(cells_csv)                                               # InSituType cluster of every cell
stopifnot("cells differ between object and InSituType table" = setequal(cells_ist$cell_id, colnames(obj)))
cl_cell <- as.character(cells_ist$insitutype_cluster[match(colnames(obj), cells_ist$cell_id)])   # object order
stopifnot(!anyNA(cl_cell))

## expression as CellTypist needs it (vignette: log1p, 10,000 per cell)
obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
X <- t(LayerData(obj, assay = "RNA", layer = "data"))                       # cells x genes (AnnData orientation)
chk <- Matrix::rowSums(expm1(X[seq_len(min(1000, nrow(X))), ]))             # first 1,000 cells
stopifnot("expression is not log1p-normalised to 10,000" = all(abs(chk - 10000) < 1))

## AnnData with the InSituType clusters
celltypist_py <- import("celltypist", convert = FALSE)                      # celltypist.annotate
anndata <- import("anndata", convert = FALSE)
pd      <- import("pandas", convert = FALSE)
obs <- pd$DataFrame(data = r_to_py(list(insitutype = cl_cell)), index = r_to_py(colnames(obj)))
var <- pd$DataFrame(index = r_to_py(rownames(obj)))
ad  <- anndata$AnnData(X = r_to_py(X)$tocsr(), obs = obs, var = var)
rm(X); invisible(gc())
message("AnnData: ", py_to_r(ad$n_obs), " cells x ", py_to_r(ad$n_vars), " genes")

per_model_mv <- list()
for (nm in names(models_use)) {
  f <- models_use[[nm]]
  n_genes <- sum(model_genes[[nm]] %in% rownames(obj))                      # model genes on the panel
  message("\n", nm, " (", f, "): ", n_genes, " of its genes are on your panel. Annotating...")
  pred <- celltypist_py$annotate(ad, model = f, majority_voting = TRUE,
                                 over_clustering = "insitutype")            # votes within InSituType clusters
  pl   <- py_to_r(pred$predicted_labels)
  conf <- as.numeric(py_to_r(pred$probability_matrix$max(axis = 1L)$to_numpy()))   # conf_score
  stopifnot("cell order changed" = identical(rownames(pl), colnames(obj)))
  
  cells <- data.table(cell_id = colnames(obj), cluster = cl_cell,
                      predicted = as.character(pl$predicted_labels),        # label of the cell on its own
                      majority_vote = as.character(pl$majority_voting),     # label of its InSituType cluster
                      conf_score = round(conf, 4))
  fwrite(cells, file.path(out_dir, paste0("celltypist_per_cell_", nm, "_", file_tag, ".csv.gz")))
  
  tab <- cells[, {
    fr <- sort(table(predicted), decreasing = TRUE) / .N                    # fraction of cells per label
    mv <- majority_vote[1]
    .(n_cells = .N,
      majority_label = mv,
      pct_cells_agreeing = round(100 * mean(predicted == mv), 1),
      median_conf_agreeing = round(median(conf_score[predicted == mv]), 3),
      top3_labels = paste0(names(fr)[1:min(3, length(fr))], " (",
                           round(100 * fr[1:min(3, length(fr))], 1), "%)", collapse = "; "))
  }, by = cluster][order(as.numeric(cluster))]
  tab[, genes_used := n_genes]
  print(tab[, .(cluster, n_cells, majority_label, pct_cells_agreeing, median_conf_agreeing)])
  per_model_mv[[nm]] <- tab
}

smry_mv <- per_model_mv[[1]][, .(cluster, n_cells)]
for (nm in names(per_model_mv)) {                                           # all models side by side
  p <- per_model_mv[[nm]]
  smry_mv[, (paste0(nm, "_label"))     := p$majority_label]
  smry_mv[, (paste0(nm, "_pct_agree")) := p$pct_cells_agreeing]
  smry_mv[, (paste0(nm, "_conf"))      := p$median_conf_agreeing]
}
print(smry_mv)
how_to_read_mv <- data.frame(item = c("Clusters", "Method", "label", "pct_agree", "conf", "top3_labels", "genes_used", "Caution"),
                             meaning = c(paste0("InSituType unsupervised clusters (", file_tag, ")"),
                                         "celltypist.annotate on log1p-normalised expression of every cell; majority_voting = TRUE with over_clustering = InSituType clusters",
                                         "majority label of the cluster (the label most of its cells received)",
                                         "% of the cluster's cells that were individually given that label (low = mixed or uncertain cluster)",
                                         "median CellTypist confidence (0-1) of those agreeing cells",
                                         "the three most frequent per-cell labels in the cluster",
                                         "model genes present on the CosMx panel (the model was trained on whole-transcriptome data)",
                                         "CellTypist always returns a label, even for cells unlike any reference type. Check labels with canonical markers; ignore cardiomyocyte/neural labels from the Heart model."))
out_mv <- file.path(out_dir, paste0("celltypist_majority_voting_", file_tag, ".xlsx"))
write.xlsx(c(list(Summary = as.data.frame(smry_mv)), lapply(per_model_mv, as.data.frame),
             list(How_to_read = how_to_read_mv)),
           out_mv, firstRow = TRUE, colWidths = "auto", overwrite = TRUE)
message("Saved: ", out_mv)
rm(obj, ad); invisible(gc())                                                # free memory


## =============================================================================
## CHUNK 6: (script 17) CellMarker 2.0 - enricher top 20/50 + counting rule
## =============================================================================
cm_file <- file.path(ref_dir, "Cell_marker_Human.xlsx")                     # downloaded by script 17
stopifnot("Cell_marker_Human.xlsx not in ref_dir - run script 17 first" = file.exists(cm_file))
cm <- as.data.table(read_excel(cm_file))
gene_col <- if ("Symbol" %in% names(cm)) "Symbol" else "marker"             # official gene symbol if present
if ("cell_type" %in% names(cm)) cm <- cm[cell_type == "Normal cell"]        # normal cells only
cm_vasc <- cm[grepl(tissue_regex, tissue_class, ignore.case = TRUE) |
                grepl(tissue_regex, tissue_type,  ignore.case = TRUE)]        # vascular + blood entries
cm_vasc <- cm_vasc[!grepl(exclude_regex, cell_name, ignore.case = TRUE)]    # drop neural cell types
t2g_cm <- unique(cm_vasc[!is.na(get(gene_col)) & get(gene_col) != "",
                         .(term = cell_name, gene = get(gene_col))])          # TERM2GENE
universe_cm <- panel                                                        # background = panel genes tested
message(length(unique(t2g_cm$term)), " CellMarker cell types; background ", length(universe_cm), " genes")

cm_text <- paste0("CellMarker 2.0 human, normal cells, tissue matching '", tissue_regex,
                  "', neural types excluded ('", exclude_regex, "')")
sets_used <- t2g_cm[gene %in% universe_cm, .(n_markers_on_panel = .N,
                                             markers_on_panel = paste(sort(gene), collapse = ", ")), by = term][order(term)]
enricher_and_save(top20, label20, file.path(out_dir, paste0("cellmarker_top20_", file_tag, ".xlsx")),
                  t2g_cm, universe_cm, cm_text, "panel genes tested by FindAllMarkers",
                  list(CellMarker_sets_used = as.data.frame(sets_used)))
enricher_and_save(top50, label50, file.path(out_dir, paste0("cellmarker_top50_", file_tag, ".xlsx")),
                  t2g_cm, universe_cm, cm_text, "panel genes tested by FindAllMarkers",
                  list(CellMarker_sets_used = as.data.frame(sets_used)))

sets_cm <- split(t2g_cm[gene %in% universe_cm, gene], t2g_cm[gene %in% universe_cm, term])   # markers per type
count_and_save(sets_cm, length(universe_cm), cm_text,
               file.path(out_dir, paste0("cellmarker_top", top_k, "_min", min_hits, "_", file_tag, ".xlsx")),
               list())


## =============================================================================
## CHUNK 7: (scripts 18A/18B) CELLxGENE Census - enricher top 20/50 + counting rule
## =============================================================================
cxg_csv <- file.path(ref_dir, paste0("CELLxGENE_reference_markers_", tissue_general, ".csv"))   # built by 18A
if (!file.exists(cxg_csv)) {                                                # not built yet: build into the InSituType folder
  cxg_csv <- file.path(out_dir, paste0("CELLxGENE_reference_markers_", tissue_general, ".csv"))
  if (!file.exists(cxg_csv)) {
    stopifnot("Python file not found - set py_file in CHUNK 1" = file.exists(py_file))
    panel_txt <- file.path(out_dir, "panel_genes_for_CELLxGENE.txt")
    writeLines(panel, panel_txt)
    source_python(py_file)                                                  # makes build_reference() available
    build_reference(panel_txt, cxg_csv, tissue_general = tissue_general, normal_only = normal_only,
                    max_cells_per_type = as.integer(max_cells), min_cells_per_type = as.integer(min_cells),
                    top_n = as.integer(ref_top_n))
  }
}
ref_cxg <- fread(cxg_csv)                                                   # cell_type, rank, gene, ...
universe_cxg <- intersect(panel, readLines(sub("\\.csv$", "_genes_in_reference.txt", cxg_csv)))   # panel genes Census has
t2g_cxg <- unique(ref_cxg[gene %in% universe_cxg, .(term = cell_type, gene)])
message(length(unique(t2g_cxg$term)), " CELLxGENE cell types; background ", length(universe_cxg), " genes")

cxg_text <- paste0("CELLxGENE Census, human ", tissue_general, if (normal_only) ", normal donors" else "",
                   "; top ", ref_top_n, " Wilcoxon markers per cell type computed on your panel genes")
enricher_and_save(top20, label20, file.path(out_dir, paste0("CELLxGENE_top20_", file_tag, ".xlsx")),
                  t2g_cxg, universe_cxg, cxg_text, "panel genes present in Census",
                  list(Reference_markers = as.data.frame(ref_cxg)))
enricher_and_save(top50, label50, file.path(out_dir, paste0("CELLxGENE_top50_", file_tag, ".xlsx")),
                  t2g_cxg, universe_cxg, cxg_text, "panel genes present in Census",
                  list(Reference_markers = as.data.frame(ref_cxg)))

sets_cxg <- split(ref_cxg[gene %in% universe_cxg, gene], ref_cxg[gene %in% universe_cxg, cell_type])
count_and_save(sets_cxg, length(universe_cxg), cxg_text,
               file.path(out_dir, paste0("CELLxGENE_top", top_k, "_min", min_hits, "_", file_tag, ".xlsx")),
               list(Reference_markers = as.data.frame(ref_cxg)))


## =============================================================================
## CHUNK 8: (scripts 19A/19B) HCA aorta (Li et al. 2020) - both label levels
## =============================================================================
for (ref_level in hca_levels) {
  hca_csv <- file.path(ref_dir, paste0("HCA_reference_markers_aorta_", ref_level, ".csv"))   # built by script 19
  if (!file.exists(hca_csv)) { message("skipped ", ref_level, ": run 19_HCA_aorta_build_reference.R first"); next }
  ref_hca <- fread(hca_csv)
  universe_hca <- intersect(readLines(sub("\\.csv$", "_genes_in_reference.txt", hca_csv)), panel)   # panel genes in the reference
  t2g_hca <- unique(ref_hca[gene %in% universe_hca, .(term = cell_type, gene)])
  message("\nHCA ", ref_level, ": ", length(unique(t2g_hca$term)), " cell types; background ", length(universe_hca), " genes")
  
  hca_text <- paste0("HCA Data Explorer project AneurysmalHumanAorticTissue (Li et al. 2020 Circulation; GSE155468), ",
                     "ascending aorta, aneurysm + control; authors' labels (", ref_level,
                     "); top Wilcoxon markers per cell type on your panel genes")
  enricher_and_save(top20, label20, file.path(out_dir, paste0("HCA_DataExplorer_top20_", ref_level, "_", file_tag, ".xlsx")),
                    t2g_hca, universe_hca, hca_text, "panel genes measured in the reference",
                    list(Reference_markers = as.data.frame(ref_hca)))
  enricher_and_save(top50, label50, file.path(out_dir, paste0("HCA_DataExplorer_top50_", ref_level, "_", file_tag, ".xlsx")),
                    t2g_hca, universe_hca, hca_text, "panel genes measured in the reference",
                    list(Reference_markers = as.data.frame(ref_hca)))
  
  sets_hca <- split(ref_hca[gene %in% universe_hca, gene], ref_hca[gene %in% universe_hca, cell_type])
  count_and_save(sets_hca, length(universe_hca), hca_text,
                 file.path(out_dir, paste0("HCA_DataExplorer_top", top_k, "_min", min_hits, "_", ref_level, "_", file_tag, ".xlsx")),
                 list(Reference_markers = as.data.frame(ref_hca)))
}
message("Done -> ", out_dir)











## =============================================================================
## 14b_Heatmap_DotPlot_top_markers.R leiden clustering
## =============================================================================
## PURPOSE
##   Heatmap + dot plot of the top N significant markers per Leiden cluster
##   (main resolution), ranked by pct.1 (as script 14), all clusters in one plot.
##   Made for N = 10, 15 and 20. Seurat vignette method:
##     top N per cluster -> ScaleData on those genes -> DoHeatmap / DotPlot.
##   Uses script 14's output (_all_significant_ranked.csv) + the 12b object.
## OUTPUT: same folder as script 14: <output_dir>/Results/<QC>_final_All_markers/
##   Heatmap_top<N>_pct1_<tag>.tiff, DotPlot_top<N>_pct1_<tag>.tiff (300 dpi)
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 1: settings - same path you gave script 14
## =============================================================================
library(data.table)   # tables
library(Seurat)       # ScaleData, DoHeatmap, DotPlot
library(ggplot2)      # theme, ggsave

clustering_path <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Clustering_PCA_UMAP_Leiden_Markers/QC03_04_05_06_07_08_09/PC20_NN250_rep3/02_SENSITIVITY_batchawarePCA_Harmony/UMAP_by_Leiden_cluster.tiff"   # <- change here
##clustering_path <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Clustering_PCA_UMAP_Leiden_Markers/QC03_04_06_07_08_09/PC20_NN250_rep3/02_SENSITIVITY_batchawarePCA_Harmony/UMAP_by_Leiden_cluster.tiff"   # <- change here

leiden_resolution <- NULL          # NULL = main resolution of 12b (as script 14)
top_ns            <- c(10, 15, 20) # genes per cluster, one heatmap + one dot plot each   <- change here
cells_per_cluster <- 500           # heatmap only: cells drawn per cluster (all 123k cells cannot be shown)   <- change here
seed              <- 42


## =============================================================================
## CHUNK 2: paths from clustering_path (identical to script 14)
## =============================================================================
branch_dir <- if (dir.exists(clustering_path)) clustering_path else dirname(clustering_path)   # file -> its folder
branch_dir <- normalizePath(branch_dir, mustWork = TRUE)                    # stop if missing
dir_seu    <- dirname(branch_dir)                                           # .../PC20_NN250_rep3
run_tag    <- basename(dir_seu)                                             # "PC20_NN250_rep3"
qc_tag     <- basename(dirname(dir_seu))                                    # QC choice
output_dir <- dirname(dirname(dirname(dir_seu)))                            # .../AAA
prefix <- if (grepl("^01_PRIMARY", basename(branch_dir))) "primary" else
  if (grepl("^02_SENSITIVITY", basename(branch_dir))) "harmony" else
    stop("branch folder must start with 01_PRIMARY or 02_SENSITIVITY")
if (is.null(leiden_resolution)) {                                           # main resolution from 12b
  params <- fread(file.path(dir_seu, paste0("parameters_used_", run_tag, ".csv")))
  leiden_resolution <- as.numeric(params[parameter == "Leiden_resolution", value])
}
file_tag <- paste0(prefix, "_", run_tag, "_res", leiden_resolution)        # "harmony_PC20_NN250_rep3_res1"
dir_mk   <- file.path(output_dir, "Results", paste0(qc_tag, "_final_All_markers"))   # script 14's folder

seurat_out <- list.files(dir_seu, pattern = paste0("^seurat_.*_SeuratFunctions_", run_tag, "\\.RDS$"),
                         full.names = TRUE)
stopifnot("expected exactly one 12b Seurat object" = length(seurat_out) == 1)
message("QC: ", qc_tag, " | ", file_tag, " | plots -> ", dir_mk)


## =============================================================================
## CHUNK 3: top N significant markers per cluster, ranked by pct.1
## =============================================================================
## read script 14's FULL FindAllMarkers table (RDS) and apply script 14's own
## filter + ranking here (the ranked CSV did not have a 'cluster' column when read)
rds_all <- file.path(dir_mk, paste0("FindAllMarkers_", file_tag, "_all_results.RDS"))
stopifnot("script 14's _all_results.RDS not found" = file.exists(rds_all))
markers <- as.data.table(readRDS(rds_all))                                  # every gene x cluster tested
print(names(markers))                                                       # should include cluster, gene, pct.1, p_val_adj
stopifnot("RDS has no 'cluster'/'gene' column" = all(c("cluster", "gene") %in% names(markers)))

ranked <- markers[p_val_adj < 0.05 & avg_log2FC > 0]                        # significant, as script 14
ranked[, cluster := as.character(cluster)]
setorder(ranked, cluster, -pct.1, -avg_log2FC)                              # rank by pct.1, ties by avg_log2FC
ranked[, rank_pct1 := seq_len(.N), by = cluster]                            # 1 = highest pct.1 in its cluster

cl_levels <- as.character(sort(unique(as.numeric(ranked$cluster))))        # 0, 1, ..., 16
top_genes <- lapply(top_ns, function(n) {
  tg <- ranked[rank_pct1 <= n][order(as.numeric(cluster), rank_pct1)]       # top n per cluster, cluster order
  fwrite(tg[, .(cluster, rank_pct1, gene, pct.1, pct.2, avg_log2FC, p_val_adj)],
         file.path(dir_mk, paste0("Top", n, "_pct1_genes_", file_tag, ".csv")))   # gene list behind each plot
  unique(tg$gene)                                                           # a gene shared by clusters shown once (first cluster)
})
names(top_genes) <- top_ns
print(sapply(top_genes, length))                                            # genes per plot


## =============================================================================
## CHUNK 4: load object, normalise, scale the plotted genes
## =============================================================================
obj <- readRDS(seurat_out)
main_col <- paste0(prefix, "_leiden_res_", leiden_resolution)              # e.g. "harmony_leiden_res_1"
stopifnot("resolution not in 12b's sweep" = main_col %in% colnames(obj[[]]))
Idents(obj) <- factor(as.character(obj[[main_col]][, 1]), levels = cl_levels)   # clusters 0..16 in order
message(nlevels(Idents(obj)), " clusters")

obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)   # as script 14
obj <- ScaleData(obj, assay = "RNA", features = unique(unlist(top_genes)))  # z-scores, on ALL cells (heatmap colour)

## heatmap cells: up to cells_per_cluster per cluster (vignette: subset(downsample = ))
set.seed(seed)
obj_hm <- subset(obj, downsample = cells_per_cluster)                       # scaled values kept from full object
message("heatmap cells: ", ncol(obj_hm))

pal <- c("#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#42D4F4", "#F032E6",
         "#9A6324", "#469990", "#BFEF45", "#800000", "#FFE119", "#000075", "#808000",
         "#FF8FB1", "#A9A9A9", "#005300", "#009FFF", "#C77CFF", "#FF6F61", "#00FFBE")   # as 12b
cl_cols <- pal[seq_along(cl_levels)]                                        # one colour per cluster


## =============================================================================
## CHUNK 5: heatmap + dot plot for each N, saved as TIFF 300 dpi
## =============================================================================
save_tiff <- function(p, f, w, h) ggsave(file.path(dir_mk, f), p, width = w, height = h, dpi = 300,
                                         device = "tiff", compression = "lzw", bg = "white", limitsize = FALSE)

for (n in top_ns) {
  g  <- top_genes[[as.character(n)]]
  ng <- length(g)
  
  ## heatmap: rows = genes (grouped by cluster), columns = cells grouped by cluster
  p_hm <- DoHeatmap(obj_hm, features = g, group.colors = cl_cols,
                    size = 4, angle = 0, hjust = 0.5) +                     # cluster labels on top
    theme(axis.text.y = element_text(size = 6)) +                           # gene names
    labs(title = paste0("Top ", n, " significant markers per cluster (ranked by pct.1)"),
         subtitle = paste0(file_tag, "; ", ncol(obj_hm), " cells shown (<= ", cells_per_cluster,
                           " per cluster); colour = scaled expression"))
  save_tiff(p_hm, paste0("Heatmap_top", n, "_pct1_", file_tag, ".tiff"),
            w = 20, h = max(12, ng * 0.11 + 3))                             # height grows with gene number
  
  ## dot plot: all cells; dot size = % cells expressing, colour = scaled average expression
  p_dot <- DotPlot(obj, features = g) + RotatedAxis() +
    theme(axis.text.x = element_text(size = 7)) +
    labs(title = paste0("Top ", n, " significant markers per cluster (ranked by pct.1)"),
         subtitle = file_tag, x = NULL, y = "Leiden cluster")
  save_tiff(p_dot, paste0("DotPlot_top", n, "_pct1_", file_tag, ".tiff"),
            w = max(20, ng * 0.16 + 4), h = 9)                              # width grows with gene number
  message("Saved top ", n, ": ", ng, " genes")
}
message("Done -> ", dir_mk)





























## =============================================================================
## 21_InSituType_heatmap_dotplot_top_pct1.R
## =============================================================================
## WHAT IT DOES
##   For the InSituType clusters of one or more runs (ist_dirs, CHUNK 1):
##   significant markers (p_val_adj < 0.05, avg_log2FC > 0) from the InSituType
##   FindAllMarkers results, ranked within each cluster by pct.1 (ties: avg_log2FC),
##   top 10 / 15 / 20 per cluster, all clusters in one plot:
##     Heatmap_top<N>_pct1_<tag>.tiff   (Seurat DoHeatmap, scaled expression)
##     DotPlot_top<N>_pct1_<tag>.tiff   (Seurat DotPlot)
##     Genes_top<N>_pct1_<tag>.csv      (the genes plotted, per cluster)
##   Saved at 300 dpi in <InSituType folder>/Heatmap_DotPlot_top_pct1/
##   Ranking is redone by pct.1 here, whatever rank_by the InSituType script used.
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 1: settings
## =============================================================================
library(data.table)   # tables
library(Seurat)       # NormalizeData, ScaleData, DoHeatmap, DotPlot
library(ggplot2)      # titles, themes, ggsave

## InSituType folders to plot (each run separately; add/remove lines)   <- change here
ist_dirs <- c(
  "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Unsupervised_clustering_folder/QC03_04_05_06_07_08_09/InSituType_PC20_NN250_rep3",
  "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Unsupervised_clustering_folder/QC03_04_06_07_08_09/InSituType_PC20_NN250_rep3"
)

top_ns            <- c(10, 15, 20)   # genes per cluster                                  <- change here
marker_padj_max   <- 0.05            # significant = p_val_adj < this (as script 14)
cells_per_cluster <- 500             # heatmap only: cells drawn per cluster (NULL = all)  <- change here
seed              <- 42

pal <- c("#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#42D4F4", "#F032E6",
         "#9A6324", "#469990", "#BFEF45", "#800000", "#FFE119", "#000075", "#808000",
         "#FF8FB1", "#A9A9A9", "#005300", "#009FFF", "#C77CFF", "#FF6F61", "#00FFBE",
         "#D4A017", "#7F7FFF", "#8B008B", "#2F4F4F", "#E9967A", "#556B2F", "#FF1493")  # as InSituType script

## 300 dpi tiff; limitsize = FALSE allows plots taller/wider than 50 inches
save_tiff <- function(p, f, w, h) ggsave(f, p, width = w, height = h, dpi = 300, device = "tiff",
                                         compression = "lzw", bg = "white", limitsize = FALSE)


## =============================================================================
## CHUNK 2: function - heatmap + dot plot for one InSituType folder
## =============================================================================
plot_top_pct1 <- function(ist_dir) {
  ## --- files of this run ---
  file_tag   <- paste0("insitutype_", sub("^InSituType_", "", basename(ist_dir)))   # "insitutype_PC20_NN250_rep3"
  rds_all    <- file.path(ist_dir, paste0("FindAllMarkers_", file_tag, "_all_results.RDS"))
  cells_csv  <- file.path(ist_dir, paste0("InSituType_cells_", file_tag, ".csv"))
  params_csv <- file.path(ist_dir, paste0("parameters_", file_tag, ".csv"))
  stopifnot("InSituType marker / cell / parameter files not found" =
              all(file.exists(c(rds_all, cells_csv, params_csv))))
  qc_tag  <- basename(dirname(ist_dir))                                     # "QC03_04_06_07_08_09"
  dir_fig <- file.path(ist_dir, "Heatmap_DotPlot_top_pct1")
  dir.create(dir_fig, showWarnings = FALSE)
  message("\n==== ", qc_tag, " | ", file_tag, " -> ", dir_fig)
  
  ## --- significant markers, ranked by pct.1 within each cluster ---
  sig <- as.data.table(readRDS(rds_all))[p_val_adj < marker_padj_max & avg_log2FC > 0]
  sig[, cluster := as.character(cluster)]
  setorder(sig, cluster, -pct.1, -avg_log2FC)                               # highest pct.1 first
  sig[, rank_pct1 := seq_len(.N), by = cluster]                             # 1, 2, 3, ... per cluster
  
  ## --- Seurat object + InSituType cluster of every cell ---
  seurat_out <- fread(params_csv)[parameter == "input_object", value]       # object InSituType used
  obj   <- readRDS(seurat_out)
  cells <- fread(cells_csv)
  cl <- cells$insitutype_cluster[match(colnames(obj), cells$cell_id)]       # cluster per cell, object order
  stopifnot("cells missing from the InSituType table" = !anyNA(cl))
  cl_levels <- as.character(sort(unique(cl)))                               # 0, 1, 2, ... (0 = largest)
  obj$insitutype_cluster <- factor(as.character(cl), levels = cl_levels)
  Idents(obj) <- "insitutype_cluster"
  stopifnot("more clusters than colours in pal" = length(cl_levels) <= length(pal))
  cl_cols <- pal[seq_along(cl_levels)]                                      # same colours as InSituType plots
  message(length(cl_levels), " clusters; ", ncol(obj), " cells")
  
  ## --- expression: log-normalise (as script 14), then scale for the heatmap ---
  obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
  set.seed(seed)
  obj_hm <- if (is.null(cells_per_cluster)) obj else subset(obj, downsample = cells_per_cluster)   # cells per cluster
  genes_max <- unique(sig[rank_pct1 <= max(top_ns), gene])                  # every gene any heatmap needs
  obj_hm <- ScaleData(obj_hm, features = genes_max)                         # z-score per gene (Seurat vignette)
  
  ## --- one heatmap + one dot plot per top N ---
  for (n in top_ns) {
    top <- sig[rank_pct1 <= n]
    top <- top[order(factor(cluster, levels = cl_levels), rank_pct1)]       # cluster 0 genes first, then 1, ...
    genes <- unique(top$gene)                                               # a gene shared by clusters is drawn once
    fwrite(top[, .(cluster, rank_pct1, gene, pct.1, pct.2, avg_log2FC, p_val_adj)],
           file.path(dir_fig, paste0("Genes_top", n, "_pct1_", file_tag, ".csv")))
    message("top ", n, ": ", nrow(top), " cluster-gene pairs, ", length(genes), " unique genes")
    
    ## heatmap: rows = genes, columns = cells grouped by cluster
    p_hm <- DoHeatmap(obj_hm, features = genes, group.by = "insitutype_cluster",
                      group.colors = cl_cols, size = 4, angle = 0) +
      theme(axis.text.y = element_text(size = 6)) +
      labs(title = paste0("InSituType clusters (", qc_tag, "): top ", n, " significant markers per cluster by pct.1"),
           subtitle = paste0(if (is.null(cells_per_cluster)) "all cells" else paste0("up to ", cells_per_cluster, " cells per cluster"),
                             "; scaled log-normalised expression"))
    save_tiff(p_hm, file.path(dir_fig, paste0("Heatmap_top", n, "_pct1_", file_tag, ".tiff")),
              w = 20, h = 4 + 0.11 * length(genes))                         # height grows with gene number
    
    ## dot plot: size = % cells expressing, colour = scaled mean expression (all cells)
    p_dot <- DotPlot(obj, features = genes, group.by = "insitutype_cluster") + RotatedAxis() +
      theme(axis.text.x = element_text(size = 7)) +
      labs(title = paste0("InSituType clusters (", qc_tag, "): top ", n, " significant markers per cluster by pct.1"),
           x = NULL, y = "InSituType cluster")
    save_tiff(p_dot, file.path(dir_fig, paste0("DotPlot_top", n, "_pct1_", file_tag, ".tiff")),
              w = 6 + 0.2 * length(genes), h = 9)                           # width grows with gene number
  }
  rm(obj, obj_hm); invisible(gc())                                          # free memory before the next run
}


## =============================================================================
## CHUNK 3: run for every folder in ist_dirs
## =============================================================================
for (d in ist_dirs) plot_top_pct1(d)
message("Done")

























## =============================================================================
## 22_InSituType_cluster_unique_genes.R
## =============================================================================
## WHAT IT DOES
##   Finds genes expressed in ONE InSituType cluster but not in the others
##   (the single dark blocks in 07_Profiles_heatmap). FindAllMarkers alone cannot
##   do this: it tests each cluster vs ALL other cells pooled, so a gene high in
##   two clusters is a "marker" of both. Here every gene is compared with EACH
##   other cluster, using three lines of evidence:
##     1. InSituType profiles (cluster mean expression), normalised per cluster,
##        -> home cluster (highest) and runner-up (2nd highest); fold home/runner-up
##     2. % of cells detecting the gene in every cluster
##        -> detected in >= min_pct_in of home cells, <= max_pct_out of cells in every other cluster
##     3. Wilcoxon test home vs runner-up cluster (Seurat FindMarkers)
##        -> still significantly higher than its closest competitor
##   Plus: gene must be a significant FindAllMarkers marker of the home cluster.
## OUTPUT (per InSituType folder, in <InSituType folder>/Cluster_unique_genes/):
##   Unique_genes_<tag>.xlsx                     Summary + one sheet per cluster + all genes
##   Heatmap_unique_genes_top<N>_<tag>.tiff      profiles of unique genes (300 dpi)
##   DotPlot_unique_genes_top<N>_<tag>.tiff      per-cell check of the same genes (300 dpi)
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## 22_InSituType_cluster_unique_genes.R
## =============================================================================
## WHAT IT DOES
##   Finds genes expressed in ONE InSituType cluster but not in the others
##   (the single dark blocks in 07_Profiles_heatmap). FindAllMarkers alone cannot
##   do this: it tests each cluster vs ALL other cells pooled, so a gene high in
##   two clusters is a "marker" of both. Here every gene is compared with EACH
##   other cluster, using three lines of evidence:
##     1. InSituType profiles (cluster mean expression), normalised per cluster,
##        -> home cluster (highest) and runner-up (2nd highest); fold home/runner-up
##     2. % of cells detecting the gene in every cluster
##        -> detected in >= min_pct_in of home cells and >= min_pct_ratio x more often
##           than in every other cluster (relative: CosMx counts per cell are low)
##     3. Wilcoxon test home vs runner-up cluster (Seurat FindMarkers)
##        -> still significantly higher than its closest competitor
##   Plus: gene must be a significant FindAllMarkers marker of the home cluster.
## OUTPUT (per InSituType folder, in <InSituType folder>/Cluster_unique_genes/):
##   Unique_genes_<tag>.xlsx                     Summary + one sheet per cluster + all genes
##   Heatmap_unique_genes_top<N>_<tag>.tiff      profiles of unique genes (300 dpi)
##   DotPlot_unique_genes_top<N>_<tag>.tiff      per-cell check of the same genes (300 dpi)
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 1: settings
## =============================================================================
library(data.table)   # tables
library(Matrix)       # sparse counts, fac2sparse
library(Seurat)       # NormalizeData, FindMarkers, DotPlot
library(ggplot2)      # ggsave
library(openxlsx)     # Excel output
library(pheatmap)     # heatmap (as plot 07)

## InSituType folders (add/remove lines)   <- change here
ist_dirs <- c(
  "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Unsupervised_clustering_folder/QC03_04_05_06_07_08_09/InSituType_PC20_NN250_rep3",
  "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Unsupervised_clustering_folder/QC03_04_06_07_08_09/InSituType_PC20_NN250_rep3"
)

## data are sparse (median 80-700 counts per cell over ~6,000 genes), so detection is
## judged RELATIVE to other clusters, not by absolute cut-offs (0.25 / 0.10 gave 0 genes)
min_pct_in    <- 0.10   # detected (count > 0) in >= 10% of home-cluster cells        <- change here
min_pct_ratio <- 3      # home detection >= 3x detection in EVERY other cluster       <- change here
min_fold      <- 2      # normalised profile: home >= 2x runner-up cluster            <- change here
qc_pct        <- 0.10   # QC only: % panel genes detected in > 10% of a cluster's cells
padj_max    <- 0.05   # significance for vs-rest and vs-runner-up tests (as script 14)
pseudo      <- 0.01   # added to profiles before the fold (avoids divide by 0)
n_plot      <- 10     # unique genes per cluster shown in heatmap / dot plot          <- change here

pal <- c("#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#42D4F4", "#F032E6",
         "#9A6324", "#469990", "#BFEF45", "#800000", "#FFE119", "#000075", "#808000",
         "#FF8FB1", "#A9A9A9", "#005300", "#009FFF", "#C77CFF", "#FF6F61", "#00FFBE",
         "#D4A017", "#7F7FFF", "#8B008B", "#2F4F4F", "#E9967A", "#556B2F", "#FF1493")  # as InSituType script


## =============================================================================
## CHUNK 2: function - unique genes for one InSituType folder
## =============================================================================
find_unique <- function(ist_dir) {
  
  ## --- files of this run (written by the InSituType script) ---
  file_tag   <- paste0("insitutype_", sub("^InSituType_", "", basename(ist_dir)))   # "insitutype_PC20_NN250_rep3"
  prof_csv   <- file.path(ist_dir, paste0("InSituType_profiles_", file_tag, ".csv"))
  rds_all    <- file.path(ist_dir, paste0("FindAllMarkers_", file_tag, "_all_results.RDS"))
  cells_csv  <- file.path(ist_dir, paste0("InSituType_cells_", file_tag, ".csv"))
  params_csv <- file.path(ist_dir, paste0("parameters_", file_tag, ".csv"))
  stopifnot("InSituType files not found" = all(file.exists(c(prof_csv, rds_all, cells_csv, params_csv))))
  qc_tag <- basename(dirname(ist_dir))                                      # "QC03_04_06_07_08_09"
  dir_u  <- file.path(ist_dir, "Cluster_unique_genes")
  dir.create(dir_u, showWarnings = FALSE)
  message("\n==== ", qc_tag, " | ", file_tag, " -> ", dir_u)
  
  ## --- 1. InSituType profiles, normalised per cluster ---
  prof <- fread(prof_csv)                                                   # gene + cluster_0, cluster_1, ...
  P <- as.matrix(prof[, -1, with = FALSE]); rownames(P) <- prof$gene        # genes x clusters
  colnames(P) <- sub("^cluster_", "", colnames(P))                          # "0", "1", ...
  cl_levels <- colnames(P)                                                  # 0 = largest cluster
  Pn <- sweep(P, 2, colSums(P), "/") * 1000                                 # per 1,000 counts: removes total-count differences between clusters
  
  ord   <- t(apply(Pn, 1, order, decreasing = TRUE))                        # cluster ranks per gene
  i_top <- ord[, 1]; i_sec <- ord[, 2]                                      # highest + 2nd highest cluster
  rows  <- seq_len(nrow(Pn))
  spec <- data.table(gene           = rownames(Pn),
                     home           = cl_levels[i_top],                     # cluster with highest expression
                     home_profile   = Pn[cbind(rows, i_top)],
                     runner_up      = cl_levels[i_sec],                     # closest competitor
                     runner_profile = Pn[cbind(rows, i_sec)])
  spec[, fold_vs_runner := (home_profile + pseudo) / (runner_profile + pseudo)]
  mx <- apply(Pn, 1, max)
  spec[, tau := ifelse(mx > 0, rowSums(1 - Pn / mx) / (ncol(Pn) - 1), NA)]   # 0 = everywhere, 1 = one cluster only
  
  ## --- 2. % cells detecting each gene, per cluster ---
  seurat_out <- fread(params_csv)[parameter == "input_object", value]       # object InSituType used
  obj   <- readRDS(seurat_out)
  cells <- fread(cells_csv)
  cl <- factor(as.character(cells$insitutype_cluster[match(colnames(obj), cells$cell_id)]), levels = cl_levels)
  stopifnot("cells missing from the InSituType table" = !anyNA(cl),
            "profile genes missing from object" = all(rownames(Pn) %in% rownames(obj)))
  obj$insitutype_cluster <- cl
  
  cnt <- LayerData(obj, assay = "RNA", layer = "counts")[rownames(Pn), ]    # genes x cells, raw counts
  det <- (cnt > 0) * 1                                                      # 1 = detected in that cell
  n_cl <- as.numeric(table(cl)[cl_levels])                                  # cells per cluster
  pct <- as.matrix(det %*% Matrix::t(Matrix::fac2sparse(cl)))               # genes x clusters, cells detecting
  pct <- sweep(pct, 2, n_cl, "/"); colnames(pct) <- cl_levels               # -> fraction of cluster's cells
  rm(cnt, det); invisible(gc())
  
  i_home <- match(spec$home, cl_levels)
  spec[, pct_home := pct[cbind(rows, i_home)]]                              # detection in home cluster
  pct_o <- pct; pct_o[cbind(rows, i_home)] <- NA                            # other clusters only
  spec[, max_pct_other := apply(pct_o, 1, max, na.rm = TRUE)]               # highest detection elsewhere
  spec[, max_pct_other_cluster := cl_levels[apply(pct_o, 1, which.max)]]
  spec[, pct_ratio := pct_home / (max_pct_other + 0.001)]                   # home vs best other cluster
  spec[, n_other_clusters_above := rowSums(pct_o > pct_home / min_pct_ratio, na.rm = TRUE)]   # clusters within 3x of home; 0 = unique
  
  ## --- 3. FindAllMarkers: significant marker of its home cluster vs rest? ---
  mk <- as.data.table(readRDS(rds_all)); mk[, cluster := as.character(cluster)]
  spec <- merge(spec, mk[, .(home = cluster, gene, avg_log2FC_vs_rest = avg_log2FC,
                             padj_vs_rest = p_val_adj)], by = c("home", "gene"), all.x = TRUE)
  spec[, marker_vs_rest := !is.na(padj_vs_rest) & padj_vs_rest < padj_max & avg_log2FC_vs_rest > 0]
  
  ## --- 4. candidates, then Wilcoxon home vs runner-up cluster ---
  cand <- spec[marker_vs_rest & pct_home >= min_pct_in & fold_vs_runner >= min_fold]
  message(nrow(cand), " candidate genes; testing each vs its runner-up cluster")
  obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)
  Idents(obj) <- "insitutype_cluster"
  pairs <- unique(cand[, .(home, runner_up)])                               # one test per home/runner-up pair
  tests <- rbindlist(lapply(seq_len(nrow(pairs)), function(i) {
    h <- pairs$home[i]; r <- pairs$runner_up[i]
    g <- cand[home == h & runner_up == r, gene]
    fm <- FindMarkers(obj, ident.1 = h, ident.2 = r, features = g, test.use = "wilcox",
                      logfc.threshold = 0, min.pct = 0)                     # test every candidate, no pre-filter
    data.table(home = h, gene = rownames(fm), avg_log2FC_vs_runner = fm$avg_log2FC,
               padj_vs_runner = fm$p_val_adj)                               # Bonferroni, as FindAllMarkers
  }))
  spec <- merge(spec, tests, by = c("home", "gene"), all.x = TRUE)
  
  ## --- 5. classify ---
  spec[, specific := marker_vs_rest & pct_home >= min_pct_in & fold_vs_runner >= min_fold &
         !is.na(padj_vs_runner) & padj_vs_runner < padj_max & avg_log2FC_vs_runner > 0]
  spec[, category := fifelse(specific & n_other_clusters_above == 0, "unique",
                             fifelse(specific & n_other_clusters_above == 1, "shared_with_1_cluster",
                                     fifelse(specific, "shared_with_2plus", "not_specific")))]
  ## how many genes survive each filter (shows which step removes them)
  funnel <- data.table(step = c("all genes", "1. marker vs rest", "2. + pct_home >= min_pct_in",
                                "3. + fold_vs_runner >= min_fold", "4. + significant vs runner-up",
                                "5. + unique (no other cluster within min_pct_ratio)"),
                       n_genes = c(nrow(spec), spec[, sum(marker_vs_rest)],
                                   spec[, sum(marker_vs_rest & pct_home >= min_pct_in)],
                                   nrow(cand), spec[, sum(specific)], spec[, sum(category == "unique")]))
  print(funnel)
  
  spec[, home := factor(home, levels = cl_levels)]
  setorder(spec, home, -fold_vs_runner, -pct_home)                          # most specific first
  spec[, (c("home_profile", "runner_profile", "fold_vs_runner", "tau", "pct_home", "max_pct_other", "pct_ratio")) :=
         lapply(.SD, round, 3), .SDcols = c("home_profile", "runner_profile", "fold_vs_runner", "tau", "pct_home", "max_pct_other", "pct_ratio")]
  uniq <- spec[category == "unique"]
  
  ## --- 6. cluster QC: clusters high for most genes (5, 3, 6 in plot 07)? ---
  md <- obj[[]]
  qc <- data.table(cluster = cl_levels, n_cells = as.integer(n_cl),
                   profile_total = round(colSums(P), 1),                    # high = high counts per cell
                   median_nCount = if ("nCount_RNA" %in% names(md)) as.numeric(tapply(md$nCount_RNA, cl, median)[cl_levels]) else NA,
                   mean_negmean = if ("negmean" %in% names(md)) signif(as.numeric(tapply(md$negmean, cl, mean)[cl_levels]), 3) else NA,   # mean: median is 0 (most cells have no neg counts)
                   pct_panel_detected_10pct = round(100 * colMeans(pct > qc_pct), 1),   # % panel genes in >10% of cells
                   n_unique = as.integer(table(factor(uniq$home, levels = cl_levels))),
                   n_shared_with_1 = as.integer(table(factor(spec[category == "shared_with_1_cluster", home], levels = cl_levels))),
                   top10_unique = sapply(cl_levels, function(k) paste(head(uniq[home == k, gene], 10), collapse = ", ")))
  print(qc[, .(cluster, n_cells, median_nCount, mean_negmean, pct_panel_detected_10pct, n_unique, n_shared_with_1)])
  
  ## --- 7. Excel ---
  cols_out <- c("home", "gene", "category", "fold_vs_runner", "runner_up", "pct_home", "max_pct_other",
                "max_pct_other_cluster", "pct_ratio", "n_other_clusters_above", "tau", "avg_log2FC_vs_rest", "padj_vs_rest",
                "avg_log2FC_vs_runner", "padj_vs_runner", "home_profile", "runner_profile")
  per_cl <- lapply(cl_levels, function(k) {                                 # unique + near-unique per cluster
    d <- spec[home == k & category != "not_specific", ..cols_out]
    if (nrow(d) == 0) data.frame(note = "no specific genes at these thresholds") else as.data.frame(d)
  })
  names(per_cl) <- paste0("cluster_", cl_levels)
  how_to_read <- data.frame(
    item = c("home / runner_up", "fold_vs_runner", "pct_home", "max_pct_other", "pct_ratio", "n_other_clusters_above",
             "tau", "padj_vs_rest", "padj_vs_runner", "unique", "shared_with_1_cluster",
             "profile_total / median_nCount / mean_negmean", "pct_panel_detected_10pct"),
    meaning = c("cluster with highest / 2nd highest InSituType profile (profiles normalised to 1,000 counts per cluster)",
                paste0("home / runner-up normalised profile (+", pseudo, "); >= ", min_fold, " required"),
                paste0("fraction of home-cluster cells with count > 0; >= ", min_pct_in, " required"),
                "highest detection fraction in any other cluster (and which cluster)",
                paste0("pct_home / max_pct_other; >= ", min_pct_ratio, " = unique"),
                paste0("other clusters detecting the gene in more than pct_home / ", min_pct_ratio, " of their cells (0 = unique)"),
                "tissue-specificity index across cluster profiles: 0 = equal everywhere, 1 = one cluster only",
                "FindAllMarkers home cluster vs all other cells (Bonferroni)",
                "FindMarkers home vs runner-up cluster (Wilcoxon, Bonferroni)",
                paste0("specific (all tests passed) and detected >= ", min_pct_ratio, "x more often than in every other cluster"),
                "specific, but one other cluster also detects it (often two states of one cell type)",
                "high values = cluster has more transcripts or background per cell; such clusters look high for most genes in plot 07",
                paste0("% of panel genes detected in > ", qc_pct, " of the cluster's cells (high = broadly expressing / possible doublets or spillover)")))
  out_xlsx <- file.path(dir_u, paste0("Unique_genes_", file_tag, ".xlsx"))
  write.xlsx(c(list(Summary = as.data.frame(qc),
                    All_unique = as.data.frame(uniq[, ..cols_out])),
               per_cl,
               list(All_genes = as.data.frame(spec[, ..cols_out]), How_to_read = how_to_read)),
             out_xlsx, firstRow = TRUE, colWidths = "auto", overwrite = TRUE)
  message("Saved: ", out_xlsx)
  
  ## --- 8. heatmap: unique genes, profiles scaled 0-1 per gene (as plot 07) ---
  hm <- uniq[, head(.SD, n_plot), by = home]                                # top n_plot per cluster, cluster order
  if (nrow(hm) == 0) { message("no unique genes - relax thresholds in CHUNK 1"); rm(obj); return(invisible(spec)) }
  mat <- Pn[hm$gene, , drop = FALSE]
  mat <- sweep(mat, 1, apply(mat, 1, max), "/")                             # 0-1 per gene
  ann_row <- data.frame(cluster = factor(as.character(hm$home), levels = cl_levels), row.names = hm$gene)
  ann_col <- list(cluster = setNames(pal[seq_along(cl_levels)], cl_levels))
  gaps <- head(cumsum(table(droplevels(ann_row$cluster))), -1)              # line between clusters
  pheatmap(mat, cluster_rows = FALSE, cluster_cols = FALSE,                 # keep block order
           color = colorRampPalette(c("white", "darkblue"))(100),
           annotation_row = ann_row, annotation_colors = ann_col, gaps_row = gaps,
           fontsize_row = 6, fontsize_col = 9, angle_col = 90,
           main = paste0(qc_tag, ": genes unique to one InSituType cluster (top ", n_plot, " per cluster)"),
           filename = file.path(dir_u, paste0("Heatmap_unique_genes_top", n_plot, "_", file_tag, ".tiff")),  # tiff = 300 dpi
           width = 8, height = 2 + 0.12 * nrow(mat))                        # height grows with gene number
  
  ## --- 9. dot plot: same genes, per cell ---
  p_dot <- DotPlot(obj, features = unique(hm$gene), group.by = "insitutype_cluster") + RotatedAxis() +
    theme(axis.text.x = element_text(size = 7)) +
    labs(title = paste0(qc_tag, ": genes unique to one InSituType cluster (top ", n_plot, " per cluster)"),
         x = NULL, y = "InSituType cluster")
  ggsave(file.path(dir_u, paste0("DotPlot_unique_genes_top", n_plot, "_", file_tag, ".tiff")), p_dot,
         width = 6 + 0.2 * nrow(hm), height = 8, dpi = 300, device = "tiff",
         compression = "lzw", bg = "white", limitsize = FALSE)
  
  rm(obj); invisible(gc())                                                  # free memory before next run
  invisible(spec)
}


## =============================================================================
## CHUNK 3: run for every folder
## =============================================================================
res_unique <- lapply(ist_dirs, find_unique)
names(res_unique) <- basename(dirname(ist_dirs))                            # by QC tag
message("Done")

























## =============================================================================
## Cell typing leiden FindAllmarkeers AI cell typist etc QC1
## =============================================================================
## =============================================================================
## 23_UMAP_Leiden_labelled_cell_types.R Cell typing leiden QC1
## =============================================================================
## WHAT IT DOES
##   Labelled versions of UMAP_by_Leiden_cluster_<file_tag>.tiff (script 14 copy):
##   same Seurat object, same UMAP, same Leiden clusters (main resolution).
##   Nothing is recomputed.
##     A. UMAP_Leiden_broad_cell_types_<file_tag>.tiff
##        cells coloured by broad cell type; C0-C16 labels on the clusters
##     B. UMAP_Leiden_cell_states_<file_tag>.tiff
##        cells coloured by cell state / sub-cell type (shades of the broad-type
##        colour); C0-C16 on the map, full names in the legend
##     Leiden_cluster_annotation_<file_tag>.csv  (cluster, broad type, state, n cells)
##   Saved at 300 dpi in the same Results folder as the original UMAP.
## Labelling follows common practice in atlas papers: broad types are split across
## the map (e.g. Fibroblast = C0, C3, C5, C11, C16), so one name per type on the
## map would sit in the wrong place -> cluster numbers on the map, names in a legend.
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 1: settings
## =============================================================================
library(data.table)   # tables
library(Seurat)       # Embeddings
library(ggplot2)      # plots
library(ggrepel)      # labels that do not overlap (installed with Seurat)

## same path as script 14 (the UMAP that was labelled)   <- change here
clustering_path <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Clustering_PCA_UMAP_Leiden_Markers/QC03_04_05_06_07_08_09/PC20_NN250_rep3/02_SENSITIVITY_batchawarePCA_Harmony/UMAP_by_Leiden_cluster.tiff"
leiden_resolution <- NULL     # NULL = main resolution from 12b (res 1)

label_style <- "number"       # plot B: "number" = C0-C16 on map (names in legend); "name" = full names on map   <- change here
pt_size  <- 0.05              # as original UMAP
pt_alpha <- 0.3               # as original UMAP
seed     <- 42

annot <- data.table(
  cluster = as.character(0:16),
  broad = c("Fibroblast", "Pericyte", "Endothelial", "Fibroblast", "VSMC", "Fibroblast", "VSMC",
            "Immune", "Pericyte", "VSMC", "Immune", "Fibroblast", "Pericyte", "VSMC", "VSMC",
            "VSMC", "Fibroblast"),
  state = c("Osteochondrogenic Fibroblast",          # C0
            "ECM-producing Pericyte",                # C1
            "Inflammatory Endothelial Cell",         # C2
            "Inflammatory Adventitial Fibroblast",   # C3
            "Inflammatory Synthetic VSMC",           # C4
            "Matrix-Remodeling Fibroblast",          # C5
            "Synthetic VSMC",                        # C6
            "T Cell",                                # C7
            "Activated Fibrogenic Pericyte",           # C8
            "Osteogenic / Calcifying VSMC",          # C9
            "Plasma Cell",                           # C10
            "Homeostatic Adventitial Fibroblast",      # C11
            "Contractile Pericyte",                  # C12
            "Fibromyocyte",                          # C13
            "Profibrotic Synthetic VSMC",            # C14
            "Contractile VSMC",                      # C15
            "Senescent Adventitial Fibroblast"))     # C16

## colours: one hue per broad type; states = shades of that hue (dark -> light)
broad_levels <- c("VSMC", "Fibroblast", "Pericyte", "Endothelial", "Immune")   # legend order
broad_cols <- c(VSMC = "#CB181D", Fibroblast = "#2171B5", Pericyte = "#F16913",
                Endothelial = "#238B45", Immune = "#6A51A3")
ramps <- list(VSMC        = c("#67000D", "#CB181D", "#FB6A4A"),
              Fibroblast  = c("#08306B", "#2171B5", "#6BAED6"),
              Pericyte    = c("#7F2704", "#F16913", "#FDAE6B"),
              Endothelial = c("#238B45"),
              Immune      = c("#3F007D", "#9E9AC8"))


## =============================================================================
## CHUNK 2: paths from clustering_path (as script 14)
## =============================================================================
branch_dir <- if (dir.exists(clustering_path)) clustering_path else dirname(clustering_path)
branch_dir <- normalizePath(branch_dir, mustWork = TRUE)                    # stop if missing
dir_seu    <- dirname(branch_dir)                                           # .../PC20_NN250_rep3
run_tag    <- basename(dir_seu)                                             # "PC20_NN250_rep3"
qc_tag     <- basename(dirname(dir_seu))                                    # "QC03_04_05_06_07_08_09"
output_dir <- dirname(dirname(dirname(dir_seu)))                            # .../AAA
prefix <- if (grepl("^01_PRIMARY", basename(branch_dir))) "primary" else
  if (grepl("^02_SENSITIVITY", basename(branch_dir))) "harmony" else
    stop("branch folder must start with 01_PRIMARY or 02_SENSITIVITY")
if (is.null(leiden_resolution)) {
  params <- fread(file.path(dir_seu, paste0("parameters_used_", run_tag, ".csv")))
  leiden_resolution <- as.numeric(params[parameter == "Leiden_resolution", value])
}
seurat_out <- list.files(dir_seu, pattern = paste0("^seurat_.*_SeuratFunctions_", run_tag, "\\.RDS$"), full.names = TRUE)
stopifnot("expected exactly one 12b Seurat object" = length(seurat_out) == 1)

dir_mk   <- file.path(output_dir, "Results", paste0(qc_tag, "_final_All_markers"))   # where the original UMAP copy is
file_tag <- paste0(prefix, "_", run_tag, "_res", leiden_resolution)        # "harmony_PC20_NN250_rep3_res1"
main_col <- paste0(prefix, "_leiden_res_", leiden_resolution)              # Leiden column
umap_name <- paste0("umap_", prefix)                                        # UMAP of the original plot
message("UMAP: ", umap_name, " | clusters: ", main_col, " | saving to ", dir_mk)


## =============================================================================
## CHUNK 3: UMAP coordinates + annotation per cell
## =============================================================================
obj <- readRDS(seurat_out)
stopifnot("UMAP not in object" = umap_name %in% Reductions(obj),
          "Leiden column not in object" = main_col %in% colnames(obj[[]]))
um <- Embeddings(obj, umap_name)                                            # same coordinates as the original plot
df <- data.table(UMAP_1 = um[, 1], UMAP_2 = um[, 2],
                 cluster = as.character(obj[[main_col]][, 1]))
stopifnot("clusters without annotation" = all(unique(df$cluster) %in% annot$cluster),
          "annotation for clusters not in data" = all(annot$cluster %in% df$cluster))
rm(obj); invisible(gc())

## state colours: shades within each broad type, ordered by type then cluster number
annot[, broad := factor(broad, levels = broad_levels)]
annot <- annot[order(broad, as.numeric(cluster))]                          # VSMC first, then Fibroblast, ...
annot[, state_col := if (.N == 1) ramps[[as.character(broad[1])]][1]
      else colorRampPalette(ramps[[as.character(broad[1])]])(.N), by = broad]
annot[, state_lab := paste0("C", cluster, "  ", state)]                     # legend text, e.g. "C9  Osteogenic / Calcifying VSMC"
annot[, n_cells := df[, .N, by = cluster][match(annot$cluster, cluster), N]]

df[, broad := factor(as.character(annot$broad[match(cluster, annot$cluster)]), levels = broad_levels)]
df[, state_lab := factor(annot$state_lab[match(cluster, annot$cluster)], levels = annot$state_lab)]
set.seed(seed); df <- df[sample(.N)]                                        # random draw order (no cluster on top)

## label positions: median UMAP position of each cluster (as Seurat LabelClusters)
cen <- df[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = cluster]
cen <- merge(cen, annot[, .(cluster, broad, state)], by = "cluster")

## legend text for broad types with cell numbers
n_broad <- df[, .N, by = broad]
broad_lab <- setNames(paste0(broad_levels, "  (n = ", format(n_broad$N[match(broad_levels, n_broad$broad)], big.mark = ","), ")"),
                      broad_levels)

fwrite(annot[, .(cluster = paste0("C", cluster), broad_cell_type = broad, cell_state = state, n_cells)],
       file.path(dir_mk, paste0("Leiden_cluster_annotation_", file_tag, ".csv")))

## shared look (clean, as in atlas papers)
theme_umap <- theme_classic(base_size = 14) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),         # UMAP units carry no meaning
        legend.text = element_text(size = 11), legend.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold"))
save_tiff <- function(p, f, w, h) ggsave(file.path(dir_mk, f), p, width = w, height = h, dpi = 300,
                                         device = "tiff", compression = "lzw", bg = "white")


## =============================================================================
## CHUNK 4: plot A - broad cell types
## =============================================================================
p_broad <- ggplot(df, aes(UMAP_1, UMAP_2, colour = broad)) +
  geom_point(size = pt_size, alpha = pt_alpha, shape = 16) +
  geom_label_repel(data = cen, aes(label = paste0("C", cluster)), colour = "black",
                   fill = alpha("white", 0.85), size = 3.8, fontface = "bold",
                   label.padding = 0.15, box.padding = 0.3, min.segment.length = 0,
                   max.overlaps = Inf, seed = seed, show.legend = FALSE) +   # cluster numbers on the map
  scale_colour_manual(values = broad_cols, labels = broad_lab, name = "Broad cell type") +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) + # big, solid legend dots
  coord_fixed() + theme_umap +
  labs(title = "Broad cell types", x = "UMAP 1", y = "UMAP 2",
       subtitle = paste0("Leiden clusters (res ", leiden_resolution, ", ", nrow(annot), " clusters), ",
                         format(nrow(df), big.mark = ","), " cells"))
save_tiff(p_broad, paste0("UMAP_Leiden_broad_cell_types_", file_tag, ".tiff"), 12, 9)


## =============================================================================
## CHUNK 5: plot B - cell states / sub-cell types
## =============================================================================
cen[, lab := if (label_style == "name") state else paste0("C", cluster)]    # text on the map
p_state <- ggplot(df, aes(UMAP_1, UMAP_2, colour = state_lab)) +
  geom_point(size = pt_size, alpha = pt_alpha, shape = 16) +
  geom_label_repel(data = cen, aes(label = lab), colour = "black",
                   fill = alpha("white", 0.85), size = if (label_style == "name") 3.2 else 3.8,
                   fontface = "bold", label.padding = 0.15, box.padding = 0.4,
                   min.segment.length = 0, max.overlaps = Inf, seed = seed, show.legend = FALSE) +
  scale_colour_manual(values = setNames(annot$state_col, annot$state_lab),
                      name = "Cell state / sub-cell type") +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5, alpha = 1))) +
  coord_fixed() + theme_umap +
  theme(legend.key.height = unit(0.55, "cm")) +                             # spacing between legend rows
  labs(title = "Cell states / sub-cell types", x = "UMAP 1", y = "UMAP 2",
       subtitle = "Colour shades grouped by broad type: VSMC (red), Fibroblast (blue), Pericyte (orange), Endothelial (green), Immune (purple)")
save_tiff(p_state, paste0("UMAP_Leiden_cell_states_", file_tag, ".tiff"), 15, 9)
message("Saved both plots + annotation table -> ", dir_mk)












## =============================================================================
## 23b_UMAP_state_label_options.R   Cell typing leiden QC1
## =============================================================================
## WHAT IT DOES
##   Five alternative ways to label the cell-state UMAP (plot B of script 23).
##   The problem with the original: labels are black-on-white "C0..C16" boxes, so
##   you cannot tell from the label which legend colour it belongs to, and the
##   shades within one broad type look alike.
##     Option 1  labels filled with the cluster's own colour        (link label <-> legend)
##     Option 2  one clearly different colour per cluster           (no look-alike shades)
##     Option 3  "C number + name" on the map, no legend            (no lookup needed)
##     Option 4  one small panel per cluster, rest of cells grey    (highlight plot)
##     Option 5  one panel per broad type, states coloured inside   (shades compared within type only)
##   Each option saves its own TIFF (300 dpi) in dir_mk, next to the script-23 plots.
## REQUIRES: script 23 CHUNKS 1-4 already run in this session (df, cen, annot,
##   theme_umap, save_tiff, file_tag, pt_size, pt_alpha, seed exist).
##   ggplot2 >= 3.5.0 (for I() inside aes; check with packageVersion("ggplot2")).
## Run chunk by chunk (Cmd+Enter). Run CHUNK 0 first, then any option you want.
## =============================================================================


## =============================================================================
## CHUNK 0: shared pieces for all options
## =============================================================================
## add legend text and the cluster's colour to the label table (cen = one row per cluster)
cen[, state_lab := annot$state_lab[match(cluster, annot$cluster)]]   # e.g. "C9  Osteogenic / Calcifying VSMC"
cen[, fill_col  := annot$state_col[match(cluster, annot$cluster)]]   # same shade as the cluster's cells

## text colour that stays readable on the fill: white on dark fills, black on light fills
text_on <- function(hex) {
  rgb <- col2rgb(hex)                                                # 3 x n matrix of red, green, blue (0-255)
  lum <- 0.299 * rgb[1, ] + 0.587 * rgb[2, ] + 0.114 * rgb[3, ]      # perceived brightness (standard luma weights)
  ifelse(lum > 150, "black", "white")
}
cen[, txt_col := text_on(fill_col)]

## label texts
cen[, lab_num  := paste0("C", cluster)]                              # "C9"
cen[, lab_full := paste0("C", cluster, "\n",                          # "C9" + name wrapped to ~16 characters per line
                         vapply(state, function(s) paste(strwrap(s, 16), collapse = "\n"), ""))]

## coloured label layer used by options 1, 2, 3, 5
## I() = use the colour values as they are (no colour scale), so the label fill
## does not clash with the colour scale of the points
lab_layer <- function(lab, fill, txt, size = 3.8) {
  geom_label_repel(data = cen,
                   aes(label = .data[[lab]], fill = I(.data[[fill]]), colour = I(.data[[txt]])),
                   size = size, fontface = "bold", label.padding = 0.15, box.padding = 0.4,
                   min.segment.length = 0, segment.colour = "black",             # line from label to cluster
                   max.overlaps = Inf, seed = seed, show.legend = FALSE)
}

## state colour scale + legend as in script 23 plot B
state_scale <- scale_colour_manual(values = setNames(annot$state_col, annot$state_lab),
                                   name = "Cell state / sub-cell type")
state_guide <- guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5, alpha = 1)))

## background cells for the panel plots (options 4, 5): this table has no
## state_lab / broad column, so ggplot draws it in every panel
bg <- df[, .(UMAP_1, UMAP_2)]


## =============================================================================
## OPTION 1: labels filled with the cluster colour
## =============================================================================
## same plot as script 23 plot B; only the label boxes change colour, so each
## "C" label matches its legend dot
p_opt1 <- ggplot(df, aes(UMAP_1, UMAP_2, colour = state_lab)) +
  geom_point(size = pt_size, alpha = pt_alpha, shape = 16) +
  lab_layer("lab_num", "fill_col", "txt_col") +
  state_scale + state_guide +
  coord_fixed() + theme_umap + theme(legend.key.height = unit(0.55, "cm")) +
  labs(title = "Cell states / sub-cell types", x = "UMAP 1", y = "UMAP 2")
save_tiff(p_opt1, paste0("UMAP_Leiden_cell_states_opt1_colouredlabels_", file_tag, ".tiff"), 15, 9)


## =============================================================================
## OPTION 2: one distinct colour per cluster (plus coloured labels)
## =============================================================================
## 17 maximally different colours (Trubetskoy "20 simple distinct colours",
## white/black/grey/beige/apricot removed); legend order stays grouped by broad type
distinct17 <- c("#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231", "#911EB4",
                "#42D4F4", "#F032E6", "#BFEF45", "#FABED4", "#469990", "#DCBEFF",
                "#9A6324", "#800000", "#AAFFC3", "#808000", "#000075")
pal2 <- setNames(distinct17[seq_len(nrow(annot))], annot$state_lab)  # one colour per legend entry
cen[, fill2 := pal2[state_lab]]                                      # label fill = new cluster colour
cen[, txt2  := text_on(fill2)]

p_opt2 <- ggplot(df, aes(UMAP_1, UMAP_2, colour = state_lab)) +
  geom_point(size = pt_size, alpha = pt_alpha, shape = 16) +
  lab_layer("lab_num", "fill2", "txt2") +
  scale_colour_manual(values = pal2, name = "Cell state / sub-cell type") + state_guide +
  coord_fixed() + theme_umap + theme(legend.key.height = unit(0.55, "cm")) +
  labs(title = "Cell states / sub-cell types", x = "UMAP 1", y = "UMAP 2")
save_tiff(p_opt2, paste0("UMAP_Leiden_cell_states_opt2_distinctcolours_", file_tag, ".tiff"), 15, 9)


## =============================================================================
## OPTION 3: cluster number + name on the map, no legend
## =============================================================================
## the reader does not need to look anything up; works if the labels do not
## crowd too much (lower size if they do)
p_opt3 <- ggplot(df, aes(UMAP_1, UMAP_2, colour = state_lab)) +
  geom_point(size = pt_size, alpha = pt_alpha, shape = 16) +
  lab_layer("lab_full", "fill_col", "txt_col", size = 2.8) +
  state_scale +
  coord_fixed() + theme_umap + theme(legend.position = "none") +    # names are on the map
  labs(title = "Cell states / sub-cell types", x = "UMAP 1", y = "UMAP 2")
save_tiff(p_opt3, paste0("UMAP_Leiden_cell_states_opt3_namesonmap_", file_tag, ".tiff"), 13, 11)


## =============================================================================
## OPTION 4: one panel per cluster (highlight plot)
## =============================================================================
## every panel shows all cells in grey and one cluster in colour; the panel
## title is the cluster name. No ambiguity, often used as a supplementary figure
p_opt4 <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(data = bg, colour = "grey85", size = pt_size, alpha = pt_alpha, shape = 16) +  # all cells
  geom_point(aes(colour = state_lab), size = pt_size * 2, alpha = 0.6, shape = 16) +       # this cluster
  scale_colour_manual(values = setNames(annot$state_col, annot$state_lab), guide = "none") +
  facet_wrap(~ state_lab, ncol = 5, labeller = label_wrap_gen(width = 25)) +              # one panel per state
  coord_fixed() + theme_umap +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 10)) +
  labs(title = "Cell states / sub-cell types (one panel per Leiden cluster)", x = "UMAP 1", y = "UMAP 2")
save_tiff(p_opt4, paste0("UMAP_Leiden_cell_states_opt4_perclusterpanels_", file_tag, ".tiff"), 20, 16)


## =============================================================================
## OPTION 5: one panel per broad type
## =============================================================================
## similar shades only have to be told apart inside their own panel; labels
## appear only in the panel of their broad type (cen has a broad column)
p_opt5 <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(data = bg, colour = "grey85", size = pt_size, alpha = pt_alpha, shape = 16) +  # all cells
  geom_point(aes(colour = state_lab), size = pt_size * 2, alpha = 0.5, shape = 16) +       # cells of this type
  lab_layer("lab_num", "fill_col", "txt_col", size = 3.2) +
  state_scale + state_guide +
  facet_wrap(~ broad, ncol = 3) +                                                          # one panel per broad type
  coord_fixed() + theme_umap +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 13),
        legend.key.height = unit(0.55, "cm")) +
  labs(title = "Cell states / sub-cell types by broad cell type", x = "UMAP 1", y = "UMAP 2")
save_tiff(p_opt5, paste0("UMAP_Leiden_cell_states_opt5_perbroadtype_", file_tag, ".tiff"), 20, 11)

message("Saved label options 1-5 -> ", dir_mk)























## =============================================================================
## 24_UMAP_label_options.R
## =============================================================================
## WHAT IT DOES
##   Four ways of labelling the SAME UMAP (umap_harmony, Leiden res 1, from 12b),
##   for both annotation levels:
##     M2 halo names         - bold names with a white outline, directly on the map, no boxes
##     M3 callouts           - names outside the map, leader lines to each cluster
##     M4 numbers + legend   - numbered circles on the map, names in the legend
##     M6 small multiples    - one mini UMAP per state, one ROW per broad cell type
##   level "state" = cell states / sub-cell types (17); level "broad" = broad types (5)
##   Labels sit on each cluster's DENSEST point (not its median): scattered clusters
##   such as C7 otherwise get their label in empty space.
## OUTPUT: <Results>/<QC>_final_All_markers/Labelled_UMAP_options/<level>_M<n>_<style>_<tag>.tiff (300 dpi)
## Packages: data.table, Seurat, ggplot2, ggrepel (installed with Seurat), MASS (base R)
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 1: settings
## =============================================================================
library(data.table)   # tables
library(Seurat)       # Embeddings
library(ggplot2)      # plots
library(ggrepel)      # non-overlapping labels

## same path as script 14 (the UMAP being labelled)   <- change here
clustering_path <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Clustering_PCA_UMAP_Leiden_Markers/QC03_04_05_06_07_08_09/PC20_NN250_rep3/02_SENSITIVITY_batchawarePCA_Harmony/UMAP_by_Leiden_cluster.tiff"
leiden_resolution <- NULL     # NULL = main resolution from 12b (res 1)

pt_size  <- 0.12              # original 0.05 looks washed out once labels are added   <- change here
pt_alpha <- 0.5               # original 0.3                                         <- change here
seed     <- 42

## annotation: one row per Leiden cluster   <- change here
annot <- data.table(
  cluster = as.character(0:16),
  broad = c("Fibroblast", "Pericyte", "Endothelial", "Fibroblast", "VSMC", "Fibroblast", "VSMC",
            "Immune", "Pericyte", "VSMC", "Immune", "Fibroblast", "Pericyte", "VSMC", "VSMC",
            "VSMC", "Fibroblast"),
  state = c("Osteochondrogenic Fibroblast",          # C0
            "ECM-producing Pericyte",                # C1
            "Inflammatory Endothelial Cell",         # C2
            "Inflammatory Adventitial Fibroblast",   # C3
            "Inflammatory Synthetic VSMC",           # C4
            "Matrix-Remodeling Fibroblast",          # C5
            "Synthetic VSMC",                        # C6
            "T Cell",                                # C7
            "Activated Fibrogenic Pericyte",           # C8
            "Osteogenic / Calcifying VSMC",          # C9
            "Plasma Cell",                           # C10
            "Homeostatic Adventitial Fibroblast",      # C11
            "Contractile Pericyte",                  # C12
            "Fibromyocyte",                          # C13
            "Profibrotic Synthetic VSMC",            # C14
            "Contractile VSMC",                      # C15
            "Senescent Adventitial Fibroblast"))     # C16

## colours: one hue per broad type; states = shades of that hue (dark -> light)
broad_levels <- c("VSMC", "Fibroblast", "Pericyte", "Endothelial", "Immune")   # legend / key order
broad_cols <- c(VSMC = "#CB181D", Fibroblast = "#2171B5", Pericyte = "#F16913",
                Endothelial = "#238B45", Immune = "#6A51A3")
ramps <- list(VSMC        = c("#67000D", "#D7301F", "#F768A1"),           # dark red -> pink
              Fibroblast  = c("#08306B", "#2171B5", "#41B6C4"),           # navy -> teal
              Pericyte    = c("#7F2704", "#EC7014", "#FEC44F"),           # brown -> yellow
              Endothelial = c("#238B45"),
              Immune      = c("#3F007D", "#B07AA1"))


## =============================================================================
## CHUNK 2: paths from clustering_path (as script 14)
## =============================================================================
branch_dir <- if (dir.exists(clustering_path)) clustering_path else dirname(clustering_path)
branch_dir <- normalizePath(branch_dir, mustWork = TRUE)                    # stop if missing
dir_seu    <- dirname(branch_dir)                                           # .../PC20_NN250_rep3
run_tag    <- basename(dir_seu)                                             # "PC20_NN250_rep3"
qc_tag     <- basename(dirname(dir_seu))                                    # "QC03_04_05_06_07_08_09"
output_dir <- dirname(dirname(dirname(dir_seu)))                            # .../AAA
prefix <- if (grepl("^01_PRIMARY", basename(branch_dir))) "primary" else
  if (grepl("^02_SENSITIVITY", basename(branch_dir))) "harmony" else
    stop("branch folder must start with 01_PRIMARY or 02_SENSITIVITY")
if (is.null(leiden_resolution)) {
  params <- fread(file.path(dir_seu, paste0("parameters_used_", run_tag, ".csv")))
  leiden_resolution <- as.numeric(params[parameter == "Leiden_resolution", value])
}
seurat_out <- list.files(dir_seu, pattern = paste0("^seurat_.*_SeuratFunctions_", run_tag, "\\.RDS$"), full.names = TRUE)
stopifnot("expected exactly one 12b Seurat object" = length(seurat_out) == 1)

file_tag  <- paste0(prefix, "_", run_tag, "_res", leiden_resolution)       # "harmony_PC20_NN250_rep3_res1"
main_col  <- paste0(prefix, "_leiden_res_", leiden_resolution)             # Leiden column
umap_name <- paste0("umap_", prefix)                                        # UMAP of the original plot
dir_opt   <- file.path(output_dir, "Results", paste0(qc_tag, "_final_All_markers"), "Labelled_UMAP_options")
dir.create(dir_opt, showWarnings = FALSE, recursive = TRUE)
message("UMAP: ", umap_name, " | clusters: ", main_col, " | saving to ", dir_opt)


## =============================================================================
## CHUNK 3: UMAP coordinates + annotation per cell
## =============================================================================
obj <- readRDS(seurat_out)
stopifnot("UMAP not in object" = umap_name %in% Reductions(obj),
          "Leiden column not in object" = main_col %in% colnames(obj[[]]))
um <- Embeddings(obj, umap_name)                                            # same coordinates as the original
df <- data.table(UMAP_1 = um[, 1], UMAP_2 = um[, 2], cluster = as.character(obj[[main_col]][, 1]))
rm(obj, um); invisible(gc())

stopifnot("clusters without annotation" = all(unique(df$cluster) %in% annot$cluster),
          "annotation for clusters not in data" = all(annot$cluster %in% df$cluster))
annot[, broad := factor(broad, levels = broad_levels)]
annot <- annot[order(broad, as.numeric(cluster))]                          # VSMC first, then Fibroblast, ...
annot[, state_col := if (.N == 1) ramps[[as.character(broad[1])]][1]
      else colorRampPalette(ramps[[as.character(broad[1])]])(.N), by = broad]
annot[, state_lab := paste0("C", cluster, "  ", state)]                     # "C9  Osteogenic / Calcifying VSMC"
annot[, n_cells := df[, .N, by = cluster][match(annot$cluster, cluster), N]]

df[, broad := factor(as.character(annot$broad[match(cluster, annot$cluster)]), levels = broad_levels)]
df[, state_lab := factor(annot$state_lab[match(cluster, annot$cluster)], levels = annot$state_lab)]
set.seed(seed); df <- df[sample(.N)]                                        # random draw order


## =============================================================================
## CHUNK 4: helpers (label anchors, theme, axis arrows, callout layout)
## =============================================================================
xr <- range(df$UMAP_1); yr <- range(df$UMAP_2)                              # data extent
h_kde <- diff(xr) / 12                                                      # smoothing for density peak

## densest point of a set of cells (where the cluster visibly is)
peak <- function(x, y) {
  if (length(x) > 5000) { i <- sample(length(x), 5000); x <- x[i]; y <- y[i] }   # speed
  d <- MASS::kde2d(x, y, h = h_kde, n = 80, lims = c(xr, yr))
  ij <- which(d$z == max(d$z), arr.ind = TRUE)[1, ]
  c(d$x[ij[1]], d$y[ij[2]])
}
set.seed(seed)
cen <- df[, { p <- peak(UMAP_1, UMAP_2); .(x = p[1], y = p[2]) }, by = cluster]   # one anchor per cluster
cen <- merge(cen, annot, by = "cluster")
cen[, num := paste0("C", cluster)]

## broad level: one label per broad type, on its largest cluster
cen_b <- cen[order(-n_cells)][, .SD[1], by = broad][, .(broad, x, y)]
cen_b[, `:=`(text = as.character(broad), col = broad_cols[as.character(broad)])]
cen_b[, members := sapply(as.character(broad), function(b) paste0("C", annot[broad == b, cluster], collapse = ", "))]

## black or white text on a coloured fill
ink <- function(hex) { v <- col2rgb(hex) / 255; ifelse(0.299 * v[1, ] + 0.587 * v[2, ] + 0.114 * v[3, ] < 0.55, "white", "black") }
wrap <- function(s, w = 16) vapply(strwrap(s, w, simplify = FALSE), paste, "", collapse = "\n")

## clean UMAP look (as in atlas papers): no box, no ticks, short axis arrows
theme_umap <- theme_void(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(size = 11),
        legend.text = element_text(size = 11), legend.title = element_text(size = 12, face = "bold"),
        plot.margin = margin(10, 10, 10, 10))
arrows <- function(len = 0.13) {
  ox <- xr[1] - 0.02 * diff(xr); oy <- yr[1] - 0.02 * diff(yr)              # just outside bottom-left
  list(annotate("segment", x = ox, y = oy, xend = ox + len * diff(xr), yend = oy,
                arrow = arrow(length = unit(2, "mm"), type = "closed"), linewidth = 0.5),
       annotate("segment", x = ox, y = oy, xend = ox, yend = oy + len * diff(yr),
                arrow = arrow(length = unit(2, "mm"), type = "closed"), linewidth = 0.5),
       annotate("text", x = ox + len * diff(xr) / 2, y = oy - 0.03 * diff(yr), label = "UMAP 1", size = 3.5),
       annotate("text", x = ox - 0.03 * diff(xr), y = oy + len * diff(yr) / 2, label = "UMAP 2", size = 3.5, angle = 90))
}
pts <- function(colour_col) geom_point(aes(colour = .data[[colour_col]]), size = pt_size, alpha = pt_alpha, shape = 16)
big_keys <- guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5, alpha = 1)))
save_tiff <- function(p, f, w, h) {
  ggsave(file.path(dir_opt, f), p, width = w, height = h, dpi = 300, device = "tiff",
         compression = "lzw", bg = "white", limitsize = FALSE)
  message("saved ", f)
}

## M3: labels stacked down the left and right of the map, in the order of their clusters
spread_y <- function(y, gap, top) {                                         # keep order, at least 'gap' apart
  n <- length(y); o <- y
  if (n > 1) for (i in 2:n) o[i] <- max(o[i], o[i - 1] + gap)
  if (o[n] > top) o <- o - (o[n] - top)
  if (n > 1) for (i in (n - 1):1) o[i] <- min(o[i], o[i + 1] - gap)
  o
}
callout <- function(lab) {
  lab <- copy(lab)
  lab[, side := ifelse(x < mean(xr), "left", "right")]
  gap <- diff(yr) / max(8, lab[, .N, by = side][, max(N)] + 1)
  pad <- 0.06 * diff(xr)
  lab[, lx := ifelse(side == "left", xr[1] - pad, xr[2] + pad)]
  setorder(lab, side, y)
  lab[, ly := spread_y(y, gap, yr[2]), by = side]
  lab[, `:=`(hj = ifelse(side == "left", 1, 0), tx = lx + ifelse(side == "left", -1, 1) * 0.015 * diff(xr))]
  lab
}


## =============================================================================
## CHUNK 5: the four chosen labelling methods (one function per method)
## =============================================================================
## lv = "state" or "broad": what the cells are coloured by and what the labels say
level_setup <- function(lv) {
  ## lab: one row per cluster; grp = what one label stands for (a state = 1 cluster, a broad type = several)
  if (lv == "state") list(col = "state_lab", cols = setNames(annot$state_col, annot$state_lab),
                          title = "Cell states / sub-cell types", legend = "Cell state / sub-cell type",
                          lab = cen[, .(grp = cluster, cl = cluster, x, y, text = state, col = state_col)])
  else list(col = "broad", cols = broad_cols, title = "Broad cell types", legend = "Broad cell type",
            lab = cen[, .(grp = as.character(broad), cl = cluster, x, y, text = as.character(broad),
                          col = broad_cols[as.character(broad)])])   # broad name on every cluster
}
subtitle <- paste0("Leiden clusters (res ", leiden_resolution, "), ", format(nrow(df), big.mark = ","), " cells")

## M2: bold names with a white halo, straight on the map (no boxes)
m2 <- function(lv) {
  s <- level_setup(lv); lab <- s$lab[, t := wrap(text)]
  p <- ggplot(df, aes(UMAP_1, UMAP_2)) + pts(s$col) +
    geom_text_repel(data = lab, aes(x, y, label = t), colour = "black", size = 3.6, fontface = "bold",
                    lineheight = 0.85, bg.color = "white", bg.r = 0.18,      # white outline around letters
                    box.padding = 0.4, min.segment.length = 0.5, segment.colour = "grey30",
                    max.overlaps = Inf, seed = seed) +
    scale_colour_manual(values = s$cols, guide = "none") + arrows() +
    coord_fixed(clip = "off") + theme_umap + labs(title = s$title, subtitle = subtitle)
  save_tiff(p, paste0(lv, "_M2_halo_names_", file_tag, ".tiff"), 11, 10)
}

## M3: callouts - names in two columns beside the map, leader lines to the clusters
##     (broad: one callout per cluster, e.g. "VSMC (C15)")
m3 <- function(lv) {
  s <- level_setup(lv)
  if (lv == "broad") s$lab[, `:=`(grp = cl, text = paste0(text, "  (C", cl, ")"))]
  lab <- callout(s$lab[, .(x = mean(x), y = mean(y), text = text[1], col = col[1]), by = grp])
  seg <- merge(s$lab[, .(grp, ax = x, ay = y)], lab[, .(grp, lx, ly)], by = "grp")   # cluster -> its label
  w_mm <- max(nchar(lab$text)) * 2.1 + 12                                   # room for the longest name
  p <- ggplot(df, aes(UMAP_1, UMAP_2)) + pts(s$col) +
    geom_segment(data = seg, aes(x = ax, y = ay, xend = lx, yend = ly), colour = "grey35", linewidth = 0.3) +
    geom_point(data = seg, aes(ax, ay), size = 1.2, colour = "black") +    # anchor on each cluster
    geom_point(data = lab, aes(lx, ly), size = 3.5, colour = lab$col) +    # colour key at the label
    geom_text(data = lab, aes(tx, ly, label = text, hjust = hj), size = 3.7, fontface = "bold") +
    scale_colour_manual(values = s$cols, guide = "none") +               # no axis arrows: they collide with the left labels
    coord_fixed(clip = "off") + theme_umap +
    theme(plot.margin = margin(10, w_mm, 10, w_mm, unit = "mm")) +
    labs(title = s$title, subtitle = subtitle)
  save_tiff(p, paste0(lv, "_M3_callouts_", file_tag, ".tiff"), 10 + 2 * w_mm / 25.4, 10)
}

## M4: numbered circles on the map, names in the legend
m4 <- function(lv) {
  s <- level_setup(lv)
  fill_col <- if (lv == "state") cen$state_col else broad_cols[as.character(cen$broad)]
  leg_lab <- if (lv == "state") annot$state_lab else                        # broad: list its clusters
    setNames(paste0(broad_levels, "  (", cen_b$members[match(broad_levels, cen_b$broad)], ")"), broad_levels)
  p <- ggplot(df, aes(UMAP_1, UMAP_2)) + pts(s$col) +
    geom_point(data = cen, aes(x, y), shape = 21, size = 7.5, fill = fill_col, colour = "white", stroke = 1) +
    geom_text(data = cen, aes(x, y, label = cluster), colour = ink(fill_col), size = 3.3, fontface = "bold") +
    scale_colour_manual(values = s$cols, labels = leg_lab, name = s$legend) + big_keys +
    arrows() + coord_fixed(clip = "off") + theme_umap + labs(title = s$title, subtitle = subtitle)
  save_tiff(p, paste0(lv, "_M4_numbers_legend_", file_tag, ".tiff"), 14, 10)
}

## M6: small multiples - each group highlighted on grey cells
##     state: one ROW per broad type (VSMC, Fibroblast, ...), its states side by side
##     broad: one row of 5 panels
m6 <- function(lv) {
  s <- level_setup(lv)
  set.seed(seed)
  bg <- df[sample(.N, min(.N, 15000)), .(UMAP_1, UMAP_2)]                   # grey background cells
  ytop <- yr[2] + 0.22 * diff(yr)                                           # room for the panel title
  if (lv == "state") {
    a <- copy(annot)[, slot := seq_len(.N), by = broad]                     # position within its broad-type row
    a[, title := paste0("C", cluster, " ", wrap(state, 24))]
    hl  <- merge(df[, .(UMAP_1, UMAP_2, cluster, state_lab)], a[, .(cluster, broad, slot)], by = "cluster")
    bgs <- a[, .(UMAP_1 = bg$UMAP_1, UMAP_2 = bg$UMAP_2), by = .(broad, slot)]   # background only in used panels
    nc <- max(a$slot); nr <- nlevels(a$broad)
    p <- ggplot(hl, aes(UMAP_1, UMAP_2)) +
      geom_point(data = bgs, colour = "grey88", size = 0.02, shape = 16) +
      geom_point(aes(colour = state_lab), size = 0.06, alpha = 0.8, shape = 16) +
      geom_text(data = a, aes(x = mean(xr), y = ytop, label = title), vjust = 1,
                size = 3.4, fontface = "bold", lineheight = 0.9) +              # panel title inside the panel
      facet_grid(broad ~ slot, switch = "y") +                              # rows = broad types
      scale_colour_manual(values = s$cols, guide = "none") +
      scale_y_continuous(limits = c(yr[1], ytop)) +
      coord_fixed() + theme_void(base_size = 12) +
      theme(strip.text.x = element_blank(),                                 # no column headers
            strip.text.y.left = element_text(face = "bold", size = 14, angle = 90, margin = margin(0, 6, 0, 0)),
            plot.title = element_text(face = "bold", size = 16),
            plot.subtitle = element_text(margin = margin(0, 0, 8, 0)),
            panel.spacing.x = unit(3, "mm"), panel.spacing.y = unit(6, "mm")) +
      labs(title = s$title, subtitle = subtitle)
    save_tiff(p, paste0(lv, "_M6_small_multiples_", file_tag, ".tiff"), 0.6 + 3 * nc, 1 + 3.3 * nr)
  } else {
    bgs <- data.table(broad = factor(rep(broad_levels, each = nrow(bg)), levels = broad_levels),
                      UMAP_1 = rep(bg$UMAP_1, length(broad_levels)), UMAP_2 = rep(bg$UMAP_2, length(broad_levels)))
    p <- ggplot(df, aes(UMAP_1, UMAP_2)) +
      geom_point(data = bgs, colour = "grey88", size = 0.02, shape = 16) +
      geom_point(aes(colour = broad), size = 0.06, alpha = 0.8, shape = 16) +
      facet_wrap(~ broad, nrow = 1) +
      scale_colour_manual(values = s$cols, guide = "none") + coord_fixed() + theme_void(base_size = 12) +
      theme(strip.text = element_text(face = "bold", size = 14, margin = margin(4, 0, 4, 0)),
            plot.title = element_text(face = "bold", size = 16),
            plot.subtitle = element_text(margin = margin(0, 0, 8, 0)), panel.spacing = unit(4, "mm")) +
      labs(title = s$title, subtitle = subtitle)
    save_tiff(p, paste0(lv, "_M6_small_multiples_", file_tag, ".tiff"), 3.2 * length(broad_levels), 4.6)
  }
}
## =============================================================================
## CHUNK 6: make the chosen options (8 files: M2, M3, M4, M6 x 2 levels)
## =============================================================================
for (lv in c("state", "broad")) {
  m2(lv); m3(lv); m4(lv); m6(lv)
}
message("Done -> ", dir_opt)



















































## =============================================================================
## Cell typing leiden FindAllmarkeers AI cell typist etc QC2 (QC03_04_06_07_08_09)
## =============================================================================
## =============================================================================
## 23_UMAP_Leiden_labelled_cell_types.R Cell typing leiden QC2
## =============================================================================
## WHAT IT DOES
##   Labelled versions of UMAP_by_Leiden_cluster_<file_tag>.tiff (script 14 copy):
##   same Seurat object, same UMAP, same Leiden clusters (main resolution).
##   Nothing is recomputed.
##     A. UMAP_Leiden_broad_cell_types_<file_tag>.tiff
##        cells coloured by broad cell type; C0-C17 labels on the clusters
##     B. UMAP_Leiden_cell_states_<file_tag>.tiff
##        cells coloured by cell state / sub-cell type (shades of the broad-type
##        colour); C0-C17 on the map, full names in the legend
##     Leiden_cluster_annotation_<file_tag>.csv  (cluster, broad type, state, n cells)
##   Saved at 300 dpi in the same Results folder as the original UMAP.
## Labelling follows common practice in atlas papers: broad types are split across
## the map (e.g. Fibroblast = C0, C3, C5, C11, C16, C17), so one name per type on the
## map would sit in the wrong place -> cluster numbers on the map, names in a legend.
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 1: settings
## =============================================================================
library(data.table)   # tables
library(Seurat)       # Embeddings
library(ggplot2)      # plots
library(ggrepel)      # labels that do not overlap (installed with Seurat)

## same path as script 14 (the UMAP that was labelled)   <- change here
clustering_path <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Clustering_PCA_UMAP_Leiden_Markers/QC03_04_06_07_08_09/PC20_NN250_rep3/02_SENSITIVITY_batchawarePCA_Harmony/UMAP_by_Leiden_cluster.tiff"
leiden_resolution <- NULL     # NULL = main resolution from 12b (res 1)

label_style <- "number"       # plot B: "number" = C0-C17 on map (names in legend); "name" = full names on map   <- change here
pt_size  <- 0.05              # as original UMAP
pt_alpha <- 0.3               # as original UMAP
seed     <- 42

annot <- data.table(
  cluster = as.character(0:17),
  broad = c("Fibroblast",                 # C0
            "Pericyte",                   # C1
            "Endothelial-like Stromal",   # C2
            "Fibroblast",                 # C3
            "VSMC",               # C4
            "Fibroblast",                 # C5
            "VSMC",               # C6
            "Aberrant / Unresolved",      # C7
            "Pericyte / Mural",           # C8
            "VSMC",               # C9
            "Immune",                     # C10
            "Fibroblast",                 # C11
            "Pericyte",                   # C12
            "VSMC",               # C13
            "VSMC",               # C14
            "VSMC",               # C15
            "Fibroblast",                 # C16
            "Fibroblast"),                # C17
  state = c("Osteochondrogenic Fibroblast",          # C0
            "ECM-producing Pericyte",                # C1
            "Inflammatory Endothelial Cell",         # C2
            "Inflammatory Adventitial Fibroblast",   # C3
            "Inflammatory Synthetic VSMC",           # C4
            "Matrix-Remodeling Fibroblast",          # C5
            "Synthetic VSMC",                        # C6
            "Unresolved / Aberrant Cluster",         # C7
            "Activated Fibrogenic Pericyte",         # C8
            "Osteogenic / Calcifying VSMC",          # C9
            "Plasma Cell",                           # C10
            "Homeostatic Fibroblast",                # C11
            "Contractile Pericyte",                  # C12
            "Fibromyocyte",                          # C13
            "Matrix-Remodeling VSMC",                # C14
            "Contractile VSMC",                      # C15
            "Senescent Fibroblast",                  # C16
            "IFN-High Senescent Fibroblast"))        # C17

## colours: one hue per broad type; states = shades of that hue (dark -> light)
broad_levels <- c("VSMC", "Fibroblast", "Pericyte", "Pericyte / Mural",
                  "Endothelial-like Stromal", "Immune", "Aberrant / Unresolved")   # legend order
broad_cols <- c("VSMC" = "#CB181D", "Fibroblast" = "#2171B5", "Pericyte" = "#F16913",
                "Pericyte / Mural" = "#8C510A", "Endothelial-like Stromal" = "#238B45",
                "Immune" = "#6A51A3", "Aberrant / Unresolved" = "#7F7F7F")
ramps <- list("VSMC"             = c("#67000D", "#CB181D", "#FB6A4A"),
              "Fibroblast"               = c("#08306B", "#2171B5", "#6BAED6"),
              "Pericyte"                 = c("#7F2704", "#F16913", "#FDAE6B"),
              "Pericyte / Mural"         = c("#8C510A"),
              "Endothelial-like Stromal" = c("#238B45"),
              "Immune"                   = c("#6A51A3"),
              "Aberrant / Unresolved"    = c("#7F7F7F"))


## =============================================================================
## CHUNK 2: paths from clustering_path (as script 14)
## =============================================================================
branch_dir <- if (dir.exists(clustering_path)) clustering_path else dirname(clustering_path)
branch_dir <- normalizePath(branch_dir, mustWork = TRUE)                    # stop if missing
dir_seu    <- dirname(branch_dir)                                           # .../PC20_NN250_rep3
run_tag    <- basename(dir_seu)                                             # "PC20_NN250_rep3"
qc_tag     <- basename(dirname(dir_seu))                                    # "QC03_04_06_07_08_09"
output_dir <- dirname(dirname(dirname(dir_seu)))                            # .../AAA
prefix <- if (grepl("^01_PRIMARY", basename(branch_dir))) "primary" else
  if (grepl("^02_SENSITIVITY", basename(branch_dir))) "harmony" else
    stop("branch folder must start with 01_PRIMARY or 02_SENSITIVITY")
if (is.null(leiden_resolution)) {
  params <- fread(file.path(dir_seu, paste0("parameters_used_", run_tag, ".csv")))
  leiden_resolution <- as.numeric(params[parameter == "Leiden_resolution", value])
}
seurat_out <- list.files(dir_seu, pattern = paste0("^seurat_.*_SeuratFunctions_", run_tag, "\\.RDS$"), full.names = TRUE)
stopifnot("expected exactly one 12b Seurat object" = length(seurat_out) == 1)

dir_mk   <- file.path(output_dir, "Results", paste0(qc_tag, "_final_All_markers"))   # where the original UMAP copy is
file_tag <- paste0(prefix, "_", run_tag, "_res", leiden_resolution)        # "harmony_PC20_NN250_rep3_res1"
main_col <- paste0(prefix, "_leiden_res_", leiden_resolution)              # Leiden column
umap_name <- paste0("umap_", prefix)                                        # UMAP of the original plot
message("UMAP: ", umap_name, " | clusters: ", main_col, " | saving to ", dir_mk)


## =============================================================================
## CHUNK 3: UMAP coordinates + annotation per cell
## =============================================================================
obj <- readRDS(seurat_out)
stopifnot("UMAP not in object" = umap_name %in% Reductions(obj),
          "Leiden column not in object" = main_col %in% colnames(obj[[]]))
um <- Embeddings(obj, umap_name)                                            # same coordinates as the original plot
df <- data.table(UMAP_1 = um[, 1], UMAP_2 = um[, 2],
                 cluster = as.character(obj[[main_col]][, 1]))
stopifnot("clusters without annotation" = all(unique(df$cluster) %in% annot$cluster),
          "annotation for clusters not in data" = all(annot$cluster %in% df$cluster))
rm(obj); invisible(gc())

## state colours: shades within each broad type, ordered by type then cluster number
annot[, broad := factor(broad, levels = broad_levels)]
annot <- annot[order(broad, as.numeric(cluster))]                          # VSMC first, then Fibroblast, ...
annot[, state_col := if (.N == 1) ramps[[as.character(broad[1])]][1]
      else colorRampPalette(ramps[[as.character(broad[1])]])(.N), by = broad]
annot[, state_lab := paste0("C", cluster, "  ", state)]                     # legend text, e.g. "C9  Osteogenic / Calcifying VSMC"
annot[, n_cells := df[, .N, by = cluster][match(annot$cluster, cluster), N]]

df[, broad := factor(as.character(annot$broad[match(cluster, annot$cluster)]), levels = broad_levels)]
df[, state_lab := factor(annot$state_lab[match(cluster, annot$cluster)], levels = annot$state_lab)]
set.seed(seed); df <- df[sample(.N)]                                        # random draw order (no cluster on top)

## label positions: median UMAP position of each cluster (as Seurat LabelClusters)
cen <- df[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = cluster]
cen <- merge(cen, annot[, .(cluster, broad, state)], by = "cluster")

## legend text for broad types with cell numbers
n_broad <- df[, .N, by = broad]
broad_lab <- setNames(paste0(broad_levels, "  (n = ", format(n_broad$N[match(broad_levels, n_broad$broad)], big.mark = ","), ")"),
                      broad_levels)

fwrite(annot[, .(cluster = paste0("C", cluster), broad_cell_type = broad, cell_state = state, n_cells)],
       file.path(dir_mk, paste0("Leiden_cluster_annotation_", file_tag, ".csv")))

## shared look (clean, as in atlas papers)
theme_umap <- theme_classic(base_size = 14) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),         # UMAP units carry no meaning
        legend.text = element_text(size = 11), legend.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold"))
save_tiff <- function(p, f, w, h) ggsave(file.path(dir_mk, f), p, width = w, height = h, dpi = 300,
                                         device = "tiff", compression = "lzw", bg = "white")


## =============================================================================
## CHUNK 4: plot A - broad cell types
## =============================================================================
p_broad <- ggplot(df, aes(UMAP_1, UMAP_2, colour = broad)) +
  geom_point(size = pt_size, alpha = pt_alpha, shape = 16) +
  geom_label_repel(data = cen, aes(label = paste0("C", cluster)), colour = "black",
                   fill = alpha("white", 0.85), size = 3.8, fontface = "bold",
                   label.padding = 0.15, box.padding = 0.3, min.segment.length = 0,
                   max.overlaps = Inf, seed = seed, show.legend = FALSE) +   # cluster numbers on the map
  scale_colour_manual(values = broad_cols, labels = broad_lab, name = "Broad cell type") +
  guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1))) + # big, solid legend dots
  coord_fixed() + theme_umap +
  labs(title = "Broad cell types", x = "UMAP 1", y = "UMAP 2",
       subtitle = paste0("Leiden clusters (res ", leiden_resolution, ", ", nrow(annot), " clusters), ",
                         format(nrow(df), big.mark = ","), " cells"))
save_tiff(p_broad, paste0("UMAP_Leiden_broad_cell_types_", file_tag, ".tiff"), 12, 9)


## =============================================================================
## CHUNK 5: plot B - cell states / sub-cell types
## =============================================================================
cen[, lab := if (label_style == "name") state else paste0("C", cluster)]    # text on the map
p_state <- ggplot(df, aes(UMAP_1, UMAP_2, colour = state_lab)) +
  geom_point(size = pt_size, alpha = pt_alpha, shape = 16) +
  geom_label_repel(data = cen, aes(label = lab), colour = "black",
                   fill = alpha("white", 0.85), size = if (label_style == "name") 3.2 else 3.8,
                   fontface = "bold", label.padding = 0.15, box.padding = 0.4,
                   min.segment.length = 0, max.overlaps = Inf, seed = seed, show.legend = FALSE) +
  scale_colour_manual(values = setNames(annot$state_col, annot$state_lab),
                      name = "Cell state / sub-cell type") +
  guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5, alpha = 1))) +
  coord_fixed() + theme_umap +
  theme(legend.key.height = unit(0.55, "cm")) +                             # spacing between legend rows
  labs(title = "Cell states / sub-cell types", x = "UMAP 1", y = "UMAP 2",
       subtitle = "Colour shades grouped by broad type: VSMC (red), Fibroblast (blue), Pericyte (orange), Pericyte / Mural (brown),\nEndothelial-like Stromal (green), Immune (purple), Aberrant / Unresolved (grey)")
save_tiff(p_state, paste0("UMAP_Leiden_cell_states_", file_tag, ".tiff"), 15, 9)
message("Saved both plots + annotation table -> ", dir_mk)












## =============================================================================
## 23b_UMAP_state_label_options.R   Cell typing leiden QC2
## =============================================================================
## WHAT IT DOES
##   Five alternative ways to label the cell-state UMAP (plot B of script 23).
##   The problem with the original: labels are black-on-white "C0..C17" boxes, so
##   you cannot tell from the label which legend colour it belongs to, and the
##   shades within one broad type look alike.
##     Option 1  labels filled with the cluster's own colour        (link label <-> legend)
##     Option 2  one clearly different colour per cluster           (no look-alike shades)
##     Option 3  "C number + name" on the map, no legend            (no lookup needed)
##     Option 4  one small panel per cluster, rest of cells grey    (highlight plot)
##     Option 5  one panel per broad type, states coloured inside   (shades compared within type only)
##   Each option saves its own TIFF (300 dpi) in dir_mk, next to the script-23 plots.
## REQUIRES: script 23 CHUNKS 1-4 already run in this session (df, cen, annot,
##   theme_umap, save_tiff, file_tag, pt_size, pt_alpha, seed exist).
##   ggplot2 >= 3.5.0 (for I() inside aes; check with packageVersion("ggplot2")).
## Run chunk by chunk (Cmd+Enter). Run CHUNK 0 first, then any option you want.
## =============================================================================


## =============================================================================
## CHUNK 0: shared pieces for all options
## =============================================================================
## add legend text and the cluster's colour to the label table (cen = one row per cluster)
cen[, state_lab := annot$state_lab[match(cluster, annot$cluster)]]   # e.g. "C9  Osteogenic / Calcifying VSMC"
cen[, fill_col  := annot$state_col[match(cluster, annot$cluster)]]   # same shade as the cluster's cells

## text colour that stays readable on the fill: white on dark fills, black on light fills
text_on <- function(hex) {
  rgb <- col2rgb(hex)                                                # 3 x n matrix of red, green, blue (0-255)
  lum <- 0.299 * rgb[1, ] + 0.587 * rgb[2, ] + 0.114 * rgb[3, ]      # perceived brightness (standard luma weights)
  ifelse(lum > 150, "black", "white")
}
cen[, txt_col := text_on(fill_col)]

## label texts
cen[, lab_num  := paste0("C", cluster)]                              # "C9"
cen[, lab_full := paste0("C", cluster, "\n",                          # "C9" + name wrapped to ~16 characters per line
                         vapply(state, function(s) paste(strwrap(s, 16), collapse = "\n"), ""))]

## coloured label layer used by options 1, 2, 3, 5
## I() = use the colour values as they are (no colour scale), so the label fill
## does not clash with the colour scale of the points
lab_layer <- function(lab, fill, txt, size = 3.8) {
  geom_label_repel(data = cen,
                   aes(label = .data[[lab]], fill = I(.data[[fill]]), colour = I(.data[[txt]])),
                   size = size, fontface = "bold", label.padding = 0.15, box.padding = 0.4,
                   min.segment.length = 0, segment.colour = "black",             # line from label to cluster
                   max.overlaps = Inf, seed = seed, show.legend = FALSE)
}

## state colour scale + legend as in script 23 plot B
state_scale <- scale_colour_manual(values = setNames(annot$state_col, annot$state_lab),
                                   name = "Cell state / sub-cell type")
state_guide <- guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5, alpha = 1)))

## background cells for the panel plots (options 4, 5): this table has no
## state_lab / broad column, so ggplot draws it in every panel
bg <- df[, .(UMAP_1, UMAP_2)]


## =============================================================================
## OPTION 1: labels filled with the cluster colour
## =============================================================================
## same plot as script 23 plot B; only the label boxes change colour, so each
## "C" label matches its legend dot
p_opt1 <- ggplot(df, aes(UMAP_1, UMAP_2, colour = state_lab)) +
  geom_point(size = pt_size, alpha = pt_alpha, shape = 16) +
  lab_layer("lab_num", "fill_col", "txt_col") +
  state_scale + state_guide +
  coord_fixed() + theme_umap + theme(legend.key.height = unit(0.55, "cm")) +
  labs(title = "Cell states / sub-cell types", x = "UMAP 1", y = "UMAP 2")
save_tiff(p_opt1, paste0("UMAP_Leiden_cell_states_opt1_colouredlabels_", file_tag, ".tiff"), 15, 9)


## =============================================================================
## OPTION 2: one distinct colour per cluster (plus coloured labels)
## =============================================================================
## 18 maximally different colours (Trubetskoy "20 simple distinct colours",
## white/black/beige/apricot removed; grey kept because there are now 18 clusters)
distinct18 <- c("#E6194B", "#3CB44B", "#FFE119", "#4363D8", "#F58231", "#911EB4",
                "#42D4F4", "#F032E6", "#BFEF45", "#FABED4", "#469990", "#DCBEFF",
                "#9A6324", "#800000", "#AAFFC3", "#808000", "#000075", "#A9A9A9")
pal2 <- setNames(distinct18[seq_len(nrow(annot))], annot$state_lab)  # one colour per legend entry
cen[, fill2 := pal2[state_lab]]                                      # label fill = new cluster colour
cen[, txt2  := text_on(fill2)]

p_opt2 <- ggplot(df, aes(UMAP_1, UMAP_2, colour = state_lab)) +
  geom_point(size = pt_size, alpha = pt_alpha, shape = 16) +
  lab_layer("lab_num", "fill2", "txt2") +
  scale_colour_manual(values = pal2, name = "Cell state / sub-cell type") + state_guide +
  coord_fixed() + theme_umap + theme(legend.key.height = unit(0.55, "cm")) +
  labs(title = "Cell states / sub-cell types", x = "UMAP 1", y = "UMAP 2")
save_tiff(p_opt2, paste0("UMAP_Leiden_cell_states_opt2_distinctcolours_", file_tag, ".tiff"), 15, 9)


## =============================================================================
## OPTION 3: cluster number + name on the map, no legend
## =============================================================================
## the reader does not need to look anything up; works if the labels do not
## crowd too much (lower size if they do)
p_opt3 <- ggplot(df, aes(UMAP_1, UMAP_2, colour = state_lab)) +
  geom_point(size = pt_size, alpha = pt_alpha, shape = 16) +
  lab_layer("lab_full", "fill_col", "txt_col", size = 2.8) +
  state_scale +
  coord_fixed() + theme_umap + theme(legend.position = "none") +    # names are on the map
  labs(title = "Cell states / sub-cell types", x = "UMAP 1", y = "UMAP 2")
save_tiff(p_opt3, paste0("UMAP_Leiden_cell_states_opt3_namesonmap_", file_tag, ".tiff"), 13, 11)


## =============================================================================
## OPTION 4: one panel per cluster (highlight plot)
## =============================================================================
## every panel shows all cells in grey and one cluster in colour; the panel
## title is the cluster name. No ambiguity, often used as a supplementary figure
p_opt4 <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(data = bg, colour = "grey85", size = pt_size, alpha = pt_alpha, shape = 16) +  # all cells
  geom_point(aes(colour = state_lab), size = pt_size * 2, alpha = 0.6, shape = 16) +       # this cluster
  scale_colour_manual(values = setNames(annot$state_col, annot$state_lab), guide = "none") +
  facet_wrap(~ state_lab, ncol = 5, labeller = label_wrap_gen(width = 25)) +              # one panel per state
  coord_fixed() + theme_umap +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 10)) +
  labs(title = "Cell states / sub-cell types (one panel per Leiden cluster)", x = "UMAP 1", y = "UMAP 2")
save_tiff(p_opt4, paste0("UMAP_Leiden_cell_states_opt4_perclusterpanels_", file_tag, ".tiff"), 20, 16)


## =============================================================================
## OPTION 5: one panel per broad type
## =============================================================================
## similar shades only have to be told apart inside their own panel; labels
## appear only in the panel of their broad type (cen has a broad column)
p_opt5 <- ggplot(df, aes(UMAP_1, UMAP_2)) +
  geom_point(data = bg, colour = "grey85", size = pt_size, alpha = pt_alpha, shape = 16) +  # all cells
  geom_point(aes(colour = state_lab), size = pt_size * 2, alpha = 0.5, shape = 16) +       # cells of this type
  lab_layer("lab_num", "fill_col", "txt_col", size = 3.2) +
  state_scale + state_guide +
  facet_wrap(~ broad, ncol = 3) +                                                          # one panel per broad type
  coord_fixed() + theme_umap +
  theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 13),
        legend.key.height = unit(0.55, "cm")) +
  labs(title = "Cell states / sub-cell types by broad cell type", x = "UMAP 1", y = "UMAP 2")
save_tiff(p_opt5, paste0("UMAP_Leiden_cell_states_opt5_perbroadtype_", file_tag, ".tiff"), 20, 11)

message("Saved label options 1-5 -> ", dir_mk)























## =============================================================================
## 24_UMAP_label_options.R
## =============================================================================
## WHAT IT DOES
##   Four ways of labelling the SAME UMAP (umap_harmony, Leiden res 1, from 12b),
##   for both annotation levels:
##     M2 halo names         - bold names with a white outline, directly on the map, no boxes
##     M3 callouts           - names outside the map, leader lines to each cluster
##     M4 numbers + legend   - numbered circles on the map, names in the legend
##     M6 small multiples    - one mini UMAP per state, one ROW per broad cell type
##   level "state" = cell states / sub-cell types (18); level "broad" = broad types (7)
##   Labels sit on each cluster's DENSEST point (not its median): scattered clusters
##   such as C7 otherwise get their label in empty space.
## OUTPUT: <Results>/<QC>_final_All_markers/Labelled_UMAP_options/<level>_M<n>_<style>_<tag>.tiff (300 dpi)
## Packages: data.table, Seurat, ggplot2, ggrepel (installed with Seurat), MASS (base R)
## Run chunk by chunk (Cmd+Enter).
## =============================================================================


## =============================================================================
## CHUNK 1: settings
## =============================================================================
library(data.table)   # tables
library(Seurat)       # Embeddings
library(ggplot2)      # plots
library(ggrepel)      # non-overlapping labels

## same path as script 14 (the UMAP being labelled)   <- change here
clustering_path <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Clustering_PCA_UMAP_Leiden_Markers/QC03_04_06_07_08_09/PC20_NN250_rep3/02_SENSITIVITY_batchawarePCA_Harmony/UMAP_by_Leiden_cluster.tiff"
leiden_resolution <- NULL     # NULL = main resolution from 12b (res 1)

pt_size  <- 0.12              # original 0.05 looks washed out once labels are added   <- change here
pt_alpha <- 0.5               # original 0.3                                         <- change here
seed     <- 42

## annotation: one row per Leiden cluster   <- change here
annot <- data.table(
  cluster = as.character(0:17),
  broad = c("Fibroblast",                 # C0
            "Pericyte",                   # C1
            "Endothelial-like Stromal",   # C2
            "Fibroblast",                 # C3
            "VSMC",               # C4
            "Fibroblast",                 # C5
            "VSMC",               # C6
            "Aberrant / Unresolved",      # C7
            "Pericyte / Mural",           # C8
            "VSMC",               # C9
            "Immune",                     # C10
            "Fibroblast",                 # C11
            "Pericyte",                   # C12
            "VSMC",               # C13
            "VSMC",               # C14
            "VSMC",               # C15
            "Fibroblast",                 # C16
            "Fibroblast"),                # C17
  state = c("Osteochondrogenic Fibroblast",          # C0
            "ECM-producing Pericyte",                # C1
            "Inflammatory Endothelial Cell",         # C2
            "Inflammatory Adventitial Fibroblast",   # C3
            "Inflammatory Synthetic VSMC",           # C4
            "Matrix-Remodeling Fibroblast",          # C5
            "Synthetic VSMC",                        # C6
            "Unresolved / Aberrant Cluster",         # C7
            "Activated Fibrogenic Pericyte",         # C8
            "Osteogenic / Calcifying VSMC",          # C9
            "Plasma Cell",                           # C10
            "Homeostatic Fibroblast",                # C11
            "Contractile Pericyte",                  # C12
            "Fibromyocyte",                          # C13
            "Matrix-Remodeling VSMC",                # C14
            "Contractile VSMC",                      # C15
            "Senescent Fibroblast",                  # C16
            "IFN-High Senescent Fibroblast"))        # C17

## colours: one hue per broad type; states = shades of that hue (dark -> light)
broad_levels <- c("VSMC", "Fibroblast", "Pericyte", "Pericyte / Mural",
                  "Endothelial-like Stromal", "Immune", "Aberrant / Unresolved")   # legend / key order
broad_cols <- c("VSMC" = "#CB181D", "Fibroblast" = "#2171B5", "Pericyte" = "#F16913",
                "Pericyte / Mural" = "#8C510A", "Endothelial-like Stromal" = "#238B45",
                "Immune" = "#6A51A3", "Aberrant / Unresolved" = "#7F7F7F")
ramps <- list("VSMC"             = c("#67000D", "#D7301F", "#F768A1"),   # dark red -> pink
              "Fibroblast"               = c("#08306B", "#2171B5", "#41B6C4"),   # navy -> teal
              "Pericyte"                 = c("#7F2704", "#EC7014", "#FEC44F"),   # brown -> yellow
              "Pericyte / Mural"         = c("#8C510A"),
              "Endothelial-like Stromal" = c("#238B45"),
              "Immune"                   = c("#6A51A3"),
              "Aberrant / Unresolved"    = c("#7F7F7F"))


## =============================================================================
## CHUNK 2: paths from clustering_path (as script 14)
## =============================================================================
branch_dir <- if (dir.exists(clustering_path)) clustering_path else dirname(clustering_path)
branch_dir <- normalizePath(branch_dir, mustWork = TRUE)                    # stop if missing
dir_seu    <- dirname(branch_dir)                                           # .../PC20_NN250_rep3
run_tag    <- basename(dir_seu)                                             # "PC20_NN250_rep3"
qc_tag     <- basename(dirname(dir_seu))                                    # "QC03_04_06_07_08_09"
output_dir <- dirname(dirname(dirname(dir_seu)))                            # .../AAA
prefix <- if (grepl("^01_PRIMARY", basename(branch_dir))) "primary" else
  if (grepl("^02_SENSITIVITY", basename(branch_dir))) "harmony" else
    stop("branch folder must start with 01_PRIMARY or 02_SENSITIVITY")
if (is.null(leiden_resolution)) {
  params <- fread(file.path(dir_seu, paste0("parameters_used_", run_tag, ".csv")))
  leiden_resolution <- as.numeric(params[parameter == "Leiden_resolution", value])
}
seurat_out <- list.files(dir_seu, pattern = paste0("^seurat_.*_SeuratFunctions_", run_tag, "\\.RDS$"), full.names = TRUE)
stopifnot("expected exactly one 12b Seurat object" = length(seurat_out) == 1)

file_tag  <- paste0(prefix, "_", run_tag, "_res", leiden_resolution)       # "harmony_PC20_NN250_rep3_res1"
main_col  <- paste0(prefix, "_leiden_res_", leiden_resolution)             # Leiden column
umap_name <- paste0("umap_", prefix)                                        # UMAP of the original plot
dir_opt   <- file.path(output_dir, "Results", paste0(qc_tag, "_final_All_markers"), "Labelled_UMAP_options")
dir.create(dir_opt, showWarnings = FALSE, recursive = TRUE)
message("UMAP: ", umap_name, " | clusters: ", main_col, " | saving to ", dir_opt)


## =============================================================================
## CHUNK 3: UMAP coordinates + annotation per cell
## =============================================================================
obj <- readRDS(seurat_out)
stopifnot("UMAP not in object" = umap_name %in% Reductions(obj),
          "Leiden column not in object" = main_col %in% colnames(obj[[]]))
um <- Embeddings(obj, umap_name)                                            # same coordinates as the original
df <- data.table(UMAP_1 = um[, 1], UMAP_2 = um[, 2], cluster = as.character(obj[[main_col]][, 1]))
rm(obj, um); invisible(gc())

stopifnot("clusters without annotation" = all(unique(df$cluster) %in% annot$cluster),
          "annotation for clusters not in data" = all(annot$cluster %in% df$cluster))
annot[, broad := factor(broad, levels = broad_levels)]
annot <- annot[order(broad, as.numeric(cluster))]                          # VSMC first, then Fibroblast, ...
annot[, state_col := if (.N == 1) ramps[[as.character(broad[1])]][1]
      else colorRampPalette(ramps[[as.character(broad[1])]])(.N), by = broad]
annot[, state_lab := paste0("C", cluster, "  ", state)]                     # "C9  Osteogenic / Calcifying VSMC"
annot[, n_cells := df[, .N, by = cluster][match(annot$cluster, cluster), N]]

df[, broad := factor(as.character(annot$broad[match(cluster, annot$cluster)]), levels = broad_levels)]
df[, state_lab := factor(annot$state_lab[match(cluster, annot$cluster)], levels = annot$state_lab)]
set.seed(seed); df <- df[sample(.N)]                                        # random draw order


## =============================================================================
## CHUNK 4: helpers (label anchors, theme, axis arrows, callout layout)
## =============================================================================
xr <- range(df$UMAP_1); yr <- range(df$UMAP_2)                              # data extent
h_kde <- diff(xr) / 12                                                      # smoothing for density peak

## densest point of a set of cells (where the cluster visibly is)
peak <- function(x, y) {
  if (length(x) > 5000) { i <- sample(length(x), 5000); x <- x[i]; y <- y[i] }   # speed
  d <- MASS::kde2d(x, y, h = h_kde, n = 80, lims = c(xr, yr))
  ij <- which(d$z == max(d$z), arr.ind = TRUE)[1, ]
  c(d$x[ij[1]], d$y[ij[2]])
}
set.seed(seed)
cen <- df[, { p <- peak(UMAP_1, UMAP_2); .(x = p[1], y = p[2]) }, by = cluster]   # one anchor per cluster
cen <- merge(cen, annot, by = "cluster")
cen[, num := paste0("C", cluster)]

## broad level: one label per broad type, on its largest cluster
cen_b <- cen[order(-n_cells)][, .SD[1], by = broad][, .(broad, x, y)]
cen_b[, `:=`(text = as.character(broad), col = broad_cols[as.character(broad)])]
cen_b[, members := sapply(as.character(broad), function(b) paste0("C", annot[broad == b, cluster], collapse = ", "))]

## black or white text on a coloured fill
ink <- function(hex) { v <- col2rgb(hex) / 255; ifelse(0.299 * v[1, ] + 0.587 * v[2, ] + 0.114 * v[3, ] < 0.55, "white", "black") }
wrap <- function(s, w = 16) vapply(strwrap(s, w, simplify = FALSE), paste, "", collapse = "\n")

## clean UMAP look (as in atlas papers): no box, no ticks, short axis arrows
theme_umap <- theme_void(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(size = 11),
        legend.text = element_text(size = 11), legend.title = element_text(size = 12, face = "bold"),
        plot.margin = margin(10, 10, 10, 10))
arrows <- function(len = 0.13) {
  ox <- xr[1] - 0.02 * diff(xr); oy <- yr[1] - 0.02 * diff(yr)              # just outside bottom-left
  list(annotate("segment", x = ox, y = oy, xend = ox + len * diff(xr), yend = oy,
                arrow = arrow(length = unit(2, "mm"), type = "closed"), linewidth = 0.5),
       annotate("segment", x = ox, y = oy, xend = ox, yend = oy + len * diff(yr),
                arrow = arrow(length = unit(2, "mm"), type = "closed"), linewidth = 0.5),
       annotate("text", x = ox + len * diff(xr) / 2, y = oy - 0.03 * diff(yr), label = "UMAP 1", size = 3.5),
       annotate("text", x = ox - 0.03 * diff(xr), y = oy + len * diff(yr) / 2, label = "UMAP 2", size = 3.5, angle = 90))
}
pts <- function(colour_col) geom_point(aes(colour = .data[[colour_col]]), size = pt_size, alpha = pt_alpha, shape = 16)
big_keys <- guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5, alpha = 1)))
save_tiff <- function(p, f, w, h) {
  ggsave(file.path(dir_opt, f), p, width = w, height = h, dpi = 300, device = "tiff",
         compression = "lzw", bg = "white", limitsize = FALSE)
  message("saved ", f)
}

## M3: labels stacked down the left and right of the map, in the order of their clusters
spread_y <- function(y, gap, top) {                                         # keep order, at least 'gap' apart
  n <- length(y); o <- y
  if (n > 1) for (i in 2:n) o[i] <- max(o[i], o[i - 1] + gap)
  if (o[n] > top) o <- o - (o[n] - top)
  if (n > 1) for (i in (n - 1):1) o[i] <- min(o[i], o[i + 1] - gap)
  o
}
callout <- function(lab) {
  lab <- copy(lab)
  lab[, side := ifelse(x < mean(xr), "left", "right")]
  gap <- diff(yr) / max(8, lab[, .N, by = side][, max(N)] + 1)
  pad <- 0.06 * diff(xr)
  lab[, lx := ifelse(side == "left", xr[1] - pad, xr[2] + pad)]
  setorder(lab, side, y)
  lab[, ly := spread_y(y, gap, yr[2]), by = side]
  lab[, `:=`(hj = ifelse(side == "left", 1, 0), tx = lx + ifelse(side == "left", -1, 1) * 0.015 * diff(xr))]
  lab
}


## =============================================================================
## CHUNK 5: the four chosen labelling methods (one function per method)
## =============================================================================
## lv = "state" or "broad": what the cells are coloured by and what the labels say
level_setup <- function(lv) {
  ## lab: one row per cluster; grp = what one label stands for (a state = 1 cluster, a broad type = several)
  if (lv == "state") list(col = "state_lab", cols = setNames(annot$state_col, annot$state_lab),
                          title = "Cell states / sub-cell types", legend = "Cell state / sub-cell type",
                          lab = cen[, .(grp = cluster, cl = cluster, x, y, text = state, col = state_col)])
  else list(col = "broad", cols = broad_cols, title = "Broad cell types", legend = "Broad cell type",
            lab = cen[, .(grp = as.character(broad), cl = cluster, x, y, text = as.character(broad),
                          col = broad_cols[as.character(broad)])])   # broad name on every cluster
}
subtitle <- paste0("Leiden clusters (res ", leiden_resolution, "), ", format(nrow(df), big.mark = ","), " cells")

## M2: bold names with a white halo, straight on the map (no boxes)
m2 <- function(lv) {
  s <- level_setup(lv); lab <- s$lab[, t := wrap(text)]
  p <- ggplot(df, aes(UMAP_1, UMAP_2)) + pts(s$col) +
    geom_text_repel(data = lab, aes(x, y, label = t), colour = "black", size = 3.6, fontface = "bold",
                    lineheight = 0.85, bg.color = "white", bg.r = 0.18,      # white outline around letters
                    box.padding = 0.4, min.segment.length = 0.5, segment.colour = "grey30",
                    max.overlaps = Inf, seed = seed) +
    scale_colour_manual(values = s$cols, guide = "none") + arrows() +
    coord_fixed(clip = "off") + theme_umap + labs(title = s$title, subtitle = subtitle)
  save_tiff(p, paste0(lv, "_M2_halo_names_", file_tag, ".tiff"), 11, 10)
}

## M3: callouts - names in two columns beside the map, leader lines to the clusters
##     (broad: one callout per cluster, e.g. "VSMC (C15)")
m3 <- function(lv) {
  s <- level_setup(lv)
  if (lv == "broad") s$lab[, `:=`(grp = cl, text = paste0(text, "  (C", cl, ")"))]
  lab <- callout(s$lab[, .(x = mean(x), y = mean(y), text = text[1], col = col[1]), by = grp])
  seg <- merge(s$lab[, .(grp, ax = x, ay = y)], lab[, .(grp, lx, ly)], by = "grp")   # cluster -> its label
  w_mm <- max(nchar(lab$text)) * 2.1 + 12                                   # room for the longest name
  p <- ggplot(df, aes(UMAP_1, UMAP_2)) + pts(s$col) +
    geom_segment(data = seg, aes(x = ax, y = ay, xend = lx, yend = ly), colour = "grey35", linewidth = 0.3) +
    geom_point(data = seg, aes(ax, ay), size = 1.2, colour = "black") +    # anchor on each cluster
    geom_point(data = lab, aes(lx, ly), size = 3.5, colour = lab$col) +    # colour key at the label
    geom_text(data = lab, aes(tx, ly, label = text, hjust = hj), size = 3.7, fontface = "bold") +
    scale_colour_manual(values = s$cols, guide = "none") +               # no axis arrows: they collide with the left labels
    coord_fixed(clip = "off") + theme_umap +
    theme(plot.margin = margin(10, w_mm, 10, w_mm, unit = "mm")) +
    labs(title = s$title, subtitle = subtitle)
  save_tiff(p, paste0(lv, "_M3_callouts_", file_tag, ".tiff"), 10 + 2 * w_mm / 25.4, 10)
}

## M4: numbered circles on the map, names in the legend
m4 <- function(lv) {
  s <- level_setup(lv)
  fill_col <- if (lv == "state") cen$state_col else broad_cols[as.character(cen$broad)]
  leg_lab <- if (lv == "state") annot$state_lab else                        # broad: list its clusters
    setNames(paste0(broad_levels, "  (", cen_b$members[match(broad_levels, cen_b$broad)], ")"), broad_levels)
  p <- ggplot(df, aes(UMAP_1, UMAP_2)) + pts(s$col) +
    geom_point(data = cen, aes(x, y), shape = 21, size = 7.5, fill = fill_col, colour = "white", stroke = 1) +
    geom_text(data = cen, aes(x, y, label = cluster), colour = ink(fill_col), size = 3.3, fontface = "bold") +
    scale_colour_manual(values = s$cols, labels = leg_lab, name = s$legend) + big_keys +
    arrows() + coord_fixed(clip = "off") + theme_umap + labs(title = s$title, subtitle = subtitle)
  save_tiff(p, paste0(lv, "_M4_numbers_legend_", file_tag, ".tiff"), 14, 10)
}

## M6: small multiples - each group highlighted on grey cells
##     state: one ROW per broad type (VSMC, Fibroblast, ...), its states side by side
##     broad: one row of 7 panels
m6 <- function(lv) {
  s <- level_setup(lv)
  set.seed(seed)
  bg <- df[sample(.N, min(.N, 15000)), .(UMAP_1, UMAP_2)]                   # grey background cells
  ytop <- yr[2] + 0.22 * diff(yr)                                           # room for the panel title
  if (lv == "state") {
    a <- copy(annot)[, slot := seq_len(.N), by = broad]                     # position within its broad-type row
    a[, title := paste0("C", cluster, " ", wrap(state, 24))]
    hl  <- merge(df[, .(UMAP_1, UMAP_2, cluster, state_lab)], a[, .(cluster, broad, slot)], by = "cluster")
    bgs <- a[, .(UMAP_1 = bg$UMAP_1, UMAP_2 = bg$UMAP_2), by = .(broad, slot)]   # background only in used panels
    nc <- max(a$slot); nr <- nlevels(a$broad)
    p <- ggplot(hl, aes(UMAP_1, UMAP_2)) +
      geom_point(data = bgs, colour = "grey88", size = 0.02, shape = 16) +
      geom_point(aes(colour = state_lab), size = 0.06, alpha = 0.8, shape = 16) +
      geom_text(data = a, aes(x = mean(xr), y = ytop, label = title), vjust = 1,
                size = 3.4, fontface = "bold", lineheight = 0.9) +              # panel title inside the panel
      facet_grid(broad ~ slot, switch = "y") +                              # rows = broad types
      scale_colour_manual(values = s$cols, guide = "none") +
      scale_y_continuous(limits = c(yr[1], ytop)) +
      coord_fixed() + theme_void(base_size = 12) +
      theme(strip.text.x = element_blank(),                                 # no column headers
            strip.text.y.left = element_text(face = "bold", size = 14, angle = 90, margin = margin(0, 6, 0, 0)),
            plot.title = element_text(face = "bold", size = 16),
            plot.subtitle = element_text(margin = margin(0, 0, 8, 0)),
            panel.spacing.x = unit(3, "mm"), panel.spacing.y = unit(6, "mm")) +
      labs(title = s$title, subtitle = subtitle)
    save_tiff(p, paste0(lv, "_M6_small_multiples_", file_tag, ".tiff"), 0.6 + 3 * nc, 1 + 3.3 * nr)
  } else {
    bgs <- data.table(broad = factor(rep(broad_levels, each = nrow(bg)), levels = broad_levels),
                      UMAP_1 = rep(bg$UMAP_1, length(broad_levels)), UMAP_2 = rep(bg$UMAP_2, length(broad_levels)))
    p <- ggplot(df, aes(UMAP_1, UMAP_2)) +
      geom_point(data = bgs, colour = "grey88", size = 0.02, shape = 16) +
      geom_point(aes(colour = broad), size = 0.06, alpha = 0.8, shape = 16) +
      facet_wrap(~ broad, nrow = 1) +
      scale_colour_manual(values = s$cols, guide = "none") + coord_fixed() + theme_void(base_size = 12) +
      theme(strip.text = element_text(face = "bold", size = 14, margin = margin(4, 0, 4, 0)),
            plot.title = element_text(face = "bold", size = 16),
            plot.subtitle = element_text(margin = margin(0, 0, 8, 0)), panel.spacing = unit(4, "mm")) +
      labs(title = s$title, subtitle = subtitle)
    save_tiff(p, paste0(lv, "_M6_small_multiples_", file_tag, ".tiff"), 3.2 * length(broad_levels), 4.6)
  }
}
## =============================================================================
## CHUNK 6: make the chosen options (8 files: M2, M3, M4, M6 x 2 levels)
## =============================================================================
for (lv in c("state", "broad")) {
  m2(lv); m3(lv); m4(lv); m6(lv)
}
message("Done -> ", dir_opt)



































































## =============================================================================
## Cell typing InSituType (unsupervised) QC1 + QC2
## =============================================================================
## =============================================================================
## 26_InSituType_celltyping_UMAPs.R   Cell typing InSituType (unsupervised) QC1 + QC2
## =============================================================================
## WHAT IT DOES
##   The same cell-typing UMAPs as scripts 23 / 23b / 24 (Leiden), now for the
##   InSituType unsupervised clusters, drawn on the InSituType log-likelihood UMAP
##   (the 01b_UMAP_of_InSituType_logliks_NN30_md0.3_sp2_rep5_* layout).
##   Run for both QC choices (runs, CHUNK 1); results -> <InSituType folder>/celltyping/
##
##   STEP 0  rebuild the 01b UMAP exactly (same log-likelihoods, settings and seed as the
##           InSituType script CHUNK 7), save its coordinates, and draw a check plot
##           (00_...) that must look like the original 01b TIFF before anything is labelled.
##   Display fixes (layout itself unchanged):
##     - far-out specks (a few hundred cells) no longer shrink the map: view window =
##       0.2-99.8% of cells (+5% margin); hidden cells are counted in the subtitle
##     - small clusters drawn on top of the large cluster 0
##     - labels on each cluster's densest point (not median)
##   PLOTS (subtype = "Final Subtype", broad = "Broad Cell Class"):
##     23A / 23B         broad / subtype UMAP, C-number labels, names in legend   (as script 23)
##     23b opt1-opt5     coloured labels, distinct colours, names on map,
##                       one panel per cluster, one panel per broad class         (as script 23b)
##     M2 halo names, M3 callouts, M4 numbers + legend, M6 small multiples       (as script 24)
##     NEW M7 outlines   density outline around each cluster/class, name on it
##     NEW M8 two-level  cells by subtype, dashed outline + name per broad class
##     NEW M9 UMAP + bar numbered UMAP + bar chart of cell numbers as the key
##     NEW M10 zoom      overview + enlarged zoom of the crowded centre, names in the zoom
## Packages: data.table, ggplot2, ggrepel, patchwork (installed with Seurat), uwot, MASS
## Run chunk by chunk (Cmd+Enter); CHUNK 4 runs everything for both QC choices.
## =============================================================================


## =============================================================================
## CHUNK 1: settings
## =============================================================================
library(data.table)   # tables
library(ggplot2)      # plots
library(ggrepel)      # non-overlapping labels
library(patchwork)    # UMAP + bar chart side by side (M9)

ist_root <- "/Users/maryamuhabilausman/Desktop/Cosmx_analysis_M/AAA/Unsupervised_clustering_folder"

## one entry per QC choice: InSituType folder + annotation (one row per InSituType cluster)   <- change here
runs <- list(
  QC1 = list(ist_dir = file.path(ist_root, "QC03_04_05_06_07_08_09", "InSituType_PC20_NN250_rep3"),
             tag_suffix = "",                          # "_refined" if refineClusters was used
             annot = data.table(
               cluster = as.character(0:13),
               broad = c("Fibroblasts",                         # 0
                         "Fibroblasts",                         # 1
                         "Plasma cells",                        # 2
                         "Macrophages",                         # 3
                         "Cytotoxic T/NK populations",          # 4
                         "Dendritic cells",                     # 5
                         "Endothelial populations",             # 6
                         "Dendritic cells",                     # 7
                         "Cytotoxic T/NK populations",          # 8
                         "Monocyte/macrophage APCs",            # 9
                         "Cytotoxic T/NK populations",          # 10
                         "Interferon-activated immune cells",   # 11
                         "Regulatory/activated T cells",        # 12
                         "Endothelial populations"),            # 13
               state = c("Osteochondrogenic Fibroblast",            # 0
                         "Inflammatory Adventitial Fibroblast",     # 1
                         "Plasma Cell",                             # 2
                         "Lipid-associated inflammatory macrophage",# 3
                         "Activated cytotoxic T/NK population",     # 4
                         "Dendritic cells",                         # 5
                         "Activated vascular endothelial cell",     # 6
                         "MHC-IIhi antigen-presenting dendritic cell", # 7
                         "GNLY+ Mature NK Cell",                    # 8
                         "Inflammatory monocyte/macrophage",        # 9
                         "TRBV7-9+ Activated Cytotoxic T Cell",     # 10
                         "Interferon-responsive immune population", # 11
                         "Activated Immunoregulatory T Cell",       # 12
                         "Vasa Vasorum Endothelial cells"))),       # 13
  QC2 = list(ist_dir = file.path(ist_root, "QC03_04_06_07_08_09", "InSituType_PC20_NN250_rep3"),
             tag_suffix = "",
             annot = data.table(
               cluster = as.character(0:13),
               broad = c("Fibroblasts",                         # 0
                         "Endothelial cells",                   # 1
                         "Monocyte/macrophage",                 # 2
                         "Plasma cells",                        # 3
                         "Dendritic cells",                     # 4
                         "Cytotoxic T/NK cells",                # 5
                         "Cytotoxic T/NK cells",                # 6
                         "Regulatory/activated T cells",        # 7
                         "Plasma cells",                        # 8
                         "Monocyte/macrophage",                 # 9
                         "Cytotoxic T/NK cells",                # 10
                         "Interferon-activated immune cells",   # 11
                         "Cytotoxic T/NK cells",                # 12
                         "Endothelial cells"),                  # 13
               state = c("Osteochondrogenic Fibroblast",            # 0
                         "Activated Vascular Endothelial Cell",     # 1
                         "Lipid-Associated Inflammatory Macrophage",# 2
                         "Plasma Cell",                             # 3
                         "Dendritic cells",                         # 4
                         "Early Activated ITGB3+ NK Cell",          # 5
                         "Inflammatory CCL3+ Effector NK Cell",     # 6
                         "Activated TIGIT+ Regulatory-like T Cell", # 7
                         "IgG-High Plasma Cell",                    # 8
                         "FOLR2+ Resident Macrophage",              # 9
                         "GNLY+ Mature NK Cell",                    # 10
                         "Interferon-Responsive Immune Cell",       # 11
                         "TRBV7-9+ Activated Cytotoxic T Cell",     # 12
                         "Vasa Vasorum Endothelial Cell")))         # 13
)

## UMAP: must equal the InSituType script CHUNK 7 (01b) settings
umap_nn  <- 30; umap_md <- 0.3; umap_sp <- 2; umap_rep <- 5
seed     <- 42
reuse_coords <- TRUE          # TRUE = use saved coordinates after the first run (identical plots every time)

## display
trim_q   <- 0.002             # view window = 0.2-99.8% of cells; 0 = show every far-out speck   <- change here
pt_size  <- 0.12              # point size   <- change here
pt_alpha <- 0.5               # point opacity

## colours: one hue per broad class (legend order = this order); subtypes = shades of that hue
broad_order <- c("Fibroblasts", "Endothelial populations", "Endothelial cells",
                 "Macrophages", "Monocyte/macrophage", "Monocyte/macrophage APCs",
                 "Dendritic cells", "Plasma cells",
                 "Cytotoxic T/NK populations", "Cytotoxic T/NK cells",
                 "Regulatory/activated T cells", "Interferon-activated immune cells")
broad_cols <- c("Fibroblasts" = "#2171B5", "Endothelial populations" = "#238B45", "Endothelial cells" = "#238B45",
                "Macrophages" = "#D94801", "Monocyte/macrophage" = "#D94801", "Monocyte/macrophage APCs" = "#8C510A",
                "Dendritic cells" = "#CE1256", "Plasma cells" = "#6A51A3",
                "Cytotoxic T/NK populations" = "#CB181D", "Cytotoxic T/NK cells" = "#CB181D",
                "Regulatory/activated T cells" = "#E6AB02", "Interferon-activated immune cells" = "#17BECF")
ramps <- list("Fibroblasts"                = c("#08306B", "#6BAED6"),          # navy -> light blue
              "Endothelial populations"    = c("#00441B", "#74C476"),          # dark -> light green
              "Endothelial cells"          = c("#00441B", "#74C476"),
              "Macrophages"                = c("#7F2704", "#FD8D3C"),          # brown -> orange
              "Monocyte/macrophage"        = c("#7F2704", "#FD8D3C"),
              "Monocyte/macrophage APCs"   = c("#8C510A"),
              "Dendritic cells"            = c("#67001F", "#DF65B0"),          # wine -> pink
              "Plasma cells"               = c("#3F007D", "#9E9AC8"),          # dark -> light purple
              "Cytotoxic T/NK populations" = c("#67000D", "#CB181D", "#FC9272"),  # dark red -> salmon
              "Cytotoxic T/NK cells"       = c("#67000D", "#CB181D", "#FC9272"),
              "Regulatory/activated T cells"      = c("#E6AB02"),
              "Interferon-activated immune cells" = c("#17BECF"))
insitutype_pal <- c("#E6194B", "#3CB44B", "#4363D8", "#F58231", "#911EB4", "#42D4F4", "#F032E6",
                    "#9A6324", "#469990", "#BFEF45", "#800000", "#FFE119", "#000075", "#808000",
                    "#FF8FB1", "#A9A9A9", "#005300", "#009FFF", "#C77CFF", "#FF6F61")       # as InSituType script


## =============================================================================
## CHUNK 2: prepare one run (rebuild 01b UMAP, annotation, colours, label anchors)
## =============================================================================
prepare_run <- function(run_name) {
  r <- runs[[run_name]]
  file_tag <- paste0("insitutype_", sub("^InSituType_", "", basename(r$ist_dir)), r$tag_suffix)
  qc_tag   <- basename(dirname(r$ist_dir))                                  # "QC03_04_06_07_08_09"
  f_res    <- file.path(r$ist_dir, paste0("InSituType_result_", file_tag, ".RDS"))
  f_cells  <- file.path(r$ist_dir, paste0("InSituType_cells_", file_tag, ".csv"))
  f_map    <- file.path(r$ist_dir, paste0("InSituType_cluster_letter_to_number_", file_tag, ".csv"))
  out      <- file.path(r$ist_dir, "celltyping")                            # new results folder
  dir.create(out, showWarnings = FALSE)
  umap_tag <- paste0("NN", umap_nn, "_md", umap_md, "_sp", umap_sp, "_rep", umap_rep)
  f_coord  <- file.path(out, paste0("UMAP_InSituType_logliks_coords_", umap_tag, "_", file_tag, ".csv"))
  message("\n==== ", run_name, " (", qc_tag, ") -> ", out)
  
  ## --- STEP 0: 01b UMAP coordinates (rebuilt as InSituType script CHUNK 7, or reused) ---
  cells <- fread(f_cells, select = c("cell_id", "insitutype_cluster"))      # object cell order, as rownames(x)
  if (reuse_coords && file.exists(f_coord)) {
    um <- fread(f_coord)
    stopifnot("saved coordinates are for other cells - delete the coords CSV" = identical(um$cell_id, cells$cell_id))
    message("UMAP coordinates reused: ", basename(f_coord))
  } else {
    stopifnot("InSituType result / letter map not found" = file.exists(f_res), file.exists(f_map))
    res <- readRDS(f_res)                                                   # clust, prob, logliks, profiles
    map <- fread(f_map, colClasses = list(character = "insitutype_letter"))
    ll  <- res$logliks[cells$cell_id, map$insitutype_letter, drop = FALSE]  # cells x clusters, object order
    rm(res); invisible(gc())
    ll <- ll - apply(ll, 1, max)                                            # relative to each cell's best cluster
    ll[ll < -50] <- -50                                                     # cap very unlikely clusters
    set.seed(seed)                                                          # same seed as 01b
    m <- uwot::umap(ll, n_neighbors = umap_nn, min_dist = umap_md, spread = umap_sp,
                    repulsion_strength = umap_rep, metric = "euclidean")
    um <- data.table(cell_id = cells$cell_id, UMAP_1 = m[, 1], UMAP_2 = m[, 2])
    fwrite(um, f_coord)                                                     # reused next time
    message("UMAP rebuilt from log-likelihoods and saved: ", basename(f_coord))
  }
  df <- data.table(UMAP_1 = um$UMAP_1, UMAP_2 = um$UMAP_2, cluster = as.character(cells$insitutype_cluster))
  
  ## --- annotation + colours ---
  a <- copy(r$annot)[, `:=`(broad = trimws(broad), state = trimws(state))]  # drop stray spaces
  stopifnot("clusters without annotation" = all(unique(df$cluster) %in% a$cluster),
            "annotation for clusters not in data" = all(a$cluster %in% df$cluster),
            "broad class missing from broad_order / broad_cols" = all(a$broad %in% broad_order))
  lv <- intersect(broad_order, a$broad)                                     # broad classes of this run, legend order
  a[, broad := factor(broad, levels = lv)]
  a <- a[order(broad, as.numeric(cluster))]
  a[, state_col := if (.N == 1) broad_cols[[as.character(broad[1])]]
    else colorRampPalette(ramps[[as.character(broad[1])]])(.N), by = broad]
  a[, state_lab := paste0("C", cluster, "  ", state)]                       # "C3  Lipid-associated ..."
  a[, n_cells := df[, .N, by = cluster][match(a$cluster, cluster), N]]
  df[, broad := factor(as.character(a$broad[match(cluster, a$cluster)]), levels = lv)]
  df[, state_lab := factor(a$state_lab[match(cluster, a$cluster)], levels = a$state_lab)]
  df[, ncl := a$n_cells[match(cluster, a$cluster)]]
  set.seed(seed); df <- df[sample(.N)]                                      # random order within a cluster
  setorder(df, -ncl)                                                        # big clusters first = small ones on top
  
  ## --- view window: far-out specks would shrink the main map ---
  q <- function(v) { r <- if (trim_q > 0) quantile(v, c(trim_q, 1 - trim_q)) else range(v); r + c(-1, 1) * 0.05 * diff(r) }
  xr <- q(df$UMAP_1); yr <- q(df$UMAP_2)
  inside <- df$UMAP_1 >= xr[1] & df$UMAP_1 <= xr[2] & df$UMAP_2 >= yr[1] & df$UMAP_2 <= yr[2]
  dfv <- df[inside]
  
  ## --- label anchor = densest point of each cluster ---
  h <- diff(xr) / 12
  peak <- function(x, y) {
    if (length(x) > 5000) { i <- sample(length(x), 5000); x <- x[i]; y <- y[i] }
    d <- MASS::kde2d(x, y, h = h, n = 100, lims = c(xr, yr))
    ij <- which(d$z == max(d$z), arr.ind = TRUE)[1, ]
    c(d$x[ij[1]], d$y[ij[2]])
  }
  set.seed(seed)
  cen <- dfv[, { p <- peak(UMAP_1, UMAP_2); .(x = p[1], y = p[2]) }, by = cluster]
  cen <- merge(cen, a, by = "cluster")
  cen[, members := sapply(as.character(broad), function(b) paste0("C", a[broad == b, cluster], collapse = ", "))]
  
  fwrite(a[, .(cluster = paste0("C", cluster), broad_cell_class = broad, final_subtype = state, n_cells, colour = state_col)],
         file.path(out, paste0("InSituType_cluster_annotation_", file_tag, ".csv")))
  n_hidden <- sum(!inside)
  list(name = run_name, qc_tag = qc_tag, file_tag = file_tag, out = out, umap_tag = umap_tag,
       df = df, dfv = dfv, annot = a, cen = cen, xr = xr, yr = yr, h = h, lv = lv,
       bcols = broad_cols[lv], n_hidden = n_hidden,
       subtitle = paste0("InSituType clusters (", qc_tag, "), ", format(nrow(df), big.mark = ","), " cells",
                         if (n_hidden > 0) paste0("; ", n_hidden, " far-out cells outside the view") else ""))
}


## =============================================================================
## CHUNK 3: plotting helpers and plot functions
## =============================================================================
save_tiff <- function(ctx, p, f, w, h) {
  ggsave(file.path(ctx$out, paste0(f, "_", ctx$file_tag, ".tiff")), p, width = w, height = h, dpi = 300,
         device = "tiff", compression = "lzw", bg = "white", limitsize = FALSE)
  message("saved ", f)
}
ink  <- function(hex) { v <- col2rgb(hex) / 255; ifelse(0.299 * v[1, ] + 0.587 * v[2, ] + 0.114 * v[3, ] < 0.55, "white", "black") }
wrap <- function(s, w = 16) vapply(strwrap(s, w, simplify = FALSE), paste, "", collapse = "\n")
big_keys <- guides(colour = guide_legend(ncol = 1, override.aes = list(size = 5, alpha = 1)))
pts <- function(col) geom_point(aes(colour = .data[[col]]), size = pt_size, alpha = pt_alpha, shape = 16)

## script-23 look (axes, no ticks) and script-24 look (no box, small axis arrows)
theme_23 <- theme_classic(base_size = 14) +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        legend.text = element_text(size = 11), legend.title = element_text(size = 12, face = "bold"),
        plot.title = element_text(face = "bold"))
theme_umap <- theme_void(base_size = 14) +
  theme(plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(size = 11),
        legend.text = element_text(size = 11), legend.title = element_text(size = 12, face = "bold"),
        plot.margin = margin(10, 10, 10, 10))
arrows <- function(ctx, len = 0.13) {
  xr <- ctx$xr; yr <- ctx$yr; ox <- xr[1]; oy <- yr[1]
  list(annotate("segment", x = ox, y = oy, xend = ox + len * diff(xr), yend = oy,
                arrow = arrow(length = unit(2, "mm"), type = "closed"), linewidth = 0.5),
       annotate("segment", x = ox, y = oy, xend = ox, yend = oy + len * diff(yr),
                arrow = arrow(length = unit(2, "mm"), type = "closed"), linewidth = 0.5),
       annotate("text", x = ox + len * diff(xr) / 2, y = oy - 0.03 * diff(yr), label = "UMAP 1", size = 3.5),
       annotate("text", x = ox - 0.03 * diff(xr), y = oy + len * diff(yr) / 2, label = "UMAP 2", size = 3.5, angle = 90))
}

## level = "subtype" (Final Subtype) or "broad" (Broad Cell Class)
level_setup <- function(ctx, lv) {
  cen <- ctx$cen
  if (lv == "subtype") list(col = "state_lab", cols = setNames(ctx$annot$state_col, ctx$annot$state_lab),
                            title = "Cell subtypes", legend = "Final subtype",
                            lab = cen[, .(grp = cluster, cl = cluster, x, y, text = state, col = state_col)])
  else list(col = "broad", cols = ctx$bcols, title = "Broad cell classes", legend = "Broad cell class",
            lab = cen[, .(grp = as.character(broad), cl = cluster, x, y, text = as.character(broad),
                          col = broad_cols[as.character(broad)])])
}

## density outline(s) of a set of cells: contour at 'frac' of its peak density
outline <- function(ctx, x, y, frac, hmul = 0.6) {
  if (length(x) > 5000) { i <- sample(length(x), 5000); x <- x[i]; y <- y[i] }
  d  <- MASS::kde2d(x, y, h = ctx$h * hmul, n = 150, lims = c(ctx$xr, ctx$yr))
  cl <- contourLines(d$x, d$y, d$z, levels = frac * max(d$z))
  rbindlist(lapply(seq_along(cl), function(k) data.table(piece = k, x = cl[[k]]$x, y = cl[[k]]$y)))
}

## ---------- STEP 0: check plot = original 01b look, all cells, InSituType colours ----------
p00_check <- function(ctx) {
  df <- ctx$df; lev <- as.character(sort(as.numeric(unique(df$cluster))))
  med <- df[, .(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2)), by = cluster]
  set.seed(seed)
  p <- ggplot(df[sample(.N)], aes(UMAP_1, UMAP_2, colour = factor(cluster, levels = lev))) +
    geom_point(size = 0.05, alpha = 0.3, shape = 16) +
    geom_label(data = med, aes(label = cluster, fill = factor(cluster, levels = lev)), colour = "black", size = 4, show.legend = FALSE) +
    scale_colour_manual(values = insitutype_pal[seq_along(lev)], aesthetics = c("colour", "fill")) +
    coord_fixed() + theme_bw(base_size = 13) +
    guides(colour = guide_legend(override.aes = list(size = 4, alpha = 1), ncol = 2)) +
    labs(title = "CHECK: rebuilt 01b UMAP - compare with the original 01b TIFF",
         subtitle = paste0(ctx$qc_tag, "; UMAP of InSituType log-likelihoods: n_neighbors ", umap_nn, ", min_dist ", umap_md,
                           ", spread ", umap_sp, ", repulsion ", umap_rep, ", seed ", seed),
         x = "UMAPIST_1", y = "UMAPIST_2", colour = "InSituType\ncluster")
  save_tiff(ctx, p, "00_CHECK_rebuilt_01b_UMAP", 11, 9)
}

## ---------- script 23 A / B ----------
p23 <- function(ctx) {
  df <- ctx$dfv; cen <- ctx$cen; a <- ctx$annot
  n_b <- df[, .N, by = broad]
  blab <- setNames(paste0(ctx$lv, "  (n = ", format(n_b$N[match(ctx$lv, n_b$broad)], big.mark = ","), ")"), ctx$lv)
  lab_c <- function(size) geom_label_repel(data = cen, aes(x, y, label = paste0("C", cluster)), colour = "black",
                                           fill = alpha("white", 0.85), size = size, fontface = "bold",
                                           label.padding = 0.15, box.padding = 0.3, min.segment.length = 0,
                                           max.overlaps = Inf, seed = seed, show.legend = FALSE)
  pA <- ggplot(df, aes(UMAP_1, UMAP_2, colour = broad)) + geom_point(size = pt_size, alpha = pt_alpha, shape = 16) +
    lab_c(3.8) + scale_colour_manual(values = ctx$bcols, labels = blab, name = "Broad cell class") + big_keys +
    coord_fixed() + theme_23 + labs(title = "Broad cell classes", subtitle = ctx$subtitle, x = "UMAP 1", y = "UMAP 2")
  save_tiff(ctx, pA, "23A_UMAP_broad_cell_classes", 13, 9)
  pB <- ggplot(df, aes(UMAP_1, UMAP_2, colour = state_lab)) + geom_point(size = pt_size, alpha = pt_alpha, shape = 16) +
    lab_c(3.8) + scale_colour_manual(values = setNames(a$state_col, a$state_lab), name = "Final subtype") + big_keys +
    coord_fixed() + theme_23 + theme(legend.key.height = unit(0.55, "cm")) +
    labs(title = "Cell subtypes", subtitle = ctx$subtitle, x = "UMAP 1", y = "UMAP 2")
  save_tiff(ctx, pB, "23B_UMAP_cell_subtypes", 15, 9)
}

## ---------- script 23b options 1-5 (labels filled with the cluster colour) ----------
## fill via scale_fill_identity; text colour = 2 layers (white on dark, black on light);
## works on every ggplot2 version (23b needed >= 3.5 for I())
col_labels <- function(dat, label, fill, size) {
  lapply(c("white", "black"), function(tc) {
    d <- dat[ink(dat[[fill]]) == tc]
    if (nrow(d) == 0) return(NULL)
    geom_label_repel(data = d, aes(x, y, label = .data[[label]], fill = .data[[fill]]), colour = tc,
                     size = size, fontface = "bold", label.padding = 0.15, box.padding = 0.4,
                     min.segment.length = 0, segment.colour = "black", max.overlaps = Inf, seed = seed,
                     show.legend = FALSE, inherit.aes = FALSE)
  })
}
p23b <- function(ctx) {
  df <- ctx$dfv; a <- ctx$annot; cen <- copy(ctx$cen)
  cen[, lab_num := paste0("C", cluster)]
  cen[, lab_full := paste0("C", cluster, "\n", wrap(state, 16))]
  sc <- scale_colour_manual(values = setNames(a$state_col, a$state_lab), name = "Final subtype")
  base <- function(p) p + coord_fixed() + theme_23 + labs(x = "UMAP 1", y = "UMAP 2", subtitle = ctx$subtitle)
  ## opt1: labels in cluster colour
  p1 <- base(ggplot(df, aes(UMAP_1, UMAP_2, colour = state_lab)) + geom_point(size = pt_size, alpha = pt_alpha, shape = 16) +
               col_labels(cen, "lab_num", "state_col", 3.8) + scale_fill_identity() + sc + big_keys +
               theme(legend.key.height = unit(0.55, "cm")) + labs(title = "Cell subtypes"))
  save_tiff(ctx, p1, "23b_opt1_coloured_labels", 15, 9)
  ## opt2: one clearly different colour per cluster
  pal2 <- setNames(insitutype_pal[seq_len(nrow(a))], a$state_lab)
  cen[, fill2 := pal2[state_lab]]
  p2 <- base(ggplot(df, aes(UMAP_1, UMAP_2, colour = state_lab)) + geom_point(size = pt_size, alpha = pt_alpha, shape = 16) +
               col_labels(cen, "lab_num", "fill2", 3.8) + scale_fill_identity() +
               scale_colour_manual(values = pal2, name = "Final subtype") + big_keys +
               theme(legend.key.height = unit(0.55, "cm")) + labs(title = "Cell subtypes"))
  save_tiff(ctx, p2, "23b_opt2_distinct_colours", 15, 9)
  ## opt3: number + name on the map, no legend
  p3 <- base(ggplot(df, aes(UMAP_1, UMAP_2, colour = state_lab)) + geom_point(size = pt_size, alpha = pt_alpha, shape = 16) +
               col_labels(cen, "lab_full", "state_col", 2.8) + scale_fill_identity() + sc +
               theme(legend.position = "none") + labs(title = "Cell subtypes"))
  save_tiff(ctx, p3, "23b_opt3_names_on_map", 13, 11)
  ## opt4: one panel per cluster
  bg <- df[, .(UMAP_1, UMAP_2)]
  p4 <- base(ggplot(df, aes(UMAP_1, UMAP_2)) +
               geom_point(data = bg, colour = "grey85", size = pt_size, alpha = pt_alpha, shape = 16) +
               geom_point(aes(colour = state_lab), size = pt_size * 2, alpha = 0.6, shape = 16) +
               scale_colour_manual(values = setNames(a$state_col, a$state_lab), guide = "none") +
               facet_wrap(~ state_lab, ncol = 5, labeller = label_wrap_gen(width = 25)) +
               theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 10)) +
               labs(title = "Cell subtypes (one panel per InSituType cluster)"))
  save_tiff(ctx, p4, "23b_opt4_per_cluster_panels", 20, 4 * ceiling(nrow(a) / 5) + 1)
  ## opt5: one panel per broad class, labels only in their own panel
  p5 <- base(ggplot(df, aes(UMAP_1, UMAP_2)) +
               geom_point(data = bg, colour = "grey85", size = pt_size, alpha = pt_alpha, shape = 16) +
               geom_point(aes(colour = state_lab), size = pt_size * 2, alpha = 0.5, shape = 16) +
               col_labels(cen, "lab_num", "state_col", 3.2) + scale_fill_identity() + sc + big_keys +
               facet_wrap(~ broad, ncol = 3) +
               theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = 13),
                     legend.key.height = unit(0.55, "cm")) +
               labs(title = "Cell subtypes by broad cell class"))
  save_tiff(ctx, p5, "23b_opt5_per_broad_class", 20, 5 * ceiling(length(ctx$lv) / 3) + 1)
}

## ---------- script 24 M2 / M3 / M4 / M6 ----------
m2 <- function(ctx, lv) {                                                   # bold names, white halo, no boxes
  s <- level_setup(ctx, lv); lab <- s$lab[, t := wrap(text)]
  p <- ggplot(ctx$dfv, aes(UMAP_1, UMAP_2)) + pts(s$col) +
    geom_text_repel(data = lab, aes(x, y, label = t), colour = "black", size = 3.6, fontface = "bold",
                    lineheight = 0.85, bg.color = "white", bg.r = 0.18, box.padding = 0.4,
                    min.segment.length = 0.5, segment.colour = "grey30", max.overlaps = Inf, seed = seed) +
    scale_colour_manual(values = s$cols, guide = "none") + arrows(ctx) +
    coord_fixed(clip = "off") + theme_umap + labs(title = s$title, subtitle = ctx$subtitle)
  save_tiff(ctx, p, paste0(lv, "_M2_halo_names"), 11, 10)
}

spread_y <- function(y, gap, top) {                                         # keep order, >= gap apart
  n <- length(y); o <- y
  if (n > 1) for (i in 2:n) o[i] <- max(o[i], o[i - 1] + gap)
  if (o[n] > top) o <- o - (o[n] - top)
  if (n > 1) for (i in (n - 1):1) o[i] <- min(o[i], o[i + 1] - gap)
  o
}
callout <- function(ctx, lab) {
  xr <- ctx$xr; yr <- ctx$yr; lab <- copy(lab)
  lab[, side := ifelse(x < mean(xr), "left", "right")]
  gap <- diff(yr) / max(8, lab[, .N, by = side][, max(N)] + 1)
  pad <- 0.06 * diff(xr)
  lab[, lx := ifelse(side == "left", xr[1] - pad, xr[2] + pad)]
  setorder(lab, side, y)
  lab[, ly := spread_y(y, gap, yr[2]), by = side]
  lab[, `:=`(hj = ifelse(side == "left", 1, 0), tx = lx + ifelse(side == "left", -1, 1) * 0.015 * diff(xr))]
  lab
}
m3 <- function(ctx, lv) {                                                   # names beside the map, leader lines
  s <- level_setup(ctx, lv)
  if (lv == "broad") s$lab[, `:=`(grp = cl, text = paste0(text, "  (C", cl, ")"))]
  lab <- callout(ctx, s$lab[, .(x = mean(x), y = mean(y), text = text[1], col = col[1]), by = grp])
  seg <- merge(s$lab[, .(grp, ax = x, ay = y)], lab[, .(grp, lx, ly)], by = "grp")
  w_mm <- max(nchar(lab$text)) * 2.1 + 12
  p <- ggplot(ctx$dfv, aes(UMAP_1, UMAP_2)) + pts(s$col) +
    geom_segment(data = seg, aes(x = ax, y = ay, xend = lx, yend = ly), colour = "grey35", linewidth = 0.3) +
    geom_point(data = seg, aes(ax, ay), size = 1.2, colour = "black") +
    geom_point(data = lab, aes(lx, ly), size = 3.5, colour = lab$col) +
    geom_text(data = lab, aes(tx, ly, label = text, hjust = hj), size = 3.7, fontface = "bold") +
    scale_colour_manual(values = s$cols, guide = "none") +
    coord_fixed(clip = "off") + theme_umap + theme(plot.margin = margin(10, w_mm, 10, w_mm, unit = "mm")) +
    labs(title = s$title, subtitle = ctx$subtitle)
  save_tiff(ctx, p, paste0(lv, "_M3_callouts"), 10 + 2 * w_mm / 25.4, 10)
}
badges <- function(cen, fill) list(
  geom_point(data = cen, aes(x, y), shape = 21, size = 7.5, fill = fill, colour = "white", stroke = 1),
  geom_text(data = cen, aes(x, y, label = cluster), colour = ink(fill), size = 3.3, fontface = "bold"))
m4 <- function(ctx, lv) {                                                   # numbered circles + legend
  s <- level_setup(ctx, lv); cen <- ctx$cen
  fill <- if (lv == "subtype") cen$state_col else broad_cols[as.character(cen$broad)]
  leg <- if (lv == "subtype") ctx$annot$state_lab else
    setNames(paste0(ctx$lv, "  (", cen$members[match(ctx$lv, cen$broad)], ")"), ctx$lv)
  p <- ggplot(ctx$dfv, aes(UMAP_1, UMAP_2)) + pts(s$col) + badges(cen, fill) +
    scale_colour_manual(values = s$cols, labels = leg, name = s$legend) + big_keys +
    arrows(ctx) + coord_fixed(clip = "off") + theme_umap + labs(title = s$title, subtitle = ctx$subtitle)
  save_tiff(ctx, p, paste0(lv, "_M4_numbers_legend"), 14, 10)
}
m6 <- function(ctx, lv) {                                                   # small multiples, one row per broad class
  s <- level_setup(ctx, lv); df <- ctx$dfv; xr <- ctx$xr; yr <- ctx$yr
  set.seed(seed); bg <- df[sample(.N, min(.N, 15000)), .(UMAP_1, UMAP_2)]
  ytop <- yr[2] + 0.22 * diff(yr)
  if (lv == "subtype") {
    a <- copy(ctx$annot)[, slot := seq_len(.N), by = broad]
    a[, title := paste0("C", cluster, " ", wrap(state, 24))]
    hl  <- merge(df[, .(UMAP_1, UMAP_2, cluster, state_lab)], a[, .(cluster, broad, slot)], by = "cluster")
    bgs <- a[, .(UMAP_1 = bg$UMAP_1, UMAP_2 = bg$UMAP_2), by = .(broad, slot)]
    a[, broad_w := factor(wrap(as.character(broad), 16), levels = wrap(ctx$lv, 16))]   # wrapped row names
    hl[, broad_w := a$broad_w[match(cluster, a$cluster)]]; bgs[, broad_w := a$broad_w[match(paste(broad, slot), paste(a$broad, a$slot))]]
    nc <- max(a$slot); nr <- length(ctx$lv)
    p <- ggplot(hl, aes(UMAP_1, UMAP_2)) +
      geom_point(data = bgs, colour = "grey88", size = 0.02, shape = 16) +
      geom_point(aes(colour = state_lab), size = 0.06, alpha = 0.8, shape = 16) +
      geom_text(data = a, aes(x = mean(xr), y = ytop, label = title), vjust = 1, size = 3.4,
                fontface = "bold", lineheight = 0.9) +
      facet_grid(broad_w ~ slot, switch = "y") +
      scale_colour_manual(values = s$cols, guide = "none") + scale_y_continuous(limits = c(yr[1], ytop)) +
      coord_fixed() + theme_void(base_size = 12) +
      theme(strip.text.x = element_blank(),
            strip.text.y.left = element_text(face = "bold", size = 12, angle = 90, margin = margin(0, 6, 0, 0)),
            plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(margin = margin(0, 0, 8, 0)),
            panel.spacing.x = unit(3, "mm"), panel.spacing.y = unit(6, "mm")) +
      labs(title = s$title, subtitle = ctx$subtitle)
    save_tiff(ctx, p, paste0(lv, "_M6_small_multiples"), 0.8 + 3 * nc, 1 + 3.3 * nr)
  } else {
    bgs <- data.table(broad = factor(rep(ctx$lv, each = nrow(bg)), levels = ctx$lv),
                      UMAP_1 = rep(bg$UMAP_1, length(ctx$lv)), UMAP_2 = rep(bg$UMAP_2, length(ctx$lv)))
    nc <- min(5, length(ctx$lv)); nr <- ceiling(length(ctx$lv) / nc)
    p <- ggplot(df, aes(UMAP_1, UMAP_2)) +
      geom_point(data = bgs, colour = "grey88", size = 0.02, shape = 16) +
      geom_point(aes(colour = broad), size = 0.06, alpha = 0.8, shape = 16) +
      facet_wrap(~ broad, ncol = nc, labeller = label_wrap_gen(width = 22)) +
      scale_colour_manual(values = s$cols, guide = "none") + coord_fixed() + theme_void(base_size = 12) +
      theme(strip.text = element_text(face = "bold", size = 13, margin = margin(4, 0, 4, 0)),
            plot.title = element_text(face = "bold", size = 16), plot.subtitle = element_text(margin = margin(0, 0, 8, 0)),
            panel.spacing = unit(4, "mm")) +
      labs(title = s$title, subtitle = ctx$subtitle)
    save_tiff(ctx, p, paste0(lv, "_M6_small_multiples"), 3.2 * nc, 1 + 3.6 * nr)
  }
}

## ---------- NEW M7: density outline around each cluster / broad class, name on it ----------
m7 <- function(ctx, lv) {
  s <- level_setup(ctx, lv); df <- ctx$dfv
  grp_col <- if (lv == "subtype") "cluster" else "broad"
  set.seed(seed)
  ol <- df[, outline(ctx, UMAP_1, UMAP_2, frac = if (lv == "subtype") 0.25 else 0.15,
                     hmul = if (lv == "subtype") 0.5 else 0.6), by = grp_col]          # contour pieces per group
  setnames(ol, grp_col, "g"); ol[, g := as.character(g)]
  ol[, col := if (lv == "subtype") ctx$annot$state_col[match(g, ctx$annot$cluster)] else broad_cols[g]]
  ol[, id := paste(g, piece)]
  lab <- if (lv == "subtype") ctx$cen[, .(x, y, t = wrap(state), col = state_col)] else {
    ctx$cen[order(-n_cells)][, .SD[1], by = broad][, .(x, y, t = wrap(as.character(broad)), col = broad_cols[as.character(broad)])] }
  p <- ggplot(df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(colour = .data[[s$col]]), size = pt_size, alpha = 0.25, shape = 16) +   # faint cells
    geom_path(data = ol, aes(x, y, group = id), colour = ol$col, linewidth = 0.7) +        # outlines
    geom_text_repel(data = lab, aes(x, y, label = t), colour = lab$col, size = 3.6, fontface = "bold",
                    lineheight = 0.85, bg.color = "white", bg.r = 0.2, box.padding = 0.5,
                    min.segment.length = 0.3, segment.colour = "grey40", max.overlaps = Inf, seed = seed) +
    scale_colour_manual(values = s$cols, guide = "none") + arrows(ctx) +
    coord_fixed(clip = "off") + theme_umap +
    labs(title = paste0(s$title, " (density outlines)"), subtitle = ctx$subtitle)
  save_tiff(ctx, p, paste0(lv, "_M7_outlines"), 11, 10)
}

## ---------- NEW M8: two levels in one plot (subtype colours + broad-class outlines) ----------
m8 <- function(ctx) {
  df <- ctx$dfv; a <- ctx$annot; cen <- ctx$cen
  set.seed(seed)
  ol <- df[, outline(ctx, UMAP_1, UMAP_2, frac = 0.12, hmul = 0.7), by = broad]
  ol[, `:=`(col = broad_cols[as.character(broad)], id = paste(broad, piece))]
  top <- ol[, .SD[which.max(y)], by = broad][, .(x, y, t = wrap(as.character(broad), 20), col)]   # name at top of outline
  p <- ggplot(df, aes(UMAP_1, UMAP_2)) + pts("state_lab") +
    geom_path(data = ol, aes(x, y, group = id), colour = ol$col, linewidth = 0.8, linetype = "22") +
    badges(cen, cen$state_col) +
    geom_text_repel(data = top, aes(x, y, label = t), colour = top$col, size = 4, fontface = "bold",
                    lineheight = 0.85, bg.color = "white", bg.r = 0.2, nudge_y = 0.04 * diff(ctx$yr),
                    min.segment.length = 0.3, segment.colour = "grey40", max.overlaps = Inf, seed = seed) +
    scale_colour_manual(values = setNames(a$state_col, a$state_lab), name = "Final subtype") + big_keys +
    arrows(ctx) + coord_fixed(clip = "off") + theme_umap +
    labs(title = "Broad cell classes (outlines) and cell subtypes (colours, numbers)", subtitle = ctx$subtitle)
  save_tiff(ctx, p, "two_level_M8_outlines_plus_subtypes", 15, 10)
}

## ---------- NEW M9: numbered UMAP + bar chart of cell numbers as the key ----------
m9 <- function(ctx, lv) {
  s <- level_setup(ctx, lv); a <- ctx$annot; cen <- ctx$cen
  fill <- if (lv == "subtype") cen$state_col else broad_cols[as.character(cen$broad)]
  p_map <- ggplot(ctx$dfv, aes(UMAP_1, UMAP_2)) + pts(s$col) + badges(cen, fill) +
    scale_colour_manual(values = s$cols, guide = "none") + arrows(ctx) + coord_fixed(clip = "off") + theme_umap
  if (lv == "subtype") {
    bars <- a[, .(lab = paste0("C", cluster, "  ", state), n = n_cells, col = state_col,
                  broad_w = factor(wrap(as.character(broad), 18), levels = wrap(ctx$lv, 18)))]
    bars[, lab := factor(lab, levels = rev(lab))]
    p_bar <- ggplot(bars, aes(n, lab, fill = lab)) + geom_col(width = 0.75) +
      facet_grid(broad_w ~ ., scales = "free_y", space = "free_y", switch = "y")
  } else {
    bars <- a[, .(n = sum(n_cells), members = paste0("C", cluster, collapse = ", ")), by = broad]
    bars[, `:=`(lab = factor(paste0(broad, "  (", members, ")"), levels = rev(paste0(broad, "  (", members, ")"))),
                col = broad_cols[as.character(broad)])]
    p_bar <- ggplot(bars, aes(n, lab, fill = lab)) + geom_col(width = 0.75)
  }
  p_bar <- p_bar + geom_text(aes(label = format(n, big.mark = ",")), hjust = -0.1, size = 3.2) +
    scale_fill_manual(values = setNames(bars$col, as.character(bars$lab)), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.3)), labels = scales::comma) +
    theme_minimal(base_size = 12) +
    theme(panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
          strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 1),
          axis.text.y = element_text(size = 10.5, colour = "black")) +
    labs(x = "Cells", y = NULL)
  p <- (p_map | p_bar) + plot_layout(widths = c(2.2, if (lv == "subtype") 1.9 else 1.3)) +
    plot_annotation(title = s$title, subtitle = ctx$subtitle,
                    theme = theme(plot.title = element_text(face = "bold", size = 16)))
  save_tiff(ctx, p, paste0(lv, "_M9_umap_plus_bar_key"), if (lv == "subtype") 18 else 16, 10)
}


## ---------- NEW M10: overview + zoom of the crowded centre (cluster 0 takes most of the map) ----------
peak_in <- function(x, y, xr, yr, h) {                                      # densest point inside a window
  if (length(x) > 5000) { i <- sample(length(x), 5000); x <- x[i]; y <- y[i] }
  d <- MASS::kde2d(x, y, h = h, n = 100, lims = c(xr, yr))
  ij <- which(d$z == max(d$z), arr.ind = TRUE)[1, ]
  c(d$x[ij[1]], d$y[ij[2]])
}
m10 <- function(ctx, lv, min_cells_zoom = 30) {
  s <- level_setup(ctx, lv); df <- ctx$dfv; a <- ctx$annot
  big <- a$cluster[which.max(a$n_cells)]                                    # the large cluster (C0)
  core <- df[cluster != big]                                                # all other clusters
  zx <- quantile(core$UMAP_1, c(0.03, 0.97)); zx <- zx + c(-1, 1) * 0.08 * diff(zx)   # zoom window
  zy <- quantile(core$UMAP_2, c(0.03, 0.97)); zy <- zy + c(-1, 1) * 0.08 * diff(zy)
  zdf <- df[UMAP_1 >= zx[1] & UMAP_1 <= zx[2] & UMAP_2 >= zy[1] & UMAP_2 <= zy[2]]
  grp <- if (lv == "subtype") "cluster" else "broad"
  set.seed(seed)
  zl <- zdf[cluster != big, if (.N >= min_cells_zoom) { p <- peak_in(UMAP_1, UMAP_2, zx, zy, diff(zx) / 14); .(x = p[1], y = p[2]) },
            by = grp]                                                       # one label per group inside the zoom
  zl[, t := if (lv == "subtype") wrap(a$state[match(cluster, a$cluster)]) else wrap(as.character(broad))]
  bl <- ctx$cen[cluster == big, .(x, y, t = if (lv == "subtype") wrap(state) else wrap(as.character(broad)))]
  halo <- function(d, size) geom_text_repel(data = d, aes(x, y, label = t), colour = "black", size = size,
                                            fontface = "bold", lineheight = 0.85, bg.color = "white", bg.r = 0.18,
                                            box.padding = 0.45, min.segment.length = 0.3, segment.colour = "grey30",
                                            max.overlaps = Inf, seed = seed)
  p_over <- ggplot(df, aes(UMAP_1, UMAP_2)) + pts(s$col) + halo(bl, 3.6) +
    annotate("rect", xmin = zx[1], xmax = zx[2], ymin = zy[1], ymax = zy[2], fill = NA, colour = "black",
             linewidth = 0.6, linetype = "22") +
    annotate("text", x = zx[2], y = zy[2], label = "zoom", hjust = 1, vjust = -0.5, size = 3.5, fontface = "italic") +
    scale_colour_manual(values = s$cols, guide = "none") + arrows(ctx) + coord_fixed(clip = "off") + theme_umap +
    labs(subtitle = "Overview")
  p_zoom <- ggplot(zdf, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(colour = .data[[s$col]]), size = pt_size * 1.8, alpha = pt_alpha, shape = 16) +   # bigger points when zoomed
    halo(zl, 3.4) + scale_colour_manual(values = s$cols, guide = "none") +
    coord_fixed(xlim = zx, ylim = zy, expand = FALSE) + theme_umap +
    theme(panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6, linetype = "22")) +
    labs(subtitle = "Zoom")
  p <- (p_over | p_zoom) + plot_layout(widths = c(1, 1.5)) +
    plot_annotation(title = s$title, subtitle = ctx$subtitle,
                    theme = theme(plot.title = element_text(face = "bold", size = 16)))
  save_tiff(ctx, p, paste0(lv, "_M10_overview_plus_zoom"), 18, 10)
}


## =============================================================================
## CHUNK 4: run everything for both QC choices
## =============================================================================
for (rn in names(runs)) {
  ctx <- prepare_run(rn)
  p00_check(ctx)                                    # STEP 0: compare with the original 01b TIFF first
  p23(ctx)                                          # 23A, 23B
  p23b(ctx)                                         # 23b options 1-5
  for (lv in c("subtype", "broad")) {
    m2(ctx, lv); m3(ctx, lv); m4(ctx, lv); m6(ctx, lv)   # script 24 methods
    m7(ctx, lv); m9(ctx, lv); m10(ctx, lv)               # new methods
  }
  m8(ctx)                                           # new two-level plot
  message("Done ", rn, " -> ", ctx$out)
}