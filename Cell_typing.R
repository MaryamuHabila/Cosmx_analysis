## =============================================================================
## cluster markers - FindAllMarkers (main resolution)
## =============================================================================
## Each Leiden cluster vs all other cells, genes HIGHER in the cluster, Wilcoxon
## (Seurat vignette: FindAllMarkers(only.pos = TRUE)). Clusters NOT recomputed.
## FindAllMarkers reads the RNA "data" layer, which still holds raw counts
## -> NormalizeData first (log1p(counts / total x 10,000)); this changes
## neither clusters nor UMAP.
## SIGNIFICANT = p_val_adj < 0.05 (Bonferroni). RANKING = pct.1 high -> low, ties by avg_log2FC.
## Faster with presto: remotes::install_github("immunogenomics/presto")
marker_padj_max <- 0.05   # significance cut-off on adjusted p-value
marker_top_n    <- 20     # genes per cluster sheet

obj <- NormalizeData(obj, assay = "RNA", normalization.method = "LogNormalize", scale.factor = 10000)

for (b in branches) {
  main_col <- paste0(b$prefix, "_leiden_res_", leiden_resolution)
  dir_mk   <- file.path(b$dir, paste0("FindAllMarkers_res", leiden_resolution))
  dir.create(dir_mk, showWarnings = FALSE)
  file_tag <- paste0(b$prefix, "_", plot_tag)                           # e.g. "primary_PC20_NN250_res1"
  
  Idents(obj) <- main_col                                               # levels already 0, 1, 2, ... (by size)
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