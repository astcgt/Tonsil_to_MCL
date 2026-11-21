#!/usr/bin/env Rscript
# 匯出 Tonsil Atlas Multiome 的 RNA counts、ATAC peak counts、metadata
# 不依賴 zellkonverter（避免 ChromatinAssay 轉換問題），直接寫 MatrixMarket。
suppressPackageStartupMessages({ library(Matrix) })

inp  <- "/raid/jupyter-alz.ee09/Tonsil_to_MCL/data/extracted/multiome/20230911_tonsil_atlas_multiome_seurat_obj.rds"
outd <- "/raid/jupyter-alz.ee09/Tonsil_to_MCL/data/processed/multiome_export"
dir.create(outd, recursive = TRUE, showWarnings = FALSE)

cat("[R] 讀 Seurat 物件 ...\n"); flush.console()
obj <- readRDS(inp)
md  <- obj@meta.data

# ---------- metadata ----------
md$cell_barcode <- rownames(md)
for (cn in colnames(md)) if (is.factor(md[[cn]])) md[[cn]] <- as.character(md[[cn]])
write.csv(md, file.path(outd, "multiome_metadata.csv"), row.names = FALSE)
cat("[R] metadata:", nrow(md), "x", ncol(md), "\n"); flush.console()

# ---------- RNA counts ----------
rna <- obj@assays$RNA
cnt <- if (!is.null(rna@counts) && length(rna@counts) > 0) rna@counts else rna@data
cat("[R] RNA counts:", paste(dim(cnt), collapse = " x "), "\n"); flush.console()
writeMM(cnt, file.path(outd, "rna_counts.mtx"))
writeLines(rownames(cnt), file.path(outd, "rna_genes.txt"))
writeLines(colnames(cnt), file.path(outd, "rna_barcodes.txt"))
cat("[R] RNA 已寫出\n"); flush.console()

# ---------- ATAC peak counts ----------
ok <- tryCatch({
  atac <- obj@assays$ATAC
  ac <- if (!is.null(atac@counts) && length(atac@counts) > 0) atac@counts else atac@data
  cat("[R] ATAC counts:", paste(dim(ac), collapse = " x "), "\n"); flush.console()
  writeMM(ac, file.path(outd, "atac_counts.mtx"))
  writeLines(rownames(ac), file.path(outd, "atac_peaks.txt"))
  writeLines(colnames(ac), file.path(outd, "atac_barcodes.txt"))
  TRUE
}, error = function(e) { cat("[R] ATAC 匯出失敗:", conditionMessage(e), "\n"); FALSE })

# ---------- 降維座標 ----------
for (r in names(obj@reductions)) {
  emb <- tryCatch(obj@reductions[[r]]@cell.embeddings, error = function(e) NULL)
  if (!is.null(emb)) {
    write.csv(as.data.frame(emb), file.path(outd, paste0("reduction_", r, ".csv")))
  }
}
cat("[R] 完成。ATAC ok =", ok, "\n")
