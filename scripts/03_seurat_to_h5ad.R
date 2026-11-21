#!/usr/bin/env Rscript
# 將 Tonsil Atlas 的 Seurat .rds 轉成 .h5ad（供 Python 端分析）
# 用法: Rscript 03_seurat_to_h5ad.R <input.rds> <output.h5ad> [assay]
suppressPackageStartupMessages({
  library(Seurat); library(SingleCellExperiment); library(zellkonverter); library(Matrix)
})
args <- commandArgs(trailingOnly = TRUE)
inp <- args[1]; outp <- args[2]
assay_use <- if (length(args) >= 3) args[3] else NA

cat("[read ]", inp, "\n"); flush.console()
obj <- readRDS(inp)
cat("[class]", class(obj)[1], "\n")
cat("[assays]", paste(Assays(obj), collapse = ", "), "\n")
cat("[dims ]", paste(dim(obj), collapse = " x "), "\n")
cat("[meta ]", paste(head(colnames(obj@meta.data), 60), collapse = ", "), "\n")
cat("[reduc]", paste(Reductions(obj), collapse = ", "), "\n")
flush.console()

if (is.na(assay_use)) {
  assay_use <- if ("RNA" %in% Assays(obj)) "RNA" else Assays(obj)[1]
}
DefaultAssay(obj) <- assay_use
cat("[use  ] assay =", assay_use, "\n")

# 只保留 RNA assay，丟掉 ChromatinAssay（zellkonverter 無法處理）
keep <- assay_use
for (a in setdiff(Assays(obj), keep)) {
  obj[[a]] <- NULL
}

# metadata 中的 list/複雜欄位會讓轉換失敗 -> 一律轉字串
md <- obj@meta.data
for (cn in colnames(md)) {
  if (!is.numeric(md[[cn]]) && !is.character(md[[cn]]) && !is.factor(md[[cn]]) && !is.logical(md[[cn]])) {
    md[[cn]] <- as.character(md[[cn]])
  }
  if (is.factor(md[[cn]])) md[[cn]] <- as.character(md[[cn]])
}
obj@meta.data <- md

sce <- as.SingleCellExperiment(obj, assay = assay_use)
# 降維座標保留
for (r in Reductions(obj)) {
  try({ reducedDim(sce, toupper(r)) <- Embeddings(obj, reduction = r) }, silent = TRUE)
}
cat("[write]", outp, "\n"); flush.console()
writeH5AD(sce, outp, X_name = "counts", compression = "gzip")
cat("[done ]", outp, " size=", file.info(outp)$size / 1e9, "GB\n")
