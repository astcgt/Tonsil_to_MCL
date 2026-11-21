#!/usr/bin/env Rscript
# 盤點 Seurat 物件結構（assay / metadata 欄位 / 註解層級 / 降維）
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
args <- commandArgs(trailingOnly = TRUE)
inp <- args[1]
cat("=== reading", inp, "===\n"); flush.console()
obj <- readRDS(inp)
cat("class :", class(obj)[1], "\n")
cat("dims  :", paste(dim(obj), collapse = " x "), "\n")
cat("assays:", paste(Assays(obj), collapse = ", "), "\n")
for (a in Assays(obj)) {
  cat("  -", a, ":", paste(dim(obj[[a]]), collapse = " x "),
      "| class:", class(obj[[a]])[1], "\n")
}
cat("reductions:", paste(Reductions(obj), collapse = ", "), "\n")
cat("n meta cols:", ncol(obj@meta.data), "\n")
cat("meta cols:\n")
print(colnames(obj@meta.data))
cat("\n=== 每個可能的註解欄位的層級數 ===\n")
for (cn in colnames(obj@meta.data)) {
  v <- obj@meta.data[[cn]]
  if (is.character(v) || is.factor(v)) {
    n <- length(unique(as.character(v)))
    if (n <= 200) cat(sprintf("%-40s n_unique=%d\n", cn, n))
  }
}
cat("\n=== head(meta) ===\n")
print(head(obj@meta.data, 3))
