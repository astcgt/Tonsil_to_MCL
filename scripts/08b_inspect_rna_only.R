#!/usr/bin/env Rscript
# 只看 metadata，不碰 ChromatinAssay（避免需要 Signac）
suppressPackageStartupMessages({ library(Matrix) })
args <- commandArgs(trailingOnly = TRUE)
obj <- readRDS(args[1])
md <- obj@meta.data
cat("n cells:", nrow(md), " n meta cols:", ncol(md), "\n")
cat("\n=== meta 欄位與層級數 ===\n")
for (cn in colnames(md)) {
  v <- md[[cn]]
  if (is.character(v) || is.factor(v)) {
    u <- unique(as.character(v)); n <- length(u)
    if (n <= 250) {
      cat(sprintf("%-38s n=%-4d  e.g. %s\n", cn, n,
                  paste(head(u, 5), collapse = " | ")))
    } else cat(sprintf("%-38s n=%-4d  (too many)\n", cn, n))
  } else {
    cat(sprintf("%-38s numeric  range=[%.3g, %.3g]\n", cn,
                suppressWarnings(min(v, na.rm = TRUE)), suppressWarnings(max(v, na.rm = TRUE))))
  }
}
cat("\n=== reductions ===\n")
print(names(obj@reductions))
cat("\n=== 若有 annotation 欄位，列出完整層級 ===\n")
for (cn in grep("annotation|cell_type|celltype|level", colnames(md), ignore.case = TRUE, value = TRUE)) {
  u <- sort(unique(as.character(md[[cn]])))
  cat("\n--", cn, "(n=", length(u), ")\n")
  print(u)
}
