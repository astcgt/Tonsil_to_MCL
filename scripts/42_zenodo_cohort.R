#!/usr/bin/env Rscript
# 步驟 42：第二個外部世代（Zenodo 18200523，5 位病人）
# 這批資料已附 cell_type 與 inferCNV 亞克隆標註，可做比 GSE303064 更乾淨的比較。
suppressPackageStartupMessages({ library(SingleCellExperiment); library(Matrix) })

f <- "/raid/jupyter-alz.ee09/Tonsil_to_MCL/data/external_mcl/zenodo18200523/MCL_subclones_SCE.rds"
out <- "/raid/jupyter-alz.ee09/Tonsil_to_MCL/results"

cat("[42] 讀取", f, "\n"); flush.console()
sce <- readRDS(f)
cat("[42] class:", class(sce)[1], " dim:", paste(dim(sce), collapse=" x "), "\n")
cat("[42] assays:", paste(assayNames(sce), collapse=", "), "\n")
cd <- as.data.frame(colData(sce))
cat("[42] colData 欄位:\n"); print(colnames(cd))
for (cn in colnames(cd)) {
  v <- cd[[cn]]
  if (is.character(v) || is.factor(v)) {
    u <- unique(as.character(v))
    if (length(u) <= 40) cat(sprintf("  %-22s n=%-3d : %s\n", cn, length(u),
                                     paste(head(u, 12), collapse=" | ")))
    else cat(sprintf("  %-22s n=%d (太多)\n", cn, length(u)))
  }
}
write.csv(cd, file.path(out, "42_zenodo_colData.csv"))
cat("[42] 已寫出 colData\n")

# 取 logcounts（若無則自行標準化）
an <- assayNames(sce)
mat <- if ("logcounts" %in% an) assay(sce, "logcounts") else {
  cat("[42] 無 logcounts，用 counts 自行標準化\n")
  cts <- assay(sce, an[1])
  lib <- Matrix::colSums(cts)
  log1p(t(t(cts) / lib * 1e4))
}
rn <- rownames(mat)
cat("[42] 基因命名範例:", paste(head(rn, 3), collapse=", "), "\n")

modules <- list(
  MHC2 = c("CIITA","HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","HLA-DQA1","HLA-DQB1","CD74","HLA-DMB"),
  MHC1 = c("HLA-A","HLA-B","HLA-C","B2M","TAP1","TAP2","NLRC5"),
  OCT  = c("POU2F2","POU2AF1","MEF2C","PAX5","EBF1"),
  AP1  = c("JUN","JUNB","FOS","FOSL2","BATF"),
  MCLm = c("CCND1","SOX11")
)
sc_scores <- sapply(modules, function(g) {
  idx <- which(rn %in% g)
  if (length(idx) < 2) return(rep(NA_real_, ncol(mat)))
  Matrix::colMeans(mat[idx, , drop = FALSE])
})
cd2 <- cbind(cd, as.data.frame(sc_scores))
write.csv(cd2, file.path(out, "42_zenodo_scores.csv"))
cat("[42] 已寫出模組分數\n")

# 若有 cell_type 欄位就直接比較
ctcol <- grep("cell_type|celltype|CellType", colnames(cd), value = TRUE)
if (length(ctcol)) {
  cc <- ctcol[1]
  cat("\n[42] 依", cc, "分組的模組平均:\n")
  agg <- aggregate(sc_scores, by = list(group = cd[[cc]]), FUN = function(x) mean(x, na.rm = TRUE))
  agg$n <- as.integer(table(cd[[cc]])[agg$group])
  print(agg[order(-agg$n), ], row.names = FALSE, digits = 3)
}
pcol <- grep("^patient$|^Patient$|^donor", colnames(cd), value = TRUE)
if (length(pcol) && length(ctcol)) {
  cat("\n[42] 病人 x 細胞型別 交叉表:\n")
  print(table(cd[[pcol[1]]], cd[[ctcol[1]]]))
}
cat("[42] DONE\n")
