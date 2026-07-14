#!/usr/bin/env Rscript
# Phase 3.1 (R) — functional COSGR test + R code-block parse checks.
# A Python-only suite would silently green-light broken R docs, so the R path is not optional.
#
#   1. Parse check: extract every ```r block from canonical/*.md and parse() it.
#   2. Functional: run COSGR (cosg) on a small Seurat object and assert it recovers markers.
#
# Run:  Rscript tests/test_r_blocks.R
# Needs: SeuratObject + COSG (remotes::install_github("genecell/COSGR")).

suppressMessages({library(SeuratObject); library(COSG)})
root <- normalizePath(file.path(dirname(sub("--file=", "",
        grep("--file=", commandArgs(FALSE), value = TRUE))), ".."))
canon <- file.path(root, "canonical")
fail <- 0

## 1. parse every ```r block ------------------------------------------------
mds <- list.files(canon, pattern = "\\.md$", recursive = TRUE, full.names = TRUE)
for (md in mds) {
  txt <- paste(readLines(md, warn = FALSE), collapse = "\n")
  blocks <- regmatches(txt, gregexpr("```r\\n(.*?)```", txt, perl = TRUE))[[1]]
  i <- 0
  for (b in blocks) {
    i <- i + 1
    code <- sub("```r\\n", "", b); code <- sub("```$", "", code)
    res <- tryCatch({ parse(text = code); TRUE },
                    error = function(e) { cat(sprintf("PARSE FAIL %s block %d: %s\n",
                                            basename(md), i, conditionMessage(e))); FALSE })
    if (!isTRUE(res)) fail <- fail + 1
  }
}
cat("R block parse checks done.\n")

## 2. functional COSGR ------------------------------------------------------
set.seed(1)
ng <- 300; nc <- 200
counts <- matrix(rpois(ng * nc, lambda = 1), nrow = ng, ncol = nc)
rownames(counts) <- paste0("Gene", seq_len(ng)); colnames(counts) <- paste0("Cell", seq_len(nc))
grp <- factor(rep(c("A", "B", "C"), length.out = nc))
counts[1:10,  grp == "A"] <- counts[1:10,  grp == "A"] + 20
counts[11:20, grp == "B"] <- counts[11:20, grp == "B"] + 20
counts[21:30, grp == "C"] <- counts[21:30, grp == "C"] + 20
counts <- as(counts, "dgCMatrix")
obj <- CreateSeuratObject(counts = counts)
libsize <- Matrix::colSums(counts)
norm <- log1p(t(t(counts) / libsize) * 1e4)
LayerData(obj, assay = "RNA", layer = "data") <- as(norm, "dgCMatrix")
Idents(obj) <- grp
res <- cosg(obj, groups = "all", assay = "RNA", slot = "data", mu = 1, n_genes_user = 10)

ok_groups <- all(c("A", "B", "C") %in% colnames(res$names))
# group A markers should be enriched for the spiked Gene1..Gene10
topA <- head(res$names$A, 10)
recovered <- sum(topA %in% paste0("Gene", 1:10))
cat(sprintf("COSGR groups OK: %s | group-A spiked markers recovered: %d/10\n",
            ok_groups, recovered))
if (!ok_groups || recovered < 5) { cat("FUNCTIONAL FAIL: COSGR did not recover markers\n"); fail <- fail + 1 }

if (fail > 0) { cat(sprintf("R TESTS FAILED: %d\n", fail)); quit(status = 1) }
cat("ALL R TESTS PASSED\n")
