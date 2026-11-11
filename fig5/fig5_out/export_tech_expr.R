
suppressPackageStartupMessages({library(Seurat); library(SeuratObject); library(Matrix)})
export_tech <- function(rds, tag, dir) {        # writes genes-x-cells .mtx + gene list + meta
  obj <- readRDS(rds)                          # RNA + ATAC assays, so Signac must be loadable
  x <- tryCatch(obj[["RNA"]]$data, error = function(e) NULL)          # Seurat v4 slot API
  if (is.null(x)) x <- tryCatch(SeuratObject::LayerData(obj, assay = "RNA", layer = "data"),
                                error = function(e) NULL)             # Seurat v5 layer API
  if (is.null(x)) x <- tryCatch(GetAssayData(obj, assay = "RNA", layer = "data"),
                                error = function(e) NULL)             # v4 accessor
  if (is.null(x)) stop("no RNA data layer found in ", rds)
  if (length(x@x) == 0) stop("RNA data layer is empty in ", rds)
  x <- as(x, "CsparseMatrix")
  mtxf <- file.path(dir, paste0(tag, "_expr.mtx.gz"))
  writeMM(x, mtxf)          # PATH, never a gzfile() connection: writeMM to a connection
  if (file.size(mtxf) < 1000) stop("writeMM produced an empty matrix file: ", mtxf)
  writeLines(rownames(x), file.path(dir, paste0(tag, "_genes.txt")))
  write.csv(obj@meta.data, file.path(dir, paste0(tag, "_meta.csv")), row.names = TRUE)
  cat(sprintf("%s: %d cells x %d genes, nnz %d\n", tag, ncol(x), nrow(x), length(x@x)))
}
args <- commandArgs(TRUE)
export_tech(args[1], "tenx",       args[3])
export_tech(args[2], "HT_harmony", args[3])
