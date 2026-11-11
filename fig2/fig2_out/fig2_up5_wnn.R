
suppressMessages({library(Seurat); library(Signac); library(harmony)})
FY  <- Sys.getenv("FIG2_DATA", Sys.getenv("DSCI_FY_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "FY")))
OUT <- Sys.getenv('FIG2_OUT', 'fig2_out')
say <- function(...) cat(sprintf(...), '\n')
plan('multicore', workers = 16); options(future.globals.maxSize = 1000 * 1024^10)

# ---- WNN_harmony.r, in full ---------------------------------------------
Process_Multiome_Integration <- function(obj, nfeatures = 3000, pc.nums = 30,
    dims_use_rna = 1:10, dims_use_atac = 2:10, res_clusters = c(0.1, 0.3, 0.5),
    var_int = 'orig.ident', reduction_use_rna = 'harmony',
    reduction_use_atac = 'harmony_lsi', run_harmony = TRUE) {
  DefaultAssay(obj) <- 'RNA'
  obj <- NormalizeData(obj); obj <- FindVariableFeatures(obj, nfeatures = nfeatures)
  obj <- ScaleData(obj); obj <- RunPCA(obj, npcs = pc.nums, reduction.name = 'pca')
  if (run_harmony) {
    obj <- RunHarmony(obj, group.by.vars = var_int, assay.use = 'RNA', reduction.use = 'pca',
                      plot.convergence = FALSE, reduction.save = reduction_use_rna)
    rna_red <- reduction_use_rna
  } else rna_red <- 'pca'
  obj <- FindNeighbors(obj, verbose = FALSE, reduction = rna_red, dims = dims_use_rna)
  obj <- FindClusters(obj, verbose = FALSE, algorithm = 3, resolution = res_clusters)
  obj <- RunUMAP(obj, reduction = rna_red, dims = dims_use_rna)              # 'umap' = panel A

  DefaultAssay(obj) <- 'ATAC'
  obj <- RunTFIDF(obj); obj <- FindTopFeatures(obj, min.cutoff = 'q50'); obj <- RunSVD(obj)
  if (run_harmony) {
    obj <- RunHarmony(obj, group.by.vars = var_int, assay.use = 'ATAC', reduction.use = 'lsi',
                      reduction.save = reduction_use_atac, project.dim = FALSE)
    atac_red <- reduction_use_atac
  } else atac_red <- 'lsi'
  obj <- FindNeighbors(obj, reduction = atac_red, dims = dims_use_atac, n.neighbors = 40)
  obj <- FindClusters(obj, verbose = FALSE, algorithm = 3, resolution = res_clusters,
                      cluster.name = 'ATAC_snn_res')
  obj <- RunUMAP(obj, reduction = atac_red, dims = dims_use_atac,
                 reduction.name = 'umap.atac', reduction.key = 'atacUMAP_')   # panel B

  obj <- FindMultiModalNeighbors(obj, reduction.list = list(rna_red, atac_red),
                                 dims.list = list(dims_use_rna, dims_use_atac),
                                 weighted.nn.name = 'wnn_harmony')
  obj <- FindClusters(obj, graph.name = 'wsnn', algorithm = 3, resolution = res_clusters,
                      verbose = FALSE)
  obj <- RunUMAP(obj, nn.name = 'wnn_harmony', reduction.name = 'wnn.umap',
                 reduction.key = 'wnnUMAP_', seed.use = 1234)                 # panel C
  return(obj)
}
# -------------------------------------------------------------------------

objf <- file.path(FY, '3k_dsci_ht_10x_gs_subset.RDS')
if (file.exists(objf)) {
  say('cached: %s (%.1f GB)', basename(objf), file.info(objf)$size / 1e9)
} else {
  obj <- readRDS(file.path(FY, '3k_dsci_ht_10x_gs_preWNN.RDS'))
  obj <- Process_Multiome_Integration(obj, pc.nums = 50, dims_use_rna = 1:50,
                                      dims_use_atac = 2:40,
                                      res_clusters = c(0.1, 0.3, 0.5, 1))
  ## gene activity from the peaks (aRNA) — used by other figures, kept here so the
  ## object is the one the panels were drawn from
  DefaultAssay(obj) <- 'ATAC'
  obj[['aRNA']] <- CreateAssayObject(counts = GeneActivity(obj))
  DefaultAssay(obj) <- 'aRNA'
  obj <- NormalizeData(obj, assay = 'aRNA', normalization.method = 'LogNormalize',
                       scale.factor = median(obj$nCount_RNA))
  saveRDS(obj, objf)
}
obj <- readRDS(objf)
say('cells: %d | assays: %s', ncol(obj), paste(Assays(obj), collapse = ', '))
say('reductions: %s', paste(Reductions(obj), collapse = ', '))
say('clusters @ RNA_snn_res.1: %d | @ wsnn_res.1: %d',
    length(levels(obj$RNA_snn_res.1)), length(levels(obj$wsnn_res.1)))
