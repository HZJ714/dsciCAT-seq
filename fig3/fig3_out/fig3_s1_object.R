
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr); library(pals)})
DATA <- Sys.getenv("FIG3_DATA", Sys.getenv("DSCI_FY_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "FY")))
OUT  <- Sys.getenv('FIG3_OUT', 'fig3_out')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
SAVE_DERIVED <- Sys.getenv('FIG3_SAVE_DERIVED', '') == '1'   # write *_from_joint.RDS

sub <- c(EL14 = 'HT/EL14_subset.RDS', FY13 = 'FY13/3k_FY13_subset.RDS', FY14 = 'FY14/3k_FY14_subset.RDS')
copclust <- c(EL14 = NA_integer_, FY13 = 29L, FY14 = 31L)   # COP/NFOL cluster of each dsciCAT library
# the rule that picks a library out of the joint object (see the text above):
rule <- list(EL14 = quote(batch == 'Cortex_EL14'),
             FY13 = quote(orig.ident == 'dsciCAT-seq-50k'),
             FY14 = quote(orig.ident == 'dsciCAT-seq-100k'))
libs <- list(); log <- character(); chk <- list()
t0 <- proc.time()[3]
ht <- readRDS(file.path(DATA, '3k_dsci_ht_10x_gs_subset.RDS')); DefaultAssay(ht) <- 'RNA'
log <- c(log, sprintf('ht_10x cells=%d loaded %.1f min', ncol(ht), (proc.time()[3] - t0) / 60))
cat(tail(log, 1), '\n')
pick <- lapply(rule, function(r) Cells(ht)[which(eval(r, ht@meta.data))])
selector <- function(nm) { r <- as.character(rule[[nm]])
                          paste0('ht$', r[2], " == '", r[3], "'") }
for (nm in names(pick)) cat(sprintf('   %s <- subset(ht, cells = Cells(ht)[which(%s)]): %d of %d cells\n',
                                    nm, selector(nm), length(pick[[nm]]), ncol(ht)))

## one pass over the three libraries: read the delivered object, check that it is
## exactly this subset of the joint object, and keep what the panels need
RNACOLS <- c('nCount_RNA', 'nFeature_RNA', 'percent.mt', 'percent.ribo', 'outlier',
             'scDblFinder.class', 'orig.ident', 'batch', 'rna_celltype')
ATACOLS <- c('nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 'FRiP', 'nucleosome_signal')
meta_long <- data.frame()
for (nm in names(sub)) {
  t0 <- proc.time()[3]
  o <- readRDS(file.path(DATA, sub[[nm]])); DefaultAssay(o) <- 'RNA'
  # the annotation of the figure: cell_type, but the COP/NFOL cluster by name
  ct <- as.character(o$cell_type)
  if (!is.na(copclust[[nm]])) ct[o$RNA_snn_res.1 == copclust[[nm]]] <- 'COP/NFOL'
  sel <- if (!is.na(copclust[[nm]])) which(o$RNA_snn_res.1 == copclust[[nm]]) else which(ct == 'OPC')
  libs[[nm]] <- list(n = ncol(o), ct = ct, cells = colnames(o),
                     rna = as.data.frame(Embeddings(o, 'umap')),          # panel S3A
                     inset_cells = colnames(o)[sel],
                     inset_counts = as.matrix(GetAssayData(o, assay = 'RNA', slot = 'counts')[, sel]),
                     gpr17 = as.numeric(GetAssayData(o, assay = 'RNA', slot = 'data')['Gpr17', ]))
  log <- c(log, sprintf('%s cells=%d types=%d inset=%d loaded %.1f min', nm, ncol(o),
                        length(unique(ct)), length(sel), (proc.time()[3] - t0) / 60))
  cat(tail(log, 1), '\n')
  jj <- pick[[nm]][match(colnames(o), pick[[nm]])]        # joint rows, in this object's order
  feq <- function(x, y) identical(as.numeric(as.matrix(x)), as.numeric(as.matrix(y)))
  g  <- intersect(rownames(o[['RNA']]),  rownames(ht[['RNA']]))[1:400]
  gp <- intersect(rownames(o[['ATAC']]), rownames(ht[['ATAC']]))[1:400]
  naeq <- function(a, b) { a <- as.character(a); b <- as.character(b)
                           mean((a == b) | (is.na(a) & is.na(b))) }   # NA == NA counts as equal
  me <- sapply(intersect(colnames(o@meta.data), colnames(ht@meta.data)), function(k)
        naeq(o@meta.data[[k]], ht@meta.data[jj, k]))
  meta_long <- rbind(meta_long, data.frame(library = nm, column = names(me), equal = as.numeric(me)))
  chk[[nm]] <- data.frame(
    library = nm, selector = selector(nm),
    cells = ncol(o), joint_cells = length(pick[[nm]]),
    same_cell_set = setequal(colnames(o), pick[[nm]]), same_cell_order = identical(colnames(o), jj),
    rna_counts_equal = feq(GetAssayData(o, assay = 'RNA',  slot = 'counts')[g, ],
                           GetAssayData(ht, assay = 'RNA',  slot = 'counts')[g, jj]),
    atac_counts_equal = feq(GetAssayData(o, assay = 'ATAC', slot = 'counts')[gp, ],
                            GetAssayData(ht, assay = 'ATAC', slot = 'counts')[gp, jj]),
    rna_genes = nrow(o[['RNA']]), rna_genes_joint = nrow(ht[['RNA']]),
    atac_peaks = nrow(o[['ATAC']]), atac_peaks_joint = nrow(ht[['ATAC']]),
    rna_qc_equal = round(mean(me[intersect(RNACOLS, names(me))]), 4),
    atac_qc_equal = round(mean(me[intersect(ATACOLS, names(me))]), 4),
    cell_type_equal = round(unname(me['cell_type']), 4),
    cluster_equal = round(unname(me['RNA_snn_res.1']), 4), stringsAsFactors = FALSE)
  cat(sprintf('   %s vs subset(ht, ...): cells %s order %s RNA matrix %s ATAC matrix %s | RNA QC meta %.4f | ATAC QC meta %.4f (%d peaks here vs %d in ht) | cell_type %.3f | RNA_snn_res.1 %.3f\n',
              nm, chk[[nm]]$same_cell_set, chk[[nm]]$same_cell_order, chk[[nm]]$rna_counts_equal,
              chk[[nm]]$atac_counts_equal, chk[[nm]]$rna_qc_equal, chk[[nm]]$atac_qc_equal,
              chk[[nm]]$atac_peaks, chk[[nm]]$atac_peaks_joint, chk[[nm]]$cell_type_equal,
              chk[[nm]]$cluster_equal))
  if (SAVE_DERIVED) {                       # the standalone subset, if you want the file
    d <- subset(ht, cells = pick[[nm]])
    saveRDS(d, file.path(OUT, paste0(nm, '_from_joint.RDS')))
    cat(sprintf('   wrote %s (%d cells, %d genes) - library-level analysis NOT included\n',
                paste0(nm, '_from_joint.RDS'), ncol(d), nrow(d)))
    rm(d); gc(verbose = FALSE)
  }
  rm(o); gc(verbose = FALSE)
}
chk <- do.call(rbind, chk)
write.csv(chk, file.path(OUT, 'fig3_subset_derivation.csv'), row.names = FALSE)
write.csv(meta_long[order(meta_long$column, meta_long$library), ],
          file.path(OUT, 'fig3_subset_meta_equal.csv'), row.names = FALSE)
log <- c(log, sprintf('derivation check written to %s', file.path(OUT, 'fig3_subset_derivation.csv')))
E <- list(atac = Embeddings(ht, 'umap.atac'), wnn = Embeddings(ht, 'wnn.umap'))
for (nm in names(libs)) {                      # same cells as the RNA panels: barcode match
  idx <- match(libs[[nm]]$cells, colnames(ht))
  libs[[nm]]$atac <- as.data.frame(E$atac[idx, , drop = FALSE])   # panel S3B
  libs[[nm]]$wnn  <- as.data.frame(E$wnn[idx,  , drop = FALSE])   # panels A-C
  log <- c(log, sprintf('%s ht_10x barcode match: %d / %d (%.1f%%)', nm, sum(!is.na(idx)),
                        length(idx), 100 * mean(!is.na(idx))))
  cat(tail(log, 1), '\n')
}
rm(ht, E); gc(verbose = FALSE)
# one colour per cell type across the three libraries
types <- sort(unique(unlist(lapply(libs, function(s) s$ct))))
cols <- unique(c(as.vector(pals::glasbey(32)), as.vector(pals::alphabet(26))))[1:length(types)]
names(cols) <- types
log <- c(log, sprintf('cell types (union of 3 libraries): %d', length(types)))
saveRDS(list(libs = libs, cols = cols, types = types), file.path(OUT, 'fig3_light.RDS'))
writeLines(log, file.path(OUT, 'fig3_s1_summary.txt'))
cat('wrote', file.path(OUT, 'fig3_light.RDS'), '\n')
