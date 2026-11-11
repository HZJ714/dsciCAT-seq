
suppressMessages({library(Seurat); library(Signac); library(rtracklayer); library(GenomeInfoDb)})
IN  <- Sys.getenv("FIG2_INPUT", Sys.getenv("FIG1_DATA", ""))
FY  <- Sys.getenv("FIG2_DATA", Sys.getenv("DSCI_FY_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "FY")))
OUT <- Sys.getenv('FIG2_OUT', 'fig2_out')
say <- function(...) cat(sprintf(...), '\n')
gtf <- file.path(Sys.getenv("DSCI_REFDATA_ROOT", ""), "mm10/scATAC_ref/genes.gtf")

## §2 needs the two objects, and the same preference rule resolves them there.
ATAC_PATH <- function(s) {
  f <- list.files(file.path(FY, s), paste0(s, '_ATAC_object.*[.]RDS$'), full.names = TRUE)
  f[order(grepl('force3k', f), decreasing = TRUE)][1]                 # the 3k-force run first
}

atacs <- list()
for (s in c('FY13', 'FY14')) {
  saved <- list.files(file.path(FY, s), paste0(s, '_ATAC_object.*[.]RDS$'), full.names = TRUE)
  saved <- saved[order(grepl('force3k', saved), decreasing = TRUE)]   # the 3k-force run first
  if (length(saved)) {
    atacs[[s]] <- readRDS(saved[1])
    say('%s: loaded %s | %d peaks x %d cells | FRiP median %.3f | TSS enrichment median %.2f',
        s, basename(saved[1]), nrow(atacs[[s]]), ncol(atacs[[s]]),
        median(atacs[[s]]$FRiP), median(atacs[[s]]$TSS.enrichment))
    next
  }
  ## built from the peak call of the same bead run
  d          <- file.path(IN, 'FY', s)
  mtx_path   <- file.path(d, 'ATAC/out/Peak/matrix.mtx')
  feat_path  <- file.path(d, 'ATAC/out/Peak/peak.bed')
  bcs_path   <- file.path(d, 'ATAC/out/Peak/barcodes.tsv')
  frag.file  <- list.files(file.path(d, 'ATAC/03.d2cfile'), '[.]fragments[.]tsv[.](gz|txt)$',
                           full.names = TRUE)[1]
  counts <- Matrix::readMM(mtx_path)
  features <- readr::read_tsv(feat_path, col_names = FALSE) %>% tidyr::unite(feature)
  barcodes <- readr::read_tsv(bcs_path, col_names = FALSE) %>% tidyr::unite(barcode)
  row.names(counts) <- features$feature; colnames(counts) <- barcodes$barcode

  gtf.gr <- if (file.exists(gtf)) rtracklayer::import(gtf) else
            readRDS(file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "mm10.gene.coords.RDS"))
  gene.coords <- gtf.gr[gtf.gr$type == 'gene']
  seqlevelsStyle(gene.coords) <- 'UCSC'
  gene.coords <- keepStandardChromosomes(gene.coords, pruning.mode = 'coarse')

  g <- StringToGRanges(rownames(counts), sep = c(':', '-'))
  counts <- counts[as.vector(seqnames(g) %in% standardChromosomes(g)), ]

  a <- CreateChromatinAssay(counts = counts, sep = c(':', '-'), fragments = frag.file,
                            min.cells = 10, annotation = gene.coords)
  obj <- CreateSeuratObject(counts[, colnames(a)], project = s)
  obj[['ATAC']] <- a
  DefaultAssay(obj) <- 'ATAC'
  obj <- TSSEnrichment(obj)
  tf <- CountFragments(frag.file); rownames(tf) <- tf$CB
  obj$fragments <- tf[colnames(obj), 'frequency_count']
  obj <- FRiP(obj, assay = 'ATAC', total.fragments = 'fragments')
  obj$loguniqueFrag <- log10(obj$fragments)
  obj$sample <- s
  saveRDS(obj, file.path(FY, s, paste0(s, '_ATAC_object.RDS')))
  atacs[[s]] <- obj
  say('%s: built from %s | %d peaks x %d cells', s, basename(frag.file), nrow(obj), ncol(obj))
}
say('checked %d objects: %s', length(atacs),
    paste(sapply(atacs, function(x) sprintf('%d peaks x %d cells', nrow(x), ncol(x))),
          collapse = ' | '))
