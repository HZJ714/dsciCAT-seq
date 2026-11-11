
suppressMessages({library(Seurat); library(Signac); library(DoubletFinder); library(dplyr)
  library(glue)})
IN  <- Sys.getenv("FIG2_INPUT", Sys.getenv("FIG1_DATA", ""))
FY  <- Sys.getenv("FIG2_DATA", Sys.getenv("DSCI_FY_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "FY")))
OUT <- Sys.getenv('FIG2_OUT', 'fig2_out')
say <- function(...) cat(sprintf(...), '\n')
set.seed(1234)
dim <- 50; minfeatures <- 300; mtegne_filter <- 10; doublets.percentage <- 0.1
atacs <- list()
for (s in c('FY13', 'FY14')) {                     # same preference rule as §1
  f <- list.files(file.path(FY, s), paste0(s, '_ATAC_object.*[.]RDS$'), full.names = TRUE)
  f <- f[order(grepl('force3k', f), decreasing = TRUE)]
  if (length(f)) atacs[[s]] <- readRDS(f[1])
}

for (s in c('FY13', 'FY14')) {
  d <- file.path(FY, s)
  fut <- file.path(d, paste0('3k_', s, '_multiome.RDS'))
  if (file.exists(fut)) {
    m <- readRDS(fut); say('%s: loaded %s | %d cells', s, basename(fut), ncol(m))
    rm(m); next
  }

  ## (a) RNA side: the QC-cut, singlet-called object of this SoupX run, or rebuild it
  rn <- list.files(d, paste0(s, '_RNA_soupx_PCA.*[.]RDS$'), full.names = TRUE)
  if (length(rn)) {
    rn <- rn[order(file.info(rn)$mtime, decreasing = TRUE)][1]     # latest run wins
    rna <- readRDS(rn)
    say('%s: RNA side loaded %s | %d cells', s, basename(rn), ncol(rna))
  } else {
    rna <- CreateSeuratObject(Read10X(file.path(d, 'SoupxMatrix'), gene.column = 1))
    rna[["percent.mt"]] <- PercentageFeatureSet(rna, pattern = "^mt-")
    om <- rna@meta.data[order(-rna@meta.data$nFeature_RNA), ]
    n95 <- as.integer(nrow(om) * 0.05)
    n95features <- as.numeric(om[n95, "nFeature_RNA"])
    rna <- subset(rna, subset = nFeature_RNA > minfeatures &
                                       nFeature_RNA < n95features &
                                       percent.mt < mtegne_filter)
    sweep.res.list <- paramSweep_v3(rna, PCs = 1:dim, sct = FALSE)
    sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
    bcmvn <- find.pK(sweep.stats)
    pK <- as.numeric(as.vector(bcmvn[bcmvn$MeanBC == max(bcmvn$MeanBC), ]$pK))
    nExp_poi <- round(doublets.percentage * ncol(rna))
    rna <- doubletFinder_v3(rna, PCs = 1:dim, pN = 0.25, pK = pK, nExp = nExp_poi,
                            reuse.pANN = FALSE, sct = FALSE)
    colnames(rna@meta.data)[ncol(rna@meta.data)] <- "doublet_info"
    write.table(rna@meta.data, file.path(d, paste0(s, "_doublets_info.txt")),
                sep = "\t", quote = FALSE)
    rna <- subset(rna, subset = doublet_info == "Singlet")
    saveRDS(rna, file.path(d, paste0(s, "_RNA_soupx_PCA.RDS")))
    say('%s: RNA side built | %d singlets', s, ncol(rna))
  }

  ## (b) keep the cells that have both modalities, then attach the peaks
  atac <- atacs[[s]]
  common <- intersect(colnames(rna), colnames(atac))
  say('%s: RNA %d | ATAC %d | intersect %d', s, ncol(rna), ncol(atac), length(common))
  rna  <- rna[, common]
  atac <- atac[, common]
  brain <- rna
  brain[["ATAC"]] <- atac[["ATAC"]]
  brain$nCount_ATAC  <- atac$nCount_ATAC
  brain$nFeature_ATAC <- atac$nFeature_ATAC
  brain$FRiP          <- atac$FRiP
  brain$TSS.enrichment <- atac$TSS.enrichment
  brain$orig.ident <- s
  saveRDS(brain, fut)
  say('%s: multiome -> %d cells', s, ncol(brain))
  rm(brain, rna, atac); invisible(gc())
}
