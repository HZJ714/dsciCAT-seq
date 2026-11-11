
suppressMessages({library(ArchR); library(Seurat); library(Signac)})
IN  <- Sys.getenv("FIG2_INPUT", Sys.getenv("FIG1_DATA", ""))
FY  <- Sys.getenv("FIG2_DATA", Sys.getenv("DSCI_FY_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "FY")))
OUT <- Sys.getenv('FIG2_OUT', 'fig2_out')
say <- function(...) cat(sprintf(...), '\n')
options(expressions = 500000); Sys.setenv(R_MAX_NUM_DLLS = 150)
plan('multicore', workers = 16); options(future.globals.maxSize = 1000 * 1024^10)
dir.create(file.path(FY, 'archr-3k_dsci_ht_10x'), showWarnings = FALSE, recursive = TRUE)
setwd(file.path(FY, 'archr-3k_dsci_ht_10x'))          # ArchR writes its logs to the cwd
addArchRGenome('mm10'); addArchRThreads(threads = 16)

prj     <- file.path(FY, 'archr-subset-copnfol')
arrows  <- file.path(FY, '3k_dsci_ht_10x_ArchR.RDS')
gs.rds  <- file.path(FY, '3k_dsci_ht_10x_gs_subset.RDS')
cand    <- readRDS(file.path(OUT, 'fig2_candidate_cells.RDS'))

## (a) arrow records: the candidate cells whose fragments ArchR can index
if (file.exists(arrows)) {
  pre <- readRDS(arrows)
  say('arrows loaded: %d cells | samples %s', length(pre$cellNames),
      paste(names(table(pre$Sample)), collapse = ', '))
} else {
  ## one fragment file per library - the delivered project's samples are
  ## FY13 FY14 EL13 EL14 EL15 EL16 10x, and sampleNames comes from these basenames
  fl <- c(list.files(file.path(IN, 'FY/FY13/ATAC/03.d2cfile'), '[.]fragments[.]tsv[.]gz$',
                     full.names = TRUE)[1],
          list.files(file.path(IN, 'FY/FY14/ATAC/03.d2cfile'), '[.]fragments[.]tsv[.]gz$',
                     full.names = TRUE)[1],
          list.files(file.path(IN, 'scCAT-seq-HT'), '[.]frag', full.names = TRUE))
  fl <- fl[file.exists(fl)]
  ArrowFiles <- createArrowFiles(inputFiles = fl, validBarcodes = cand$cell,
                                 sampleNames = sub('[.]frag.*', '', basename(fl)),
                                 addTileMat = TRUE, addGeneScoreMat = TRUE,
                                 threads = 16, force = TRUE, cleanTmp = TRUE)
  pre <- ArchRProject(ArrowFiles = ArrowFiles,
                      outputDirectory = file.path(FY, 'archr-3k_dsci_ht_10x'),
                      copyArrows = FALSE)
  saveRDS(pre, arrows)
}
stopifnot(length(pre$cellNames) == 83135)

## (b) ArchR's own default filter, measured on the arrows
tss <- as.numeric(pre@cellColData$TSSEnrichment)
nfr <- as.numeric(pre@cellColData$nFrags)
keep <- !(tss < 4 | nfr < 1000)
say('arrows: %d cells | dropped by the default filter: %d | survivors: %d',
    length(pre$cellNames), sum(!keep), sum(keep))
say('below minFrags = 1000: %d | below minTSS = 4: %d', sum(nfr < 1000), sum(tss < 4))
print(table(pre$Sample[keep]))
stopifnot(sum(keep) == 70643)

## (c) the subsetted project - the object §9 draws its tracks from
if (file.exists(file.path(prj, 'Save-ArchR-Project.rds'))) {
  proj <- loadArchRProject(prj, showLogo = FALSE)
  say('project loaded: %s | %d cells | matrices: %s', basename(prj), length(proj$cellNames),
      paste(getAvailableMatrices(proj), collapse = ', '))
} else {
  proj <- subsetArchRProject(pre, cells = pre$cellNames[keep], outputDirectory = prj,
                             dropCells = TRUE, force = TRUE)
  proj <- addIterativeLSI(proj, useMatrix = 'TileMatrix', name = 'IterativeLSI', iterations = 2,
           clusterParams = list(resolution = c(0.1, 0.3, 0.5), sampleCells = 30000, n.start = 10),
           varFeatures = 25000, dimsToUse = 2:40, force = TRUE)
  proj <- addGroupCoverages(proj, groupby = 'Sample', force = TRUE)
  proj <- addPeakCalls(proj, ID = 'Conservative', type = 'Quantile', quantile = 0.99, force = TRUE)
  proj <- addPeakMatrix(proj)
  proj <- addImputeWeights(proj, useLogScaling = TRUE, parallel = 8)
  proj <- saveArchRProject(proj, outputDirectory = prj, load = FALSE)
}

## (d) the gene-score matrix: the third modality, and the fingerprint of this step
Matrix <- getMatrixFromProject(proj, useMatrix = 'GeneScoreMatrix', threads = 16)
geneDF <- ArchR:::.getFeatureDF(getArrowFiles(proj), 'GeneScoreMatrix')
rownames(Matrix) <- geneDF[, 'name']
Matrix <- assays(Matrix)[['GeneScoreMatrix']]
colnames(Matrix) <- sub('^.*#', '', colnames(Matrix))
say('GeneScoreMatrix: %d genes x %d cells | every cell sums to 1e4: %s',
    nrow(Matrix), ncol(Matrix), all(round(colSums(Matrix)) == 1e4))
stopifnot(nrow(Matrix) == 24333, ncol(Matrix) == 70643)

## (e) the pre-WNN object, built only if §5 still needs it
out.rds <- file.path(FY, '3k_dsci_ht_10x_gs_preWNN.RDS')
if (file.exists(out.rds)) {
  say('preWNN object present: %s (%.1f GB)', basename(out.rds), file.info(out.rds)$size / 1e9)
} else if (file.exists(gs.rds)) {
  say('preWNN object not rebuilt: %s already exists', basename(gs.rds))
} else {
  merged <- readRDS(file.path(FY, '3k_dsci_ht_10x_merged.RDS'))
  obj <- merged[, intersect(colnames(merged), colnames(Matrix))]
  stopifnot(ncol(obj) == 70643)
  obj[['gene_score']] <- CreateAssayObject(counts = Matrix[, colnames(obj), drop = FALSE])
  obj$nCount_gene_score <- colSums(Matrix[, colnames(obj), drop = FALSE])
  print(table(obj$orig.ident))
  rm(Matrix, merged); invisible(gc())
  saveRDS(obj, out.rds)
}
