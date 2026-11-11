
suppressPackageStartupMessages({library(Seurat)})
DATA <- Sys.getenv("FIGS3_DATA", Sys.getenv("DSCI_FY_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "FY"))); OUT <- Sys.getenv('FIGS3_OUT', 'figS3_out')
MK   <- c('Gpr17', 'Bmp4', 'C1ql1', 'Tmem144', 'Ptprz1', 'Plp1', 'Mbp')
RES  <- c(0.5, 1.0, 1.5, 2.0); PRIMARY <- 'RNA_snn_res.1'
NSEED <- as.integer(Sys.getenv('NSEED', '20'))
loadk <- function(k) CreateSeuratObject(counts = readRDS(file.path(OUT, paste0('S3_lineage_', k, '.rds'))),
                                        project = k, assay = 'RNA')
recly <- function(o) { set.seed(42)                       # same seed for every condition
  o <- NormalizeData(o, verbose = FALSE)
  o <- suppressWarnings(FindVariableFeatures(o, selection.method = 'vst', nfeatures = 3000, verbose = FALSE))
  # no ScaleData here on purpose: the cached assay brings the whole-brain scale.data,
  # RunPCA then uses the genes the two HVG lists share - exactly what the delivered
  # run did. The count is printed so the reader can see how small that space is.
  cat(sprintf('   %s cells: %d vst genes, %d of them scaled
', ncol(o),
              length(VariableFeatures(o)),
              sum(VariableFeatures(o) %in% rownames(GetAssayData(o, assay = 'RNA',
                                                                slot = 'scale.data')))))
  o <- RunPCA(o, npcs = 30, verbose = FALSE)
  o <- FindNeighbors(o, reduction = 'pca', dims = 1:10, graph.name = 'RNA_snn', verbose = FALSE)
  FindClusters(o, algorithm = 3, resolution = RES, verbose = FALSE) }
tab <- function(o, cond, seed = NA) {
  gs <- intersect(MK, rownames(o))
  ex <- as.matrix(GetAssayData(o, assay = 'RNA', slot = 'data')[gs, , drop = FALSE]); storage.mode(ex) <- 'double'
  out <- NULL
  for (cn in paste0('RNA_snn_res.', RES)) {
    cl <- as.character(o@meta.data[[cn]]); lv <- sort(unique(cl))
    W <- sapply(lv, function(k) as.numeric(cl == k)); n <- colSums(W)
    # column-wise division: (ex %*% W) / n would recycle element-wise and is WRONG
    m <- sweep(ex %*% W, 2, n, '/'); p <- 100 * sweep((ex > 0) %*% W, 2, n, '/')
    d <- data.frame(cond = cond, seed = ifelse(is.na(seed), -1, seed), res = cn, cl = lv,
                    n = as.integer(n), nlin = sum(n), stringsAsFactors = FALSE)
    for (k in MK) d[[paste0('mean_', k)]] <- as.numeric(m[k, ])
    for (k in MK) d[[paste0('pct_', k)]]  <- as.numeric(p[k, ])
    d$cop <- with(d, mean_Bmp4 >= 0.25 & mean_Gpr17 >= 0.5 & pct_Gpr17 >= 30 &
                          mean_C1ql1 <= 0.2 & mean_Ptprz1 <= 1.5 & n >= 10)
    d$celllevel_copfrac <- 100 * mean((ex['Gpr17', ] > 0) & (ex['Bmp4', ] > 0) &
                                      (ex['C1ql1', ] == 0) & (ex['Ptprz1', ] < 1))
    out <- rbind(out, d) }
  stopifnot(all(out$pct_Gpr17 <= 100 + 1e-6, na.rm = TRUE))
  out }
subN <- function(o, n, seed) { set.seed(seed)
  # materialise the sample: subset(o, cells = sample(...)) lets Seurat reuse the
  # promise, and every seed then draws the SAME cells (this is what turned 20
  # independent draws into one and changed the published x/20 numbers)
  sel <- sample(colnames(o), min(n, ncol(o))); subset(o, cells = sel) }
ALL <- NULL; t0 <- proc.time()[3]
lin <- lapply(c('FY13', 'FY14', 'HTp', 'EL14'), loadk); names(lin) <- c('FY13', 'FY14', 'HTp', 'EL14')
for (k in names(lin)) { d <- tab(recly(lin[[k]]), paste0(k, '_full_lineage')); ALL <- rbind(ALL, d)
  cat(sprintf('%-22s lineage %5d | COP cluster at %s: %s\n', paste0(k, '_full_lineage'), d$nlin[1],
              PRIMARY, ifelse(any(d$cop & d$res == PRIMARY), 'YES', 'no'))) }
N4029 <- round(ncol(lin$FY13) * 4029 / 17093)   # lineage cells left when 50k is thinned to 4,029 nuclei
cat(sprintf('(dsciCAT whole brain at 4,029 nuclei = lineage %d cells)\n', N4029))
SPECS <- list(c('FY13_wholebrain_to_4029n', 'FY13', N4029),
              c('FY13_lineage_to_HTpooledN', 'FY13', ncol(lin$HTp)),
              c('FY13_lineage_to_EL14N', 'FY13', ncol(lin$EL14)),
              c('HTpooled_lineage_to_EL14N', 'HTp', ncol(lin$EL14)))
for (sd in seq_len(NSEED)) {
  for (sp in SPECS) {
    d <- tab(recly(subN(lin[[sp[2]]], as.integer(sp[3]), sd)), sp[1], sd); ALL <- rbind(ALL, d) }
  cat('.', ''); flush(stdout()) }
cat('\n')
write.csv(ALL, file.path(OUT, 'S3_u13_clusters.csv'), row.names = FALSE)
# panel G's table is the four full conditions at the primary resolution
mt <- ALL[grepl('_full_lineage', ALL$cond) & ALL$res == PRIMARY, ]
mt$group <- sub('_lineage', '', mt$cond)
write.csv(mt, file.path(OUT, 'S3_marker_table.csv'), row.names = FALSE)
# the strict depth: thin the WHOLE BRAIN to 4,029 nuclei first, then take the lineage
fy <- readRDS(file.path(DATA, 'FY13/3k_FY13_subset.RDS')); DefaultAssay(fy) <- 'RNA'
S <- NULL
for (sd in seq_len(NSEED)) { set.seed(sd)
  sel <- sample(colnames(fy), min(4029, ncol(fy))); w <- subset(fy, cells = sel)
  keep <- as.character(w$cell_type) %in% c('Oligo', 'OPC')
  w2 <- subset(w, cells = colnames(w)[keep])
  d <- tab(recly(w2), 'FY13_strict_4029nuclei', sd); S <- rbind(S, d)
  cat(sprintf('strict seed=%2d lineage %d cells | COP cluster %s\n', sd, sum(d$n),
              any(d$cop & d$res == PRIMARY))) }
write.csv(S, file.path(OUT, 'S3_strict4029.csv'), row.names = FALSE)
cat(sprintf('%.0f s | %d condition rows, %d marker rows, %d strict rows\n',
            (proc.time()[3] - t0), nrow(ALL), nrow(mt), nrow(S)))
