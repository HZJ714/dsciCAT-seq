
suppressPackageStartupMessages({library(Seurat); library(harmony); library(parallel)})
DATA <- Sys.getenv("S2_DATA", Sys.getenv("DSCI_FY_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "FY")))
OUT  <- Sys.getenv('S2_OUT', 'figS2_out')
FY   <- file.path(DATA, '3k_dsci')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- scores ---------------------------------------------------------------
# aricode is the reference (NMI default variant = 'max'). It is not always
# installed and cannot be fetched offline, so the embedded implementation below
# is used as the fallback; it was verified against the package on simulated
# partitions (ARI diff = 0, NMI max-variant diff ~ 1e-16).
USE_ARICODE <- requireNamespace('aricode', quietly = TRUE)
cat('scoring backend:', if (USE_ARICODE) 'aricode package' else 'embedded (max variant)', '\n')

emb_nmi <- function(x, y) {
  ax  <- table(x, y); n <- sum(ax)
  pxy <- ax / n
  px  <- rowSums(ax) / n; py <- colSums(ax) / n
  ent <- function(p) -sum(p[p > 0] * log2(p[p > 0]))
  mi  <- sum(pxy[pxy > 0] * log2(pxy[pxy > 0] /
            as.vector(outer(px, py, '*'))[pxy > 0]))
  mi / max(ent(px), ent(py))
}
emb_ari <- function(x, y) {
  tab <- table(x, y); n <- sum(tab)
  a <- rowSums(tab); b <- colSums(tab)
  C2 <- function(v) sum(v * (v - 1) / 2)
  obs <- C2(as.vector(tab)); expd <- C2(a) * C2(b) / C2(n)
  maxd <- 0.5 * (C2(a) + C2(b))
  (obs - expd) / (maxd - expd)
}
score <- function(a, b) {
  if (!USE_ARICODE) return(list(NMI = emb_nmi(a, b), ARI = emb_ari(a, b)))
  lev <- union(as.character(unique(a)), as.character(unique(b)))   # aricode wants
  ia  <- match(as.character(a), lev); ib <- match(as.character(b), lev)  # integer codes
  return(list(NMI = aricode::NMI(ia, ib), ARI = aricode::ARI(ia, ib)))
}

# ---- per-type metrics from the contingency table --------------------------
per_type_metrics <- function(truth, clust, min_cluster_size = 5) {
  ct <- table(clust, truth)
  if (ncol(ct) == 0) return(NULL)
  plurality <- colnames(ct)[apply(ct, 1, which.max)]
  names(plurality) <- rownames(ct)
  pred <- plurality[as.character(clust)]
  types <- sort(unique(truth))
  do.call(rbind, lapply(types, function(t) {
    tp <- sum(pred == t & truth == t); fp <- sum(pred == t & truth != t)
    fn <- sum(pred != t & truth == t)
    det <- any(apply(ct, 1, function(r) colnames(ct)[which.max(r)] == t) &
                 rowSums(ct) >= min_cluster_size)
    rec  <- ifelse(tp + fn > 0, tp / (tp + fn), NA)
    prec <- ifelse(tp + fp > 0, tp / (tp + fp), NA)
    f1   <- ifelse(is.na(rec) | is.na(prec) | rec + prec == 0, NA,
                   2 * rec * prec / (rec + prec))
    data.frame(type = t, n_true = tp + fn, detected = as.integer(det),
               recall = rec, precision = prec, F1 = f1, stringsAsFactors = FALSE)
  }))
}

# ---- one down-sample = one full re-analysis -------------------------------
Process_RNA <- function(proj, nfeatures = 3000, pc.nums = 50) {
  DefaultAssay(proj) <- 'RNA'
  proj <- NormalizeData(proj)
  proj <- FindVariableFeatures(proj, nfeatures = nfeatures)
  proj <- ScaleData(proj, verbose = FALSE)
  proj <- RunPCA(proj, npcs = pc.nums, reduction.name = 'pca', verbose = FALSE)
  proj
}
run_one <- function(seu, df_full, target_cell_count, seed) {
  set.seed(seed)
  cells <- sample(colnames(seu), size = target_cell_count, replace = FALSE)
  obj <- Process_RNA(subset(seu, cells = cells))
  obj <- RunHarmony(obj, group.by.vars = 'orig.ident', assay.use = 'RNA',
                    reduction.use = 'pca', plot_convergence = FALSE)
  obj <- FindNeighbors(obj, reduction = 'harmony', dims = 1:50, verbose = FALSE)
  obj <- FindClusters(obj, algorithm = 3, resolution = 1, verbose = FALSE)
  clust <- as.character(obj$RNA_snn_res.1)
  truth <- as.character(obj$cell_type)

  ct <- table(clust, truth)
  auto <- colnames(ct)[apply(ct, 1, which.max)]; names(auto) <- rownames(ct)
  sub <- data.frame(cellname = colnames(obj), Cell_type = as.character(auto[clust]))
  both <- merge(df_full, sub, by = 'cellname', all = FALSE)
  stopifnot(nrow(both) == length(clust))
  s <- score(both$Cell_type.x, both$Cell_type.y)
  pt <- per_type_metrics(truth, clust)
  list(seed = seed, NMI = s$NMI, ARI = s$ARI, pt = pt)
}

# ---- the 8 x 20 loop, cached ---------------------------------------------
LEVELS  <- c('3k','5k','10k','15k','20k','25k','30k','40k')
# parallelism per level, as used by the run that produced the delivered CSVs:
# memory-bound, not core-bound (each worker holds its own copy of the object)
NCORES  <- c('3k'=20,'5k'=16,'10k'=12,'15k'=10,'20k'=8,'25k'=6,'30k'=5,'40k'=4)
if (nzchar(Sys.getenv('S2_NCORES'))) NCORES[] <- as.integer(Sys.getenv('S2_NCORES'))
if (nzchar(Sys.getenv('S2_LEVELS')))
  LEVELS <- intersect(LEVELS, strsplit(Sys.getenv('S2_LEVELS'), ',')[[1]])
seeds  <- if (nzchar(Sys.getenv('S2_SEEDS')))
            as.integer(strsplit(Sys.getenv('S2_SEEDS'), ',')[[1]]) else seq(1, 200, by = 10)
FORCE  <- Sys.getenv('S2_FORCE') == '1'

todo <- LEVELS[sapply(LEVELS, function(lv) FORCE || !(
  file.exists(file.path(OUT, paste0(lv, '_ARI_NMI.csv'))) &&
  file.exists(file.path(OUT, paste0(lv, '_per_seed_per_type.csv')))))]
cat('levels to run:', if (length(todo)) paste(todo, collapse = ',') else 'none (cache complete)',
    '| seeds:', length(seeds), '\n')

if (length(todo)) {
  seu <- readRDS(file.path(FY, '3k_dsci_full_annotCOP.RDS'))
  for (a in setdiff(Assays(seu), 'RNA')) seu[[a]] <- NULL     # RAM
  gc()
  df_full <- data.frame(cellname = rownames(seu@meta.data),
                        Cell_type = as.character(seu@meta.data$cell_type))
  cat('truth object:', ncol(seu), 'cells |', length(unique(df_full$Cell_type)),
      'types | COP/NFOL:', sum(df_full$Cell_type == 'COP/NFOL'), '\n')

  for (lv in todo) {
    n <- as.integer(gsub('k', '000', lv)); t0 <- Sys.time()
    res <- mclapply(seeds, function(s) tryCatch(run_one(seu, df_full, n, s),
                       error = function(e) { cat('[ERROR]', lv, s, conditionMessage(e), '\n'); NULL }),
                    mc.cores = min(NCORES[[lv]], length(seeds)), mc.preschedule = FALSE)
    # mclapply silently drops a worker that dies (OOM / fork instability) --
    # those seeds come back as NULL from the tryCatch OR as missing elements.
    # Retry any seed that did not deliver, one at a time, in the parent process.
    retry <- seeds[!vapply(seq_along(seeds), function(k)
              !is.null(res[[k]]) && !is.null(res[[k]]$seed), logical(1))]
    if (length(retry)) {
      cat('[retry]', length(retry), 'seeds:', paste(retry, collapse = ','), '\n')
      for (s in retry) {
        res[[match(s, seeds)]] <- tryCatch(run_one(seu, df_full, n, s),
                     error = function(e) { cat('[ERROR]', lv, s, conditionMessage(e), '\n'); NULL })
      }
    }
    ok  <- Filter(Negate(is.null), res)
    met <- do.call(rbind, lapply(ok, function(o) data.frame(Seed = o$seed, NMI = o$NMI, ARI = o$ARI)))
    pt  <- do.call(rbind, lapply(ok, function(o) { p <- o$pt; p$Seed <- o$seed; p$Level <- lv; p }))
    write.csv(met, file.path(OUT, paste0(lv, '_ARI_NMI.csv')), row.names = FALSE)
    write.csv(pt,  file.path(OUT, paste0(lv, '_per_seed_per_type.csv')), row.names = FALSE)
    cat(sprintf('[%s] %d/%d seeds in %.1f min | NMI %.3f+-%.3f | ARI %.3f+-%.3f | rows %d\n',
                lv, length(ok), length(seeds),
                as.numeric(difftime(Sys.time(), t0, units = 'mins')),
                mean(met$NMI), sd(met$NMI), mean(met$ARI), sd(met$ARI), nrow(pt)))
  }
}
