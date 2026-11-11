
suppressPackageStartupMessages({library(Seurat)})
DATA <- Sys.getenv("S2_DATA", Sys.getenv("DSCI_FY_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "FY")))
OUT  <- Sys.getenv('S2_OUT', 'figS2_out')
FY   <- file.path(DATA, '3k_dsci')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

suppressPackageStartupMessages({library(Seurat)})
DATA <- Sys.getenv("S2_DATA", Sys.getenv("DSCI_FY_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "FY")))
OUT  <- Sys.getenv('S2_OUT', 'figS2_out')
FY   <- file.path(DATA, '3k_dsci')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ---- the rule, derived from the ORIGINAL object (cell_type has 36 types,     --------
#      no COP/NFOL yet), with the compact marker table as a cross-check only.   --------
#      annotCOP_marker_df.rds is a POST-split export (its cell_type already     --------
#      contains the 269 COP/NFOL) -- using it as the derivation input silently  --------
#      yields 0 splits and "passes" because the composition is unchanged.       --------
seu <- readRDS(file.path(FY, '3k_dsci_full.RDS'))
mdf <- readRDS(file.path(FY, 'annotCOP_marker_df.rds'))
stopifnot(ncol(seu) == nrow(mdf))                       # same cell set
rna <- GetAssayData(seu, assay = 'RNA', slot = 'data')
g <- as.numeric(rna['Gpr17', ]); b <- as.numeric(rna['Bmp4', ])
old_oligo <- as.character(seu$cell_type) == 'Oligo'
cop <- old_oligo & (g > 0 | b > 0)
new <- as.character(seu$cell_type); new[cop] <- 'COP/NFOL'
cat('cells:', ncol(seu), '| old Oligo:', sum(old_oligo),
    '| split to COP/NFOL:', sum(cop),
    sprintf('(%.2f%%)', 100 * mean(cop)), '
')
cat('types after split:', length(unique(new)), '
')
cat('COP/NFOL by RNA cluster:',
    paste(names(table(seu$RNA_snn_res.1[cop])), table(seu$RNA_snn_res.1[cop]),
          sep = '=', collapse = ', '), '
')
cat(sprintf('within old Oligo: Gpr17>0 %.1f%% | Bmp4>0 %.1f%% | both %.1f%%
',
            100 * mean(g[old_oligo] > 0), 100 * mean(b[old_oligo] > 0),
            100 * mean(g[old_oligo] > 0 & b[old_oligo] > 0)))
# cross-check against the post-split marker export: same cells, same labels
mk <- match(mdf$cell, colnames(seu))
stopifnot(identical(as.character(mdf$cell_type[mk]), new))
cat('marker-df cross-check: identical labels for all', nrow(mdf), 'cells
')
write.csv(data.frame(cell = colnames(seu), cell_type = new,
                     rna_cluster = as.character(seu$RNA_snn_res.1)),
          file.path(OUT, 'truth_annotCOP.csv'), row.names = FALSE)

# ---- composition of the truth object (panel D's x order comes from this) ---
comp <- sort(table(new), decreasing = TRUE)
write.csv(data.frame(type = names(comp), Freq = as.integer(comp)),
          file.path(OUT, 'truth_composition.csv'), row.names = FALSE)
fdel <- file.path(FY, '3k_dsci_full_annotCOP_meancelltype.csv')
if (!file.exists(fdel)) stop('delivered composition table missing: ', fdel)
deliv <- read.csv(fdel)
stopifnot(nrow(deliv) == length(comp), sum(deliv$Freq) == sum(comp),
          identical(sort(deliv$type), sort(names(comp))))
cat('composition matches the delivered table (37 types,', sum(comp), 'cells)
')

# ---- optionally rebuild the RDS the experiment runs on --------------------
if (Sys.getenv('S2_BUILD_TRUTH') == '1') {
  idx <- which(cop)                       # derived above, no second 5.4 GB read
  ct <- as.character(seu$cell_type); ct[idx] <- 'COP/NFOL'
  seu$cell_type <- factor(ct)
  if (!is.null(seu$cell_type_wnn)) {
    cw <- as.character(seu$cell_type_wnn)
    ow <- cw == 'Oligo'
    cw[ow & (g > 0 | b > 0)] <- 'COP/NFOL'
    seu$cell_type_wnn <- factor(cw)
  }
  saveRDS(seu, file.path(OUT, '3k_dsci_full_annotCOP.RDS'))
  cat('wrote', file.path(OUT, '3k_dsci_full_annotCOP.RDS'),
      '| COP/NFOL:', sum(ct == 'COP/NFOL'), '
')
}

# ---- composition of the truth object (panel D's x order comes from this) ---
comp <- sort(table(new), decreasing = TRUE)
write.csv(data.frame(type = names(comp), Freq = as.integer(comp)),
          file.path(OUT, 'truth_composition.csv'), row.names = FALSE)
fdel <- file.path(FY, '3k_dsci_full_annotCOP_meancelltype.csv')
if (!file.exists(fdel)) stop('delivered composition table missing: ', fdel)
deliv <- read.csv(fdel)
stopifnot(nrow(deliv) == length(comp), sum(deliv$Freq) == sum(comp),
          identical(sort(deliv$type), sort(names(comp))))
cat('composition matches the delivered table (37 types,', sum(comp), 'cells)\n')

# ---- optionally rebuild the RDS the experiment runs on --------------------
if (Sys.getenv('S2_BUILD_TRUTH') == '1') {
  seu <- readRDS(file.path(FY, '3k_dsci_full.RDS'))
  stopifnot(ncol(seu) == nrow(mdf))
  rna <- GetAssayData(seu, assay = 'RNA', slot = 'data')
  g <- as.numeric(rna['Gpr17', ]); b <- as.numeric(rna['Bmp4', ])
  idx <- which(as.character(seu$cell_type) == 'Oligo' & (g > 0 | b > 0))
  ct <- as.character(seu$cell_type); ct[idx] <- 'COP/NFOL'
  seu$cell_type <- factor(ct)
  if (!is.null(seu$cell_type_wnn)) {
    cw <- as.character(seu$cell_type_wnn)
    ow <- cw == 'Oligo'
    cw[ow & (g > 0 | b > 0)] <- 'COP/NFOL'
    seu$cell_type_wnn <- factor(cw)
  }
  saveRDS(seu, file.path(OUT, '3k_dsci_full_annotCOP.RDS'))
  cat('wrote', file.path(OUT, '3k_dsci_full_annotCOP.RDS'),
      '| COP/NFOL:', sum(ct == 'COP/NFOL'), '\n')
}
