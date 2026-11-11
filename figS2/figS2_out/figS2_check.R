
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr)})
DATA <- Sys.getenv("S2_DATA", Sys.getenv("DSCI_FY_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "FY")))
OUT  <- Sys.getenv('S2_OUT', 'figS2_out')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
need <- c('3k_dsci_ht_10x_gs_subset.RDS', '3k_dsci/3k_dsci_full.RDS',
          '3k_dsci/annotCOP_marker_df.rds', '3k_dsci/3k_dsci_full_annotCOP.RDS',
          '3k_dsci/3k_dsci_full_annotCOP_meancelltype.csv')
for (f in need) cat(sprintf('%-52s %s  %6.2f GB\n', f,
    file.exists(file.path(DATA, f)), ifelse(file.exists(file.path(DATA, f)),
    file.info(file.path(DATA, f))$size / 1e9, NA)))
lv <- c('3k','5k','10k','15k','20k','25k','30k','40k')
cat('\n§4 cache present (skipped when both CSVs of a level exist):\n')
print(data.frame(level = lv,
                 ari_nmi = file.exists(file.path(DATA, '3k_dsci',
                                                 paste0('3k_dsci_down', lv, '_ARI_NMI.csv'))),
                 per_type = file.exists(file.path(DATA, '3k_dsci',
                                                  paste0(lv, '_per_seed_per_type.csv')))))
cat('\noptional packages:',
    paste(sapply(c('aricode', 'ggtext'), function(p)
        sprintf('%s=%s', p, requireNamespace(p, quietly = TRUE))), collapse = '  '), '\n')
