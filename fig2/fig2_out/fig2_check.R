
suppressPackageStartupMessages({library(Seurat); library(Signac); library(ggplot2)
  library(dplyr); library(tidyr); library(pals)})
IN  <- Sys.getenv("FIG2_INPUT", Sys.getenv("FIG1_DATA", ""))
FY  <- Sys.getenv("FIG2_DATA", Sys.getenv("DSCI_FY_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "FY")))
OUT <- Sys.getenv('FIG2_OUT', 'fig2_out')
p <- c(file.path(FY, 'FY13/SoupxMatrix', 'matrix.mtx.gz'),
       file.path(FY, 'FY14/SoupxMatrix', 'matrix.mtx.gz'),
       list.files(file.path(FY, 'FY13'), 'RNA_soupx_PCA.*[.]RDS$', full.names = TRUE)[1],
       list.files(file.path(FY, 'FY14'), 'RNA_soupx_PCA.*[.]RDS$', full.names = TRUE)[1],
       file.path(IN, 'FY/FY13/ATAC/out/Peak/matrix.mtx'),
       file.path(IN, 'FY/FY14/ATAC/out/Peak/matrix.mtx'),
       file.path(FY, 'FY13/FY13_ATAC_object.RDS'),
       file.path(FY, 'FY14/FY14_ATAC_object_force3k.RDS'),
       file.path(IN, 'scCAT-seq-HT/Mbrain_process.rds'),
       file.path(FY, 'archr-subset-copnfol/Save-ArchR-Project.rds'),
       file.path(FY, '3k_dsci_ht_10x_gs_subset.RDS'))
for (f in p) cat(sprintf('%-62s %s\n', sub(IN, '<IN>', sub(FY, '<DATA>', f, fixed = TRUE),
                              fixed = TRUE), if (file.exists(f)) 'present' else 'MISSING'))
cat('R', R.version.string, '| Seurat', as.character(packageVersion('Seurat')),
    '| Signac', as.character(packageVersion('Signac')), '| ArchR',
    tryCatch(as.character(packageVersion('ArchR')), error = function(e) 'NO'),
    '| harmony', tryCatch(as.character(packageVersion('harmony')), error = function(e) 'NO'),
    '| pals', as.character(packageVersion('pals')), '\n')
