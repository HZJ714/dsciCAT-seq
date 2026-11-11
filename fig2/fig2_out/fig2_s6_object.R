
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr)
  library(tidyr); library(pals)})
DATA <- Sys.getenv("FIG2_DATA", Sys.getenv("DSCI_FY_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "FY")))
OUT  <- Sys.getenv('FIG2_OUT', 'fig2_out')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
rds <- file.path(DATA, '3k_dsci_ht_10x_gs_subset.RDS')
t0 <- proc.time()[3]
obj <- readRDS(rds)
cat('loaded', rds, 'in', round((proc.time()[3] - t0) / 60, 1), 'min
')

# ---- the final annotation (used by every panel) ---------------------------
obj$cell_type_updated <- as.character(obj$cell_type)
obj$cell_type_updated[obj$RNA_snn_res.1 == 34] <- 'COP/NFOL'   # 280 cells
obj$technique <- ifelse(obj$orig.ident %in% c('dsciCAT-seq-50k', 'dsciCAT-seq-100k'),
                        'dsciCAT-seq', as.character(obj$orig.ident))

# fixed 33-type order = track/dotplot order; numbering 1=AST ... 33=Ex_L6b
dot_order <- c('AST','Oligo','COP/NFOL','OPC','MG','EC','VLMC','Ependyma',
               'In_D1','In_D2','In_HAP1','In_LAMP5','In_LHX6','In_SST','In_VIP',
               'In_PVALB','In_TH','Neuron_MB','Ex_LHA','Ex_HY-MM','Ex_DG','Ex_TH',
               'Ex_MH-LH','Ex_SUB_Pros','Ex_CA1_Pros','Ex_CAR3','Ex_L2/3_IT',
               'Ex_L2/3_IT_PIR','Ex_L4/5_IT','Ex_L5/6_NP','Ex_L6_IT','Ex_L6_CT','Ex_L6b')
num_map <- setNames(seq_along(dot_order), dot_order)

# shared palette = global_colors_full (unique(c(glasbey(32), alphabet(26)))[1:42])
all_types_full <- sort(unique(c(
  'AST','COP/NFOL','EC','Ependyma','Ex_CA1_Pros','Ex_CA3','Ex_CAR3','Ex_CR',
  'Ex_DG','Ex_HY-MM','Ex_IC','Ex_L2/3_IT','Ex_L2/3_IT_PIR','Ex_L4','Ex_L4/5_IT',
  'Ex_L5_ET','Ex_L5_IT','Ex_L5_NP','Ex_L5/6_NP','Ex_L6_CT','Ex_L6_IT','Ex_L6b',
  'Ex_LHA','Ex_MH-LH','Ex_SUB_Pros','Ex_TH','In_D1','In_D2','In_D1-D2','In_HAP1',
  'In_LAMP5','In_LHX6','In_MEIS2','In_PVALB','In_SST','In_TH','In_VIP','MG',
  'Neuron_MB','Oligo','OPC','VLMC')))
global_colors_full <- unique(c(as.vector(pals::glasbey(32)), as.vector(pals::alphabet(26))))[1:length(all_types_full)]
names(global_colors_full) <- all_types_full
cols <- global_colors_full[dot_order]

# ---- measured numbers -----------------------------------------------------
cat('total cells:', ncol(obj), '
')
print(table(obj$technique)); print(table(obj$orig.ident))
cat('cell types (cell_type_updated):', length(unique(obj$cell_type_updated)), '
')
cat('COP/NFOL:', sum(obj$cell_type_updated == 'COP/NFOL'),
    sprintf('(%.1f%% of %d)', 100 * mean(obj$cell_type_updated == 'COP/NFOL'), ncol(obj)), '
')
stopifnot(all(dot_order %in% unique(obj$cell_type_updated)))

# ---- small derived cache for panels A-D -----------------------------------
red <- list(RNA = Embeddings(obj, 'umap'),
            ATAC = Embeddings(obj, 'umap.atac'),
            WNN  = Embeddings(obj, 'wnn.umap'))
m48 <- c('Qk','Slc1a2','Slc1a3','Mbp','Plp1','Gpr17','Pdgfra','Tnr','P2ry12',
         'Hexb','Flt1','Pecam1','Eya2','Foxd1','Dnah12','Adgrv1','Rbfox3','Gad1',
         'Gad2','Drd1','Drd2','Hap1','Lamp5','Lhx6','Sst','Vip','Pvalb','Th',
         'Slc17a6','Rmst','Ndst4','Foxb1','Prox1','Prkcd','Tac2','Pou4f1','Hcn1',
         'Fibcd1','Spink8','Car3','Slc17a7','Cux2','Igfn1','Rorb','Tshz2','Il1rapl2',
         'Tle4','Shisa6')
DefaultAssay(obj) <- 'RNA'
ex <- as.matrix(GetAssayData(obj, slot = 'data')[m48, , drop = FALSE])
m48_present <- intersect(m48, rownames(ex))
cat('marker genes in RNA assay:', length(m48_present), '/', length(m48), '
')
saveRDS(list(meta = obj@meta.data[, c('orig.ident', 'technique', 'cell_type_updated')],
             red = red, ex = ex, m48 = m48_present,
             dot_order = dot_order, num_map = num_map, cols = cols),
        file.path(OUT, 'gs_subset_light.RDS'))
writeLines(c(sprintf('total=%d', ncol(obj)),
             sprintf('tech=%s', paste(names(table(obj$technique)), table(obj$technique),
                                      sep = ':', collapse = ',')),
             sprintf('types=%d', length(unique(obj$cell_type_updated))),
             sprintf('copnfol=%d', sum(obj$cell_type_updated == 'COP/NFOL')),
             sprintf('markers=%d', length(m48_present))),
           file.path(OUT, 'summary.txt'))
cat('wrote gs_subset_light.RDS + summary.txt
')
