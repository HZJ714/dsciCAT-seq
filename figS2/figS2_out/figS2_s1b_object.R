
suppressPackageStartupMessages({library(Seurat); library(pals)})
DATA <- Sys.getenv("S2_DATA", Sys.getenv("DSCI_FY_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "FY")))
OUT  <- Sys.getenv('S2_OUT', 'figS2_out')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
t0 <- proc.time()[3]
obj <- readRDS(file.path(DATA, '3k_dsci_ht_10x_gs_subset.RDS'))
cat('loaded in', round((proc.time()[3] - t0) / 60, 1), 'min |', ncol(obj), 'cells\n')

md <- obj@meta.data
md$cell_type_updated <- as.character(md$cell_type)
md$cell_type_updated[md$RNA_snn_res.1 == 34] <- 'COP/NFOL'      # 280 cells
md$technique <- ifelse(md$orig.ident %in% c('dsciCAT-seq-50k', 'dsciCAT-seq-100k'),
                       'dsciCAT-seq', as.character(md$orig.ident))
cat('libraries:\n'); print(table(md$orig.ident))
cat('techniques:\n'); print(table(md$technique))
cat('types:', length(unique(md$cell_type_updated)),
    '| COP/NFOL:', sum(md$cell_type_updated == 'COP/NFOL'),
    sprintf('(%.2f%%)', 100 * mean(md$cell_type_updated == 'COP/NFOL')), '\n')
stopifnot(ncol(obj) == 70643, length(unique(md$cell_type_updated)) == 33)

# fixed 33-type order + numbering + shared palette (= Fig 2's, = track/dot plot)
dot_order <- c('AST','Oligo','COP/NFOL','OPC','MG','EC','VLMC','Ependyma',
               'In_D1','In_D2','In_HAP1','In_LAMP5','In_LHX6','In_SST','In_VIP',
               'In_PVALB','In_TH','Neuron_MB','Ex_LHA','Ex_HY-MM','Ex_DG','Ex_TH',
               'Ex_MH-LH','Ex_SUB_Pros','Ex_CA1_Pros','Ex_CAR3','Ex_L2/3_IT',
               'Ex_L2/3_IT_PIR','Ex_L4/5_IT','Ex_L5/6_NP','Ex_L6_IT','Ex_L6_CT','Ex_L6b')
num_map <- setNames(seq_along(dot_order), dot_order)
all_types_full <- sort(unique(c(
  'AST','COP/NFOL','EC','Ependyma','Ex_CA1_Pros','Ex_CA3','Ex_CAR3','Ex_CR',
  'Ex_DG','Ex_HY-MM','Ex_IC','Ex_L2/3_IT','Ex_L2/3_IT_PIR','Ex_L4','Ex_L4/5_IT',
  'Ex_L5_ET','Ex_L5_IT','Ex_L5_NP','Ex_L5/6_NP','Ex_L6_CT','Ex_L6_IT','Ex_L6b',
  'Ex_LHA','Ex_MH-LH','Ex_SUB_Pros','Ex_TH','In_D1','In_D2','In_D1-D2','In_HAP1',
  'In_LAMP5','In_LHX6','In_MEIS2','In_PVALB','In_SST','In_TH','In_VIP','MG',
  'Neuron_MB','Oligo','OPC','VLMC')))
global_colors_full <- unique(c(as.vector(pals::glasbey(32)),
                               as.vector(pals::alphabet(26))))[1:length(all_types_full)]
names(global_colors_full) <- all_types_full
stopifnot(all(dot_order %in% unique(md$cell_type_updated)))

saveRDS(list(wnn = SeuratObject::Embeddings(obj, 'wnn.umap'), meta = md,
             dot_order = dot_order, num_map = num_map,
             cols = global_colors_full[dot_order]),
        file.path(OUT, 's2_light.RDS'))
cat('wrote s2_light.RDS\n')
