
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr)
  library(tidyr); library(pals)})
OUT  <- Sys.getenv('FIG2_OUT', 'fig2_out')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
S <- readRDS(file.path(OUT, 'gs_subset_light.RDS'))
meta <- S$meta; dot_order <- S$dot_order; num_map <- S$num_map; cols <- S$cols

# number ink colour: dark type colour -> white digits, light -> black
lum <- apply(col2rgb(cols), 2, function(v) 0.299*v[1] + 0.587*v[2] + 0.114*v[3])
num_col <- ifelse(lum < 0.5, 'white', 'black'); names(num_col) <- dot_order
leg_lab <- setNames(paste0(num_map[dot_order], ' ', dot_order), dot_order)

reductions <- list(rna = 'RNA', atac = 'ATAC', wnn = 'WNN')
for (nm in names(reductions)) {
  e <- as.data.frame(S$red[[reductions[[nm]]]]); colnames(e) <- c('UMAP_1', 'UMAP_2')
  e$cell_type <- factor(meta$cell_type_updated, levels = dot_order)
  ctr <- e %>% group_by(cell_type) %>%
    summarise(UMAP_1 = mean(UMAP_1), UMAP_2 = mean(UMAP_2), .groups = 'drop')
  ctr$num <- num_map[as.character(ctr$cell_type)]; ctr$txt_col <- num_col[as.character(ctr$cell_type)]
  p <- ggplot(e, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(color = cell_type), size = 0.06, alpha = 0.5) +
    scale_color_manual(values = cols, labels = leg_lab, name = 'Cell Type (n = number)') +
    geom_text(data = ctr, aes(UMAP_1, UMAP_2, label = num), color = 'white', size = 6.6/2.845) +
    geom_text(data = ctr[ctr$txt_col == 'black', ], aes(UMAP_1, UMAP_2, label = num),
              color = 'black', size = 6/2.845) +
    geom_text(data = ctr[ctr$txt_col == 'white', ], aes(UMAP_1, UMAP_2, label = num),
              color = 'white', size = 6/2.845) +
    guides(color = guide_legend(ncol = 2, keywidth = 0.18, keyheight = 0.18,
                                override.aes = list(size = 1.6, alpha = 1))) +
    theme_bw(base_size = 9) +
    theme(legend.position = 'right', legend.text = element_text(size = 6),
          legend.title = element_text(size = 7), legend.key = element_blank(),
          legend.spacing.y = unit(0.02, 'cm'), panel.grid = element_blank())
  ggsave(file.path(OUT, sprintf('fig2_umap_%s_celltype.pdf', nm)), p, width = 18, height = 12, units = 'cm')
  ggsave(file.path(OUT, sprintf('fig2_umap_%s_celltype.png', nm)), p, width = 18, height = 12, units = 'cm', dpi = 150)
  cat('saved fig2_umap_%s_celltype.pdf
', nm)
}
cat('order:', paste(paste0(num_map, names(num_map)), collapse = ' '), '
')
