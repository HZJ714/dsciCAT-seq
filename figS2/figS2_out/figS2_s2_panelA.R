
suppressPackageStartupMessages({library(ggplot2)})
OUT <- Sys.getenv('S2_OUT', 'figS2_out')
S   <- readRDS(file.path(OUT, 's2_light.RDS'))
co  <- as.data.frame(S$wnn); colnames(co) <- c('UMAP_1', 'UMAP_2')
batch_order  <- c('10x', 'HT-scCAT-seq-HT', 'dsciCAT-seq-50k', 'dsciCAT-seq-100k')
batch_colors <- c('10x' = '#E64B35FF', 'HT-scCAT-seq-HT' = '#3C5488FF',
                  'dsciCAT-seq-50k' = '#00A087FF', 'dsciCAT-seq-100k' = '#4DBBD5FF')
co$Batch <- factor(as.character(S$meta$orig.ident), levels = batch_order)
stopifnot(!any(is.na(co$Batch)), nrow(co) == 70643)
cat('cells per facet:\n'); print(table(co$Batch))

p_a <- ggplot() +
  geom_point(data = co[, c('UMAP_1', 'UMAP_2')], aes(UMAP_1, UMAP_2),
             color = '#E0E0E0', size = 0.12, alpha = 0.35, stroke = 0) +
  geom_point(data = co, aes(UMAP_1, UMAP_2, color = Batch),
             size = 0.15, alpha = 0.8, stroke = 0) +
  scale_color_manual(values = batch_colors, breaks = batch_order, name = NULL) +
  facet_wrap(~Batch, ncol = 4) +
  theme_void(base_size = 7) +
  theme(strip.text = element_text(size = 6.5),
        panel.spacing = unit(0.05, 'cm'), legend.position = 'bottom',
        legend.text = element_text(size = 6),
        legend.key.height = unit(0.16, 'cm'), legend.key.width = unit(0.16, 'cm'))
ggsave(file.path(OUT, 'S2a_batch_shadow_facets.pdf'), p_a,
       width = 16.97, height = 4.97, units = 'cm', device = grDevices::cairo_pdf)
ggsave(file.path(OUT, 'S2a_batch_shadow_facets.png'), p_a,
       width = 16.97, height = 4.97, units = 'cm', dpi = 150)
cat('saved S2a_batch_shadow_facets.pdf\n')
