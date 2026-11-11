
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr)})
OUT <- Sys.getenv('FIG3_OUT', 'fig3_out')
S   <- readRDS(file.path(OUT, 'fig3_light.RDS')); libs <- S$libs
# common colour limit: max Gpr17 over the inset cells of all three libraries, rounded up
max_gpr <- ceiling(max(sapply(names(libs), function(k) max(libs[[k]]$gpr17[match(libs[[k]]$inset_cells,
                                                                                 libs[[k]]$cells)]))) * 10) / 10
cat('common Gpr17 scale: 0 -', max_gpr, '\n')
ps <- list(); info <- character()
for (nm in names(libs)) {
  s <- libs[[nm]]
  sub <- CreateSeuratObject(s$inset_counts, project = nm); DefaultAssay(sub) <- 'RNA'
  sub <- NormalizeData(sub, verbose = FALSE)
  sub <- suppressWarnings(FindVariableFeatures(sub, nfeatures = min(3000, nrow(sub)), verbose = FALSE))
  sub <- ScaleData(sub, verbose = FALSE)
  sub <- RunPCA(sub, npcs = 30, verbose = FALSE)
  sub <- FindNeighbors(sub, reduction = 'pca', dims = 1:10, verbose = FALSE)
  sub <- FindClusters(sub, algorithm = 3, resolution = c(0.1, 0.3, 0.5), verbose = FALSE)
  sub <- RunUMAP(sub, reduction = 'pca', dims = 1:10, seed.use = 42, verbose = FALSE)
  tb <- sort(table(sub$RNA_snn_res.0.5), decreasing = TRUE)
  info <- c(info, sprintf('%s n=%d -> %d cluster(s) at res 0.5 (%s), Gpr17 mean %.2f', nm,
           ncol(sub), length(tb), paste(tb, collapse = ','),
           mean(GetAssayData(sub, assay = 'RNA', slot = 'data')['Gpr17', ])))
  df <- data.frame(x = Embeddings(sub, 'umap')[, 1], y = Embeddings(sub, 'umap')[, 2],
                   gpr = as.numeric(GetAssayData(sub, assay = 'RNA', slot = 'data')['Gpr17', ]))
  ps[[nm]] <- ggplot(df, aes(x, y)) +
    geom_point(aes(color = gpr), size = 1.6, alpha = 0.95) +
    scale_color_gradientn(colours = c('grey80', '#ffcccc', '#ff0000', '#b30000'),
                          limits = c(0, max_gpr), oob = scales::squish, name = 'Gpr17 (log-norm)') +
    labs(title = sprintf('%s (n = %d, re-clustered)', nm, nrow(df))) +
    coord_fixed() + theme_void() +
    theme(plot.title = element_text(size = 6, hjust = 0.5),
          legend.key.size = unit(2.2, 'mm'), legend.text = element_text(size = 5),
          legend.title = element_text(size = 5.5))
}
writeLines(info, file.path(OUT, 'fig3_insets_summary.txt')); cat(paste(info, collapse = '\n'), '\n')
fig <- patchwork::wrap_plots(ps, ncol = 3)
ggsave(file.path(OUT, 'fig3bc_gpr17_insets.pdf'), fig, width = 15, height = 4.2, units = 'cm')
ggsave(file.path(OUT, 'fig3bc_gpr17_insets.png'), fig, width = 15, height = 4.2, units = 'cm', dpi = 300)
cat('wrote fig3bc_gpr17_insets.pdf\n')
