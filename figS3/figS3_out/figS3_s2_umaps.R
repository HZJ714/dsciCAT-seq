
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(ggrepel); library(grid)})
OUT <- Sys.getenv('FIGS3_OUT', 'figS3_out')
S   <- readRDS(file.path(OUT, 'S3_light.RDS')); libs <- S$libs; cols <- S$cols
pt  <- c(EL14 = 0.07, FY13 = 0.055, FY14 = 0.05)
ttl <- c(EL14 = 'HT-scCAT-seq (n = 4,029)', FY13 = 'dsciCAT-seq-50k (n = 17,093)',
         FY14 = 'dsciCAT-seq-100k (n = 30,350)')
norm2 <- function(coords, span = 2.4) {
  x <- coords[, 1]; y <- coords[, 2]; s <- max(diff(range(x)), diff(range(y)))
  data.frame(x = (x - mean(range(x))) / s * span, y = (y - mean(range(y))) / s * span)
}
panel <- function(df, t) {
  lab <- df %>% group_by(cell_type) %>% summarise(x = median(x), y = median(y), .groups = 'drop')
  ggplot(df, aes(x, y)) +
    geom_point(aes(color = cell_type), size = df$pt[1], stroke = 0, alpha = 0.95) +
    scale_color_manual(values = cols, guide = 'none') +
    geom_text_repel(data = lab, aes(x, y, label = cell_type, color = cell_type),
                    size = 5 / 2.845276, box.padding = unit(0.25, 'lines'),
                    point.padding = unit(0.2, 'lines'), min.segment.length = 0.1,
                    segment.size = 0.3, segment.alpha = 0.7, max.overlaps = Inf,
                    force = 4, force_pull = 1.2, bg.color = 'white', bg.r = 0.15, seed = 42) +
    coord_fixed() + xlim(-1.7, 1.7) + ylim(-1.7, 1.7) + theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0)) +
    annotate('text', x = -1.68, y = 1.66, label = t, size = 5 / 2.845276, hjust = 0)
}
ps <- list()
for (emb in c('rna', 'atac')) for (nm in names(libs)) {
  df <- norm2(libs[[nm]][[emb]]); df$cell_type <- libs[[nm]]$ct; df$pt <- pt[[nm]]
  ps[[length(ps) + 1]] <- panel(df, paste(ttl[[nm]], toupper(emb)))
}
draw <- function(g, x, y, w, h) { pushViewport(viewport(x = x, y = y, width = w, height = h,
  just = c('left', 'bottom'), default.units = 'cm')); grid.draw(g); popViewport() }
for (dev in c('pdf', 'png')) {
  if (dev == 'pdf') pdf(file.path(OUT, 'figS3ab_umaps.pdf'), width = 15 / 2.54, height = 10 / 2.54)
  else png(file.path(OUT, 'figS3ab_umaps.png'), width = 15, height = 10, units = 'cm', res = 300)
  grid.newpage()
  for (i in seq_along(ps)) draw(ggplotGrob(ps[[i]]), ((i - 1) %% 3) * 5,
                                (1 - (i - 1) %/% 3) * 5, 5, 5)
  dev.off()
}
cat('wrote figS3ab_umaps.pdf + .png (2 rows x 3 columns, 5 x 5 cm each)\n')
