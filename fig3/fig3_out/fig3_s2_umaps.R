
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(dplyr)
  library(tidyr); library(ggrepel); library(grid)})
OUT <- Sys.getenv('FIG3_OUT', 'fig3_out')
S   <- readRDS(file.path(OUT, 'fig3_light.RDS')); libs <- S$libs; cols <- S$cols
titles <- c(EL14 = 'HT-scCAT-seq (n = 4,029) WNN', FY13 = 'dsciCAT-seq-50k (n = 17,093) WNN',
            FY14 = 'dsciCAT-seq-100k (n = 30,350) WNN')
pt    <- c(EL14 = 0.07, FY13 = 0.055, FY14 = 0.05)

# coordinates are drawn on a square canvas: normalise each library to the same span
norm2 <- function(coords, span = 2.4) {
  x <- coords[, 1]; y <- coords[, 2]; s <- max(diff(range(x)), diff(range(y)))
  data.frame(x = (x - mean(range(x))) / s * span, y = (y - mean(range(y))) / s * span)
}
make_panel <- function(df, letter, ttl, mark = NULL) {
  lab <- df %>% group_by(cell_type) %>% summarise(x = median(x), y = median(y), .groups = 'drop')
  ggplot(df, aes(x, y)) +
    geom_point(aes(color = cell_type), size = pt[[df$lib[1]]], stroke = 0, alpha = 0.95) +
    scale_color_manual(values = cols, guide = 'none') +
    geom_text_repel(data = lab, aes(x, y, label = cell_type, color = cell_type),
                    size = 5 / 2.845276, fontface = 'bold', box.padding = unit(0.25, 'lines'),
                    point.padding = unit(0.2, 'lines'), min.segment.length = 0.1,
                    segment.size = 0.3, segment.alpha = 0.7, max.overlaps = Inf,
                    force = 4, force_pull = 1.2, bg.color = 'white', bg.r = 0.15, seed = 42) +
    coord_fixed() + xlim(-1.7, 1.7) + ylim(-1.7, 1.7) + theme_void() +
    theme(plot.margin = margin(0, 0, 0, 0)) +
    annotate('text', x = -1.6, y = 1.62, label = letter, size = 6 / 2.845276, hjust = 0) +
    annotate('text', x = -1.42, y = 1.62, label = ttl, size = 5 / 2.845276, hjust = 0) +
    if (!is.null(mark)) geom_point(data = mark, aes(x, y), inherit.aes = FALSE, shape = 21,
                                   colour = 'red', fill = NA, size = 1.6, linewidth = 0.25)
}
ps <- list(); hull <- NULL
for (i in seq_along(libs)) {
  nm <- names(libs)[i]
  df <- norm2(libs[[nm]]$wnn); df$cell_type <- libs[[nm]]$ct; df$lib <- nm
  g <- df[df$cell_type == 'COP/NFOL', ]            # outline of the COP/NFOL cluster
  k <- if (nrow(g) > 3) grDevices::chull(g$x, g$y) else NULL
  if (!is.null(k)) hull <- rbind(hull, data.frame(lib = nm, x = g$x[k], y = g$y[k],
                                                  cell = libs[[nm]]$cells[df$cell_type == 'COP/NFOL'][k],
                                                  n = nrow(g)))
  ps[[i]] <- make_panel(df, LETTERS[i], titles[[nm]], if (is.null(k)) NULL else g)
}
write.csv(hull, file.path(OUT, 'fig3_copnfol_hull.csv'), row.names = FALSE)
cat('COP/NFOL cells: FY13', sum(libs$FY13$ct == 'COP/NFOL'), '| FY14',
    sum(libs$FY14$ct == 'COP/NFOL'), '| EL14', sum(libs$EL14$ct == 'COP/NFOL'), '\n')
draw <- function(g, x, y, w, h) { pushViewport(viewport(x = x, y = y, width = w, height = h,
  just = c('left', 'bottom'), default.units = 'cm')); grid.draw(g); popViewport() }
pdf(file.path(OUT, 'fig3abc_wnn_umaps.pdf'), width = 15 / 2.54, height = 5 / 2.54)
grid.newpage()
for (i in 1:3) draw(ggplotGrob(ps[[i]]), (i - 1) * 5, 0, 5, 5)
dev.off()
png(file.path(OUT, 'fig3abc_wnn_umaps.png'), width = 15, height = 5, units = 'cm', res = 300)
grid.newpage()
for (i in 1:3) draw(ggplotGrob(ps[[i]]), (i - 1) * 5, 0, 5, 5)
dev.off()
cat('wrote fig3abc_wnn_umaps.pdf\n')
