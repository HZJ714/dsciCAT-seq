
suppressPackageStartupMessages({library(Seurat); library(ggplot2); library(grid)})
OUT  <- Sys.getenv('FIG2_OUT', 'fig2_out')
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
S <- readRDS(file.path(OUT, 'gs_subset_light.RDS'))
obj <- CreateSeuratObject(S$ex); obj@meta.data$cell_type_updated <- S$meta$cell_type_updated
Idents(obj) <- 'cell_type_updated'
present <- S$m48
lv <- levels(factor(S$meta$cell_type_updated, levels = S$dot_order))

# Panel geometry below was measured off the delivered D panel (its source,
# 'Figure S2-dot-track.pdf', is on this server; its own build script is not).
# dot.scale 16 (the value in the analysis notebook) draws 12 mm dots here, which
# is what the delivered panel does NOT look like: its largest dot is 1.37 mm.
PANEL_W <- 8.90   # cm -> 48 columns at 0.186 cm pitch (delivered: 0.186)
PANEL_H <- 7.48   # cm -> 33 rows at 0.227 cm pitch     (delivered: 0.227)
DOT_SCALE <- 1.5  # mm at 100% expressed (delivered largest dot: 1.37 mm)

p <- DotPlot(obj, group.by = 'cell_type_updated', features = present,
             scale = TRUE, cluster.idents = FALSE, col.min = 0, col.max = 3,
             dot.scale = DOT_SCALE) +
  theme_classic(base_size = 6) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, face = 'italic', size = 6),
        axis.text.y = element_text(size = 6),
        axis.title  = element_text(size = 8, face = 'bold'),
        legend.key.width = unit(0.22, 'cm'), legend.key.height = unit(0.30, 'cm'),
        legend.key.size  = unit(0.22, 'cm'), legend.spacing = unit(0, 'cm'),
        legend.position = 'top', legend.box = 'horizontal', legend.box.just = 'left',
        legend.box.spacing = unit(0.25, 'cm'),
        legend.margin = margin(0, 0, 0, 0), legend.box.margin = margin(0, 0, 0, 0),
        plot.title = element_blank(),
        plot.margin = unit(c(0.1, 0.7, 0.15, 0.1), 'cm')) +
  scale_color_gradientn(colours = c(rgb(252/255, 247/255, 243/255),
                                    rgb(253/255, 200/255, 180/255),
                                    rgb(229/255, 51/255, 38/255),
                                    rgb(102/255, 0/255, 13/255))) +
  # Seurat builds these two guides itself; guide-level themes are the only ones
  # that survive the complete theme it attaches (cowplot's legend.text = 8 pt).
  guides(size = guide_legend(title = 'Percent Expressed', order = 1,
                            title.theme = element_text(size = 7.6),
                            label.theme = element_text(size = 7.6)),
         colour = guide_colorbar(title = 'Average Expression', order = 2,
                                 barwidth = unit(1.9, 'cm'), barheight = unit(0.28, 'cm'),
                                 title.theme = element_text(size = 7.6,
                                             margin = margin(l = 0.35, unit = 'cm')),
                                 label.theme = element_text(size = 7.6))) +
  labs(x = 'Gene expression', y = 'Celltypes') +
  scale_y_discrete(limits = rev(lv)) + scale_x_discrete(limits = present)

# Fix the panel box to the measured size, then let the page shrink to the content.
g <- ggplotGrob(p)
pc <- which(sapply(g$grobs, function(x) grepl('^panel', x$name)))
g$widths[[g$layout$l[pc]]]  <- unit(PANEL_W, 'cm')
g$heights[[g$layout$t[pc]]] <- unit(PANEL_H, 'cm')
w <- as.numeric(convertWidth(sum(g$widths), 'cm'))
h <- as.numeric(convertHeight(sum(g$heights), 'cm'))
cat(sprintf('page %.2f x %.2f cm
', w, h))
pdf(file.path(OUT, 'fig2_full_celltype_dotplot.pdf'), width = w/2.54, height = h/2.54)
grid.newpage(); grid.draw(g); dev.off()
png(file.path(OUT, 'fig2_full_celltype_dotplot.png'), width = w/2.54, height = h/2.54,
    units = 'in', res = 200)
grid.newpage(); grid.draw(g); dev.off()
cat('dotplot genes:', length(present), '| types:', length(lv), '
')
cat(paste(present, collapse = ' '), '
')
