
suppressPackageStartupMessages({library(Seurat); library(Signac); library(Matrix)
  library(readr); library(magrittr); library(ggplot2); library(patchwork)
  library(GenomicRanges); library(EnsDb.Hsapiens.v86)})
DATA   <- Sys.getenv("FIG1_DATA", "")
OUT    <- Sys.getenv("FIG1_OUT", "fig1_out")
REGION <- "chr13-91340000-91370000"   # hg38, Mir17hg. NIH3T3: "chr11-69090000-69130000" (mm10)
libs <- list(                       # library -> (fragments.tsv.gz, per-species peak matrix dir)
  "dsciCAT-seq-50k"  = c("HG/HG12/sub/human_fragments.tsv.gz",  "HG/HG12/sub/ATAC/Human"),
  "dsciCAT-seq-100k" = c("HG/HG14/sub/human_fragments.tsv.gz",  "HG/HG14/sub/ATAC/Human"),
  "HT-scCAT-seq"     = c("scCAT-seq-HT/IL9/ATAC-sub/human_fragments.tsv.gz",
                         "scCAT-seq-HT/IL9/ATAC-sub/Human"))
col <- c("dsciCAT-seq-50k" = "#E15759", "dsciCAT-seq-100k" = "#59A14F", "HT-scCAT-seq" = "#4E79A7")

ann <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)   # gene models for the annotation row
seqlevels(ann) <- paste0("chr", seqlevels(ann)); genome(ann) <- "hg38"

load_one <- function(id, frags, mtx) {                  # peak matrix -> ChromatinAssay object
  f <- as.character(read_tsv(file.path(mtx, "features.tsv"), col_names = FALSE, show_col_types = FALSE)[[1]])
  b <- as.character(read_tsv(file.path(mtx, "barcodes.tsv"), col_names = FALSE, show_col_types = FALSE)[[1]])
  X <- readMM(file.path(mtx, "matrix.mtx"))
  if (nrow(X) == length(b)) X <- t(X)                   # stored cells x peaks -> peaks x cells
  pk  <- GenomicRanges::GRanges(f, sep = c(":", "-"))   # names are chr1:100000-100200
  sel <- which(width(pk) > 100 & width(pk) < 5000)      # the peak-width filter of the original
  ca <- CreateChromatinAssay(X[sel, , drop = FALSE] %>% set_rownames(f[sel]) %>% set_colnames(b),
        sep = c(":", "-"), fragments = frags, min.cells = 10, annotation = ann)
  o <- CreateSeuratObject(ca, assay = "ATAC"); o$orig.ident <- id
  DefaultAssay(o) <- "ATAC"; Idents(o) <- "orig.ident"
  cat(id, ":", ncol(o), "cells,", nrow(o), "peaks
")
  o
}

# The three libraries are plotted separately, not merged: CoveragePlot/TilePlot read
# the fragments file, not the assay, so merging would only pay for a shared peak set
# (BindMatrices over ~780k peak names) that no row of this panel uses.
# Signac draws its row titles ("Normalized signal (range 0 - 130)", "Genes") as
# annotate() text and puts the real axes on inner plots of a patchwork tree, so the
# only way to get a stackable row is to walk the tree and blank them. The
# patchwork branch has to be tested first: a patchwork is also a ggplot.
kill <- theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
              axis.title.x = element_blank(), axis.text.y = element_blank(),
              axis.ticks.y = element_blank(), axis.title.y = element_blank(),
              strip.text = element_blank(), plot.margin = margin(1, 1, 0, 1),
              legend.position = "none")
clean <- function(x) {
  if (inherits(x, "patchwork")) {
    for (i in seq_along(x)) x[[i]] <- tryCatch(clean(x[[i]]), error = function(e) x[[i]])
    return(x)
  }
  if (!inherits(x, "ggplot")) return(x)
  keep <- vapply(seq_along(x$layers), function(i) {
    d <- tryCatch(ggplot2::layer_data(x, i), error = function(e) NULL)
    if (is.null(d)) return(TRUE)
    bad <- FALSE
    for (cc in d) if (is.character(cc) && any(grepl("Normalized signal|^Genes$", cc))) bad <- TRUE
    !bad
  }, TRUE)
  x$layers <- x$layers[keep]
  x + kill
}
# CoveragePlot maps the coverage area to fill = group and then installs its own
# fill scale, so a scale_fill_manual() added here is silently replaced on the row
# that also carries the gene models. aes_params wins over the mapping, so the
# colour is set on the GeomArea layer directly.
setfill <- function(x, colour) {
  if (inherits(x, "patchwork")) {
    for (i in seq_along(x)) x[[i]] <- tryCatch(setfill(x[[i]], colour), error = function(e) x[[i]])
    return(x)
  }
  if (!inherits(x, "ggplot")) return(x)
  for (i in seq_along(x$layers))
    if (inherits(x$layers[[i]]$geom, "GeomArea")) x$layers[[i]]$aes_params$fill <- colour
  x
}
last <- tail(names(libs), 1)
cov <- list(); til <- list()
for (id in names(libs)) {
  o <- load_one(id, file.path(DATA, libs[[id]][1]), file.path(DATA, libs[[id]][2]))
  # the gene models ride on the last coverage row (annotation = TRUE)
  cov[[id]] <- setfill(clean(CoveragePlot(object = o, region = REGION, annotation = (id == last),
                                          peaks = FALSE, extend.upstream = 0,
                                          extend.downstream = 0)), col[id])
  til[[id]] <- clean(TilePlot(object = o, region = REGION, tile.cells = 300, order.by = "top") +
                     scale_fill_gradientn(colours = c("white", "blue", "darkblue"),
                                          limits = c(0, 2), breaks = c(0, 1, 2)))
}
til[[last]] <- til[[last]] + theme(legend.position = "right")   # one shared key

fig <- wrap_plots(c(cov, til), ncol = 1, heights = c(2, 2, 2.8, 2, 2, 2)) +
  plot_annotation(title = "chr13:91,340,000-91,370,000 (hg38, MIR17HG)",
                  theme = theme(plot.title = element_text(size = 6.5, hjust = 0,
                                                          margin = margin(0, 0, 2, 0))))
ggsave(file.path(OUT, "fig1c_atac_tracks_HEK293T.pdf"), fig, width = 9, height = 13, units = "cm")
ggsave(file.path(OUT, "fig1c_atac_tracks_HEK293T.png"), fig, width = 9, height = 13, units = "cm", dpi = 300)
cat("wrote", file.path(OUT, "fig1c_atac_tracks_HEK293T.pdf"), "
")
