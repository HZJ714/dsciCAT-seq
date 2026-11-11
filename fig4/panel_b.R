DATA <- Sys.getenv("DATA"); OUT <- Sys.getenv("OUT")

suppressPackageStartupMessages({library(Seurat); library(Signac); library(scMEGA)
                                library(ArchR); library(SummarizedExperiment)
                                library(ComplexHeatmap); library(readxl)})
LayerData <- function(object, assay=NULL, layer="counts", cells=NULL) {
  mat <- GetAssayData(object, assay=assay, slot=layer)
  if (!is.null(cells)) mat <- mat[, cells]; mat
}
JoinLayers <- function(object, ...) object
obj.t  <- readRDS(file.path(DATA, "scmega_obj_target.rds"))
df.p2g <- read.csv(file.path(DATA, "scmega_p2g.csv"))
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
# ---- smooth both assays over the same pseudotime axis (100 bins) ----
trajRNA  <- GetTrajectory(obj.t, assay="RNA", trajectory.name="Trajectory",
                          groupEvery=1, slot="data", smoothWindow=7, log2Norm=TRUE)
trajATAC <- GetTrajectory(obj.t, assay="ATAC", trajectory.name="Trajectory",
                          groupEvery=1, slot="data", smoothWindow=7, log2Norm=TRUE)
trajATAC <- trajATAC[df.p2g$peak, ]
trajRNA  <- trajRNA[df.p2g$gene, ]
# ---- row order: combined z-score matrix, TrajectoryHeatmap's row dendrogram ----
trajCombined <- trajATAC
assay(trajCombined, withDimnames=FALSE) <-
  t(apply(assay(trajRNA), 1, scale)) + t(apply(assay(trajATAC), 1, scale))
combinedMat <- TrajectoryHeatmap(trajCombined, returnMatrix=TRUE, varCutOff=0)
rowOrder <- match(rownames(combinedMat), rownames(trajATAC))
# ---- left heatmap: chromatin accessibility (solarExtra, z clipped at +-2) ----
ht1 <- TrajectoryHeatmap(trajATAC, pal=paletteContinuous(set="solarExtra"),
                         varCutOff=0, maxFeatures=nrow(trajATAC), rowOrder=rowOrder,
                         limits=c(-2, 2), labelRows=FALSE, labelTop=0,
                         name="Chromatin accessibility")
# ---- right heatmap: gene expression, same row order, z clipped at +-2 ----
mat2 <- assay(trajRNA)[rowOrder, , drop=FALSE]
mat2 <- sweep(mat2 - rowMeans(mat2), 1, matrixStats::rowSds(mat2), `/`)
mat2[mat2 > 2] <- 2; mat2[mat2 < -2] <- -2
# ---- 14 stage markers (anti-overlap: DFS over occurrences, min row gap 200) ----
mk <- intersect(c("Bmp4","Mbp","Plp1","Mag","Mog","Gpr17","Olig1","Olig2",
                  "Tcf7l2","Car2","Opalin","Mal","Mobp","Apod"), df.p2g$gene)
mk <- mk[mk %in% rownames(trajRNA)]
genes2 <- rownames(trajRNA)[rowOrder]
occ <- lapply(mk, function(g) which(genes2 == g)); names(occ) <- mk
first_pos <- sapply(occ, `[`, 1); target <- first_pos
target["Mbp"] <- 14792; target["Mag"] <- 14094; target["Plp1"] <- 11107
g_order <- mk[order(first_pos)]
MIN_GAP <- 200; chosen_pos <- integer(); chosen_g <- character()
dfs <- function(i) {
  if (i > length(g_order)) return(TRUE)
  g <- g_order[i]
  for (cp in occ[[g]][order(abs(occ[[g]] - target[g]))]) {
    if (all(abs(cp - chosen_pos) >= MIN_GAP)) {
      chosen_pos <<- c(chosen_pos, cp); chosen_g <<- c(chosen_g, g)
      if (dfs(i + 1)) return(TRUE)
      chosen_pos <<- chosen_pos[-length(chosen_pos)]
      chosen_g <<- chosen_g[-length(chosen_g)]
    }
  }
  FALSE
}
dfs(1)
lbl <- rep("", length(genes2)); lbl[chosen_pos] <- chosen_g
# ---- colourbars + pseudotime annotation (blueYellow, same ramp as panel A) ----
cols2 <- colorRamp2(seq(-2, 2, length.out=100),
                    paletteContinuous(set="horizonExtra", n=100))
pseudo_cols <- paletteContinuous(set="blueYellow", n=100)
pseudo_fun  <- colorRamp2(seq(0, 100, length.out=100), pseudo_cols)
pseudo_lab  <- rep("", 100); pseudo_lab[c(1,26,51,76,100)] <- c("0","25","50","75","100")
bottom_anno <- HeatmapAnnotation(
  "Pseudotime" = anno_simple((0:99) + 0.5, col=pseudo_fun, border=FALSE,
                             height=unit(0.8, "cm")),
  "tick" = anno_text(pseudo_lab, location=0.5,
                     gp=gpar(fontsize=9, col="grey30")),
  which="column", annotation_name_gp=gpar(fontsize=11),
  show_legend=c("Pseudotime"=FALSE),
  show_annotation_name=c("Pseudotime"=TRUE, "tick"=FALSE))
ht2 <- Heatmap(mat2, name="Gene expression", col=cols2,
               cluster_rows=FALSE, cluster_columns=FALSE, show_column_names=FALSE,
               row_labels=lbl, row_names_gp=gpar(fontsize=12, fontface="bold.italic"),
               row_names_side="right",
               heatmap_legend_param=list(title="Gene expression"),
               bottom_annotation=bottom_anno)
ht1@bottom_annotation <- bottom_anno
# vertical 05f-equivalent (12x16 in), identical to make_05_heatmaps.R
pdf(file.path(OUT, "fig4b_05f.pdf"), width=12, height=16)
draw(ht1 + ht2); dev.off()
# marker positions (row indices in the heatmap row order) for panel layout
write.csv(data.frame(gene=chosen_g, pos=chosen_pos),
          file.path(OUT, "fig4b_labels.csv"), row.names=FALSE)
cat("B done\n")
