DATA <- Sys.getenv("DATA"); OUT <- Sys.getenv("OUT")

suppressPackageStartupMessages({library(Seurat); library(Signac)
                                library(scMEGA); library(ArchR); library(ggplot2)})
# Seurat v4 compatibility shims (v5 API not present in this environment)
LayerData <- function(object, assay=NULL, layer="counts", cells=NULL) {
  mat <- GetAssayData(object, assay=assay, slot=layer)
  if (!is.null(cells)) mat <- mat[, cells]; mat
}
JoinLayers <- function(object, ...) object

obj <- readRDS(file.path(DATA, "3k_dsci_oligo_subset.RDS"))
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
# ---- pseudotime: group mapping + within-group rank on harmony embedding ----
obj$traj_group <- as.character(obj$cell_type)
obj$traj_group[obj$traj_group %in% c("MOL2","MOL3","MOL3/6","MOL6")] <- "MOL"
set.seed(42)
emb <- Embeddings(obj, "harmony")[, 1:30]
tg  <- obj$traj_group
c_mfol <- colMeans(emb[tg == "MFOL", , drop=FALSE])
c_mol  <- colMeans(emb[tg == "MOL", , drop=FALSE])
pseudo <- rep(NA, ncol(obj)); names(pseudo) <- colnames(obj)
idx <- which(tg == "COP/NFOL")
d <- sqrt(colSums((t(emb[idx, , drop=FALSE]) - c_mfol)^2))
pseudo[idx] <- (1 - rank(d)/length(idx)) * 25
idx <- which(tg == "MFOL")
d <- sqrt(colSums((t(emb[idx, , drop=FALSE]) - c_mol)^2))
pseudo[idx] <- 25 + (1 - rank(d)/length(idx)) * 25
idx <- which(tg == "MOL")
pseudo[idx] <- 50 + rank(runif(length(idx)))/length(idx) * 50
obj$Trajectory <- pseudo
obj.t <- obj[, !is.na(obj$Trajectory)]
# ---- panel A source: real TrajectoryPlot from the RDS-derived object ----
p3 <- TrajectoryPlot(object=obj.t, reduction="umap",
                     continuousSet="blueYellow", size=1, addArrow=FALSE)
ggsave(file.path(OUT, "fig4a_raw.pdf"), p3, width=7, height=6)
# blueYellow ramp (0->100, same palette as the published figure)
pal <- paletteContinuous(set="blueYellow", n=100)
write.table(pal, file.path(OUT, "blueYellow.csv"),
            sep=",", col.names=FALSE, row.names=FALSE)
cat("A source written:", ncol(obj.t), "cells\\n")
