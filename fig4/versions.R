DATA <- Sys.getenv("DATA"); OUT <- Sys.getenv("OUT")

ip <- installed.packages()
for (p in c("Seurat","Signac","scMEGA","ArchR","ComplexHeatmap","SummarizedExperiment",
            "chromVAR","ggplot2","dplyr","matrixStats","readxl")) {
  if (p %in% rownames(ip)) cat(sprintf("  %-22s %s\n", p, ip[p, "Version"]))
}
cat("  R ", R.version.string, "\n", sep="")
