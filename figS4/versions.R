DATA <- Sys.getenv("DATA"); OUT <- Sys.getenv("OUT")

ip <- installed.packages()
for (p in c("scMEGA","ArchR","Seurat","Signac")) {
  if (p %in% rownames(ip)) cat(sprintf("  %-22s %s\n", p, ip[p, "Version"]))
}
cat("  R ", R.version.string, "\n", sep="")
