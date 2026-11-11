DATA <- Sys.getenv("DATA"); OUT <- Sys.getenv("OUT")

suppressPackageStartupMessages({library(scMEGA); library(ArchR)})
pal <- paletteContinuous(set="blueYellow", n=100)
write.table(pal, file.path(OUT, "blueYellow.csv"),
            sep=",", col.names=FALSE, row.names=FALSE)
cat("ramp rows:", length(pal), "\n")
