
suppressMessages({library(Seurat); library(Signac)})
IN  <- Sys.getenv("FIG2_INPUT", Sys.getenv("FIG1_DATA", ""))
FY  <- Sys.getenv("FIG2_DATA", Sys.getenv("DSCI_FY_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "FY")))
OUT <- Sys.getenv('FIG2_OUT', 'fig2_out')
say <- function(...) cat(sprintf(...), '\n')

## What §4 and the legend need from §3 is one thing: the candidate cell list, one
## row per cell, with the library it came from. Assembling the *object* behind
## that list - merge() over 83,242 cells x ~1e5 peaks - is the most expensive step
## in this notebook and feeds nothing else, so it is built only when the chain
## really has to be re-run; the file.exists() rule below says so out loud.
libs <- list(FY13 = file.path(FY, 'FY13', '3k_FY13_multiome.RDS'),
             FY14 = file.path(FY, 'FY14', '3k_FY14_multiome.RDS'))
dsci <- do.call(rbind, lapply(names(libs), function(s)
  data.frame(cell = colnames(readRDS(libs[[s]])), sample = s)))
auth <- readRDS(file.path(IN, 'scCAT-seq-HT/Mbrain_process.rds'))
acell <- data.frame(cell = colnames(auth),
                    sample = ifelse(grepl('^Cortex_EL', colnames(auth)), 'HT-scCAT-seq-HT', '10x'))
cand <- rbind(dsci, acell)
say('candidates: %d cells = %d dsciCAT-seq + %d HT/10x',
    nrow(cand), nrow(dsci), nrow(acell))
print(table(cand$sample))
stopifnot(nrow(cand) == 83242, !anyDuplicated(cand$cell))
saveRDS(cand, file.path(OUT, 'fig2_candidate_cells.RDS'))
rm(dsci, auth); invisible(gc())

## the merged object, built only if something downstream still needs it
out.rds <- file.path(FY, '3k_dsci_ht_10x_merged.RDS')
gs.rds  <- file.path(FY, '3k_dsci_ht_10x_gs_subset.RDS')
if (file.exists(out.rds)) {
  say('merged object present: %s (%.1f GB)', basename(out.rds), file.info(out.rds)$size / 1e9)
} else if (file.exists(gs.rds)) {
  say('merged object not rebuilt: %s already exists and nothing downstream reads the intermediate',
      basename(gs.rds))
} else {
  ## merge() refuses objects that carry reductions, graphs or neighbours and
  ## silently mixes them when it accepts, so every slot is cleared first. Cell
  ## names must come through unchanged - they already start with FY13_ / FY14_ and
  ## §4 matches them against the fragment files - so no add.cell.ids.
  obj <- lapply(libs, function(f) {
    o <- readRDS(f); o@graphs <- list(); o@neighbors <- list(); o@reductions <- list(); o })
  multi <- merge(obj[[1]], y = obj[2], add.cell.ids = NULL)
  say('dsciCAT-seq merged: %d cells', ncol(multi))
  auth <- readRDS(file.path(IN, 'scCAT-seq-HT/Mbrain_process.rds'))
  auth@graphs <- list(); auth@neighbors <- list(); auth@reductions <- list()
  auth$orig.ident <- ifelse(grepl('^Cortex_EL', colnames(auth)), 'HT-scCAT-seq-HT', '10x')
  DefaultAssay(multi) <- 'RNA'; DefaultAssay(auth) <- 'RNA'
  merged <- merge(auth, y = multi, add.cell.ids = NULL)
  stopifnot(ncol(merged) == 83242, !anyDuplicated(colnames(merged)))
  saveRDS(merged, out.rds)
}
