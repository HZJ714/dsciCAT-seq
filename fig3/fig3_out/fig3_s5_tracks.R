
suppressPackageStartupMessages({library(ArchR); library(grid)})
addArchRGenome('mm10'); addArchRThreads(threads = 8)
OLIGO <- Sys.getenv("FIG3_OLIGO", file.path(Sys.getenv("DSCI_FY_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "FY")), "oligo"))
OUT   <- Sys.getenv('FIG3_OUT', 'fig3_out')
proj  <- loadArchRProject(file.path(OLIGO, '3kdsci/archr-subset'), showLogo = FALSE)
ST6 <- c('COP/NFOL', 'MFOL', 'MOL2', 'MOL3', 'MOL3/6', 'MOL6')
stopifnot(all(ST6 %in% unique(proj$cell_type)))
cat('project cells:', nCells(proj), '| groups:',
    paste(names(table(proj$cell_type)), table(proj$cell_type), sep = '=', collapse = ' '), '\n')
p <- plotBrowserTrack(ArchRProj = proj, groupBy = 'cell_type',
                      geneSymbol = c('Bmp4', 'Synpr', 'Anln', 'Nrxn3', 'Grm3'),
                      upstream = 50000, downstream = 50000,
                      plotSummary = c('bulkTrack', 'geneTrack'),
                      features = NULL, loops = NULL, verbose = FALSE)
pdf(file.path(OUT, 'fig3f_archr_tracks.pdf'), width = 10, height = 14)
for (i in seq_along(p)) { grid.newpage(); grid.draw(p[[i]]) }
dev.off()
for (i in seq_along(p)) {
  png(file.path(OUT, sprintf('fig3f_archr_track_%d.png', i)),
      width = 10, height = 14, units = 'in', res = 150)
  grid.newpage(); grid.draw(p[[i]]); dev.off()
}
cat('wrote fig3f_archr_tracks.pdf with', length(p), 'pages\n')
