
cat('R', R.version.string, '\n')
pkgs <- c('Seurat','SeuratObject','harmony','ggplot2','dplyr','tidyr','stringr',
          'pals','parallel','aricode','ggtext','Matrix','Rcpp')
for (p in pkgs) cat(sprintf('%-14s %s\n', p,
    tryCatch(as.character(packageVersion(p)), error = function(e) 'NOT INSTALLED')))
