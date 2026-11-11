
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)})
OUT  <- Sys.getenv('S2_OUT', 'figS2_out')
FY   <- file.path(Sys.getenv("S2_DATA", Sys.getenv("DSCI_FY_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "FY"))), '3k_dsci')
LEVELS <- c('3k','5k','10k','15k','20k','25k','30k','40k')
LEVEL_N <- as.integer(gsub('k', '000', LEVELS))
read_lv <- function(lv) {
  f <- file.path(OUT, paste0(lv, '_ARI_NMI.csv'))
  if (!file.exists(f)) f <- file.path(FY, paste0('3k_dsci_down', lv, '_ARI_NMI.csv'))
  read.csv(f)
}
dN <- do.call(rbind, lapply(seq_along(LEVELS), function(i) {
  x <- read_lv(LEVELS[i])
  data.frame(Level = LEVELS[i], n = LEVEL_N[i], NMI = x$NMI, ARI = x$ARI) }))
cat('seeds per level:', paste(table(dN$Level), collapse = ' '), '\n')
s <- dN %>% group_by(Level, n) %>%
  summarise(nmi_m = mean(NMI), nmi_sd = sd(NMI), ari_m = mean(ARI), ari_sd = sd(ARI),
            .groups = 'drop') %>%
  pivot_longer(c(nmi_m, ari_m), names_to = 'stat', values_to = 'm')
s$sd    <- ifelse(s$stat == 'nmi_m', s$nmi_sd, s$ari_sd)
s$label <- ifelse(s$stat == 'nmi_m', 'NMI', 'ARI')
print(s %>% group_by(Level) %>% summarise(NMI = mean(m[label == 'NMI']),
                                          ARI = mean(m[label == 'ARI']), .groups = 'drop'))

p_c <- ggplot(s, aes(x = n, y = m, color = label)) +
  geom_errorbar(aes(ymin = m - sd, ymax = m + sd), linewidth = 0.25, width = 0.03) +
  geom_line(linewidth = 0.5) + geom_point(size = 1.1) +
  scale_x_log10(breaks = LEVEL_N, labels = LEVELS) +
  scale_color_manual(values = c(NMI = '#E41A1C', ARI = '#377EB8')) +
  labs(x = 'Downsample throughput (cells)', y = 'Score') +
  guides(color = guide_legend(title = NULL, keywidth = 0.25, keyheight = 0.25)) +
  theme_bw(base_size = 7) +
  theme(axis.text = element_text(size = 6.5), axis.text.x = element_text(size = 6.5,
        angle = 45, hjust = 1), axis.title = element_text(size = 7),
        panel.grid = element_blank(), panel.border = element_rect(linewidth = 0.3, fill = NA),
        legend.position = 'right', legend.text = element_text(size = 6.5),
        legend.key.size = unit(0.28, 'cm'),
        plot.margin = margin(1.5, 1.5, 1.5, 1.5, unit = 'mm'))
ggsave(file.path(OUT, 'S2c_NMI_ARI_vs_throughput.pdf'), p_c,
       width = 8.8, height = 6.55, units = 'cm', device = grDevices::cairo_pdf)
ggsave(file.path(OUT, 'S2c_NMI_ARI_vs_throughput.png'), p_c,
       width = 8.8, height = 6.55, units = 'cm', dpi = 150)
cat('saved S2c_NMI_ARI_vs_throughput.pdf\n')
