
suppressPackageStartupMessages({library(ggplot2); library(dplyr); library(tidyr)
  library(stringr)})
OUT  <- Sys.getenv('S2_OUT', 'figS2_out')
FY   <- file.path(Sys.getenv("S2_DATA", Sys.getenv("DSCI_FY_ROOT", file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "FY"))), '3k_dsci')
LEVELS <- c('3k','5k','10k','15k','20k','25k','30k','40k')
read_lv <- function(lv, what) {
  f <- file.path(OUT, paste0(lv, what))
  if (!file.exists(f)) f <- file.path(FY, paste0(if (grepl('ARI', what)) '3k_dsci_down' else '', lv, what))
  read.csv(f)
}
d <- do.call(rbind, lapply(LEVELS, function(lv) {
  x <- read_lv(lv, '_per_seed_per_type.csv'); x$Level <- lv; x }))
cat('rows:', nrow(d), '| types:', length(unique(d$type)), '| levels:',
    length(unique(d$Level)), '\n')
cat('rows per level:', paste(sapply(LEVELS, function(lv) sum(d$Level == lv)),
                             collapse = ' '), '\n')

# the composition §1a wrote, else the delivered copy of it
fc <- file.path(OUT, 'truth_composition.csv')
if (!file.exists(fc)) fc <- file.path(FY, '3k_dsci_full_annotCOP_meancelltype.csv')
comp <- read.csv(fc)
comp$ab <- comp$Freq / sum(comp$Freq)
d <- merge(d, comp[, c('type', 'ab')], by = 'type', all.x = TRUE)

# ---- paired test, recomputed here (not read from a table) -----------------
a <- d %>% filter(Level == '3k')  %>% select(Seed, type, recall)
b <- d %>% filter(Level == '40k') %>% select(Seed, type, recall)
m <- merge(a, b, by = c('Seed', 'type'), suffixes = c('_3k', '_40k'))
res <- m %>% group_by(type) %>%
  summarise(diff = mean(recall_40k - recall_3k),
            pairs = n(),
            p = tryCatch(suppressWarnings(
                   wilcox.test(recall_3k, recall_40k, paired = TRUE)$p.value),
                   error = function(e) NA_real_), .groups = 'drop') %>%
  mutate(padj = p.adjust(p, method = 'BH'))
write.csv(res, file.path(OUT, 'FigD_paired_recall_3k_vs_40k.csv'), row.names = FALSE)
sig <- res$type[which(res$padj < 0.05 & res$diff >= 0.5)]   # which(): NA-safe
cat('marked types (BH-FDR < 0.05 & delta >= 0.5):', length(sig), '\n')
cat(sort(sig), collapse = ', ', '\n')
cat('mean recall by level:', paste(sprintf('%s %.3f', LEVELS,
    sapply(LEVELS, function(lv) mean(d$recall[d$Level == lv], na.rm = TRUE))),
    collapse = '  '), '\n')
cat('COP/NFOL recall by level:', paste(sprintf('%s %.2f', LEVELS,
    sapply(LEVELS, function(lv) mean(d$recall[d$Level == lv & d$type == 'COP/NFOL'],
                                      na.rm = TRUE))), collapse = '  '), '\n')

# ---- the heatmap ----------------------------------------------------------
dd <- d %>% group_by(type, Level) %>%
  summarise(m = mean(recall, na.rm = TRUE), ab = first(ab), .groups = 'drop')
ord <- dd %>% arrange(ab) %>% pull(type) %>% unique()
dd$type  <- factor(dd$type, levels = ord)
dd$Level <- factor(dd$Level, levels = LEVELS)
wrapped <- str_wrap(ord, width = 10)
have_ggtext <- requireNamespace('ggtext', quietly = TRUE)
if (have_ggtext) {
  lab <- ifelse(ord %in% sig, sprintf("<span style='color:#B2182B'>%s</span>", wrapped),
                sprintf("<span style='color:grey25'>%s</span>", wrapped))
} else {
  cat('[note] ggtext absent -> significant labels get a "+" marker instead of red ink\n')
  lab <- ifelse(ord %in% sig, paste0(wrapped, ' +'), wrapped)
}
p_d <- ggplot(dd, aes(x = type, y = Level, fill = m)) +
  geom_tile(color = 'white', linewidth = 0.2) +
  scale_fill_gradient2(low = '#2166AC', mid = '#F7F7F7', high = '#B2182B',
                       midpoint = 0.5, na.value = 'grey92', name = 'Mean\nrecall') +
  scale_x_discrete(labels = setNames(lab, ord)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 7) +
  theme(axis.text.x = if (have_ggtext) ggtext::element_markdown(angle = 90, hjust = 0.5,
            vjust = 0.5, size = 6, lineheight = 0.85) else
            element_text(angle = 90, hjust = 0.5, vjust = 0.5, size = 6,
                         lineheight = 0.85, color = 'grey25'),
        axis.text.y = element_text(size = 6.5), axis.ticks = element_blank(),
        panel.border = element_blank(), panel.grid = element_blank(),
        legend.position = 'right', legend.title = element_text(size = 6.5),
        legend.text = element_text(size = 6), legend.key.size = unit(0.28, 'cm'),
        plot.margin = margin(1.5, 1.5, 1.5, 1.5, unit = 'mm'))
ggsave(file.path(OUT, 'S2d_pertype_recall_heatmap.pdf'), p_d,
       width = 16.97, height = 4.90, units = 'cm', device = grDevices::cairo_pdf)
ggsave(file.path(OUT, 'S2d_pertype_recall_heatmap.png'), p_d,
       width = 16.97, height = 4.90, units = 'cm', dpi = 150)
cat('saved S2d_pertype_recall_heatmap.pdf\n')
