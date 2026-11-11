
suppressPackageStartupMessages({library(ggplot2); library(dplyr)})
OUT <- Sys.getenv('S2_OUT', 'figS2_out')
S   <- readRDS(file.path(OUT, 's2_light.RDS'))
md  <- S$meta; dot_order <- S$dot_order; num_map <- S$num_map; cols <- S$cols
leg_lab <- setNames(paste0(num_map[dot_order], ' ', dot_order), dot_order)

prop <- md %>% count(technique, cell_type_updated) %>%
  group_by(technique) %>% mutate(proportion = n / sum(n) * 100) %>% ungroup()
tech_levels <- c('10x', 'dsciCAT-seq', 'HT-scCAT-seq-HT')
prop$technique <- factor(prop$technique, levels = tech_levels)

cat('per-technology totals:\n'); print(prop %>% group_by(technique) %>%
  summarise(cells = sum(n), .groups = 'drop'))
cat('\ntypes the legend quotes (%, by technology):\n')
q <- rbind('In_D1', 'AST', 'Oligo', 'COP/NFOL', 'OPC', 'Ex_L2/3_IT', 'Ex_L6_CT')
w <- prop %>% filter(cell_type_updated %in% q) %>%
  select(technique, cell_type_updated, proportion) %>%
  tidyr::pivot_wider(names_from = technique, values_from = proportion)
print(as.data.frame(
  cbind(cell_type_updated = w$cell_type_updated,
        round(w[, tech_levels, drop = FALSE], 1)), row.names = NULL))

p_b <- ggplot(prop, aes(x = technique, y = proportion,
                        fill = factor(cell_type_updated, levels = dot_order))) +
  geom_col(width = 0.68, colour = NA) +
  scale_x_discrete(limits = tech_levels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = cols, breaks = dot_order, labels = leg_lab,
                    name = 'Cell type (n = number)') +
  guides(fill = guide_legend(ncol = 2, keywidth = 0.20, keyheight = 0.20)) +
  labs(x = NULL, y = 'Cell type proportion (%)') +
  theme_bw(base_size = 7) +
  theme(axis.text = element_text(size = 6.5), axis.text.x = element_text(size = 6.5,
        angle = 30, hjust = 1, vjust = 1), panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.3, fill = NA),
        legend.position = 'right', legend.text = element_text(size = 6),
        legend.title = element_text(size = 6.5), legend.key = element_blank(),
        legend.background = element_blank(),
        plot.margin = margin(1.5, 1.5, 1.5, 1.5, unit = 'mm'))
ggsave(file.path(OUT, 'S2b_proportion_by_technique.pdf'), p_b,
       width = 8.8, height = 6.55, units = 'cm', device = grDevices::cairo_pdf)
ggsave(file.path(OUT, 'S2b_proportion_by_technique.png'), p_b,
       width = 8.8, height = 6.55, units = 'cm', dpi = 150)
cat('\nsaved S2b_proportion_by_technique.pdf\n')
