# ---------------------------------------------------------------------------
# Invasion assay — figures from invasion_conditions.csv (invasion_count.py aggregate).
# Style: theme_bw, CTRL grey, treated in colour, shape = biological replicate.
# Usage: Rscript plot_invasion.R <invasion_conditions.csv> <out_dir> [treated_colour]
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
csv  <- args[1]
out  <- args[2]
col  <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "darkred"

x <- read_csv(csv, show_col_types = FALSE) %>%
  mutate(is_ctrl   = toupper(condition) == "CTRL",
         condition = factor(condition, levels = unique(condition[order(!is_ctrl)])),
         line      = factor(line, levels = unique(line)),
         replicate = factor(replicate))

treated <- setdiff(levels(x$condition), levels(x$condition)[toupper(levels(x$condition)) == "CTRL"])
fill_val <- setNames(c("grey75", rep(col, length(treated))),
                     c(levels(x$condition)[toupper(levels(x$condition)) == "CTRL"], treated))

med <- function(d, v) d %>% group_by(line, condition) %>%
  summarise(m = mean({{ v }}, na.rm = TRUE),
            se = sd({{ v }}, na.rm = TRUE) / sqrt(sum(!is.na({{ v }}))), .groups = "drop")

# --- 1. invaded cells per field, one point per replicate
a <- med(x, cells_mean)
p_area <- ggplot(a, aes(condition, m, fill = condition)) +
  geom_col(colour = "black", width = 0.65) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.18) +
  geom_point(data = x, aes(condition, cells_mean, shape = replicate),
             inherit.aes = FALSE, size = 2, position = position_jitter(width = 0.12, height = 0)) +
  scale_fill_manual(values = fill_val) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Invaded cells per field", shape = NULL) +
  facet_wrap(~ line, nrow = 1) +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1)) +
  guides(fill = "none")

# --- 2. fold change vs the CTRL of the same replicate
f <- med(x, fold_change)
p_fc <- ggplot(f, aes(condition, m, fill = condition)) +
  geom_col(colour = "black", width = 0.65) +
  geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.18) +
  geom_point(data = x, aes(condition, fold_change, shape = replicate),
             inherit.aes = FALSE, size = 2, position = position_jitter(width = 0.12, height = 0)) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  scale_fill_manual(values = fill_val) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Invaded cells, fold change vs CTRL", shape = NULL) +
  facet_wrap(~ line, nrow = 1) +
  theme_bw() +
  theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1)) +
  guides(fill = "none")

# --- 3. spread of the fields within each condition (quality control)
fld <- file.path(dirname(csv), "invasion_fields.csv")
if (file.exists(fld)) {
  g <- read_csv(fld, show_col_types = FALSE) %>%
    filter(used) %>%
    mutate(condition = factor(condition, levels = levels(x$condition)),
           line = factor(line, levels = levels(x$line)), replicate = factor(replicate))
  p_fields <- ggplot(g, aes(condition, cells, fill = condition)) +
    geom_boxplot(outlier.shape = NA, width = 0.6) +
    geom_point(aes(shape = replicate), size = 1.4, alpha = .7,
               position = position_jitter(width = 0.14, height = 0)) +
    scale_fill_manual(values = fill_val) +
    labs(x = NULL, y = "Invaded cells, single field", shape = NULL) +
    facet_wrap(~ line, nrow = 1) +
    theme_bw() +
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 20, hjust = 1)) +
    guides(fill = "none")
  ggsave(file.path(out, "invasion_fields.png"), p_fields, width = 9, height = 4.2, dpi = 300)
}

ggsave(file.path(out, "invasion_cells.png"), p_area, width = 8, height = 4.2, dpi = 300)
ggsave(file.path(out, "invasion_fc.png"),   p_fc,   width = 8, height = 4.2, dpi = 300)
cat("wrote invasion_cells.png · invasion_fc.png · invasion_fields.png to", out, "\n")
