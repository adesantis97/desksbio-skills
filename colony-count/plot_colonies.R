# ---------------------------------------------------------------------------
# Colony formation — figures from counts.csv (written by colony_count.py apply).
# Style: theme_bw, CTRL grey, treated in colour, dashed line at 1.
# Usage: Rscript plot_colonies.R <counts.csv> <out_dir> [treated_colour]
# ---------------------------------------------------------------------------
suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
csv  <- args[1]
out  <- args[2]
col  <- if (length(args) >= 3 && nzchar(args[3])) args[3] else "darkred"

x <- read_csv(csv, show_col_types = FALSE) %>%
  mutate(dose  = if ("dose" %in% names(.)) as.character(dose) else NA_character_,
         group = if_else(is.na(dose), line, paste(line, dose)),
         group = factor(group, levels = unique(group)),
         is_ctrl = toupper(condition) == "CTRL",
         condition = factor(condition, levels = unique(condition[order(!is_ctrl)])))

treated <- levels(x$condition)[!toupper(levels(x$condition)) == "CTRL"]
fill_val <- setNames(c("grey75", rep(col, length(treated))),
                     c(levels(x$condition)[toupper(levels(x$condition)) == "CTRL"], treated))
lab_pct  <- scales::percent_format(accuracy = 1)

# --- 1. raw values: colonies (n) and mean size (px), CTRL vs treated
raw <- x %>%
  select(group, condition, `Colonies (n)` = colonies, `Mean colony size (px)` = mean_size_px) %>%
  pivot_longer(-c(group, condition), names_to = "metric", values_to = "value") %>%
  mutate(metric = factor(metric, levels = c("Colonies (n)", "Mean colony size (px)")))

p_raw <- ggplot(raw, aes(x = condition, y = value, fill = condition)) +
  geom_col(colour = "black", width = 0.65) +
  geom_text(aes(label = round(value, 1)), vjust = -0.4, size = 3) +
  scale_fill_manual(values = fill_val) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = NULL) +
  facet_grid(metric ~ group, scales = "free_y", switch = "y") +
  theme_bw() +
  theme(legend.position = "none", strip.placement = "outside",
        axis.text.x = element_text(angle = 20, hjust = 1))

# --- 1b. colony number only
p_n <- ggplot(filter(raw, metric == "Colonies (n)"), aes(x = condition, y = value, fill = condition)) +
  geom_col(colour = "black", width = 0.65) +
  geom_text(aes(label = value), vjust = -0.4, size = 3) +
  scale_fill_manual(values = fill_val) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "Colonies (n)") +
  facet_wrap(~ group, nrow = 1) +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))

# --- 1c. size of every colony (colonies.csv): box + points
cols <- read_csv(file.path(dirname(csv), "colonies.csv"), show_col_types = FALSE) %>%
  mutate(dose  = if ("dose" %in% names(.)) as.character(dose) else NA_character_,
         group = factor(if_else(is.na(dose), line, paste(line, dose)), levels = levels(x$group)),
         condition = factor(condition, levels = levels(x$condition)))

p_size <- ggplot(cols, aes(x = condition, y = size_px, fill = condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.6) +
  geom_point(size = 0.8, alpha = 0.6, position = position_jitter(width = 0.15, seed = 1)) +
  scale_fill_manual(values = fill_val) +
  scale_y_log10() +
  labs(x = NULL, y = "Colony size (px, log)") +
  facet_wrap(~ group, nrow = 1, drop = FALSE) +
  theme_bw() +
  theme(legend.position = "none", axis.text.x = element_text(angle = 20, hjust = 1))

# --- 2. fold change vs CTRL: number and mean size
fc <- x %>% filter(!is_ctrl) %>%
  select(group, Colonies = fold_change, `Mean size` = size_fold_change) %>%
  pivot_longer(-group, names_to = "metric", values_to = "fc") %>%
  mutate(metric = factor(metric, levels = c("Colonies", "Mean size")))

p_fc <- ggplot(fc, aes(x = group, y = fc)) +
  geom_col(fill = col, colour = "black", width = 0.55) +
  geom_text(aes(label = lab_pct(fc)), vjust = -0.4, size = 3) +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  scale_y_continuous(labels = lab_pct, limits = c(0, max(1.25, max(fc$fc, na.rm = TRUE) * 1.1)),
                     expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = paste0(treated[1], " / CTRL")) +
  facet_wrap(~ metric) +
  theme_bw()

w <- 4 + 3.5 * nlevels(x$group)
ggsave(file.path(out, "colonies_raw.png"), p_raw, width = w, height = 11, units = "cm", dpi = 300)
ggsave(file.path(out, "colonies_n.png"),   p_n,   width = w, height = 8,  units = "cm", dpi = 300)
ggsave(file.path(out, "colonies_size.png"), p_size, width = w, height = 8, units = "cm", dpi = 300)
ggsave(file.path(out, "colonies_fc.png"),  p_fc,  width = w + 3, height = 9, units = "cm", dpi = 300)
cat("done:", file.path(out, c("colonies_n.png", "colonies_size.png", "colonies_raw.png", "colonies_fc.png")), sep = "\n  ")
