library(gcplyr)
library(lubridate)
library(ggplot2)
library(dplyr)
library(readxl)
library(stringr)
library(forcats)
library(emmeans)
library(pheatmap)
library(tidyverse)
library(cowplot)


source("utilities.R")


######### A loop to generate stats and all growth curves in the master folder
# Define the main directory to search for folders
main_directory <- "/Tetresistance/Tet/RealData"

# Get a list of all directories in the main directory
all_dirs <- list.dirs(main_directory, recursive = TRUE)

# Loop through each directory
for (dir in all_dirs) {
  # Get all Excel files in the current directory

  all_csv_files <- list.files(path = dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(all_csv_files) == 0) {
    message(paste("No .csv files found in directory:", dir))
    next # Skip to the next iteration of the loop
  }
  # all_stats_files <- list.files(path = dir, pattern = "\\_stats.csv$", full.names = TRUE)
  if (grepl("figures_stats$", basename(dir))) {
    message(paste("Skipping directory ending in figures_stats:", dir))
    next
  }
  # Filter out files that contain "media", "strain", or "deepwell"
  datafiles <- all_csv_files[!grepl("media|strain|deepwell|stats|_growthnormalizedtomin|stickers|concentration", all_csv_files, ignore.case = TRUE)]
  # datafiles <- list.files(path = dir, pattern = "^(?!.*(media|strain|deepwell)).*\\.csv$", full.names = TRUE)
  mediafiles <- list.files(path = dir, pattern = "media.*\\.csv$", full.names = TRUE)
  strainfiles <- list.files(path = dir, pattern = "strains.*\\.csv$", full.names = TRUE)
  concentrationfiles <- list.files(path = dir, pattern = "concentration.*\\.csv$", full.names = TRUE)

  imported_widedata_file <- read_wides(files = datafiles, startrow = 1)

  imported_wides_now_tidy_file <- trans_wide_to_tidy(
    wides = imported_widedata_file,
    id_cols = c("file", "T° 600", "Time")
  )

  my_design_file <- import_blockdesigns(files = c(mediafiles, strainfiles, concentrationfiles), block_names = c("Media", "Strain_letters", "Concentration"), startrow = 2)
  ex_dat_mrg_file <- merge_dfs(imported_wides_now_tidy_file, my_design_file)

  ex_dat_mrg_file$Time <- time_length(hms(ex_dat_mrg_file$Time), unit = "hour")

  ex_dat_mrg_file$Well <- factor(ex_dat_mrg_file$Well, levels = paste0(rep(LETTERS[1:8], each = 12), 1:12))

  p0 <- ggplot(data = ex_dat_mrg_file, aes(x = Time, y = Measurements)) +
    geom_line() +
    facet_wrap(~Well, nrow = 8, ncol = 12)

  plot_file <- file.path(dir, "/figures_stats/", paste0(Sys.Date(), tools::file_path_sans_ext(basename(datafiles[1])), "_growth_plot.png"))

  ex_dat_mrg_file <- mutate(group_by(ex_dat_mrg_file, Well, Media, Strain_letters), deriv_percap5 = calc_deriv(x = Time, y = Measurements, percapita = TRUE, blank = (first_minima(y = Measurements, x = Time, return = "y")), window_width_n = 5, trans_y = "log"))


  ex_dat_mrg_sum_file_minmax <-
    summarize(group_by(ex_dat_mrg_file, Well, Media, Strain_letters, Concentration),
      max_dens = max_gc(Measurements, na.rm = TRUE),
      max_time = extr_val(Time, which_max_gc(Measurements)),
      min_dens = first_minima(y = Measurements, x = Time, return = "y"),
      min_time = first_minima(y = Measurements, x = Time, return = "x"),
      # lag_time = lag_time(y = Measurements, x = Time, deriv = deriv_percap5, y0 = min_dens, blank = min_dens),
      # deriv_percap5 = calc_deriv(x = Time, y = Measurements, percapita = TRUE, blank = min_dens, window_width_n = 5, trans_y = "log"),
      max_percap = max_gc(deriv_percap5, na.rm = TRUE), max_percap_time = extr_val(Time, which_max_gc(deriv_percap5)), doub_time = doubling_time(y = max_percap),
      auc = auc(x = Time, y = Measurements, blank = min(Measurements))
    )

  ex_dat_mrg_sum_file_minmax$maxminusmin <- ex_dat_mrg_sum_file_minmax$max_dens - ex_dat_mrg_sum_file_minmax$min_dens

  ex_dat_mrg_sum_file_minmax$filepath <- file.path(dir, "/figures_stats/", paste0(tools::file_path_sans_ext(basename(datafiles[1])), "_stats.csv"))

  output_file <- file.path(dir, "/figures_stats/", paste0(tools::file_path_sans_ext(basename(datafiles[1])), "_stats.csv"))

  write.csv(ex_dat_mrg_sum_file_minmax, output_file, row.names = FALSE)
}


stats_194_201_X_WT <- read.csv("/Tetresistance/Tet/RealData/250815_MMM_194_201_X_WT/figures_stats/250815_MMM_194_201_X_WT_stats.csv")
stats_213_216_X_WT <- read.csv("/Tetresistance/Tet/RealData/250816_MMM_213_216_X_WT/figures_stats/250816_MMM_213_216_X_WT_stats.csv")
stats_194_201_X_tet38 <- read.csv("/Tetresistance/Tet/RealData/250817_MMM_194_201_X_Tet38KO/figures_stats/250817_MMM_194_201_X_Tet38KO_stats.csv")
stats_213_216_X_tet38 <- read.csv("/Tetresistance/Tet/RealData/250819_MMM_213_216_tet38KO/figures_stats/250819_MMM_213_216_tet38KO_stats.csv")


# Stack then summary
df_all <- bind_rows(
  stats_194_201_X_WT %>% mutate(source = "stats_194_201_X_WT"),
  stats_213_216_X_WT %>% mutate(source = "stats_213_216_X_WT"),
  stats_194_201_X_tet38 %>% mutate(source = "stats_194_201_X_tet38"),
  stats_213_216_X_tet38 %>% mutate(source = "stats_213_216_X_tet38")
)

# Now apply your summary pipeline once
df_summary_all <- df_all %>%
  group_by(Media, Concentration, Strain_letters) %>%
  summarise(
    mean_AUC = mean(auc, na.rm = TRUE),
    sd_AUC   = sd(auc, na.rm = TRUE),
    n        = n(),
    .groups  = "drop"
  ) %>%
  mutate(
    se_AUC = sd_AUC / sqrt(n)
  ) %>%
  group_by(Media, Strain_letters) %>%
  mutate(
    percent_max_AUC = 100 * mean_AUC / max(mean_AUC, na.rm = TRUE),
    percent_se      = 100 * se_AUC / max(mean_AUC, na.rm = TRUE)
  )


summarywt <- df_summary_all %>%
  filter(Strain_letters %in% c("JE2"))

# Plot with error bars
p1 <- ggplot(summarywt, aes(
  x = as.factor(Concentration),
  y = percent_max_AUC,
  color = Media,
  shape = Media,
  group = Media
)) +
  geom_line() +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = percent_max_AUC - percent_se,
      ymax = percent_max_AUC + percent_se
    ),
    width = 0.1
  ) +
  scale_color_manual(values = c(
    "MMM" = "black",
    "MMM_194" = "#9F00BB",
    "MMM_201" = "#6A0DAD",
    "MMM_213" = "#4B0082",
    "MMM_216" = "#6A5ACD"
  )) +
  scale_shape_manual(values = c(
    "MMM" = 16, # filled circle
    "MMM_194" = 17, # triangle
    "MMM_201" = 15, # square
    "MMM_213" = 18, # diamond
    "MMM_216" = 8
  )) +
  labs(x = "Tet (ug/mL)", y = "% of max AUC in WT") +
  theme_pub()


plot_file <- paste0("Tetresistance/Analysis/Figures/250924/", Sys.Date(), "tetresistanceWT.pdf")


summarytet38 <- df_summary_all %>%
  filter(Strain_letters %in% c("tet38"))

# Plot with error bars
p2 <- ggplot(summarytet38, aes(
  x = as.factor(Concentration),
  y = percent_max_AUC,
  color = Media,
  shape = Media,
  group = Media
)) +
  geom_line() +
  geom_point(size = 3) +
  geom_errorbar(
    aes(
      ymin = percent_max_AUC - percent_se,
      ymax = percent_max_AUC + percent_se
    ),
    width = 0.1
  ) +
  scale_color_manual(values = c(
    "MMM" = "black",
    "MMM_194" = "#9F00BB",
    "MMM_201" = "#6A0DAD",
    "MMM_213" = "#4B0082",
    "MMM_216" = "#6A5ACD"
  )) +
  scale_shape_manual(values = c(
    "MMM" = 16, # filled circle
    "MMM_194" = 17, # triangle
    "MMM_201" = 15, # square
    "MMM_213" = 18, # diamond
    "MMM_216" = 8
  )) +
  labs(x = "Tet (ug/mL)", y = "% of max AUC in WT") +
  theme_pub()


# Statistics
group_means <- df_all %>%
  group_by(Media, Strain_letters, Concentration) %>%
  summarise(mean_AUC = mean(auc, na.rm = TRUE), .groups = "drop")

# Step 2: find max mean_AUC per (Media, Strain_letters)
max_means <- group_means %>%
  group_by(Media, Strain_letters) %>%
  summarise(max_mean_AUC = max(mean_AUC, na.rm = TRUE), .groups = "drop")

# Step 3: merge back into replicate-level data
df_all_replicates <- df_all %>%
  left_join(max_means, by = c("Media", "Strain_letters")) %>%
  mutate(percent_max_AUC = 100 * auc / max_mean_AUC)


df_all_replicatesWT <- df_all_replicates %>%
  filter(Strain_letters %in% c("JE2"))


lm_resultWT <- lm(percent_max_AUC ~ Media * as.factor(Concentration),
  data = df_all_replicatesWT
)

anova_resultWT <- anova(lm_resultWT)
print(anova_resultWT)

# Post-hoc pairwise strain comparisons (within each concentration)
emmWT <- emmeans(lm_resultWT, ~ Media | Concentration)
pairwise_strainsWT <- contrast(emmWT, method = "pairwise", adjust = "tukey")
print(pairwise_strainsWT)


df_all_replicatestet38 <- df_all_replicates %>%
  filter(Strain_letters %in% c("tet38"))


lm_resulttet38 <- lm(percent_max_AUC ~ Media * as.factor(Concentration),
  data = df_all_replicatestet38
)

anova_resulttet38 <- anova(lm_resulttet38)
print(anova_resulttet38)

# Post-hoc pairwise strain comparisons (within each concentration)
emmtet38 <- emmeans(lm_resulttet38, ~ Media | Concentration)
pairwise_strainstet38 <- contrast(emmtet38, method = "pairwise", adjust = "tukey")
print(pairwise_strainstet38)
