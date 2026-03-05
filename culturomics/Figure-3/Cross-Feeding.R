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
main_directory <- "Crossfeeding/WT_NanA"

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
  datafiles <- all_csv_files[!grepl("media|strain|deepwell|stats|_growthnormalizedtomin|stickers", all_csv_files, ignore.case = TRUE)]
  # datafiles <- list.files(path = dir, pattern = "^(?!.*(media|strain|deepwell)).*\\.csv$", full.names = TRUE)
  mediafiles <- list.files(path = dir, pattern = "media.*\\.csv$", full.names = TRUE)
  strainfiles <- list.files(path = dir, pattern = "strains.*\\.csv$", full.names = TRUE)

  imported_widedata_file <- read_wides(files = datafiles, startrow = 1)

  imported_wides_now_tidy_file <- trans_wide_to_tidy(
    wides = imported_widedata_file,
    id_cols = c("file", "T° 600", "Time")
  )

  my_design_file <- import_blockdesigns(files = c(mediafiles, strainfiles), block_names = c("Media", "Strain_letters"), startrow = 2)
  ex_dat_mrg_file <- merge_dfs(imported_wides_now_tidy_file, my_design_file)

  ex_dat_mrg_file$Time <- time_length(hms(ex_dat_mrg_file$Time), unit = "hour")

  ex_dat_mrg_file$Well <- factor(ex_dat_mrg_file$Well, levels = paste0(rep(LETTERS[1:8], each = 12), 1:12))

  p0 <- ggplot(data = ex_dat_mrg_file, aes(x = Time, y = Measurements)) +
    geom_line() +
    facet_wrap(~Well, nrow = 8, ncol = 12)

  plot_file <- file.path(dir, "/figures_stats/", paste0(Sys.Date(), tools::file_path_sans_ext(basename(datafiles[1])), "_growth_plot.png"))
  ex_dat_mrg_file <- mutate(group_by(ex_dat_mrg_file, Well, Media, Strain_letters), deriv_percap5 = calc_deriv(x = Time, y = Measurements, percapita = TRUE, blank = (first_minima(y = Measurements, x = Time, return = "y")), window_width_n = 5, trans_y = "log"))


  ex_dat_mrg_sum_file_minmax <-
    summarize(group_by(ex_dat_mrg_file, Well, Media, Strain_letters),
      max_dens = max_gc(Measurements, na.rm = TRUE),
      max_time = extr_val(Time, which_max_gc(Measurements)),
      min_dens = first_minima(y = Measurements, x = Time, return = "y"),
      min_time = first_minima(y = Measurements, x = Time, return = "x"),
      max_percap = max_gc(deriv_percap5, na.rm = TRUE), max_percap_time = extr_val(Time, which_max_gc(deriv_percap5)), doub_time = doubling_time(y = max_percap),
      auc = auc(x = Time, y = Measurements, blank = min_dens)
    )

  ex_dat_mrg_sum_file_minmax$maxminusmin <- ex_dat_mrg_sum_file_minmax$max_dens - ex_dat_mrg_sum_file_minmax$min_dens

  ex_dat_mrg_sum_file_minmax$filepath <- file.path(dir, "/figures_stats/", paste0(tools::file_path_sans_ext(basename(datafiles[1])), "_stats.csv"))

  output_file <- file.path(dir, "/figures_stats/", paste0(tools::file_path_sans_ext(basename(datafiles[1])), "_stats.csv"))
}


###### Making figures with the merged stats file

# Noting samples that are obviously contaminated or plate reader noise manually
# I believe this was potentially a batch of bad filters, we had people with filters failing around the entire lab around this week
# 250205_211_207_201_154_X_WT_growth_plot.png shows that the entire 5 column had contamination. Also A9, B9, C9, G9, H9 are contaminated
# 250129_211_207_194_X_WT_growth_plot.png shows that the 5 column had some contamination. Also it looks like D3, E3, and E1, had some severe plate reader noise
# 250209_123_154_204_208_213_220_X_WT_growth_plot.png has the entire column 6 contaminated
# 250210_128_185_194_206_216_217_X_WT_growth_plot.png has A5,B5,C5, entire 12 column, and G5, H5, D6,E6,F6


##################### Remaking all figures, but with normalizing to lowest value######################################
######### A loop to generate growth curves normalized to min density for the purpose of adjusted AUC calculation in the master folder
# Define the main directory to search for folders
main_directory <- "~/Library/CloudStorage/Box-Box/2024_Dome2.0/Crossfeeding/WT_NanA"

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
  if (grepl("figures_stats$", basename(dir))) {
    message(paste("Skipping directory ending in figures_stats:", dir))
    next
  }
  # Filter out files that contain "media", "strain", or "deepwell"
  datafiles <- all_csv_files[!grepl("media|strain|deepwell|stats|_growthnormalizedtomin|stickers", all_csv_files, ignore.case = TRUE)]
  statsfiles <- list.files(path = paste0(dir, "/figures_stats/"), pattern = "stats.*\\.csv$", full.names = TRUE)

  # Read the first spreadsheet where column names are well names
  first_sheet <- read.csv(statsfiles)

  # Read the second spreadsheet where 'Well' and 'min_dens' are columns
  second_sheet <- read.csv(datafiles)

  # Convert the second sheet to a named vector for easier matching
  min_dens_values <- setNames(first_sheet$min_dens, first_sheet$Well)

  # Subtract the matched 'min_dens' values from the first spreadsheet
  result_sheet <- second_sheet
  for (well in names(result_sheet)) {
    if (well %in% names(min_dens_values)) {
      result_sheet[[well]] <- result_sheet[[well]] - min_dens_values[well]
    }
  }
  output_file <- file.path(dir, "/figures_stats/", paste0(tools::file_path_sans_ext(basename(datafiles[1])), "_growthnormalizedtomin.csv"))
}


#### Remaking stats but with the adjusted growth curves

######### A loop to generate stats on all growth curves in the master folder
# Define the main directory to search for folders
main_directory <- "~/Library/CloudStorage/Box-Box/2024_Dome2.0/Crossfeeding/WT_NanA"

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
  if (grepl("figures_stats$", basename(dir))) {
    message(paste("Skipping directory ending in figures_stats:", dir))
    next
  }
  # Filter out files that contain "media", "strain", or "deepwell"
  datafiles <- list.files(path = paste0(dir, "/figures_stats/"), pattern = "growthnormalizedtomin.*\\.csv$", full.names = TRUE)
  mediafiles <- list.files(path = dir, pattern = "media.*\\.csv$", full.names = TRUE)
  strainfiles <- list.files(path = dir, pattern = "strains.*\\.csv$", full.names = TRUE)

  imported_widedata_file <- read_wides(files = datafiles, startrow = 1)

  imported_wides_now_tidy_file <- trans_wide_to_tidy(
    wides = imported_widedata_file,
    id_cols = c("file", "T..600", "Time")
  )

  my_design_file <- import_blockdesigns(files = c(mediafiles, strainfiles), block_names = c("Media", "Strain_letters"), startrow = 2)
  ex_dat_mrg_file <- merge_dfs(imported_wides_now_tidy_file, my_design_file)

  ex_dat_mrg_file$Time <- time_length(hms(ex_dat_mrg_file$Time), unit = "hour")

  ex_dat_mrg_file$Well <- factor(ex_dat_mrg_file$Well, levels = paste0(rep(LETTERS[1:8], each = 12), 1:12))

  p0 <- ggplot(data = ex_dat_mrg_file, aes(x = Time, y = Measurements)) +
    geom_line() +
    facet_wrap(~Well, nrow = 8, ncol = 12)

  plot_file <- file.path(dir, "/figures_stats/", paste0(Sys.Date(), tools::file_path_sans_ext(basename(datafiles[1])), "_growth_plot_adjusted.png"))
  next
}


####### To merge stats files
master_dir <- "Crossfeeding/WT_NanA"
#
# # Get all subdirectories within the master directory
all_subdirs <- list.dirs(path = master_dir, full.names = TRUE, recursive = TRUE)
all_subdirs <- grep("figures_stats", all_subdirs, value = TRUE)
#
# # Initialize an empty list to store dataframes
merged_data <- list()
#
# # Loop through each subdirectory to find and merge '_stats' files
for (subdir in all_subdirs) {
  #
  #   # List all files in the current subdirectory that contain '_stats' in their names
  stats_files <- list.files(path = subdir, pattern = "_stats.*\\.csv$", full.names = TRUE)
  #
  #   # If no '_stats' files are found, skip this directory
  if (length(stats_files) == 0) {
    message(paste("No '_stats' files found in directory:", subdir))
    next
  }
  #
  #   # Read each '_stats' file and add it to the list
  for (file in stats_files) {
    file_data <- read.csv(file)
    merged_data[[length(merged_data) + 1]] <- file_data
  }
}
#
# # If there are no '_stats' files to merge, exit
if (length(merged_data) == 0) {
  stop("No '_stats' files found in any of the subdirectories.")
}

# Merge all dataframes in the list by columns
final_merged_data <- Reduce(function(x, y) merge(x, y, all = TRUE), merged_data)

# Now going to clean up this so that I get rid of the above mentioned bad data


# 250205_211_207_201_154_X_WT_growth_plot.png shows that the entire 5 column had contamination. Also A9, B9, C9, G9, H9 are contaminated
my_values <- c("A5", "B5", "C5", "D5", "E5", "F5", "G5", "H5", "A9", "B9", "C9", "G9", "H9")
df_filtered_1 <- final_merged_data %>% filter(!(str_detect(filepath, "250205_211_207_201_154_X_WT") & Well %in% my_values))

# 250129_211_207_194_X_WT_growth_plot.png shows that the 5 column had some contamination. Also it looks like D3, E3, and E1, had some severe plate reader noise
my_values <- c("A5", "B5", "C5", "D5", "E5", "F5", "G5", "H5", "D3", "E3", "E1")
df_filtered_2 <- df_filtered_1 %>% filter(!(str_detect(filepath, "250129_211_207_194_X_WT") & Well %in% my_values))


# 250209_123_154_204_208_213_220_X_WT_growth_plot.png has the entire column 6 contaminated
my_values <- c("A6", "B6", "C6", "D6", "E6", "F6", "G6", "H6")
df_filtered_3 <- df_filtered_2 %>% filter(!(str_detect(filepath, "250209_123_154_204_208_213_220_X_WT") & Well %in% my_values))

# 250210_128_185_194_206_216_217_X_WT_growth_plot.png has A5,B5,C5, entire 12 column, and G5, H5, D6,E6,F6

my_values <- c("A12", "B12", "C12", "D12", "E12", "F12", "G12", "H12", "A5", "B5", "C5", "G5", "H5", "D6", "E6", "F6")
cleanedstatsdata <- df_filtered_3 %>% filter(!(str_detect(filepath, "250210_128_185_194_206_216_217_X_WT") & Well %in% my_values))

taxonomy <- read.csv("~/Library/CloudStorage/Box-Box/2024_Dome2.0/Culturing/IsolateAnalysis/Tree/FinalGTDBTK/250917_taxonomymap.csv")

taxonomy$Media <- sub("^MMM_0", "", taxonomy$Media)
taxonomy$Media <- sub("^0", "", taxonomy$Media)

cleanedstatsdata$Media <- sub("^MMM_", "", cleanedstatsdata$Media)

stats_taxon <- merge(cleanedstatsdata, taxonomy, by = "Media", all.x = TRUE)

stats_taxon <- stats_taxon %>%
  filter(!is.na(Species) | Media %in% c("MMM", "MMMG"))

stats_taxon <- stats_taxon %>%
  mutate(category = case_when(
    Media == "Empty" | Strain_letters == "Empty" ~ "Empty",
    !is.na(Species) ~ Species,
    TRUE ~ Media
  ))

stats_taxon$category <- paste0(toupper(substr(stats_taxon$category, 1, 1)), substr(stats_taxon$category, 2, nchar(stats_taxon$category)))


stats_taxon <- stats_taxon %>%
  mutate(
    Genus = word(Species, 1), # first word of Species
    genus_color = ifelse(Genus %in% names(genus.colors),
      genus.colors[Genus],
      "black"
    )
  )

P1 <- stats_taxon %>%
  filter(Strain_letters != "NanA") %>%
  filter(category != "Empty" & category != "MMG" & category != "MM") %>%
  mutate(category = fct_reorder(category, auc, .fun = median, .desc = FALSE)) %>%
  ggplot(aes(x = category, y = auc, fill = genus_color)) +
  coord_flip() +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = .5, alpha = 0.7, aes(color = "#899499")) +
  scale_fill_identity() +
  scale_color_identity() +
  labs(
    title = "AUC by Media/Supernatant",
    x = "Media/Supernatant source",
    y = "Area Under Curve (AUC)"
  ) +
  theme_pub()

print(P1)

tmp.fname <- paste0("~/Library/CloudStorage/Box-Box/2024_Dome2.0/Crossfeeding/Analysis/250918Figures/", Sys.Date(), "AUCbymedia.pdf")

P2 <- stats_taxon %>%
  filter(Strain_letters != "NanA") %>%
  # filter(!(category != "Empty" & Strain_letters == "Empty")) %>%
  mutate(category = fct_reorder(category, maxminusmin, .fun = median, .desc = FALSE)) %>%
  ggplot(aes(x = category, y = maxminusmin)) +
  coord_flip() +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = .5, alpha = 0.7) +
  theme_minimal() +
  labs(
    title = "Maxminusmin by Media/Supernatant",
    x = "Media/Supernatant source",
    y = "Maxminusmin"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))

plot_file <- paste0("Crossfeeding/Analysis/250730Figures/", Sys.Date(), "Maxminusminbymedia.png")

P3 <- stats_taxon %>%
  filter(Strain_letters != "NanA") %>%
  filter(category != "Empty") %>%
  filter(category != "MM") %>%
  # filter(!(category != "Empty" & Strain_letters == "Empty")) %>%
  mutate(category = fct_reorder(category, max_percap, .fun = median, .desc = FALSE)) %>%
  ggplot(aes(x = category, y = max_percap)) +
  coord_flip() +
  geom_jitter(width = 0.2, size = .5, alpha = 0.7) +
  geom_boxplot(outlier.shape = NA) +
  theme_minimal() +
  labs(
    title = "Max percap growth rate by Media/Supernatant",
    x = "Media/Supernatant source",
    y = "Max percap growth rate"
  ) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 1))


plot_file <- paste0("Crossfeeding/Analysis/250730Figures/", Sys.Date(), "Growthratebymedia.png")

# Statistics
ref_level <- "MMM"
stats_taxon$category <- factor(stats_taxon$category)
stats_taxon$category <- relevel(stats_taxon$category, ref = ref_level)

# Outcome variables
outcomes <- c("auc", "maxminusmin", "max_percap")

# Initialize results list
all_results <- list()

for (outcome in outcomes) {
  # Fit model
  formula <- as.formula(paste0(outcome, " ~ category"))
  model <- lm(formula, data = stats_taxon)

  # Get emmeans
  em <- emmeans(model, "category")
  ref_index <- which(levels(stats_taxon$category) == ref_level)

  # Contrast: treatment vs control
  cont <- contrast(em, method = "trt.vs.ctrl", ref = ref_index)
  cont_summary <- summary(cont, adjust = "BH") %>%
    as.data.frame() %>%
    mutate(outcome = outcome)

  # Calculate N for each category
  n_by_category <- stats_taxon %>%
    filter(!is.na(.data[[outcome]])) %>%
    count(category, name = "n")

  # Add N to contrast summary
  # Extract category names from contrast column
  cont_summary <- cont_summary %>%
    mutate(
      ref_category = ref_level,
      comparison_category = gsub(paste0(" - ", ref_level), "", contrast)
    ) %>%
    left_join(n_by_category %>% rename(n_ref = n),
      by = c("ref_category" = "category")
    ) %>%
    left_join(n_by_category %>% rename(n_comparison = n),
      by = c("comparison_category" = "category")
    ) %>%
    mutate(N_total = paste0("ref=", n_ref, ", comp=", n_comparison))

  all_results[[outcome]] <- cont_summary
}

# Combine all into one data frame
final_results <- bind_rows(all_results)

# Optional: select and rename columns for clarity
final_results_clean <- final_results %>%
  select(outcome, contrast, N_total, estimate, SE, df, t.ratio, p.value, adj.p.value = p.value) %>%
  arrange(outcome, adj.p.value)

print(final_results_clean)

output_file <- file.path("Fig3/3BC/251213_aureusstatscrossfeed_Nadded.csv")

# Remake AUC figure with significance asterisks

final_results_clean$category <- sub("-.*", "", final_results_clean$contrast)

sig_results <- final_results_clean %>%
  filter(outcome == "auc") %>%
  mutate(
    category = trimws(as.character(category)),
    sig_stars = case_when(
      adj.p.value < 0.001 ~ "***",
      adj.p.value < 0.01 ~ "**",
      adj.p.value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    sig_label = ifelse(sig_stars == "",
      "",
      paste0("bold('", sig_stars, "')*'\n'*bold(' ", signif(adj.p.value, 2), "')")
    )
  ) %>%
  select(category, sig_label)

# 2. Constant y-position above all boxes
global_y <- max(stats_taxon$auc, na.rm = TRUE) * 1.15

sig_positions <- stats_taxon %>%
  mutate(category = trimws(as.character(category))) %>%
  distinct(category) %>%
  left_join(sig_results, by = "category") %>%
  mutate(y_pos = global_y)


sig_positions <- sig_positions %>% filter(category != "Empty" & category != "MMG" & category != "MM")

p4 <- P1 +
  geom_text(
    data = sig_positions %>% filter(sig_label != ""),
    aes(x = category, y = y_pos, label = sig_label),
    inherit.aes = FALSE,
    parse = TRUE, # lets bold() render
    angle = 0, # upright text
    vjust = 0, # above the axis
    size = 4.5
  ) +
  theme(
    axis.text.y = element_text(size = 16, face = "bold") # bigger, bold category labels
  )


tmp.fname <- paste0("Crossfeeding/Analysis/250918Figures/", Sys.Date(), "AUCbymedia.pdf")

######################################
################### COGNATE####################
########################################

main_directory <- "Crossfeeding/CognateCombinations"

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
  datafiles <- all_csv_files[!grepl("media|strain|deepwell|stats|_growthnormalizedtomin|stickers", all_csv_files, ignore.case = TRUE)]
  mediafiles <- list.files(path = dir, pattern = "media.*\\.csv$", full.names = TRUE)
  strainfiles <- list.files(path = dir, pattern = "strains.*\\.csv$", full.names = TRUE)

  imported_widedata_file <- read_wides(files = datafiles, startrow = 1)

  imported_wides_now_tidy_file <- trans_wide_to_tidy(
    wides = imported_widedata_file,
    id_cols = c("file", "T° 600", "Time")
  )

  my_design_file <- import_blockdesigns(files = c(mediafiles, strainfiles), block_names = c("Media", "Strain_letters"), startrow = 2)
  ex_dat_mrg_file <- merge_dfs(imported_wides_now_tidy_file, my_design_file)

  ex_dat_mrg_file$Time <- time_length(hms(ex_dat_mrg_file$Time), unit = "hour")

  ex_dat_mrg_file$Well <- factor(ex_dat_mrg_file$Well, levels = paste0(rep(LETTERS[1:8], each = 12), 1:12))

  p0 <- ggplot(data = ex_dat_mrg_file, aes(x = Time, y = Measurements)) +
    geom_line() +
    facet_wrap(~Well, nrow = 8, ncol = 12)

  plot_file <- file.path(dir, "/figures_stats/", paste0(Sys.Date(), tools::file_path_sans_ext(basename(datafiles[1])), "_growth_plot.png"))

  ex_dat_mrg_file <- mutate(group_by(ex_dat_mrg_file, Well, Media, Strain_letters), deriv_percap5 = calc_deriv(x = Time, y = Measurements, percapita = TRUE, blank = (first_minima(y = Measurements, x = Time, return = "y")), window_width_n = 5, trans_y = "log"))

  ex_dat_mrg_sum_file_minmax <-
    summarize(group_by(ex_dat_mrg_file, Well, Media, Strain_letters),
      max_dens = max_gc(Measurements, na.rm = TRUE),
      max_time = extr_val(Time, which_max_gc(Measurements)),
      min_dens = first_minima(y = Measurements, x = Time, return = "y"),
      min_time = first_minima(y = Measurements, x = Time, return = "x"),
      max_percap = max_gc(deriv_percap5, na.rm = TRUE), max_percap_time = extr_val(Time, which_max_gc(deriv_percap5)), doub_time = doubling_time(y = max_percap),
      auc = auc(x = Time, y = Measurements, blank = min_dens)
    )

  ex_dat_mrg_sum_file_minmax$maxminusmin <- ex_dat_mrg_sum_file_minmax$max_dens - ex_dat_mrg_sum_file_minmax$min_dens

  ex_dat_mrg_sum_file_minmax$filepath <- file.path(dir, "/figures_stats/", paste0(tools::file_path_sans_ext(basename(datafiles[1])), "_stats.csv"))

  output_file <- file.path(dir, "/figures_stats/", paste0(tools::file_path_sans_ext(basename(datafiles[1])), "_stats.csv"))

  write.csv(ex_dat_mrg_sum_file_minmax, output_file, row.names = FALSE)
}

# Catching bad wells/contamination

# 2025-07-31250123_201_203_204_205_206_207_208_X_0043_growth_plot.png D9, D11 severe noise
# 2025-07-31250307_120_123_128_X_55_56___194_X_60___220_X_77_growth_plot.png entire 6 column
# 2025-07-31250309_128_X_55_56___135_137_X_84_growth_plot entire 5 and 11 columns
# 2025-07-31250310_193_X_567_568_569_570_growth_plot.png H2 has noise

####### To merge stats files
master_dir <- "~/Library/CloudStorage/Box-Box/2024_Dome2.0/Crossfeeding/CognateCombinations"
#
# # Get all subdirectories within the master directory
all_subdirs <- list.dirs(path = master_dir, full.names = TRUE, recursive = TRUE)
all_subdirs <- grep("figures_stats", all_subdirs, value = TRUE)
#
# # Initialize an empty list to store dataframes
merged_data <- list()
#
# # Loop through each subdirectory to find and merge '_stats' files
for (subdir in all_subdirs) {
  #
  #   # List all files in the current subdirectory that contain '_stats' in their names
  stats_files <- list.files(path = subdir, pattern = "_stats.*\\.csv$", full.names = TRUE)
  #
  #   # If no '_stats' files are found, skip this directory
  if (length(stats_files) == 0) {
    message(paste("No '_stats' files found in directory:", subdir))
    next
  }
  #
  #   # Read each '_stats' file and add it to the list
  for (file in stats_files) {
    file_data <- read.csv(file)
    merged_data[[length(merged_data) + 1]] <- file_data
  }
}
#
# # If there are no '_stats' files to merge, exit
if (length(merged_data) == 0) {
  stop("No '_stats' files found in any of the subdirectories.")
}


# Merge all dataframes in the list by columns
final_merged_data <- Reduce(function(x, y) merge(x, y, all = TRUE), merged_data)


# 2025-07-31250123_201_203_204_205_206_207_208_X_0043_growth_plot.png D9, D11 severe noise

my_values <- c("D9", "D11")
df_filtered_1 <- final_merged_data %>% filter(!(str_detect(filepath, "250123_201_203_204_205_206_207_208_X_0043") & Well %in% my_values))


# 2025-07-31250307_120_123_128_X_55_56___194_X_60___220_X_77_growth_plot.png entire 6 column
my_values <- c("A6", "B6", "C6", "D6", "E6", "F6", "G6", "H6")
df_filtered_2 <- df_filtered_1 %>% filter(!(str_detect(filepath, "250307_120_123_128_X_55_56___194_X_60___220_X_77") & Well %in% my_values))


# 2025-07-31250309_128_X_55_56___135_137_X_84_growth_plot entire 5 and 11 columns
my_values <- c("A11", "B11", "C11", "D11", "E11", "F11", "G11", "H11", "A5", "B5", "C5", "D5", "E5", "F5", "G5", "H5")
df_filtered_3 <- df_filtered_2 %>% filter(!(str_detect(filepath, "250309_128_X_55_56___135_137_X_84") & Well %in% my_values))

# 2025-07-31250310_193_X_567_568_569_570_growth_plot.png H2 has noise

my_values <- c("H2")
cleanedstatsdata <- df_filtered_3 %>% filter(!(str_detect(filepath, "250310_193_X_567_568_569_570") & Well %in% my_values))


taxonomy <- read.csv("FinalGTDBTK/250917_taxonomymap.csv")

taxonomy$Species <- paste0(toupper(substr(taxonomy$Species, 1, 1)), substr(taxonomy$Species, 2, nchar(taxonomy$Species)))

taxonomy$Media <- sub("^MMM_0", "", taxonomy$Media)
taxonomy$Media <- sub("^0", "", taxonomy$Media)

cleanedstatsdata$Media <- sub("^MMM_", "", cleanedstatsdata$Media)
cleanedstatsdata$Strain_letters <- sub("^DOME_", "", cleanedstatsdata$Strain_letters)
cleanedstatsdata$Strain_letters <- sub("^0", "", cleanedstatsdata$Strain_letters)
cleanedstatsdata$Strain_letters <- sub("^0", "", cleanedstatsdata$Strain_letters)

names(taxonomy)[names(taxonomy) == "Media"] <- "Number"

stats_taxon <- merge(cleanedstatsdata, taxonomy, by.x = "Media", by.y = "Number", all.x = TRUE)

names(stats_taxon)[names(stats_taxon) == "Species"] <- "media_species"

stats_taxon <- merge(stats_taxon, taxonomy, by.x = "Strain_letters", by.y = "Number", all.x = TRUE)
names(stats_taxon)[names(stats_taxon) == "Species"] <- "Strain_species"


stats_taxon <- stats_taxon %>%
  mutate(mediacategory = case_when(
    Media == "Empty" | Strain_letters == "Empty" ~ "Empty",
    !is.na(media_species) ~ media_species,
    TRUE ~ Media
  ))

stats_taxon <- stats_taxon %>%
  mutate(straincategory = case_when(
    Media == "Empty" | Strain_letters == "Empty" ~ "Empty",
    !is.na(Strain_species) ~ Strain_species,
    TRUE ~ Media
  ))

unique_combos <- unique(stats_taxon[c("media_species", "Strain_species")])
unique_combos <- na.omit(unique_combos)


stats_taxon_filtered <- stats_taxon %>%
  filter(straincategory != "Bacillus licheniformis")


ref_level <- "MMM"
outcomes <- c("auc", "maxminusmin", "max_percap")

stats_taxon_filtered2 <- stats_taxon_filtered %>%
  filter(straincategory != "Empty")

stats_taxon_filtered2 <- stats_taxon_filtered2 %>%
  filter(!is.na(media_species) | mediacategory %in% c("MMM", "MMMG"))

stats_taxon_filtered2 <- stats_taxon_filtered2 %>%
  filter(!is.na(Strain_species))

# Initialize results list
all_results <- list()

# Loop through each unique straincategory
for (strain in unique(stats_taxon_filtered2$straincategory)) {
  # Subset data for current strain
  strain_data <- stats_taxon_filtered2 %>%
    filter(straincategory == strain)

  # Ensure category is a factor and relevel
  strain_data$mediacategory <- relevel(factor(strain_data$mediacategory), ref = ref_level)

  # Loop through outcomes
  for (outcome in outcomes) {
    # Define formula and fit model
    formula <- as.formula(paste0(outcome, " ~ mediacategory"))
    model <- lm(formula, data = strain_data)

    # Compute emmeans and contrasts
    em <- emmeans(model, "mediacategory")
    ref_index <- which(levels(strain_data$mediacategory) == ref_level)

    cont <- contrast(em, method = "trt.vs.ctrl", ref = ref_index)

    cont_summary <- summary(cont, adjust = "BH") %>%
      as.data.frame() %>%
      mutate(
        outcome = outcome,
        straincategory = strain
      )

    # Calculate N for each mediacategory within this strain and outcome
    n_by_category <- strain_data %>%
      filter(!is.na(.data[[outcome]])) %>%
      count(mediacategory, name = "n")

    # Add N to contrast summary
    cont_summary <- cont_summary %>%
      mutate(
        ref_category = ref_level,
        comparison_category = gsub(paste0(" - ", ref_level), "", contrast)
      ) %>%
      left_join(n_by_category %>% rename(n_ref = n),
        by = c("ref_category" = "mediacategory")
      ) %>%
      left_join(n_by_category %>% rename(n_comparison = n),
        by = c("comparison_category" = "mediacategory")
      ) %>%
      mutate(N_total = paste0("ref=", n_ref, ", comp=", n_comparison))

    # Store result
    all_results[[paste(strain, outcome, sep = "_")]] <- cont_summary
  }
}

# Combine all into one data frame
final_results <- bind_rows(all_results)

# Optional: clean and arrange
final_results_clean <- final_results %>%
  select(straincategory, outcome, contrast, N_total, estimate, SE, df, t.ratio, p.value, adj.p.value = p.value) %>%
  arrange(straincategory, outcome, adj.p.value)

print(final_results_clean)

output_file <- file.path("Fig3/3BC/251213_cognatestatscrossfeed_Nadded.csv")

sig_results <- final_results_clean %>%
  filter(
    outcome == "auc",
    adj.p.value < 0.05,
    estimate > 0
  ) %>%
  mutate(mediacategory = gsub(" - MMM", "", contrast)) %>%
  filter(mediacategory != "MMMG") %>% # Exclude MMMG
  select(straincategory, mediacategory)
mmm_controls <- sig_results %>%
  distinct(straincategory) %>%
  mutate(mediacategory = "MMM")
sig_combos_with_ref <- bind_rows(sig_results, mmm_controls)

sig_auc_data <- stats_taxon %>%
  filter(straincategory != "Empty") %>%
  semi_join(sig_combos_with_ref, by = c("straincategory", "mediacategory")) %>%
  mutate(
    combo = paste(straincategory, mediacategory, sep = " | "),
    straincategory = factor(straincategory) # ensure proper ordering
  )


sig_auc_data <- sig_auc_data %>%
  mutate(
    media_genus = word(media_species, 1), # get genus from species
    genus_color = ifelse(media_genus %in% names(genus.colors),
      genus.colors[media_genus],
      "black"
    )
  )

# Step 4: Plot with straincategory-based clustering (faceting)
P_sig_auc <- sig_auc_data %>%
  mutate(mediacategory = fct_reorder(mediacategory, auc, .fun = median, .desc = FALSE)) %>%
  ggplot(aes(x = mediacategory, y = auc, fill = genus_color)) +
  coord_flip() +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 0.5, alpha = 0.7, aes(color = "#899499")) +
  facet_wrap(~straincategory, scales = "free_y", ncol = 1) +
  scale_fill_identity() + # trust the hex codes directly
  scale_color_identity() +
  labs(
    title = "Significant AUC Increases by Strain in Cognate Combinations",
    x = "Media Condition/Supernatant",
    y = "Area Under Curve (AUC)"
  ) +
  theme_pub()

plot_file <- paste0("Crossfeeding/Analysis/250730Figures/", Sys.Date(), "cognateauc.png")


# 1. Create 'category' from contrast
final_results_clean <- final_results_clean %>%
  mutate(category = sub("-.*", "", contrast))

# 2. Build significance labels
sig_results <- final_results_clean %>%
  filter(outcome == "auc", estimate > 0) %>%
  mutate(
    mediacategory = trimws(as.character(category)),
    sig_stars = case_when(
      adj.p.value < 0.001 ~ "***",
      adj.p.value < 0.01 ~ "**",
      adj.p.value < 0.05 ~ "*",
      TRUE ~ ""
    ),
    sig_label = ifelse(sig_stars == "",
      NA_character_,
      paste0("bold('", sig_stars, " ", signif(adj.p.value, 2), "')")
    )
  ) %>%
  select(straincategory, mediacategory, sig_label) %>%
  filter(!mediacategory %in% c("Empty", "MMG", "MM", "MMMG")) %>%
  filter(!is.na(sig_label))

# 3. Compute a single global y-position above all boxplots
global_y <- max(sig_auc_data$auc, na.rm = TRUE) * 1.15

sig_results <- sig_results %>%
  mutate(y_pos = global_y)

# 4. Add to faceted plot
p5 <- P_sig_auc +
  geom_text(
    data = sig_results,
    aes(x = mediacategory, y = y_pos, label = sig_label),
    inherit.aes = FALSE,
    parse = TRUE,
    angle = 0,
    vjust = 0,
    size = 5.5
  )

p5 <- p5 +
  # Expand space between axis limits and data
  scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) + # 5% below, 15% above
  scale_x_discrete(expand = expansion(add = 1)) + # adds space on sides
  theme(
    panel.spacing = unit(1, "cm"), # gap between facets
    panel.background = element_rect(
      fill = "white",
      color = "grey80",
      size = 1
    )
  ) +
  theme(
    axis.text.y = element_text(size = 16, face = "bold") # bigger, bold category labels
  )
tmp.fname <- paste0("~/Library/CloudStorage/Box-Box/2024_Dome2.0/Crossfeeding/Analysis/250918Figures/", Sys.Date(), "sigcognateAUCbymedia.pdf")


all_auc_results <- final_results_clean %>%
  filter(outcome == "auc") %>%
  mutate(
    mediacategory = gsub(" - MMM", "", contrast),
    significance = case_when(
      adj.p.value < 0.05 & estimate > 0 ~ "Upregulated (sig)",
      adj.p.value < 0.05 & estimate < 0 ~ "Downregulated (sig)",
      TRUE ~ "Not Significant"
    )
  )
plot_data_tested <- stats_taxon %>%
  filter(straincategory != "Empty") %>%
  inner_join(all_auc_results %>% select(straincategory, mediacategory, significance),
    by = c("straincategory", "mediacategory")
  ) %>%
  mutate(combo = paste(straincategory, mediacategory, sep = " | "))

# MMM data per strain (not in model results but should be plotted)
mmm_data <- stats_taxon %>%
  filter(mediacategory == "MMM", straincategory != "Empty") %>%
  mutate(
    significance = "MMM (reference)",
    combo = paste(straincategory, mediacategory, sep = " | ")
  )

# Combine both
plot_data_all <- bind_rows(plot_data_tested, mmm_data)


P_all_auc <- plot_data_all %>%
  mutate(combo = fct_reorder(combo, auc, .fun = median)) %>%
  ggplot(aes(x = combo, y = auc, color = significance)) +
  coord_flip() +
  geom_boxplot(outlier.shape = NA, aes(fill = significance), alpha = 0.2, color = NA) +
  geom_jitter(width = 0.2, size = 0.6, alpha = 0.7) +
  scale_color_manual(values = c(
    "Upregulated (sig)" = "#1b9e77",
    "Downregulated (sig)" = "#d95f02",
    "Not Significant" = "grey50"
  )) +
  scale_fill_manual(values = c(
    "Upregulated (sig)" = "#1b9e77",
    "Downregulated (sig)" = "#d95f02",
    "Not Significant" = "grey80"
  )) +
  facet_wrap(~straincategory, scales = "free_y", ncol = 1) +
  theme_classic() +
  labs(
    title = "AUC by Strain and Media — All Tested Comparisons",
    x = "Strain | Media",
    y = "Area Under Curve (AUC)",
    color = "Significance",
    fill = "Significance"
  ) +
  theme(axis.text.y = element_text(size = 7))
