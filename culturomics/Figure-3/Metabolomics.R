# Load required libraries
library(gcplyr)
library(lubridate)
library(ggplot2)
library(dplyr)
library(readxl)
library(purrr)
library(stringr)
library(tidyr)
library(forcats)
library(tidyverse)
library(broom)
library(pheatmap)
library(xcms)
library(magrittr)
library(Spectra)
library(jrtools)
library(MSnbase)
source("utilities.R")
library(cowplot)
library(eulerr)
library(ggforce)
library(data.table)
library(parallel)
library(igraph)
library(doParallel)
library(foreach)
# library(extrafont)

# Set working directory
setwd("Metabolomics/HigResRun")

#######################################
# PART 1: DATA LOADING AND PREPROCESSING
#######################################

# Load sample mapping and raw metabolomics data
samplemap <- read_excel("idx_testrun_samples.xlsx")
rawdata <- read_excel("250324_mucins_HILIC_neg_untargeted_compounds.xlsx")

# Simplify column names by extracting sample identifiers
new_colnames <- sapply(colnames(rawdata), function(colname) {
  match <- sub(
    pattern = "^Area: 250324_mucins_HILIC_neg_untargeted_(.*?)\\.raw.*$",
    replacement = "\\1",
    x = colname
  )
  # Only change if the match actually occurred (i.e., colname was transformed)
  if (match != colname) {
    return(match)
  } else {
    return(colname)
  }
})

# Assign new column names
colnames(rawdata) <- new_colnames

# Fix specific column name issue (MR093 has timestamp appended)
colnames(rawdata)[colnames(rawdata) == "MR093_20250325120021"] <- "MR093"

# Keep only the first 102 columns (area measurements)
rawdata <- rawdata[, 1:102]

# Rename columns based on samplemap
for (i in seq_len(nrow(samplemap))) {
  sample_id <- samplemap$sample[i]
  replacement_name <- samplemap$name[i]

  # Find which column(s) match the sample ID
  match_index <- which(colnames(rawdata) == sample_id)

  # Replace the column name if a match is found
  if (length(match_index) > 0) {
    colnames(rawdata)[match_index] <- replacement_name
  }
}

# Log10 transform all area measurements (columns 20 onwards)
rawdata[, 20:ncol(rawdata)] <- log10(rawdata[, 20:ncol(rawdata)])

# Strip off 'g' suffix from column names (blue/green media designation)
colnames(rawdata) <- sub("g$", "", colnames(rawdata))

#######################################
# PART 2: CALCULATE REPLICATE AVERAGES
#######################################

all_cols <- colnames(rawdata)

# Strip off the last character if it's A, B, or C (technical replicates)
base_names <- sub("[ABC]$", "", all_cols)

# Find unique base names that had A/B/C variants
unique_bases <- unique(base_names[duplicated(base_names)])

# Initialize result dataframe with original columns
averaged_data <- rawdata

# Loop over each base name and calculate row means for replicates
for (base in unique_bases) {
  # Identify columns like glucoseA/B/C
  matching_cols <- grep(paste0("^", base, "[ABC]$"), all_cols, value = TRUE)

  if (length(matching_cols) > 1) {
    # Compute row means and add as "base_average"
    averaged_data[[paste0(base, "_average")]] <- rowMeans(rawdata[, matching_cols], na.rm = TRUE)
  }
}

avg_cols <- colnames(averaged_data)

# Extract base names for *_Pr_average and *_PoT_average columns
pr_bases <- sub("Pr_average$", "", grep("Pr_average$", avg_cols, value = TRUE))
pot_bases <- sub("PoT_average$", "", grep("PoT_average$", avg_cols, value = TRUE))

# Identify common base names that have both Pr (Pre) and PoT (Post) averages
filter_bases <- intersect(pr_bases, pot_bases)

# Confirm MMMPr_average exists (control medium)
if (!"MMMPr_average" %in% avg_cols) stop("MMMPr_average column not found.")

#######################################
# PART 3: ANOVA AND TUKEY-KRAMER TESTS
#######################################

results_list <- list()

# Loop through each base to identify significantly consumed metabolites
# Uses one-way ANOVA with Tukey-Kramer post-hoc tests and BH FDR correction
for (base in filter_bases) {
  message("Processing base: ", base)

  basePr_col <- paste0(base, "Pr_average")
  basePoT_col <- paste0(base, "PoT_average")
  mmmpr_col <- "MMMPr_average"

  # Sanity check that columns exist
  if (!all(c(basePr_col, basePoT_col, mmmpr_col) %in% colnames(averaged_data))) next

  # Step 1: Magnitude filter (> 0.5 log units over MMMPr)
  # This filters for metabolites that increase during bacterial growth
  keep_rows <- averaged_data[[basePr_col]] > averaged_data[[mmmpr_col]] + 0.5

  filtered_IDs <- averaged_data$ID[keep_rows]

  if (length(filtered_IDs) == 0) next # skip if nothing passes filter

  # Step 2: Set up sample columns for the three conditions
  basePr_reps <- grep(paste0("^", base, "Pr[ABC]$"), colnames(rawdata), value = TRUE)
  basePoT_reps <- grep(paste0("^", base, "PoT[ABC]$"), colnames(rawdata), value = TRUE)
  MMMPr_reps <- grep("^MMMPr[ABC]?$", colnames(rawdata), value = TRUE)

  if (length(basePr_reps) < 2 | length(basePoT_reps) < 2 | length(MMMPr_reps) < 2) next

  all_samples <- c(basePr_reps, basePoT_reps, MMMPr_reps)

  # Storage for ANOVA results
  anova_pvals <- numeric(length(filtered_IDs))
  tukey_results <- vector("list", length(filtered_IDs))
  mean_diffs <- matrix(NA, nrow = length(filtered_IDs), ncol = 3)
  colnames(mean_diffs) <- c("BasePr_vs_MMMPr", "BasePr_vs_BasePoT", "MMMPr_vs_BasePoT")

  # Run ANOVA for each filtered feature
  for (j in seq_along(filtered_IDs)) {
    metab_id <- filtered_IDs[j]
    row_index <- which(rawdata$ID == metab_id)

    y <- as.numeric(rawdata[row_index, all_samples])
    cond <- c(
      rep("BasePr", length(basePr_reps)),
      rep("BasePoT", length(basePoT_reps)),
      rep("MMMPr", length(MMMPr_reps))
    )

    df <- data.frame(
      LogArea = y,
      Condition = factor(cond, levels = c("BasePr", "MMMPr", "BasePoT"))
    )

    # Fit one-way ANOVA
    anova_model <- try(aov(LogArea ~ Condition, data = df), silent = TRUE)

    if (inherits(anova_model, "try-error")) next

    # Extract ANOVA p-value
    anova_summary <- summary(anova_model)
    anova_pvals[j] <- anova_summary[[1]]$`Pr(>F)`[1]

    # Calculate mean differences for reporting
    means <- tapply(df$LogArea, df$Condition, mean, na.rm = TRUE)
    mean_diffs[j, "BasePr_vs_MMMPr"] <- means["BasePr"] - means["MMMPr"]
    mean_diffs[j, "BasePr_vs_BasePoT"] <- means["BasePr"] - means["BasePoT"]
    mean_diffs[j, "MMMPr_vs_BasePoT"] <- means["MMMPr"] - means["BasePoT"]
  }

  # Apply BH correction to ANOVA p-values
  adj_anova_pvals <- p.adjust(anova_pvals, method = "BH")

  # Find significant features (FDR < 0.05)
  sig_anova_rows <- which(adj_anova_pvals < 0.05)

  if (length(sig_anova_rows) == 0) next # skip if no significant ANOVA results

  # Run Tukey-Kramer post-hoc tests for significant features
  tukey_sig <- matrix("", nrow = length(sig_anova_rows), ncol = 3)
  colnames(tukey_sig) <- c("BasePr_vs_MMMPr_sig", "BasePr_vs_BasePoT_sig", "MMMPr_vs_BasePoT_sig")

  tukey_pvals <- matrix(NA, nrow = length(sig_anova_rows), ncol = 3)
  colnames(tukey_pvals) <- c("BasePr_vs_MMMPr_p", "BasePr_vs_BasePoT_p", "MMMPr_vs_BasePoT_p")

  for (k in seq_along(sig_anova_rows)) {
    j <- sig_anova_rows[k]
    metab_id <- filtered_IDs[j]
    row_index <- which(rawdata$ID == metab_id)

    y <- as.numeric(rawdata[row_index, all_samples])
    cond <- c(
      rep("BasePr", length(basePr_reps)),
      rep("BasePoT", length(basePoT_reps)),
      rep("MMMPr", length(MMMPr_reps))
    )

    df <- data.frame(
      LogArea = y,
      Condition = factor(cond, levels = c("BasePr", "MMMPr", "BasePoT"))
    )

    anova_model <- aov(LogArea ~ Condition, data = df)
    tukey_test <- TukeyHSD(anova_model)

    tukey_df <- as.data.frame(tukey_test$Condition)

    # Extract p-values for each pairwise comparison
    tukey_pvals[k, "BasePr_vs_MMMPr_p"] <- tukey_df["MMMPr-BasePr", "p adj"]
    tukey_pvals[k, "BasePr_vs_BasePoT_p"] <- tukey_df["BasePoT-BasePr", "p adj"]
    tukey_pvals[k, "MMMPr_vs_BasePoT_p"] <- tukey_df["BasePoT-MMMPr", "p adj"]

    # Mark significance (* if p < 0.05, ns otherwise)
    tukey_sig[k, "BasePr_vs_MMMPr_sig"] <- ifelse(tukey_df["MMMPr-BasePr", "p adj"] < 0.05, "*", "ns")
    tukey_sig[k, "BasePr_vs_BasePoT_sig"] <- ifelse(tukey_df["BasePoT-BasePr", "p adj"] < 0.05, "*", "ns")
    tukey_sig[k, "MMMPr_vs_BasePoT_sig"] <- ifelse(tukey_df["BasePoT-MMMPr", "p adj"] < 0.05, "*", "ns")
  }

  # Output filtered significant hits
  final_IDs <- filtered_IDs[sig_anova_rows]
  filtered_df <- rawdata[rawdata$ID %in% final_IDs, c("ID", all_samples)]

  filtered_df$ANOVA_FDR <- adj_anova_pvals[sig_anova_rows]
  filtered_df$diff_BasePr_vs_MMMPr <- mean_diffs[sig_anova_rows, "BasePr_vs_MMMPr"]
  filtered_df$diff_BasePr_vs_BasePoT <- mean_diffs[sig_anova_rows, "BasePr_vs_BasePoT"]
  filtered_df$diff_MMMPr_vs_BasePoT <- mean_diffs[sig_anova_rows, "MMMPr_vs_BasePoT"]

  filtered_df <- cbind(filtered_df, tukey_pvals, tukey_sig)

  results_list[[base]] <- filtered_df
  assign(paste0(base, "_consumed_filtered"), filtered_df)

  message(sprintf("  Found %d significant features (ANOVA FDR < 0.05)", nrow(filtered_df)))
}

#######################################
# PART 4: CLASSIFY METABOLITES BY PATTERN
#######################################

all_results <- list()
plot_list <- list()
summ_list <- list()

# Load taxonomy data
taxonomy <- read.csv("~/Library/CloudStorage/Box-Box/2024_Dome2.0/Culturing/IsolateAnalysis/Tree/FinalGTDBTK/250917_taxonomymap.csv")

taxonomy <- taxonomy %>%
  mutate(
    BaseMatch = sub("^MMM_0*", "", Media) # remove leading zeros and MMM_ prefix
  )

# Define species of interest
species_to_plot <- c(
  "Arthrobacter_B koreensis", "Cellulosimicrobium funkei",
  "Nocardiopsis alba", "Shouchella clausii"
)

species_bases <- taxonomy %>%
  filter(Species %in% species_to_plot) %>%
  pull(BaseMatch)

# Loop only over bases that belong to species_of_interest AND have significant ANOVA results
for (base in intersect(filter_bases, species_bases)) {
  # Get species name for this base
  species <- taxonomy$Species[taxonomy$BaseMatch == base][1]
  if (is.na(species) | length(species) == 0) species <- base # fallback
  message("Processing base: ", base, " (Species: ", species, ")")

  # Check if this base has results from ANOVA
  consumed_df_name <- paste0(base, "_consumed_filtered")
  if (!exists(consumed_df_name)) {
    message("  No ANOVA results found for ", base)
    next
  }

  consumed_df <- get(consumed_df_name)

  if (nrow(consumed_df) == 0) {
    message("  No significant features for ", base)
    next
  }

  basePr_reps <- grep(paste0("^", base, "Pr[ABC]$"), colnames(rawdata), value = TRUE)
  basePoT_reps <- grep(paste0("^", base, "PoT[ABC]$"), colnames(rawdata), value = TRUE)
  MMMPr_reps <- grep("^MMMPr[ABC]?$", colnames(rawdata), value = TRUE)

  all_samples <- c(basePr_reps, basePoT_reps, MMMPr_reps)

  # Classify metabolites based on Tukey results
  # IMPORTANT: Increase+Decrease and IncreaseOnly are MUTUALLY EXCLUSIVE in this classification
  # - Increase+Decrease: BasePr significantly higher than BOTH MMMPr AND BasePoT
  #   (metabolite increases during growth, then decreases after)
  # - IncreaseOnly: BasePr significantly higher than MMMPr, but NOT significantly different from BasePoT
  #   (metabolite increases during growth and stays elevated)
  res_base <- consumed_df %>%
    mutate(Group = case_when(
      BasePr_vs_MMMPr_sig == "*" & diff_BasePr_vs_MMMPr >= 0 &
        BasePr_vs_BasePoT_sig == "*" & diff_BasePr_vs_BasePoT >= 0 ~ "Increase+Decrease",
      BasePr_vs_MMMPr_sig == "*" & diff_BasePr_vs_MMMPr >= 0 ~ "IncreaseOnly",
      TRUE ~ NA_character_
    )) %>%
    filter(!is.na(Group)) %>%
    select(
      ID, Group, ANOVA_FDR,
      diff_BasePr_vs_MMMPr, BasePr_vs_MMMPr_p, BasePr_vs_MMMPr_sig,
      diff_BasePr_vs_BasePoT, BasePr_vs_BasePoT_p, BasePr_vs_BasePoT_sig
    )

  if (nrow(res_base) == 0) {
    message("  No metabolites meeting group criteria for ", base)
    next
  }

  all_results[[base]] <- res_base

  # Build plot dataframe from replicate columns
  plot_df <- rawdata %>%
    filter(ID %in% res_base$ID) %>%
    select(ID, all_of(all_samples)) %>%
    pivot_longer(-ID, names_to = "Sample", values_to = "LogArea") %>%
    mutate(
      Time = case_when(
        grepl("MMMPr", Sample) ~ "MMM",
        grepl("Pr", Sample) & !grepl("MMMPr", Sample) ~ "Pre",
        grepl("PoT", Sample) ~ "Post",
        TRUE ~ NA_character_
      ),
      Base = base
    ) %>%
    left_join(res_base %>% select(ID, Group), by = "ID") %>%
    mutate(Time = factor(Time, levels = c("MMM", "Pre", "Post")))

  # Summarize replicates per metabolite × time (mean ± SE)
  summ_df <- plot_df %>%
    group_by(ID, Base, Group, Time) %>%
    summarise(
      meanLog = mean(LogArea, na.rm = TRUE),
      seLog = ifelse(sum(!is.na(LogArea)) > 1,
        sd(LogArea, na.rm = TRUE) / sqrt(sum(!is.na(LogArea))),
        NA_real_
      ),
      nRep = sum(!is.na(LogArea)),
      .groups = "drop"
    ) %>%
    group_by(ID, Base) %>%
    mutate(
      Pre_mean = if (any(Time == "Pre")) meanLog[Time == "Pre"] else NA_real_,
      FoldChange = ifelse(!is.na(Pre_mean), meanLog - Pre_mean, NA_real_)
    ) %>%
    ungroup()

  # Store per-base summarized metabolite data
  summ_list[[base]] <- summ_df
}

#######################################
# PART 5: CREATE SPAGHETTI PLOTS
#######################################

# Bind all summaries collected in the loop
all_summ <- bind_rows(summ_list, .id = "BaseLoop")

# Join taxonomy to bring in Species names
all_summ <- all_summ %>%
  left_join(taxonomy %>% select(Species, BaseMatch),
    by = c("Base" = "BaseMatch")
  )

all_summ_out <- merge(all_summ, rawdata, by = "ID")

output_file <- file.path("/Users/gdlab/Library/CloudStorage/Box-Box/2024_Dome2.0/Manuscript/FinalizedFiguresAndScripts/Fig3/3D/251129_metabolomicstatscrossfeed.csv")

# write.csv(all_summ_out, output_file, row.names = FALSE)

# Function to generate spaghetti plot for a species
plot_species_spaghetti <- function(species_of_interest, df) {
  # Filter for the species
  df_species <- df %>% filter(Species == species_of_interest)

  if (nrow(df_species) == 0) {
    message("No data for species: ", species_of_interest)
    return(NULL)
  }

  # Split by Group
  group_list <- split(df_species, df_species$Group)

  plot_list <- lapply(names(group_list), function(gr) {
    df_group <- group_list[[gr]]

    # Calculate n per Base (for legend labels)
    n_per_base <- df_group %>%
      group_by(Base) %>%
      summarise(n = n_distinct(ID), .groups = "drop")

    # Merge n into plotting df
    df_group <- df_group %>%
      left_join(n_per_base, by = "Base") %>%
      mutate(
        Base_n = paste0(Base, " (n=", n, ")"),
        ID_base = paste(ID, Base, sep = "_")
      ) # unique grouping per base+ID

    # Make the plot
    p <- ggplot(df_group, aes(x = Time, y = FoldChange, color = Base_n)) +
      geom_line(aes(group = ID_base), alpha = 0.025, size = 0.6) +
      geom_point(aes(group = ID_base), size = 1, alpha = 0.6) +
      stat_summary(fun = mean, geom = "line", aes(group = Base_n), size = 1.2, alpha = 0.9) +
      stat_summary(fun = mean, geom = "point", aes(group = Base_n), size = 2, alpha = 0.9) +
      theme_minimal(base_size = 14) +
      labs(
        title = paste0(species_of_interest, " — ", gr),
        y = "FoldChange vs Pre",
        color = "Base (n)"
      ) +
      scale_x_discrete(limits = c("MMM", "Pre", "Post"))

    return(p)
  })

  names(plot_list) <- names(group_list)
  return(plot_list)
}

# Generate plots for all species
species_list <- c(
  "Arthrobacter_B koreensis",
  "Cellulosimicrobium funkei",
  "Nocardiopsis alba",
  "Shouchella clausii"
)

plot_list_species <- lapply(species_list, function(sp) plot_species_spaghetti(sp, all_summ))
names(plot_list_species) <- species_list

# Save plots
# plot_dir <- "/Users/gdlab/Library/CloudStorage/Box-Box/2024_Dome2.0/Metabolomics/HigResRun/251030_JRedits_metabolomics/Figures"

# for (sp in names(plot_list_species)) {
#   group_plots <- plot_list_species[[sp]]
#   if (is.null(group_plots)) next
#
#   for (gr in names(group_plots)) {
#     safe_sp <- gsub("[^a-zA-Z0-9_-]", "_", sp)
#     safe_gr <- gsub("[^a-zA-Z0-9_-]", "_", gr)
#     plot_file <- file.path(
#       plot_dir,
#       paste0(Sys.Date(), "_", safe_sp, "_", safe_gr, "_ANOVA_plot.pdf")
#     )
#
#     cairo_pdf(plot_file, width = 9, height = 6)
#     print(group_plots[[gr]])
#     dev.off()
#   }
# }

#######################################
# PART 6: DEFINE SPECIES-LEVEL METABOLITE LISTS
# This is the critical step that ensures consistency across all figures
#######################################

base_to_species <- taxonomy %>%
  mutate(BaseMatch = sub("^MMM_0*", "", Media)) %>%
  select(BaseMatch, Species) %>%
  distinct()

species_increase_only <- list()
species_increase_decrease <- list()

for (base in names(all_results)) {
  result_df <- all_results[[base]]

  species <- base_to_species$Species[base_to_species$BaseMatch == base][1]

  if (is.na(species) || length(species) == 0) next

  # Get IncreaseOnly IDs (ONLY those labeled "IncreaseOnly", NOT including Increase+Decrease)
  increase_only_ids <- result_df %>%
    filter(Group == "IncreaseOnly") %>%
    pull(ID)

  # Get Increase+Decrease IDs
  increase_decrease_ids <- result_df %>%
    filter(Group == "Increase+Decrease") %>%
    pull(ID)

  # Add to species lists
  if (length(increase_only_ids) > 0) {
    if (species %in% names(species_increase_only)) {
      species_increase_only[[species]] <- c(species_increase_only[[species]], increase_only_ids)
    } else {
      species_increase_only[[species]] <- increase_only_ids
    }
  }

  if (length(increase_decrease_ids) > 0) {
    if (species %in% names(species_increase_decrease)) {
      species_increase_decrease[[species]] <- c(species_increase_decrease[[species]], increase_decrease_ids)
    } else {
      species_increase_decrease[[species]] <- increase_decrease_ids
    }
  }
}

# Get unique IDs per species
species_increase_only <- lapply(species_increase_only, unique)
species_increase_decrease <- lapply(species_increase_decrease, unique)

# Print summary
message("\n=== Metabolite Counts by Group ===")
for (sp in names(species_increase_only)) {
  message(sp, ":")
  message("  IncreaseOnly (ONLY): ", length(species_increase_only[[sp]]), " metabolites")
  if (sp %in% names(species_increase_decrease)) {
    message("  Increase+Decrease: ", length(species_increase_decrease[[sp]]), " metabolites")
    message(
      "  Total (Combined): ",
      length(species_increase_only[[sp]]) + length(species_increase_decrease[[sp]])
    )
  }
}

#######################################
# PART 7: CREATE HEATMAP (INCREASE+DECREASE ONLY)
#######################################

# Extract only Increase+Decrease IDs from all results
increase_decrease_ids <- c()

for (base in names(all_results)) {
  result_df <- all_results[[base]]
  ids <- result_df %>%
    filter(Group == "Increase+Decrease") %>%
    pull(ID)
  increase_decrease_ids <- c(increase_decrease_ids, ids)
}

# Get unique IDs
unique_increase_decrease_ids <- unique(increase_decrease_ids)
message("Found ", length(unique_increase_decrease_ids), " unique Increase+Decrease metabolites")

# Subset rawdata to get only Increase+Decrease rows
increase_decrease_rawdata <- rawdata %>%
  filter(ID %in% unique_increase_decrease_ids)

# Average replicates for Pr, PoT, and MMMPr
all_cols <- colnames(increase_decrease_rawdata)

pr_cols <- grep("Pr[A-Z]$", all_cols, value = TRUE)
pot_cols <- grep("PoT[A-Z]$", all_cols, value = TRUE)
mmmpr_cols <- grep("^MMMPr[A-Z]$", all_cols, value = TRUE)

get_base <- function(colname, suffix) {
  sub(paste0(suffix, "[A-Z]$"), "", colname)
}
pr_bases <- unique(get_base(pr_cols, "Pr"))
pot_bases <- unique(get_base(pot_cols, "PoT"))
common_bases <- intersect(pr_bases, pot_bases)

# Average Pr and PoT replicates per sample base
averaged_matrix <- lapply(common_bases, function(base) {
  pr_group <- grep(paste0("^", base, "Pr[A-Z]$"), all_cols, value = TRUE)
  pot_group <- grep(paste0("^", base, "PoT[A-Z]$"), all_cols, value = TRUE)

  pr_mean <- rowMeans(increase_decrease_rawdata[, pr_group, drop = FALSE], na.rm = TRUE)
  pot_mean <- rowMeans(increase_decrease_rawdata[, pot_group, drop = FALSE], na.rm = TRUE)

  df <- data.frame(pr_mean, pot_mean)
  colnames(df) <- c(paste0(base, "_Pr"), paste0(base, "_PoT"))
  return(df)
})
heatmap_matrix <- do.call(cbind, averaged_matrix)

# Average MMMPr replicates into one column
if (length(mmmpr_cols) > 0) {
  mmmpr_avg <- rowMeans(increase_decrease_rawdata[, mmmpr_cols, drop = FALSE], na.rm = TRUE)
  heatmap_matrix <- cbind(heatmap_matrix, MMMPr = mmmpr_avg)
}

# Assign rownames
rownames(heatmap_matrix) <- increase_decrease_rawdata$ID

# Build sample metadata: Species and Sample Type
heatmap_colnames <- colnames(heatmap_matrix)
sample_bases <- str_extract(heatmap_colnames, "^[^_]+")

sample_map <- data.frame(
  ColName = heatmap_colnames,
  SampleBase = sample_bases,
  stringsAsFactors = FALSE
)

# Join with taxonomy
taxonomy <- taxonomy %>%
  mutate(
    SampleID = str_remove(Media, "^MMM_"),
    SampleID = str_replace(SampleID, "^0+", "")
  )

sample_map <- sample_map %>%
  left_join(taxonomy, by = c("SampleBase" = "SampleID"))

# Label sample types
sample_map$Type <- case_when(
  str_detect(sample_map$ColName, "_Pr$") ~ "Pr",
  str_detect(sample_map$ColName, "_PoT$") ~ "PoT",
  str_detect(sample_map$ColName, "^MMMPr") ~ "MMMPr",
  TRUE ~ "Other"
)

sample_map$Species[is.na(sample_map$Species)] <- "MMM"
sample_map$Type[is.na(sample_map$Type)] <- "MMMPr"

sample_map$Type <- factor(sample_map$Type, levels = c("Pr", "PoT", "MMMPr"))
sample_map$Species <- factor(sample_map$Species)

# Order columns: Species → Pr → PoT → MMMPr
sample_map_ordered <- sample_map %>%
  arrange(Species, Type)

ordered_cols <- sample_map_ordered$ColName
mmm_cols <- grep("^MMM", ordered_cols, value = TRUE)
other_cols <- setdiff(ordered_cols, mmm_cols)
ordered_cols <- c(other_cols, mmm_cols)

heatmap_matrix_ordered <- heatmap_matrix[, ordered_cols]

# Prepare annotation and colors
sample_map <- sample_map %>%
  mutate(Genus = word(Species, 1))

annotation_col <- data.frame(
  Type = sample_map$Type,
  Genus = sample_map$Genus
)
rownames(annotation_col) <- sample_map$ColName

annotation_col$Genus <- factor(annotation_col$Genus, levels = names(genus.colors))

data_genera <- unique(annotation_col$Genus)
data_genera[data_genera == "MMM"] <- "MMM_black"

used_genera <- intersect(names(genus.colors), data_genera)
missing_genera <- setdiff(data_genera, names(genus.colors))

genus.colors <- c(
  genus.colors[used_genera],
  setNames(rep("black", length(missing_genera)), missing_genera)
)

annotation_col$Genus[annotation_col$Genus == "MMM"] <- "MMM"

ann_colors <- list(
  Genus = genus.colors,
  Type = c("Pr" = "steelblue", "PoT" = "firebrick", "MMMPr" = "darkgrey")
)

# Scale data and plot heatmap
scaled_matrix <- t(scale(t(heatmap_matrix_ordered)))

# Remove MMMPr column
scaled_matrix <- subset(scaled_matrix, select = -c(MMMPr))

# Remove contaminated genomes
scaled_matrix <- scaled_matrix[, !(colnames(scaled_matrix) %in% c("185_Pr", "185_PoT", "220_Pr", "220_PoT"))]

# Create heatmap
p1 <- pheatmap(scaled_matrix,
  cluster_cols = FALSE,
  cluster_rows = TRUE,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  fontsize_row = 6,
  main = "Increase+Decrease metabolites (consumed then released)"
)

# Save heatmap
plot_file <- file.path(plot_dir, paste0(Sys.Date(), "_heatmap.pdf"))
# ggsave(filename = plot_file, plot = p1, width = 9, height = 9)

#######################################
# PART 8: CREATE EULER DIAGRAMS
#######################################

# Get colors from genus.colors
genera_decrease <- sapply(names(species_increase_decrease), function(x) word(x, 1))
venn_colors_decrease <- genus.colors[genera_decrease]

# Create Increase+Decrease Euler diagram
fit_decrease <- euler(species_increase_decrease)

message("\n=== Increase+Decrease Euler Diagram Fit ===")
message("Stress: ", round(fit_decrease$stress, 4))
# Stress: 5e-04
message("DiagError: ", round(fit_decrease$diagError, 4))
# DiagError: 0.0184

# Plot Increase+Decrease Euler
# pdf(file.path(plot_dir, paste0(Sys.Date(), "_IncreaseDecrease_Euler_proportional.pdf")),
#     width = 10, height = 8)
# plot(fit_decrease,
#      quantities = TRUE,
#      fills = venn_colors_decrease,
#      alpha = 0.6,
#      lty = 1,
#      lwd = 2,
#      fontsize = 10,
#      labels = list(font = 2, cex = 1.2),
#      main = "Percentage of total S. aureus consumed metabolites attributed to each species of MD",
#      cex = 1.5)
# dev.off()

# Display in viewer
# plot(fit_decrease,
#      quantities = TRUE,
#      fills = venn_colors_decrease,
#      alpha = 0.6,
#      lty = 1,
#      lwd = 2,
#      fontsize = 10,
#      labels = list(font = 2, cex = 1.2),
#      main = "Percentage of total S. aureus consumed metabolites attributed to each species of MD",
#      cex = 1.5)

# Print overlap statistics for Increase+Decrease
message("\n=== Increase+Decrease Overlap Statistics ===")
all_ids_decrease <- unique(unlist(species_increase_decrease))
message("Total unique metabolites: ", length(all_ids_decrease))

if (length(species_increase_decrease) >= 2) {
  shared_all <- Reduce(intersect, lapply(species_increase_decrease, function(x) x))
  message("Shared by ALL species: ", length(shared_all))

  # Pairwise overlaps
  species_names <- names(species_increase_decrease)
  for (i in 1:(length(species_names) - 1)) {
    for (j in (i + 1):length(species_names)) {
      sp1 <- species_names[i]
      sp2 <- species_names[j]
      overlap <- length(intersect(species_increase_decrease[[sp1]], species_increase_decrease[[sp2]]))
      message(sp1, " & ", sp2, ": ", overlap, " shared")
    }
  }
}


# --- Increase+Decrease Euler diagram ---

# Extract genus colors
genera_decrease <- sapply(names(species_increase_decrease), function(x) word(x, 1))
venn_colors_decrease <- genus.colors[genera_decrease]

# Fit Euler diagram
fit_decrease <- euler(species_increase_decrease)

# Compute total metabolites and percentage labels
all_ids_decrease <- unique(unlist(species_increase_decrease))
total_metabolites <- length(all_ids_decrease)
quantity_labels_decrease <- paste0(round(fit_decrease$original.values / total_metabolites * 100, 1), "%")

# Display in viewer
# plot(fit_decrease,
#      quantities = quantity_labels_decrease,
#      fills = venn_colors_decrease,
#      alpha = 0.6,
#      lty = 1,
#      lwd = 2,
#      fontsize = 10,
#      labels = list(font = 2, cex = 1.2),
#      main = "Percentage of total S. aureus consumed metabolites attributed to each species of MD",
#      cex = 1.5)

# # Save PDF
# pdf(file.path(plot_dir, paste0(Sys.Date(), "_IncreaseDecrease_Euler_percentage.pdf")),
#     width = 10, height = 8)
# plot(fit_decrease,
#      quantities = quantity_labels_decrease,
#      fills = venn_colors_decrease,
#      alpha = 0.6,
#      lty = 1,
#      lwd = 2,
#      fontsize = 10,
#      labels = list(font = 2, cex = 1.2),
#      main = "Percentage of total S. aureus consumed metabolites attributed to each species of MD",
#      cex = 1.5)
# dev.off()


# Create IncreaseOnly Euler diagram
genera_only <- sapply(names(species_increase_only), function(x) word(x, 1))
venn_colors_only <- genus.colors[genera_only]

fit_only <- euler(species_increase_only)

message("\n=== IncreaseOnly Euler Diagram Fit ===")
message("Stress: ", round(fit_only$stress, 4))
# Stress: 0.0114
message("DiagError: ", round(fit_only$diagError, 4))
# DiagError: 0.0374


# Plot IncreaseOnly Euler
# pdf(file.path(plot_dir, paste0(Sys.Date(), "_IncreaseOnly_Euler_proportional.pdf")),
#     width = 10, height = 8)
# plot(fit_only,
#      quantities = TRUE,
#      fills = venn_colors_only,
#      alpha = 0.6,
#      lty = 1,
#      lwd = 2,
#      fontsize = 10,
#      labels = list(font = 2, cex = 1.2),
#      main = "Metabolites Released by Mucin Degraders (Increase Only)",
#      cex = 1.5)
# #dev.off()

# Display in viewer
# plot(fit_only,
#      quantities = TRUE,
#      fills = venn_colors_only,
#      alpha = 0.6,
#      lty = 1,
#      lwd = 2,
#      fontsize = 10,
#      labels = list(font = 2, cex = 1.2),
#      main = "Metabolites Released by Mucin Degraders (Increase Only)",
#      cex = 1.5)

# Print overlap statistics for IncreaseOnly
message("\n=== IncreaseOnly Overlap Statistics ===")
all_ids_only <- unique(unlist(species_increase_only))
message("Total unique metabolites: ", length(all_ids_only))

if (length(species_increase_only) >= 2) {
  shared_all <- Reduce(intersect, lapply(species_increase_only, function(x) x))
  message("Shared by ALL species: ", length(shared_all))

  # Pairwise overlaps
  species_names <- names(species_increase_only)
  for (i in 1:(length(species_names) - 1)) {
    for (j in (i + 1):length(species_names)) {
      sp1 <- species_names[i]
      sp2 <- species_names[j]
      overlap <- length(intersect(species_increase_only[[sp1]], species_increase_only[[sp2]]))
      message(sp1, " & ", sp2, ": ", overlap, " shared")
    }
  }
}


# --- IncreaseOnly Euler diagram ---

# Extract genus colors
genera_only <- sapply(names(species_increase_only), function(x) word(x, 1))
venn_colors_only <- genus.colors[genera_only]

# Fit Euler diagram
fit_only <- euler(species_increase_only)

# Compute total metabolites and percentage labels
all_ids_only <- unique(unlist(species_increase_only))
total_metabolites_only <- length(all_ids_only)
quantity_labels_only <- paste0(round(fit_only$original.values / total_metabolites_only * 100, 1), "%")

# Display in viewer
# plot(fit_only,
#      quantities = quantity_labels_only,
#      fills = venn_colors_only,
#      alpha = 0.6,
#      lty = 1,
#      lwd = 2,
#      fontsize = 10,
#      labels = list(font = 2, cex = 1.2),
#      main = "Metabolites Released by Mucin Degraders (Increase Only)",
#      cex = 1.5)

# # Save PDF
# pdf(file.path(plot_dir, paste0(Sys.Date(), "_IncreaseOnly_Euler_percentage.pdf")),
#     width = 10, height = 8)
# plot(fit_only,
#      quantities = quantity_labels_only,
#      fills = venn_colors_only,
#      alpha = 0.6,
#      lty = 1,
#      lwd = 2,
#      fontsize = 10,
#      labels = list(font = 2, cex = 1.2),
#      main = "Metabolites Released by Mucin Degraders (Increase Only)",
#      cex = 1.5)
# dev.off()


#######################################
# PART 9: CREATE STACKED BARPLOTS
#######################################

# Calculate the counts per species
ratio_data <- data.frame()

# Get all species that appear in either list
all_species <- union(names(species_increase_only), names(species_increase_decrease))

for (sp in all_species) {
  # Get counts for each group
  n_increase_only <- if (sp %in% names(species_increase_only)) {
    length(species_increase_only[[sp]])
  } else {
    0
  }

  n_increase_decrease <- if (sp %in% names(species_increase_decrease)) {
    length(species_increase_decrease[[sp]])
  } else {
    0
  }

  # Total is the sum of both groups (they are mutually exclusive)
  n_total <- n_increase_only + n_increase_decrease

  genus <- word(sp, 1)

  ratio_data <- rbind(ratio_data, data.frame(
    Species = sp,
    Genus = genus,
    IncreaseDecrease = n_increase_decrease,
    IncreaseOnly = n_increase_only,
    Total = n_total,
    Proportion_IncreaseDecrease = if (n_total > 0) n_increase_decrease / n_total else 0,
    stringsAsFactors = FALSE
  ))
}

# Reshape for stacked bar chart
ratio_long <- ratio_data %>%
  pivot_longer(
    cols = c(IncreaseOnly, IncreaseDecrease),
    names_to = "Category",
    values_to = "Count"
  ) %>%
  mutate(Category = factor(Category,
    levels = c("IncreaseOnly", "IncreaseDecrease"),
    labels = c("Increase Only", "Increase + Decrease")
  ))

# Create stacked bar chart (absolute counts)
p_stacked <- ggplot(ratio_long, aes(
  x = reorder(Species, -Proportion_IncreaseDecrease),
  y = Count,
  fill = Category
)) +
  geom_bar(stat = "identity", color = "black", size = 0.3) +
  geom_text(
    data = ratio_data,
    aes(
      x = Species, y = Total,
      label = paste0(round(Proportion_IncreaseDecrease * 100, 1), "%")
    ),
    vjust = -0.5, size = 4, fontface = "bold", inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c(
    "Increase Only" = "lightgray",
    "Increase + Decrease" = "steelblue"
  )) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "Species",
    y = "Number of Metabolites",
    title = "Metabolites Released by Mucin Degraders",
    subtitle = "Proportion that Increase then Decrease (consumed by S. aureus)",
    fill = ""
  ) +
  ylim(0, max(ratio_data$Total) * 1.1)

print(p_stacked)

# Save absolute counts version
# ggsave(file.path(plot_dir, paste0(Sys.Date(), "_Stacked_Bar_IncreaseDecrease_Ratio.pdf")),
#        plot = p_stacked, width = 10, height = 7)

# Create 100% stacked bar chart (proportions)
p_proportion <- ggplot(ratio_long, aes(
  x = reorder(Species, -Proportion_IncreaseDecrease),
  y = Count,
  fill = Category
)) +
  geom_bar(stat = "identity", position = "fill", color = "black", size = 0.3) +
  geom_text(
    data = ratio_data,
    aes(
      x = Species, y = 0.5,
      label = paste0(
        IncreaseDecrease, " / ", Total, "\n",
        round(Proportion_IncreaseDecrease * 100, 1), "%"
      )
    ),
    size = 4, fontface = "bold", inherit.aes = FALSE
  ) +
  scale_fill_manual(values = c(
    "Increase Only" = "lightgray",
    "Increase + Decrease" = "steelblue"
  )) +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  ) +
  labs(
    x = "Species",
    y = "Proportion",
    title = "Proportion of each MDs' released metabolites that are consumed by S. aureus",
    subtitle = "(Consumed by S. aureus)",
    fill = ""
  )

print(p_proportion)

# Save proportion version
# ggsave(file.path(plot_dir, paste0(Sys.Date(), "_Proportion_Bar_IncreaseDecrease_Ratio.pdf")),
#        plot = p_proportion, width = 10, height = 7)

# Print summary table
message("\n=== Increase+Decrease Ratio Summary ===")
print(ratio_data %>%
  select(Species, Total, IncreaseDecrease, IncreaseOnly, Proportion_IncreaseDecrease) %>%
  arrange(desc(Proportion_IncreaseDecrease)))

#######################################
# ANALYSIS COMPLETE
#######################################

message("\n=== Analysis Complete ===")
message("All figures have been saved to: ", plot_dir)
message("\nFigure summary:")
message("1. Spaghetti plots: One per species per group")
message("2. Heatmap: Increase+Decrease metabolites only")
message("3. Euler diagrams: Separate for IncreaseOnly and Increase+Decrease")
message("4. Stacked bar charts: Absolute counts and proportions")


#######################################
# PART 10: PULL SPECTRA FOR SIRIUS
#######################################


cf_ids <- all_summ %>%
  filter(Species == "Cellulosimicrobium funkei") %>%
  pull(ID) %>%
  unique()

# Step 2: Filter rawdata using those IDs
rawdata_cf <- rawdata %>%
  filter(ID %in% cf_ids)


# Step 7: Keep only those with real MS2 (not "No MS2") and select top 50
cfsig_withms2 <- rawdata_cf %>%
  filter(MS2 != "No MS2")

# Now from this list, I'm going to pull MS2 spectra to search

# Run these two lines together for some reason
setwd("mucins_IDX/mucins/data/")
ms2mzxml <- readMSData(files = paste0("aquireX_data/", list.files("aquireX_data", pattern = "*mzXML")), mode = "onDisk", verbose = TRUE)

myms2 <- filterMsLevel(ms2mzxml, 2) # this grabs only MS2 spectra, not the MS1 spectra (which are used to see which peaks are available to grab MS2 of)


# 3. Create output directory
# dir.create("sirius_export", showWarnings = FALSE)

# 4. Function to export filtered spectra to MGF (optimized for your Orbitrap data)
export_spectra_to_mgf <- function(ms2_obj, feature_name, output_dir = "sirius_export") {
  if (length(ms2_obj) == 0) {
    warning(paste("No spectra found for", feature_name))
    return(NULL)
  }

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  # CRITICAL: With AcquireX you have multiple spectra per feature
  # Strategy: Select best quality spectrum

  # Calculate quality metrics for each spectrum
  quality_metrics <- lapply(seq_along(ms2_obj), function(i) {
    spec <- ms2_obj[[i]]
    mzs <- mz(spec)
    ints <- intensity(spec)

    list(
      idx = i,
      tic = sum(ints),
      num_peaks = length(mzs),
      base_peak_intensity = max(ints),
      high_quality_peaks = sum(ints > (max(ints) * 0.01)) # 1% threshold for Orbitrap
    )
  })

  # Select best spectrum (highest TIC, most informative)
  scores <- sapply(quality_metrics, function(m) {
    tic_score <- m$tic / max(sapply(quality_metrics, function(x) x$tic))
    peak_score <- m$num_peaks / max(sapply(quality_metrics, function(x) x$num_peaks))
    quality_score <- m$high_quality_peaks / max(sapply(quality_metrics, function(x) x$high_quality_peaks))

    0.4 * tic_score + 0.3 * peak_score + 0.3 * quality_score
  })

  best_idx <- which.max(scores)
  spec <- ms2_obj[[best_idx]]
  best_metrics <- quality_metrics[[best_idx]]

  prec_mz <- precursorMz(spec)
  prec_int <- precursorIntensity(spec)
  rt_sec <- rtime(spec)
  mzv <- mz(spec)
  inten <- intensity(spec)

  # ---- CRITICAL NOISE FILTERING FOR ORBITRAP ----
  # Orbitrap produces very clean data but can have low-intensity noise
  # Filter peaks below 0.5% of base peak (adjust if needed)
  keep <- inten > 0.005 * max(inten)
  mzv <- mzv[keep]
  inten <- inten[keep]

  # Additional filter: remove peaks very close to precursor (in-source fragments)
  # SIRIUS expects only MS2 fragments, not isotopes of precursor
  precursor_filter <- abs(mzv - prec_mz) > 0.5 # Remove peaks within 0.5 Da of precursor
  mzv <- mzv[precursor_filter]
  inten <- inten[precursor_filter]

  # Normalize intensities (optional but can help SIRIUS)
  inten <- (inten / max(inten)) * 1000

  # Create individual MGF file for this feature
  mgf_file <- file.path(output_dir, paste0(feature_name, ".mgf"))

  con <- file(mgf_file, "w")
  writeLines(c(
    "BEGIN IONS",
    paste0("TITLE=", feature_name),
    paste0("PEPMASS=", format(prec_mz, nsmall = 6)),
    paste0("RTINSECONDS=", format(rt_sec, nsmall = 2)),
    # CRITICAL: Negative mode for HILIC
    "CHARGE=1-",
    # IMPORTANT: Specify adduct (common in negative HILIC)
    # SIRIUS needs this for accurate formula prediction
    "ADDUCT=[M-H]-", # Most common; may also see [M+FA-H]-, [M+Cl]-
    # Optional: Add instrument type hint
    "MSLEVEL=2",
    "INSTRUMENT=Orbitrap",
    paste(
      format(mzv, nsmall = 6),
      format(inten, scientific = FALSE),
      sep = "\t"
    ),
    "END IONS"
  ), con)
  close(con)

  return(list(
    feature = feature_name,
    file = mgf_file,
    precursor_mz = prec_mz,
    precursor_intensity = prec_int,
    rt = rt_sec,
    num_peaks = length(mzv),
    num_peaks_before_filter = best_metrics$num_peaks,
    tic = best_metrics$tic,
    base_peak_int = best_metrics$base_peak_intensity,
    high_quality_peaks = best_metrics$high_quality_peaks,
    quality_score = scores[best_idx],
    num_spectra_found = length(ms2_obj)
  ))
}

# 5. Alternative: Export all to single MGF file (also works well)
# export_all_to_single_mgf <- function(spectra_list, output_file = "sirius_export/all_features.mgf") {
#
#   con <- file(output_file, "w")
#
#   for (i in seq_along(spectra_list)) {
#     feature_data <- spectra_list[[i]]
#     if (is.null(feature_data)) next
#
#     spec <- feature_data$spectrum
#
#     writeLines(c(
#       "BEGIN IONS",
#       paste0("TITLE=", feature_data$feature_name),
#       paste0("PEPMASS=", format(precursorMz(spec), nsmall = 4)),
#       paste0("RTINSECONDS=", format(rtime(spec), nsmall = 2)),
#       "CHARGE=1-",
#       paste(
#         format(mz(spec), nsmall = 4),
#         format(intensity(spec), scientific = FALSE),
#         sep = "\t"
#       ),
#       "END IONS",
#       ""
#     ), con)
#   }
#
#   close(con)
#   cat("Exported", length(spectra_list), "spectra to", output_file, "\n")
# }

# 6. Loop through your features and export each
results_list <- list()

for (i in seq_len(nrow(cfsig_withms2))) {
  mz_val <- cfsig_withms2$`m/z`[i]
  rt_center_min <- cfsig_withms2$`RT [min]`[i]

  # Create a meaningful feature name
  feature_name <- paste0(
    "Feature_",
    sprintf("%04d", i),
    "_mz", format(mz_val, nsmall = 4),
    "_rt", format(rt_center_min, nsmall = 2)
  )

  # Convert center RT to seconds
  rt_center_sec <- rt_center_min * 60
  rt_window_sec <- 5 # ±5 seconds
  rt_range_sec <- c(rt_center_sec - rt_window_sec, rt_center_sec + rt_window_sec)

  # Filter MS2 spectra
  possibleMS2 <- myms2 %>%
    filterPrecursorMz(mz_val, ppm = 5) %>%
    filterRt(rt_range_sec)

  # Export this feature
  result <- export_spectra_to_mgf(possibleMS2, feature_name)

  if (!is.null(result)) {
    results_list[[i]] <- result
    cat(sprintf(
      "Exported %s: %d spectra found, %d peaks in best spectrum\n",
      feature_name, result$num_spectra_found, result$num_peaks
    ))
  }
}

# 7. Create summary table with quality metrics
results_df <- do.call(rbind, lapply(results_list, function(x) {
  if (is.null(x)) {
    return(NULL)
  }
  data.frame(
    feature = x$feature,
    precursor_mz = x$precursor_mz,
    precursor_intensity = x$precursor_intensity,
    rt_seconds = x$rt,
    num_peaks = x$num_peaks,
    tic = x$tic,
    base_peak_int = x$base_peak_int,
    high_quality_peaks = x$high_quality_peaks,
    quality_score = x$quality_score,
    num_spectra = x$num_spectra_found,
    mgf_file = x$file,
    stringsAsFactors = FALSE
  )
}))

# Save summary
# write.csv(results_df, "sirius_export/251205_export_summary.csv", row.names = FALSE)

#######################################
# PART 11: PULL SPECTRA FOR MZCLOUD
#######################################

# Process
# Take all_summ, and merge in Name column from cfsig_withms2, and the MS2 column
filtered_df_cf <- all_summ %>%
  filter(
    Species == "Cellulosimicrobium funkei",
    Group == "Increase+Decrease"
  )

# 2. Merge in columns from cfsig_withms2 by shared ID column
merged_df_cf <- filtered_df_cf %>%
  left_join(cfsig_withms2, by = "ID")


# Then subset to the top 25 or 50 named unique ID numbers with MS2 spectra, make a table with these predictions as the figure/supplement, or can make a pie chart based on how the numbers look

# 3. Rank from most negative to least negative FoldChange
ranked_df_cf <- merged_df_cf %>%
  arrange(FoldChange)

# 4. Keep only rows with non-NA Name and take top 50 UNIQUE IDs
top50_df_cf <- ranked_df_cf %>%
  filter(!is.na(Name)) %>% # must have Name
  distinct(ID, .keep_all = TRUE) %>% # ensure unique IDs
  distinct(Name, .keep_all = TRUE) %>%
  slice_head(n = 50)

# write.csv(top50_df_cf, "/Users/gdlab/Library/CloudStorage/Box-Box/2024_Dome2.0/Metabolomics/HigResRun/251030_JRedits_metabolomics/251205MZCloudValidation/251205_Top50MZCLOUDvalidationexport_summary.csv", row.names = FALSE)


# Then use Jrtools to export spectra to folders named by the predicted compound
setwd("mucins_IDX/mucins/data/")
ms2mzxml <- readMSData(files = paste0("aquireX_data/", list.files("aquireX_data", pattern = "*mzXML")), mode = "onDisk", verbose = TRUE)

myms2 <- filterMsLevel(ms2mzxml, 2) # this grabs only MS2 spectra, not the MS1 spectra (which are used to see which peaks are available to grab MS2 of)


for (i in seq_len(nrow(top50_df_cf))) {
  mz_val <- top50_df_cf$`m/z`[i]
  rt_center_min <- top50_df_cf$`RT [min]`[i]

  # Convert center RT to seconds
  rt_center_sec <- rt_center_min * 60
  rt_window_sec <- 5 # ±5 seconds
  rt_range_sec <- c(rt_center_sec - rt_window_sec, rt_center_sec + rt_window_sec)

  # Filter MS2 spectra
  possibleMS2 <- myms2 %>%
    filterPrecursorMz(mz_val, ppm = 5) %>%
    filterRt(rt_range_sec)

  if (length(possibleMS2) > 0) {
    for (j in seq_along(possibleMS2)) {
      plot(possibleMS2[[j]])

      # Define output path
      dir_path <- paste0("peak_", top50_df_cf$ID[i], "_", top50_df_cf$Name[i], "_", round(mz_val, 2), "_", round(rt_center_min, 2), "/spectra_", j)
      file_path <- paste0(dir_path, "/_spectra_", j, ".csv")

      # Create directory if it doesn't exist
      if (!dir.exists(dir_path)) {
        dir.create(dir_path, recursive = TRUE, showWarnings = FALSE)
      }

      # Export spectrum
      jrtools::export_idx_spectrum(
        possibleMS2[[j]],
        file_path
      )
    }
  } else {
    message("No MS2 spectra found for m/z = ", mz_val, " at RT ~", rt_center_min, " min")
  }
}


# Then manually use MZcloud, and note if it matches the prediction and if cosine is above 70, screenshot the match


#######################################
# PART 12: FLASH ENTROPY
#######################################


# ============================================================================
# STEP 1: Extract and prepare your MS2 spectra (ALL 1200+ features that are significant in C.funkei)
# ============================================================================

cat("Extracting MS2 spectra from all features...\n")

# Collect all MS2 spectra
all_ms2_spectra <- list()
spectrum_metadata <- list()

# Get total for progress tracking
n_features <- nrow(cfsig_withms2) # Rename your dataframe if needed!

for (i in seq_len(n_features)) {
  if (i %% 100 == 0) cat("  Processing feature", i, "/", n_features, "\n")

  mz_val <- cfsig_withms2$`m/z`[i]
  rt_center_min <- cfsig_withms2$`RT [min]`[i]

  # Convert center RT to seconds
  rt_center_sec <- rt_center_min * 60
  rt_window_sec <- 5 # ±5 seconds
  rt_range_sec <- c(rt_center_sec - rt_window_sec, rt_center_sec + rt_window_sec)

  # Filter MS2 spectra
  possibleMS2 <- myms2 %>%
    filterPrecursorMz(mz_val, ppm = 5) %>%
    filterRt(rt_range_sec)

  # Store if we found spectra
  if (length(possibleMS2) > 0) {
    # Extract all spectra at once as a list
    spec_list <- spectra(possibleMS2)

    for (j in seq_along(spec_list)) {
      spec_id <- paste0("feat_", i, "_spec_", j)

      # Extract m/z and intensity from the Spectrum1 object
      # These are stored in slots, not accessed by functions
      mz_vals <- spec_list[[j]]@mz
      int_vals <- spec_list[[j]]@intensity

      # Create peaks matrix
      peaks <- cbind(mz = mz_vals, intensity = int_vals)

      # Store spectrum (need mz and intensity columns)
      if (nrow(peaks) > 0) {
        all_ms2_spectra[[spec_id]] <- peaks

        # Store metadata
        spectrum_metadata[[spec_id]] <- tibble(
          spectrum_id = spec_id,
          feature_idx = i,
          precursor_mz = mz_val,
          precursor_rt = rt_center_min,
          n_peaks = nrow(peaks)
        )
      }
    }
  }
}

# Combine metadata
spectrum_meta_df <- bind_rows(spectrum_metadata)

cat(
  "Collected", length(all_ms2_spectra), "MS2 spectra from",
  length(unique(spectrum_meta_df$feature_idx)), "features\n"
)


# Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.

# John modified code

## STEP 2: Clean and filter spectra

# Spectra with only 1-2 peaks are pretty commonplace, although for similarity calculations I might require at least 2.


# Remove spectra with too few peaks
min_peaks <- 2
good_spectra <- all_ms2_spectra[spectrum_meta_df$n_peaks >= min_peaks]
good_meta <- spectrum_meta_df %>% filter(n_peaks >= min_peaks)

# normalize spectral intensities so sum of each spectrum is 0.5
good_spectra <- lapply(X = good_spectra, FUN = function(sp) {
  sp[, "intensity"] <- sp[, "intensity"] / (sum(sp[, "intensity"]) * 2)
  return(sp)
})

cat("After filtering: ", length(good_spectra), "spectra with >=", min_peaks, "peaks\n")


## STEP 3: Create flash entropy search library


cat("\nBuilding spectral library...\n")

# Combine all spectra into a single matrix for fast searching
n_frags <- sapply(good_spectra, nrow)

# Create library matrix with spectrum ID column
clean_spectra_frags <- do.call(rbind, good_spectra)
idx <- rep(seq_along(good_spectra), times = n_frags)
clean_spectra_frags <- cbind(Spectrum = idx, clean_spectra_frags)

# Sort by fragment mass for efficient searching
clean_spectra_frags <- clean_spectra_frags[order(clean_spectra_frags[, "mz"]), ]

cat("Created spectral library with", nrow(clean_spectra_frags), "total fragments\n")


## STEP 4: PARALLEL all-vs-all entropy similarity calculation


# Note this introduces some RNG, so we set a seed
set.seed(42)

n_spec <- length(good_spectra)
cat("\n========================================\n")
cat("PARALLELIZED SIMILARITY CALCULATION\n")
cat("========================================\n")
cat("Total spectra:", n_spec, "\n")
cat("Total comparisons:", format(n_spec^2, big.mark = ","), "\n\n")

# Set up parallel backend
n_cores <- detectCores() - 2 # Leave one core free
cat("Using", n_cores, "CPU cores\n\n")

cl <- makeCluster(n_cores)
clusterSetRNGStream(cl, 123)
registerDoParallel(cl)

# Export necessary objects to workers
# clusterExport(cl, c("good_spectra", "clean_spectra_frags"), envir = environment())

# Load required packages on workers
# clusterEvalQ(cl, {
#   library(jrtools)  # or library(msentropy) if using that instead
# })

start_time <- Sys.time()

# PARALLEL LOOP - each worker processes different spectra
cat("Starting parallel calculation...\n")

similarity_results <- foreach(
  i = 1:n_spec,
  .combine = "cbind",
  .packages = c("jrtools"),
  .verbose = FALSE
) %dopar% {
  # Get query spectrum
  sp <- good_spectra[[i]]

  # Do similarity search
  sp_sim <- jrtools::flash_entropy_search(
    fragment_library = clean_spectra_frags,
    query_spectrum = sp,
    ms2_tol_ppm = 10
  )

  # Create result vector
  sim_vec <- numeric(n_spec)
  sim_vec[sp_sim$Spectrum] <- sp_sim$Similarity

  sim_vec
}

# Stop cluster
stopCluster(cl)

# Convert to matrix
full_es_mat <- as.matrix(similarity_results)

# Make symmetric and remove diagonal
full_es_mat[lower.tri(full_es_mat, diag = TRUE)] <- 0

elapsed_total <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
cat("\n✓ Similarity calculation complete in", round(elapsed_total, 1), "minutes\n")
cat("  (", round(n_spec^2 / elapsed_total / 60, 0), "comparisons/second)\n\n")


## STEP 5: Analyze similarity distribution


# Plot similarity distribution
par(mfrow = c(1, 2))
hist(full_es_mat[full_es_mat > 0],
  breaks = 50,
  main = "Entropy Similarity Distribution\n(All matches > 0)",
  xlab = "Entropy Similarity", col = "lightblue"
)

hist(full_es_mat[full_es_mat > 0.5],
  breaks = 50,
  main = "High Similarity Matches\n(> 0.5)",
  xlab = "Entropy Similarity", col = "coral"
)
par(mfrow = c(1, 1))

# Summary statistics
cat("\nSimilarity Statistics:\n")
cat("  Non-zero matches:", sum(full_es_mat > 0), "\n")
cat("  Mean:", round(mean(full_es_mat[full_es_mat > 0]), 3), "\n")
cat("  Median:", round(median(full_es_mat[full_es_mat > 0]), 3), "\n")
cat("  75th percentile:", round(quantile(full_es_mat[full_es_mat > 0], 0.75), 3), "\n")
cat("  95th percentile:", round(quantile(full_es_mat[full_es_mat > 0], 0.95), 3), "\n")
cat("  Max:", round(max(full_es_mat), 3), "\n\n")


# ============================================================================
# STEP 6: Filter by similarity threshold and build network
# ============================================================================

# Set threshold (somewhat arbitrary, but seems conservative to go a bit beyond the 95th percentile up to 0.6)
similarity_threshold <- 0.6

cat("Filtering with threshold =", similarity_threshold, "\n")

# Get match indices
es_match_idx <- which(full_es_mat > similarity_threshold, arr.ind = TRUE)
es_match_idx <- es_match_idx[es_match_idx[, "row"] != es_match_idx[, "col"], ]
es_match_sim <- full_es_mat[es_match_idx]

cat("Found", nrow(es_match_idx), "spectral matches above threshold\n")
cat(
  "This represents",
  length(unique(c(es_match_idx[, 1], es_match_idx[, 2]))),
  "unique spectra in the network\n\n"
)

# Create vertices (nodes) with metadata - SMART LABELING
# First, identify top 500 ranked features
top_ranked_ids <- ranked_df_cf %>%
  slice_head(n = 500) %>%
  pull(ID) # Adjust column name if your ID column is named differently

# Add ranking information and feature names
es_match_vertices <- tibble(
  spectrum_idx = seq_along(good_spectra),
  spectrum_id = names(good_spectra)
) %>%
  left_join(good_meta, by = "spectrum_id") %>%
  # Add feature data including ID for ranking lookup
  left_join(
    cfsig_withms2 %>%
      mutate(feature_idx = row_number()) %>%
      select(feature_idx, ID, feature_name = Name, feature_mz = `m/z`, feature_rt = `RT [min]`),
    by = "feature_idx"
  ) %>%
  # Add ranking information
  mutate(
    in_top_500 = ID %in% top_ranked_ids,
    rank = match(ID, top_ranked_ids), # Gets rank position (NA if not in top 500)
    # Short label for node
    label = paste0(
      "F", feature_idx, "_S",
      gsub(".*spec_(\\d+)", "\\1", spectrum_id)
    ),
    # Full label with name for detailed view
    full_label = paste0(
      ifelse(!is.na(feature_name) & feature_name != "", feature_name, paste0("Feature ", feature_idx)),
      "\n",
      round(precursor_mz, 2), " m/z @ ",
      round(precursor_rt, 2), " min",
      ifelse(in_top_500, paste0("\nRank: ", rank), "")
    ),
    # SMART DISPLAY LABEL: Only show name if it's in top 500 AND has a name
    display_label = case_when(
      in_top_500 & !is.na(feature_name) & feature_name != "" ~ feature_name,
      in_top_500 ~ paste0("F", feature_idx),
      TRUE ~ NA_character_ # Don't label if not in top 500
    )
  )

cat("\nLabeling statistics:\n")
cat("  Features in top 500 ranked:", sum(unique(es_match_vertices$feature_idx) %in%
  match(top_ranked_ids, cfsig_withms2$ID), na.rm = TRUE), "\n")
cat("  Of those, have names:", sum(!is.na(es_match_vertices$display_label)), "\n")

# Only keep vertices that are in the network
vertices_in_network <- unique(c(es_match_idx[, "row"], es_match_idx[, "col"]))
es_match_vertices <- es_match_vertices %>%
  filter(spectrum_idx %in% vertices_in_network)

# Create edges with similarity scores
es_match_edges <- tibble(
  from = es_match_idx[, "row"],
  to = es_match_idx[, "col"],
  entropy_similarity = es_match_sim
)

# Build igraph object
g_spectral_network <- graph_from_data_frame(
  d = es_match_edges,
  directed = FALSE,
  vertices = es_match_vertices
)


# ============================================================================
# STEP 7: Detect clusters
# ============================================================================

cat("Analyzing network structure...\n")

# Find connected components (separate clusters)
clusters <- components(g_spectral_network)

# Add cluster membership to vertices
V(g_spectral_network)$cluster <- clusters$membership

cat("Network contains", clusters$no, "connected components\n")

# Show distribution of cluster sizes
cluster_size_table <- table(clusters$membership)
cat("\nCluster size distribution:\n")
cat("  1 spectrum:", sum(cluster_size_table == 1), "clusters\n")
cat("  2-5 spectra:", sum(cluster_size_table >= 2 & cluster_size_table <= 5), "clusters\n")
cat("  6-10 spectra:", sum(cluster_size_table >= 6 & cluster_size_table <= 10), "clusters\n")
cat("  11-20 spectra:", sum(cluster_size_table >= 11 & cluster_size_table <= 20), "clusters\n")
cat("  >20 spectra:", sum(cluster_size_table > 20), "clusters\n")

# Show largest clusters
largest_clusters <- sort(cluster_size_table, decreasing = TRUE)[1:min(10, length(cluster_size_table))]
cat("\nTop 10 largest clusters:\n")
print(largest_clusters)

# Detect communities within the network
cat("\nDetecting communities within network...\n")
communities_louvain <- cluster_louvain(g_spectral_network)
communities_walktrap <- cluster_walktrap(g_spectral_network)

cat("Louvain algorithm:", length(communities_louvain), "communities\n")
cat("Walktrap algorithm:", length(communities_walktrap), "communities\n\n")


# ============================================================================
# STEP 8: Visualize clusters across (GGRAPH VERSION)
# ============================================================================


# Note that this may produce a clustering figure that will have nodes on different positions of the PDF. However, I believe that the actual nodes and edges look like they are the same, especially using the seed setting up above, that should deal with the RNG generated by parallelization. But in the end this methodology really is a screen to look for features that look like other features that are consumed (as defined by our above stats methods).


cat("\n========================================\n")
cat("GENERATING COMPREHENSIVE VISUALIZATIONS\n")
cat("========================================\n")

# Load required libraries
library(ggraph)
library(tidygraph)
library(ggplot2)
library(ggrepel)
library(igraph)
library(dplyr)

# Set graphics options to use Cairo (better font handling)
options(bitmapType = "cairo")

# Check if Cairo is available
if (!capabilities("cairo")) {
  warning("Cairo graphics not available. Font rendering may have issues.")
}

# Analyze cluster size distribution
cluster_sizes <- table(clusters$membership)
sorted_clusters <- sort(cluster_sizes, decreasing = TRUE)

cat("Total clusters:", clusters$no, "\n")
cat("Size distribution:\n")
cat("  Singletons (1 node):", sum(cluster_sizes == 1), "\n")
cat("  Small (2-10 nodes):", sum(cluster_sizes >= 2 & cluster_sizes <= 10), "\n")
cat("  Medium (11-50 nodes):", sum(cluster_sizes >= 11 & cluster_sizes <= 50), "\n")
cat("  Large (51-200 nodes):", sum(cluster_sizes >= 51 & cluster_sizes <= 200), "\n")
cat("  Very large (>200 nodes):", sum(cluster_sizes > 200), "\n\n")

# Helper function to add unique labels to any graph
add_unique_labels <- function(g) {
  vertices_df <- as_data_frame(g, what = "vertices")

  # Identify which nodes should get labels (one per unique name)
  label_map <- vertices_df %>%
    filter(!is.na(display_label)) %>%
    group_by(display_label) %>%
    slice_max(order_by = n_peaks, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    pull(name)

  # Create unique label vector
  V(g)$unique_label <- ifelse(V(g)$name %in% label_map,
    V(g)$display_label,
    NA_character_
  )
  return(g)
}

# Helper function to create ggraph plot
create_cluster_ggraph <- function(g, title, subtitle = NULL) {
  # Add unique labels
  g <- add_unique_labels(g)

  # Compute layout
  layout_matrix <- layout_with_fr(g, niter = 500)

  # Convert to tidygraph
  g_tidy <- as_tbl_graph(g)

  # Create plot with red nodes, transparency, and white borders
  p <- ggraph(g_tidy, layout = layout_matrix) +
    geom_edge_link(aes(width = entropy_similarity, alpha = entropy_similarity),
      color = "gray50"
    ) +
    scale_edge_width(range = c(0.3, 2)) +
    scale_edge_alpha(range = c(0.1, 0.5)) +
    geom_node_point(aes(size = 4),
      color = "red",
      fill = "red",
      alpha = 0.4,
      stroke = 0.8,
      shape = 21,
      colour = "white"
    ) +
    scale_size_continuous(range = c(2, 8)) +
    geom_node_text(
      aes(
        label = unique_label,
        filter = !is.na(unique_label)
      ),
      repel = TRUE,
      size = 3.5,
      max.overlaps = Inf,
      box.padding = 0.5,
      point.padding = 0.5,
      segment.color = "gray70",
      segment.size = 0.3,
      force = 2
    ) +
    theme_graph() +
    labs(
      title = title,
      subtitle = subtitle
    ) +
    theme(
      legend.position = "right",
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 10)
    )

  return(p)
}

# Prepare graph aesthetics
V(g_spectral_network)$size <- ifelse(
  V(g_spectral_network)$in_top_500,
  scales::rescale(V(g_spectral_network)$n_peaks, to = c(5, 12)),
  scales::rescale(V(g_spectral_network)$n_peaks, to = c(3, 8))
)

# ===========================================================================
# PLOT: Full network overview (up to 1000 nodes)
# ===========================================================================

cat("Creating full network overview...\n")

if (vcount(g_spectral_network) <= 1000) {
  g_full <- g_spectral_network
} else {
  cat("Network too large (", vcount(g_spectral_network), " nodes) - showing top components\n")
  top_components <- order(clusters$csize, decreasing = TRUE)[1:min(5, clusters$no)]
  g_full <- induced_subgraph(
    g_spectral_network,
    V(g_spectral_network)[clusters$membership %in% top_components]
  )
}

p_full <- create_cluster_ggraph(
  g_full,
  title = "Full Spectral Similarity Network",
  subtitle = paste0(
    vcount(g_full), " spectra in ",
    length(unique(V(g_full)$cluster)), " clusters"
  )
)

ggsave(
  filename = file.path(
    "FlashClust",
    paste0(Sys.Date(), "_01_FullNetwork_Overview.pdf")
  ),
  plot = p_full,
  width = 18,
  height = 14,
  units = "in",
  device = cairo_pdf
)

cat("✓ Saved full network overview\n")


nodes <- as_data_frame(g_spectral_network, what = "vertices")
write.csv(nodes, "251209_nodes.csv", row.names = FALSE)

edges <- as_data_frame(g_spectral_network, what = "edges")
write.csv(edges, "251209_edges.csv", row.names = FALSE)

write.csv(
  data.frame(
    id = names(clusters$membership),
    cluster = clusters$membership
  ),
  "251209_clusters.csv",
  row.names = FALSE
)
