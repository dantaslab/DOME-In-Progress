library(readxl)
library(dplyr)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(NatParksPalettes)
library(RColorBrewer)
library(phytools)
library(ggtree)
library(pheatmap)
library(tidytree)
library(ggnewscale)
library(lubridate)
library(kableExtra)
library(ggpubr)
library(qgraph)
library(phangorn)
library(purrr)
library(colorspace)

culture <- read_excel("250603_CulturingExtractionMasterlist_V9_SP.xlsx")
swab <- read.csv("240304_NasalMasterList_V1_SP.csv")
culture$SampleNameLower <- tolower(culture$`Sample Name`)
swab$SampleIDLower <- tolower(swab$Sample.ID)

# Remove all asterisks from swab$SampleIDLower
swab$SampleIDLower <- gsub("\\*", "", swab$SampleIDLower)

# Merge using the cleaned, lowercase columns
cultureswab <- merge(
  culture,
  swab,
  by.x = "SampleNameLower",
  by.y = "SampleIDLower"
)

# (Optional) Remove the helper columns
culture$SampleNameLower <- NULL
swab$SampleIDLower <- NULL


cultureswab <- data.frame(cultureswab, row.names = cultureswab$DomeNumber)


cultureswab <- cultureswab %>%
  mutate(
    color_group = case_when(
      grepl("Cow", Subject, ignore.case = TRUE) ~ "Cow",
      grepl("Worker", Subject, ignore.case = TRUE) ~ "Farmer",
      grepl("Control", Subject, ignore.case = TRUE) ~ "Non-Farmer",
      TRUE ~ "Other"
    ),
    color_code = case_when(
      color_group == "Cow" ~ "#873e23",
      color_group == "Farmer" ~ "#18678d",
      color_group == "Non-Farmer" ~ "#626262",
      TRUE ~ "#bbbbbb" # fallback color
    )
  )

subject_df <- cultureswab[, "color_group", drop = F]
subject_color_df <- cultureswab[, "color_code", drop = F]

# Confirm this is the final taxonomy table
taxonomycsv <- read.csv("gtdbtk.bac120.summary.csv", header = F)

parse_gtdb_taxonomy <- function(tax_string) {
  tax_split <- str_split(tax_string, ";")[[1]]
  tax_list <- setNames(rep(NA_character_, 7), c("domain", "phylum", "class", "order", "family", "genus", "species"))
  for (entry in tax_split) {
    level <- substr(entry, 1, 1)
    val <- sub("^[a-z]__+", "", entry)
    if (level == "d") tax_list["domain"] <- val
    if (level == "p") tax_list["phylum"] <- val
    if (level == "c") tax_list["class"] <- val
    if (level == "o") tax_list["order"] <- val
    if (level == "f") tax_list["family"] <- val
    if (level == "g") tax_list["genus"] <- val
    if (level == "s") tax_list["species"] <- val
  }
  return(tax_list)
}

names(taxonomycsv)[names(taxonomycsv) == "V2"] <- "taxonomy"

taxonomycsv <- data.frame(taxonomycsv, row.names = taxonomycsv$V1)

parsed_mat <- do.call(rbind, lapply(taxonomycsv$taxonomy, parse_gtdb_taxonomy))
parsedtax <- cbind(taxonomycsv, parsed_mat)

family_df <- parsedtax[, "family", drop = F]
genus_df <- parsedtax[, "genus", drop = F]
species_df <- parsedtax[, "species", drop = F]


species_list <- unique(species_df$species)

species_df_processing <- data.frame(
  species = species_list,
  genus = sub(" .*", "", species_list), # extract genus before underscore
  stringsAsFactors = FALSE
)

# Step 2: Assign each genus a base hue
genera <- sort(unique(species_df_processing$genus)) # alphabetize
genus_hues <- seq(15, 375, length.out = length(genera) + 1)[-1]
genus_base <- setNames(genus_hues, genera)

priority_genera <- c(
  "Cellulosimicrobium", "Brevibacterium",
  "Micrococcus", "Bacillus", "Arthrobacter_b", "Shouchella", "Kocuria",
  "Nocardiopsis", "Staphylococcus"
)

# all genera (alphabetical)
genera <- sort(unique(species_df_processing$genus))

# evenly spaced hues across 360 degrees
genus_hues_all <- seq(15, 375, length.out = length(genera) + 1)[-1]

# assign priority genera first, with maximally distinct hues
n_priority <- length(priority_genera)
priority_hues <- seq(15, 375, length.out = n_priority + 1)[-1]

# start mapping with priority genera
genus_base <- setNames(rep(NA, length(genera)), genera)
genus_base[priority_genera] <- priority_hues

# assign remaining genera to leftover hues
remaining_genera <- setdiff(genera, priority_genera)
remaining_hues <- setdiff(genus_hues_all, priority_hues)

genus_base[remaining_genera] <- remaining_hues


# Step 3: Assign species colors by tweaking lightness within each genus
species_colors <- species_df_processing %>%
  group_by(genus) %>%
  mutate(
    n_species = n(),
    species_idx = row_number(),
    color = hcl(
      h = genus_base[genus],
      c = 60,
      l = seq(40, 80, length.out = first(n_species)) # ensure scalar
    )[species_idx] # pick the right color for each row
  ) %>%
  ungroup()

species_colors <- species_colors %>%
  filter(!is.na(color))

# Step 4: Compute genus midpoint color as average of its species
genus_colors <- species_colors %>%
  group_by(genus) %>%
  summarise(
    color = {
      rgb_mat <- hex2RGB(color)@coords # convert hex → RGB matrix (rows = colors, cols = R,G,B)
      avg_rgb <- colMeans(rgb_mat) # average each channel separately
      hex(RGB(avg_rgb[1], avg_rgb[2], avg_rgb[3])) # convert back to hex
    }
  )


cat("species.colors <- c(",
  paste0('"', species_colors$species, '" = "', species_colors$color, '"', collapse = ", \n "),
  ")",
  sep = ""
)

cat("genus.colors <- c(",
  paste0('"', genus_colors$genus, '" = "', genus_colors$color, '"', collapse = ", \n"),
  ")",
  sep = ""
)


ggplot(species_colors, aes(x = genus, y = species, fill = color)) +
  geom_tile(color = "white") +
  scale_fill_identity() +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.title = element_blank()
  ) +
  ggtitle("Species-level colors")

ggplot(genus_colors, aes(x = genus, y = 1, fill = color)) +
  geom_tile(color = "white", height = 0.9) +
  scale_fill_identity() +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_blank(),
    axis.title = element_blank(),
    axis.ticks.y = element_blank()
  ) +
  ggtitle("Genus-level colors")


# Hard changing colors on the genera, since we have a few problem genera that don't pop out nicely

genus_colors <- read.csv("250917colorsFINALIZED.csv")


species_df_colors_process <- species_colors %>%
  transmute(
    Name = species,
    Hex = color,
    Category = "species"
  )

# genus dataframe
genus_df_colors_process <- genus_colors %>%
  transmute(
    Name = genus,
    Hex = color,
    Category = "genus"
  )

# combine
combined_colors <- bind_rows(species_df_colors_process, genus_df_colors_process)

combined_colors <- combined_colors %>%
  bind_rows(tibble(Name = "Cow", Hex = "#B03A1C", Category = "subject")) %>% # vivid orange
  bind_rows(tibble(Name = "Farmer", Hex = "#0F7FBF", Category = "subject")) %>% # strong blue
  bind_rows(tibble(Name = "Non-Farmer", Hex = "#3F3F3F", Category = "subject")) # true black


colors <- combined_colors


# Plotting the trees!!
genus_colors_tree <- colors$Hex[colors$Category == "genus"]
names(genus_colors_tree) <- colors$Name[colors$Category == "genus"]

# mucin tree first

tree_mucin <- read.tree("250617mucin-isolates.tree") # tips: 436
tree_mucin <- midpoint(tree_mucin)

treeplot <- ggtree(tree_mucin) %<+% genus_df
circle_treeplot_mucin <- ggtree(tree_mucin, layout = "fan") %<+% genus_df


circle_genus_tree_mucin <- gheatmap(circle_treeplot_mucin, genus_df,
  colnames_angle = 85,
  colnames_offset_y = -0.1, # keep labels close
  offset = 0.01, # very small gap
  width = 0.1, # thinner ring
  font.size = 5
) +
  scale_fill_manual(
    values = genus_colors_tree,
    breaks = names(genus_colors_tree)
  ) +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  ) + # remove extra white space
  labs(fill = "Bacterial host")


circle_genus_tree_mucin2 <- circle_genus_tree_mucin + new_scale_fill()


# adding abx group ring
subject_colors_tree <- colors$Hex[colors$Category == "subject"]
names(subject_colors_tree) <- colors$Name[colors$Category == "subject"]


circle_genus_subject_tree_mucin <- gheatmap(circle_genus_tree_mucin2, subject_df,
  colnames_angle = 85,
  colnames_offset_y = -0.2, # keep labels close
  offset = 0.1, # small spacing from previous ring
  width = 0.08, # slimmer ring
  font.size = 5
) +
  scale_fill_manual(
    values = subject_colors_tree,
    breaks = names(subject_colors_tree)
  ) +
  theme(
    legend.position = "none", # hide legend
    plot.margin = margin(0, 0, 0, 0)
  ) + # minimize extra whitespace
  labs(fill = "Infant antibiotic exposure")

circle_genus_subject_tree_mucin2 <- circle_genus_subject_tree_mucin + new_scale_fill()

# plot_file <- file.path('./figures/', paste0(Sys.Date(), "abx_genus_phage.pdf"))
# ggsave(filename = plot_file, plot = circle_abx_mash_nicu_tree2, width = 18, height = 18)


# Staph tree next

species_colors_tree <- colors$Hex[colors$Category == "species"]
names(species_colors_tree) <- colors$Name[colors$Category == "species"]

# As of
tree_staph <- read.tree("250617msa-isolates.tree") # tips: 436
tree_staph <- midpoint(tree_staph)

treeplot <- ggtree(tree_staph) %<+% species_df
circle_treeplot_staph <- ggtree(tree_staph, layout = "fan") %<+% species_df


circle_species_tree_staph <- gheatmap(circle_treeplot_staph, species_df,
  colnames_angle = 85,
  colnames_offset_y = -0.1, # keep labels close
  offset = 0.01, # very small gap
  width = 0.1, # thinner ring
  font.size = 5
) +
  # geom_treescale(x = .3) +
  scale_fill_manual(
    values = species_colors_tree,
    breaks = names(species_colors_tree)
  ) +
  theme(
    legend.position = "none",
    plot.margin = margin(0, 0, 0, 0)
  ) + # remove extra white space
  labs(fill = "Bacterial host")


circle_species_tree_staph2 <- circle_species_tree_staph + new_scale_fill()


# adding abx group ring
subject_colors_tree <- colors$Hex[colors$Category == "subject"]
names(subject_colors_tree) <- colors$Name[colors$Category == "subject"]


circle_species_subject_tree_staph <- gheatmap(circle_species_tree_staph2, subject_df,
  colnames_angle = 85,
  colnames_offset_y = -0.2, # keep labels close
  offset = 0.1, # small spacing from previous ring
  width = 0.08, # slimmer ring
  font.size = 5
) +
  scale_fill_manual(
    values = subject_colors_tree,
    breaks = names(subject_colors_tree)
  ) +
  theme(
    legend.position = "none", # hide legend
    plot.margin = margin(0, 0, 0, 0)
  ) + # minimize extra whitespace
  labs(fill = "Infant antibiotic exposure")

circle_species_subject_tree_staph2 <- circle_species_subject_tree_staph + new_scale_fill()


# Alternative approach, put all the isolates on the same tree, and simply add a new ring to show culturing source

media_df <- cultureswab[, "PlateType", drop = F]


colors <- colors %>%
  bind_rows(
    tibble(
      Name = "M",
      Hex = "#000000",
      Category = "media"
    )
  )

colors <- colors %>%
  bind_rows(
    tibble(
      Name = "S",
      Hex = "#FFFFFF",
      Category = "media"
    )
  )


media_colors_tree <- colors$Hex[colors$Category == "media"]
names(media_colors_tree) <- colors$Name[colors$Category == "media"]

tree_all <- read.tree("250917_allisolates_rootedtoacinetobacter_gtdbtk.bac120.decorated.tree")

treeplot <- ggtree(tree_all) %<+% species_df

tree_all <- midpoint(tree_all)

circle_treeplot_all <- ggtree(tree_all, layout = "fan") %<+% species_df


circle_media_tree_all <- gheatmap(circle_treeplot_all, media_df,
  colnames_angle = 85,
  colnames_offset_y = -0.1, # keep labels close
  offset = 0.0, # very small gap
  width = 0.1, # thinner ring
  font.size = 5, color = "black"
) +
  geom_treescale() +
  scale_fill_manual(
    values = media_colors_tree,
    breaks = names(media_colors_tree)
  ) +
  theme(
    legend.position = "right",
    plot.margin = margin(0, 0, 0, 0)
  ) + # remove extra white space
  labs(fill = "Culturing media")


circle_media_tree_all2 <- circle_media_tree_all + new_scale_fill()

circle_media_subject_tree_all <- gheatmap(circle_media_tree_all2, subject_df,
  colnames_angle = 85,
  colnames_offset_y = -0.2, # keep labels close
  offset = 0.1, # small spacing from previous ring
  width = 0.1, # slimmer ring
  font.size = 5, color = NA
) +
  scale_fill_manual(
    values = subject_colors_tree,
    breaks = names(subject_colors_tree)
  ) +
  theme(
    legend.position = "right", # hide legend
    plot.margin = margin(0, 0, 0, 0)
  ) + # minimize extra whitespace
  labs(fill = "Subject type")

circle_media_subject_tree_all2 <- circle_media_subject_tree_all + new_scale_fill()


circle_media_subject_genus_tree_all <- gheatmap(circle_media_subject_tree_all2, genus_df,
  colnames_angle = 85,
  colnames_offset_y = -0.2, # keep labels close
  offset = 0.2, # small spacing from previous ring
  width = 0.1, # slimmer ring
  font.size = 5, color = NA
) +
  scale_fill_manual(
    values = genus_colors_tree,
    breaks = names(genus_colors_tree)
  ) +
  theme(
    legend.position = "right", # hide legend
    plot.margin = margin(0, 0, 0, 0)
  ) + # minimize extra whitespace
  labs(fill = "Genus")

circle_media_subject_genus_tree_all2 <- circle_media_subject_genus_tree_all + new_scale_fill()


# --- Step 1: Extract tree data ---
# --- Step 1: Extract tree tips ---
tree_data <- circle_media_subject_genus_tree_all2$data
tip_data <- tree_data %>% filter(isTip)

# --- Step 2: Initial first tip for each genus ---
genus_labels <- tip_data %>%
  left_join(genus_df %>% mutate(tip = rownames(genus_df)), by = c("label" = "tip")) %>%
  group_by(genus) %>%
  arrange(angle) %>%
  slice_max(angle) %>%
  ungroup()

# --- Step 3: Cheat to avoid overlap ---
min_spacing <- 5 # degrees along the circle
genus_labels <- genus_labels %>% arrange(angle)

for (i in 2:nrow(genus_labels)) {
  angle_diff <- genus_labels$angle[i] - genus_labels$angle[i - 1]
  if (angle_diff < min_spacing) {
    this_genus <- genus_labels$genus[i]
    tips_of_genus <- tip_data %>%
      filter(label %in% (genus_df %>% filter(genus == this_genus) %>% rownames())) %>%
      arrange(angle)
    # pick second tip if available
    if (nrow(tips_of_genus) > 1) {
      genus_labels$x[i] <- tips_of_genus$x[2]
      genus_labels$y[i] <- tips_of_genus$y[2]
      genus_labels$angle[i] <- tips_of_genus$angle[2]
    }
  }
}

# --- Step 4: Radial offset beyond heatmaps ---
heatmap_offsets <- c(0.0, 0.1, 0.2)
heatmap_widths <- c(0.1, 0.1, 0.1)
radial_offset <- max(heatmap_offsets + heatmap_widths) + 1.1
genus_labels$x <- radial_offset

# --- Step 5: Adjust label_angle for readability ---
genus_labels$label_angle <- genus_labels$angle
genus_labels$label_angle <- ifelse(
  genus_labels$label_angle > 90 & genus_labels$label_angle < 270,
  genus_labels$label_angle + 180, # flip left side
  genus_labels$label_angle
)

# --- Step 6: Plot ---
final_plot <- circle_media_subject_genus_tree_all2 +
  geom_text(
    data = genus_labels,
    aes(x = x, y = y, label = genus, angle = label_angle),
    size = 3,
    color = "black",
    hjust = 0.5,
    vjust = 0.5
  )

final_plot


plot_file <- file.path("IsolateAnalysis/Tree/250917_figures/", paste0(Sys.Date(), "circle_media_subject_genus_labeled.pdf"))
# ggsave(filename = plot_file, plot = final_plot, width = 18, height = 18)

writeoutmetadata <- species_df %>%
  rownames_to_column("rowname") %>%
  left_join(cultureswab %>% rownames_to_column("rowname"),
    by = "rowname"
  ) %>%
  column_to_rownames("rowname")

output_file <- file.path("251129_isolatetreemetadata.csv")

# ---------------------------
# PCA-like MDS from tree + stats comparing distances
# ---------------------------

library(ape)

library(reshape2)

# 1) Cophenetic distances from the tree (use tree_all which you created above)
#    If you want a different tree, change to read.tree(...) accordingly.
dist_mat <- cophenetic.phylo(tree_all) # matrix of pairwise patristic distances

# 2) Classical MDS (2D coordinates)
mds_res <- cmdscale(as.dist(dist_mat), k = 2, eig = TRUE) # k = 2 for 2D
mds_df <- data.frame(
  label = rownames(mds_res$points),
  Dim1 = mds_res$points[, 1],
  Dim2 = mds_res$points[, 2],
  stringsAsFactors = FALSE
)

# 3) Prepare metadata for merging
#    Attempt to use cultureswab rownames (you set cultureswab rownames to DomeNumber earlier)
#    Make a small metadata frame: label must match tip names in the tree
meta <- cultureswab %>%
  mutate(label = row.names(cultureswab)) %>%
  select(label, color_group, color_code)

# Quick check: which mds labels are present in meta
missing_in_meta <- setdiff(mds_df$label, meta$label)
if (length(missing_in_meta) > 0) {
  warning(
    "Some tree tip labels do not match cultureswab labels. Number missing: ", length(missing_in_meta),
    ". First few missing: ", paste0(head(missing_in_meta, 10), collapse = ", ")
  )
}
# Keep only matching tips to avoid merge issues
mds_df <- inner_join(mds_df, meta, by = "label")

# 4) Plot MDS colored by subject type (Cow / Farmer / Non-Farmer)
#    Use scale_color_manual so groups get the correct hex codes.
group_colors <- meta %>%
  distinct(color_group, color_code) %>%
  tidyr::drop_na()
col_map <- setNames(group_colors$color_code, group_colors$color_group)

p_mds <- ggplot(mds_df, aes(x = Dim1, y = Dim2, color = color_group)) +
  geom_point(size = 3) +
  stat_ellipse(aes(fill = color_group), geom = "polygon", alpha = 0.12, show.legend = FALSE) +
  scale_color_manual(values = col_map) +
  scale_fill_manual(values = col_map) +
  theme_minimal(base_size = 14) +
  labs(
    title = "MDS (cmdscale) of tree cophenetic distances",
    x = "MDS1", y = "MDS2", color = "Subject type"
  )
print(p_mds)

# 5) Extract pairwise distances (long format) and attach group labels
dist_mat <- as.matrix(dist_mat)
dist_long <- melt(dist_mat,
  varnames = c("Sample1", "Sample2"), value.name = "TreeDist",
  stringsAsFactors = FALSE
)

# remove self comparisons and duplicates (keep Sample1 < Sample2 lexicographically)
dist_long$Sample1 <- as.character(dist_long$Sample1)
dist_long$Sample2 <- as.character(dist_long$Sample2)
dist_long <- dist_long[dist_long$Sample1 != dist_long$Sample2, ]
# To avoid duplicates, keep only one of each unordered pair:
dist_long <- dist_long[!duplicated(t(apply(dist_long[, c("Sample1", "Sample2")], 1, function(x) paste(sort(x), collapse = "__")))), ]

# Merge group info for both ends
meta_small <- meta # label, color_group, color_code
colnames(meta_small)[colnames(meta_small) == "color_group"] <- "Group"
dist_long <- left_join(dist_long, meta_small, by = c("Sample1" = "label"))
colnames(dist_long)[ncol(dist_long)] <- "Group1"
dist_long <- left_join(dist_long, meta_small, by = c("Sample2" = "label"))
colnames(dist_long)[ncol(dist_long)] <- "Group2"

# Remove any pairs where group info is missing
dist_long <- dist_long %>% filter(!is.na(Group1) & !is.na(Group2))

# 6) Build two sets: Cow–Farmer and Cow–Non-Farmer (both directions included)
cow_farmer <- dist_long %>%
  filter((Group.x == "Cow" & Group.y == "Farmer") | (Group.x == "Farmer" & Group.y == "Cow")) %>%
  pull(TreeDist)

cow_nonfarmer <- dist_long %>%
  filter((Group.x == "Cow" & Group.y == "Non-Farmer") | (Group.x == "Non-Farmer" & Group.y == "Cow")) %>%
  pull(TreeDist)

# Quick summaries
cat("Cow-Farmer distances: n =", length(cow_farmer), "mean =", mean(cow_farmer), "median =", median(cow_farmer), "\n")
cat("Cow-NonFarmer distances: n =", length(cow_nonfarmer), "mean =", mean(cow_nonfarmer), "median =", median(cow_nonfarmer), "\n")

# 7) Statistical test: compare distributions (non-parametric Wilcoxon rank-sum)
#    Use wilcox.test (Mann-Whitney). If very small sample sizes, interpret cautiously.
if (length(cow_farmer) >= 2 && length(cow_nonfarmer) >= 2) {
  wilcox_res <- wilcox.test(cow_farmer, cow_nonfarmer, alternative = "two.sided")
  print(wilcox_res)
} else {
  warning("Not enough pairwise distances in one of the groups to run Wilcoxon test.")
}

# 8) Boxplot showing the two distance distributions side-by-side
plot_df <- data.frame(
  Distance = c(cow_farmer, cow_nonfarmer),
  Comparison = c(rep("Cow-Farmer", length(cow_farmer)), rep("Cow-NonFarmer", length(cow_nonfarmer)))
)

p_dist_box <- ggplot(plot_df, aes(x = Comparison, y = Distance, fill = Comparison)) +
  geom_boxplot() +
  geom_jitter(width = 0.15, size = 1, alpha = 0.6) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Tree distances: Cow vs Farmer and Cow vs Non-Farmer",
    y = "Patristic (cophenetic) distance"
  ) +
  scale_fill_brewer(palette = "Set2")
print(p_dist_box)

# Print p-value if computed
if (exists("wilcox_res")) {
  cat("Wilcoxon p-value:", wilcox_res$p.value, "\n")
}
