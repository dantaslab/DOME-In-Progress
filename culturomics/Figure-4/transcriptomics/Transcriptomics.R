library(DESeq2)
library(purrr)
library(rtracklayer)
library(dplyr)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)
library(vegan)
library(VennDiagram)
library(tidyverse)
library(DESeq2)
library(broom)
library(car)
library(clusterProfiler)
library(stringr)
library(httr)
library(ggrepel)
source("utilities.R")
library(cowplot)

# Set to parent folder containing all sample folders
setwd("Transcriptomics/Analysis/250626run/250627_d03_featurecounts")
# List all subfolders that match your sample naming pattern
sample_dirs <- list.dirs(recursive = FALSE)
sample_dirs <- sample_dirs[grepl("^JE2_", basename(sample_dirs))]

# Build full paths to the sample.txt files
sample_files <- file.path(sample_dirs, "sample.txt")
names(sample_files) <- basename(sample_dirs) # name by sample


# Function to extract the counts
read_counts <- function(file, sample_name) {
  df <- read.delim(file, comment.char = "#")

  # Inspect the available columns if unsure
  if (!"Geneid" %in% names(df)) stop(paste("Geneid column not found in:", file))

  # Count column is usually the last column (or second to last if there's a "Length" column)
  count_col <- setdiff(names(df), c("Geneid", "Chr", "Start", "End", "Strand", "Length"))

  if (length(count_col) != 1) stop(paste("Unexpected number of count columns in:", file))

  counts <- df[, c("Geneid", count_col)]
  colnames(counts)[2] <- sample_name
  return(counts)
}

# Read in all sample files
count_list <- mapply(read_counts, sample_files, names(sample_files), SIMPLIFY = FALSE)

# Merge all into a single count matrix
counts_merged <- purrr::reduce(count_list, full_join, by = "Geneid")
rownames(counts_merged) <- counts_merged$Geneid
count_matrix <- as.matrix(counts_merged[, -1])
mode(count_matrix) <- "integer"

sample_info <- tibble(
  sample = colnames(count_matrix)
) %>%
  separate(sample, into = c("Strain", "Media", "Time", "Replicate"), sep = "_", remove = FALSE) %>%
  mutate(condition = paste(Media, Time, sep = "_")) %>%
  column_to_rownames("sample")
# confirm rownames equal column names for deseq

sample_info <- sample_info %>%
  mutate(Media_category = ifelse(grepl("^[0-9]+$", Media),
    "Crossfeed",
    Media
  ))

# The Media_category column will be used towards the end, when I start looking at higher level comparisons for figures that will be better in the main text

all(rownames(sample_info) == colnames(count_matrix))


dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ Media + Time + Media:Time
)

dds$Media <- relevel(factor(dds$Media), ref = "MMM") # Control media
dds$Time <- factor(dds$Time, levels = c("T5", "T9")) # Optional

dds <- DESeq(dds)

vsd <- vst(dds, blind = FALSE)
rld <- rlog(dds, blind = FALSE)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep = "-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
  clustering_distance_rows = sampleDists,
  clustering_distance_cols = sampleDists,
  col = colors
)
# Save these to the supplement, maybe add in confidence interval circles
pcacondition <- plotPCA(vsd, intgroup = c("condition"))
plot_file <- paste0("~/Library/CloudStorage/Box-Box/2024_Dome2.0/Transcriptomics/Analysis/250626run/Figures/250922/", Sys.Date(), "PCAcondition.pdf")

pcatime <- plotPCA(vsd, intgroup = c("Time")) # Time clearly has an effect, going to do
plot_file <- paste0("~/Library/CloudStorage/Box-Box/2024_Dome2.0/Transcriptomics/Analysis/250626run/Figures/250922/", Sys.Date(), "PCAtime.pdf")

pcamedia <- plotPCA(vsd, intgroup = c("Media")) # MMMG clusters out really strongly, I think it's something to think about
plot_file <- paste0("~/Library/CloudStorage/Box-Box/2024_Dome2.0/Transcriptomics/Analysis/250626run/Figures/250922/", Sys.Date(), "PCAMedia.pdf")

keep_genes <- rowSums(count_matrix >= 10) >= 2
filtered_counts <- count_matrix[keep_genes, ]

dds <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = sample_info,
  design = ~ Media + Time + Media:Time
)

# Ensure correct reference levels
dds$Media <- relevel(dds$Media, ref = "MMM") # Media reference = MMM
dds$Time <- relevel(dds$Time, ref = "T5") # Time reference = T5

# Run DESeq (Wald test for individual contrasts)
dds <- DESeq(dds)

#----------------------------
# 2. Wald tests
#----------------------------


# (a) Differences between each Media vs MMM at T5 (reference timepoint)
res_194_T5 <- results(dds, contrast = c("Media", "194", "MMM"))
res_201_T5 <- results(dds, contrast = c("Media", "201", "MMM"))
res_213_T5 <- results(dds, contrast = c("Media", "213", "MMM"))
res_216_T5 <- results(dds, contrast = c("Media", "216", "MMM"))
res_MMMG_T5 <- results(dds, contrast = c("Media", "MMMG", "MMM"))

# (b) Differences between each Media vs MMM at T9
# Combine: main effect of media + interaction term

res_194_T9 <- results(dds, contrast = list(c("Media_194_vs_MMM", "Media194.TimeT9")))
res_201_T9 <- results(dds, contrast = list(c("Media_201_vs_MMM", "Media201.TimeT9")))
res_213_T9 <- results(dds, contrast = list(c("Media_213_vs_MMM", "Media213.TimeT9")))
res_216_T9 <- results(dds, contrast = list(c("Media_216_vs_MMM", "Media216.TimeT9")))
res_MMMG_T9 <- results(dds, contrast = list(c("Media_MMMG_vs_MMM", "MediaMMMG.TimeT9")))

# Optional: check the names of coefficients to avoid typos
resultsNames(dds)

#----------------------------
# 3. LRT: Is there an overall Media effect (agnostic of Time)?
#----------------------------

dds_lrt <- DESeq(dds,
  test = "LRT",
  reduced = ~Time
)

res_LRT_media <- results(dds_lrt)

# This tells you which genes have a significant average difference between media types.


#----------------------------
# 4. plotting
#----------------------------

par(mfrow = c(3, 4)) # 3 rows, 4 cols for compact layout

plotMA(res_194_T5, main = "194 vs MMM (T5)", ylim = c(-5, 5))
plotMA(res_201_T5, main = "201 vs MMM (T5)", ylim = c(-5, 5))
plotMA(res_213_T5, main = "213 vs MMM (T5)", ylim = c(-5, 5))
plotMA(res_216_T5, main = "216 vs MMM (T5)", ylim = c(-5, 5))
plotMA(res_MMMG_T5, main = "MMMG vs MMM (T5)", ylim = c(-5, 5))

plotMA(res_194_T9, main = "194 vs MMM (T9)", ylim = c(-5, 5))
plotMA(res_201_T9, main = "201 vs MMM (T9)", ylim = c(-5, 5))
plotMA(res_213_T9, main = "213 vs MMM (T9)", ylim = c(-5, 5))
plotMA(res_216_T9, main = "216 vs MMM (T9)", ylim = c(-5, 5))
plotMA(res_MMMG_T9, main = "MMMG vs MMM (T9)", ylim = c(-5, 5))

plotMA(res_LRT_media, main = "Media Effect (LRT)", ylim = c(-5, 5))


#----------------------------
# 5. Extract significant results
#----------------------------

sig_194_T5 <- subset(res_194_T5, padj < 0.05)
sig_201_T5 <- subset(res_201_T5, padj < 0.05)
sig_213_T5 <- subset(res_213_T5, padj < 0.05)
sig_216_T5 <- subset(res_216_T5, padj < 0.05)
sig_MMMG_T5 <- subset(res_MMMG_T5, padj < 0.05)

sig_194_T9 <- subset(res_194_T9, padj < 0.05)
sig_201_T9 <- subset(res_201_T9, padj < 0.05)
sig_213_T9 <- subset(res_213_T9, padj < 0.05)
sig_216_T9 <- subset(res_216_T9, padj < 0.05)
sig_MMMG_T9 <- subset(res_MMMG_T9, padj < 0.05)

sig_LRT_media <- subset(res_LRT_media, padj < 0.05)


############
############ Making the same figures at a higher level (Crossfeed vs MMM), and adding in annotations for main figure text
############

dds_cross <- DESeqDataSetFromMatrix(
  countData = filtered_counts,
  colData = sample_info,
  design = ~ Media_category + Time + Media_category:Time
)

# Ensure correct reference levels
dds_cross$Media_category <- relevel(dds_cross$Media_category, ref = "MMM") # Media reference = MMM
dds_cross$Time <- relevel(dds_cross$Time, ref = "T5") # Time reference = T5

# Run DESeq (Wald test for individual contrasts)
dds_cross <- DESeq(dds_cross)

resultsNames(dds_cross)

#----------------------------
# 2. Wald tests
#----------------------------
gtf_file <- "Transcriptomics/Analysis/250626run/genomic.gtf"

# Read the GTF file
gtf <- import(gtf_file)


# Filter to only 'gene' type entries
gene_annot <- gtf[gtf$type == "gene"]

# Extract data frame with gene_id and gene_name
gene_df <- as.data.frame(mcols(gene_annot)) %>%
  select(gene_id, gene) %>%
  distinct()


# (a) Differences between each Media_category vs MMM at T5 (reference timepoint)
res_cross_T5 <- results(dds_cross, contrast = c("Media_category", "Crossfeed", "MMM"))

res_cross_T5_df <- as.data.frame(res_cross_T5)
res_cross_T5_df$gene_id <- rownames(res_cross_T5_df)

# Join with gene name annotations
res_cross_T5_df_annotated <- left_join(res_cross_T5_df, gene_df, by = c("gene_id" = "gene_id"))


res_cross_T9 <- results(dds_cross, contrast = list(c("Media_category_Crossfeed_vs_MMM", "Media_categoryCrossfeed.TimeT9")))

res_cross_T9_df <- as.data.frame(res_cross_T9)
res_cross_T9_df$gene_id <- rownames(res_cross_T9_df)

# Join with gene name annotations
res_cross_T9_df_annotated <- left_join(res_cross_T9_df, gene_df, by = c("gene_id" = "gene_id"))


output_file <- file.path("Fig4/Main_T5/251129_transcriptomicsT5statscrossfeed.csv")

output_file <- file.path("Fig4/Supp_T9/251129_transcriptomicsT9statscrossfeed.csv")


# Modified volcanoplot function


modified_volcano_plot <- function(df, title, highlight_genes = NULL, top_n = 10, cap_quantile = 0.99) {
  # Handle cases where padj is NA or 0
  df$padj[is.na(df$padj) | df$padj == 0] <- 1e-300
  df$logP <- -log10(df$padj)

  # Significance logic
  df$sig <- ifelse(df$padj < 0.05 & abs(df$log2FoldChange) > 1, "yes", "no")

  # Dynamic cap for plotting
  y_max <- ceiling(quantile(df$logP, cap_quantile, na.rm = TRUE))
  df$logP_plot <- pmin(df$logP, y_max) # capped for points

  # Flag points exceeding the cap
  df$above_cap <- df$logP > y_max

  # Determine top N genes by absolute log2FoldChange
  top_genes <- head(df$gene[order(-abs(df$log2FoldChange))], top_n)
  label_genes <- if (!is.null(highlight_genes)) unique(c(top_genes, highlight_genes)) else top_genes
  df$label <- ifelse(df$gene %in% label_genes, df$gene, "")

  # Plot
  ggplot(df, aes(x = log2FoldChange, y = logP_plot, color = sig)) +
    # regular points
    geom_point(alpha = 0.4, size = .8, data = df[!df$above_cap, ]) +
    # triangle markers for points above cap
    geom_point(shape = 24, alpha = 0.4, fill = "#9F00BB", size = 1, data = df[df$above_cap, ]) +
    # labels use true logP
    geom_text_repel(aes(label = label, y = logP, x = log2FoldChange), size = 8, max.overlaps = Inf, segment.color = "grey50", segment.size = .5, segment.alpha = .7) +
    scale_color_manual(values = c("no" = "grey", "yes" = "#9F00BB")) +
    theme_minimal(base_size = 14) +
    labs(title = title, x = "Log2 Fold Change", y = expression(-log[10](padj))) +
    coord_cartesian(ylim = c(0, y_max)) +
    theme(legend.position = "none", plot.title = element_text(hjust = 0.5)) +
    theme_pub() +
    theme(
      axis.text.y = element_text(size = 16, face = "bold") # bigger, bold category labels
    ) +
    theme(
      axis.text.x = element_text(size = 16, face = "bold") # bigger, bold category labels
    )
}

highlight <- c("tet(38)", "spa", "hld", "coa")
p1 <- modified_volcano_plot(res_cross_T5_df_annotated, "C. funkei supernatant vs MMM at T5", highlight_genes = highlight, top_n = 20)

plot_file <- paste0("~/Library/CloudStorage/Box-Box/2024_Dome2.0/Transcriptomics/Analysis/250626run/Figures/250922/", Sys.Date(), "crossfeedvsMMMT5.pdf")


p2 <- modified_volcano_plot(res_cross_T9_df_annotated, "C. funkei supernatant vs MMM at T9", highlight_genes = highlight, top_n = 20)

plot_file <- paste0("~/Library/CloudStorage/Box-Box/2024_Dome2.0/Transcriptomics/Analysis/250626run/Figures/250922/", Sys.Date(), "crossfeedvsMMMT9.pdf")

sig_cross_T5 <- subset(res_cross_T5, padj < 0.05)
sig_genes_T5 <- rownames(sig_cross_T5)
bg_genes_T5 <- rownames(res_cross_T5)

sig_cross_T9 <- subset(res_cross_T9, padj < 0.05)
sig_genes_T9 <- rownames(sig_cross_T9)
bg_genes_T9 <- rownames(res_cross_T9)


# -------------------------
# Setup mapping tables
# -------------------------
gtf_file <- "Transcriptomics/Analysis/250626run/genomic.gtf"
gtf <- rtracklayer::import(gtf_file)

gtf_map <- mcols(gtf) %>%
  as.data.frame() %>%
  select(gene_id, protein_id) %>%
  filter(!is.na(protein_id)) %>%
  distinct()

eggnog <- read_tsv(
  "Transcriptomics/Analysis/250626run/250718_webEGGNOGrun/MM_pj7nwwpq.emapper.annotations.tsv",
  comment = "#"
)

ko_map <- eggnog %>%
  select(query, KEGG_ko) %>%
  filter(!is.na(KEGG_ko)) %>%
  separate_rows(KEGG_ko, sep = ",") %>%
  rename(protein_id = query, ko = KEGG_ko) %>%
  mutate(ko = sub("ko:", "", ko)) %>%
  distinct()

ko2pathway <- read.delim("https://rest.kegg.jp/link/pathway/ko",
  header = FALSE, sep = "\t",
  col.names = c("ko", "pathway")
) %>%
  mutate(
    ko = sub("ko:", "", ko),
    pathway = sub("path:", "", pathway)
  )

pathway_names <- read.delim("https://rest.kegg.jp/list/pathway/ko",
  header = FALSE, sep = "\t",
  col.names = c("pathway", "description")
) %>%
  mutate(pathway = sub("path:", "", pathway))

protein2pathway <- inner_join(ko_map, ko2pathway, by = "ko") %>%
  inner_join(pathway_names, by = "pathway") %>%
  distinct()

# -------------------------
# Function: map sig genes → KEGG pathways
# -------------------------
map_to_pathways <- function(sig_genes, gtf_map, protein2pathway) {
  sig_proteins <- gtf_map %>%
    filter(gene_id %in% sig_genes) %>%
    pull(protein_id)

  df <- protein2pathway %>%
    filter(protein_id %in% sig_proteins) %>%
    select(description, protein_id) %>%
    distinct()

  if (nrow(df) == 0) {
    warning("No pathway mapping found for this set of genes")
    return(tibble(description = character(), protein_count = integer()))
  }

  df <- dplyr::count(df, description, sort = TRUE, name = "protein_count")

  return(df)
}

# -------------------------
# Example for T5
# -------------------------
sig_cross_T5 <- subset(res_cross_T5, padj < 0.05)
sig_genes_T5 <- rownames(sig_cross_T5)

pathway_counts_T5 <- map_to_pathways(sig_genes_T5, gtf_map, protein2pathway)

sig_cross_T9 <- subset(res_cross_T9, padj < 0.05)
sig_genes_T9 <- rownames(sig_cross_T9)

pathway_counts_T9 <- map_to_pathways(sig_genes_T9, gtf_map, protein2pathway)

# Modules

gtf_file <- "Transcriptomics/Analysis/250626run/genomic.gtf"
gtf <- rtracklayer::import(gtf_file)

gtf_map <- mcols(gtf) %>%
  as.data.frame() %>%
  select(gene_id, protein_id) %>%
  filter(!is.na(protein_id)) %>%
  distinct()

# -------------------------
# Step 2: Identify significant genes from DESeq2
# -------------------------
# Example: res_cross_T5
sig_cross_T5 <- subset(res_cross_T5, padj < 0.05)
sig_genes_T5 <- rownames(sig_cross_T5)

# Map sig genes → proteins
sig_proteins_T5 <- gtf_map %>%
  filter(gene_id %in% sig_genes_T5) %>%
  pull(protein_id)


sig_cross_T9 <- subset(res_cross_T9, padj < 0.05)
sig_genes_T9 <- rownames(sig_cross_T9)

# Map sig genes → proteins
sig_proteins_T9 <- gtf_map %>%
  filter(gene_id %in% sig_genes_T9) %>%
  pull(protein_id)


# -------------------------
# Step 3: Load EggNOG mapping (proteins → KOs)
# -------------------------
eggnog <- read_tsv(
  "Transcriptomics/Analysis/250626run/250718_webEGGNOGrun/MM_pj7nwwpq.emapper.annotations.tsv",
  comment = "#"
)

ko_map <- eggnog %>%
  select(query, KEGG_ko) %>%
  filter(!is.na(KEGG_ko)) %>%
  separate_rows(KEGG_ko, sep = ",") %>%
  rename(protein_id = query, ko = KEGG_ko) %>%
  mutate(ko = sub("ko:", "", ko)) %>%
  distinct()

# -------------------------
# Step 4: Load module mapping (KO → module) and module descriptions
# -------------------------
module_ko <- read.delim(
  "https://rest.kegg.jp/link/ko/module",
  header = FALSE,
  sep = "\t",
  col.names = c("module", "ko"),
  stringsAsFactors = FALSE
)

module_desc <- read.delim(
  "https://rest.kegg.jp/list/module",
  header = FALSE,
  sep = "\t",
  col.names = c("module", "description"),
  stringsAsFactors = FALSE
)

module_ko_clean <- module_ko %>%
  mutate(
    ko = sub("^ko:", "", ko), # remove 'ko:' from KOs
    module = sub("^md:", "", module) # remove 'md:' from modules
  )

module_desc_clean <- module_desc %>%
  mutate(
    module = sub("^md:", "", module) # remove 'md:' from module IDs
  )


# -------------------------
# Step 5: Map proteins → KOs → modules
# -------------------------
map_to_modules <- function(sig_proteins, ko_map, module_ko, module_desc) {
  # Proteins → KOs
  protein_ko <- ko_map %>%
    filter(protein_id %in% sig_proteins) %>%
    select(protein_id, ko) %>%
    distinct()

  if (nrow(protein_ko) == 0) {
    warning("No KOs found for these proteins")
    return(tibble(module = character(), description = character(), protein_count = integer()))
  }

  # KOs → modules (with cleaned IDs)
  df <- protein_ko %>%
    inner_join(module_ko_clean, by = "ko") %>%
    inner_join(module_desc_clean, by = "module") %>%
    distinct()

  if (nrow(df) == 0) {
    warning("No modules found for these proteins")
    return(tibble(module = character(), description = character(), protein_count = integer()))
  }

  # Count proteins per module
  df_counts <- dplyr::count(df, module, description, sort = TRUE, name = "protein_count")

  return(df_counts)
}
module_counts_T5 <- map_to_modules(sig_proteins_T5, ko_map, module_ko, module_desc)
module_counts_T9 <- map_to_modules(sig_proteins_T9, ko_map, module_ko, module_desc)


# -------------------------
# Function to create annotated gene labels
# -------------------------
create_gene_labels <- function(res_df_annotated, gtf_map, protein2pathway, module_counts) {
  # Start with gene_id and gene columns
  labels_df <- res_df_annotated %>%
    select(gene_id, gene, log2FoldChange, padj) %>%
    mutate(
      # Use 'gene' column if available, otherwise trim B7H15_ from gene_id
      display_name = ifelse(!is.na(gene) & gene != "",
        gene,
        sub("^B7H15_", "", gene_id)
      )
    )

  # Map gene_id -> protein_id
  labels_df <- labels_df %>%
    left_join(gtf_map %>% select(gene_id, protein_id), by = "gene_id")

  # Get pathway annotations (trim at comma or semicolon)
  pathway_annot <- protein2pathway %>%
    select(protein_id, description) %>%
    distinct() %>%
    mutate(pathway_short = str_extract(description, "^[^,;]+")) %>%
    group_by(protein_id) %>%
    slice(1) %>% # Take first pathway if multiple
    ungroup() %>%
    select(protein_id, pathway_short)

  # Get module annotations (trim at comma or semicolon)
  module_annot <- module_counts %>%
    select(module, description) %>%
    distinct() %>%
    mutate(module_short = str_extract(description, "^[^,;]+"))

  # Map protein_id to modules via ko_map and module_ko_clean
  protein_modules <- ko_map %>%
    inner_join(module_ko_clean, by = "ko") %>%
    left_join(module_annot, by = "module") %>%
    group_by(protein_id) %>%
    slice(1) %>% # Take first module if multiple
    ungroup() %>%
    select(protein_id, module_short)

  # Join pathway and module info
  labels_df <- labels_df %>%
    left_join(pathway_annot, by = "protein_id") %>%
    left_join(protein_modules, by = "protein_id")

  # Create final label: gene | pathway (or module if no pathway)
  labels_df <- labels_df %>%
    mutate(
      annotation = case_when(
        !is.na(pathway_short) ~ paste0(" | ", pathway_short),
        !is.na(module_short) ~ paste0(" | ", module_short),
        TRUE ~ ""
      ),
      final_label = paste0(display_name, annotation)
    )

  return(labels_df)
}

# -------------------------
# Create heatmaps for T5 and T9
# -------------------------

# Get top 30 genes by absolute log2FoldChange for T5
top30_T5 <- res_cross_T5_df_annotated %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(30)

# Create labels
labels_T5 <- create_gene_labels(top30_T5, gtf_map, protein2pathway, module_counts_T5)

# Get the gene_ids in order
top_genes_T5 <- labels_T5$gene_id

# Filter samples to only Crossfeed and MMM
samples_to_include <- sample_info %>%
  filter(Media_category %in% c("Crossfeed", "MMM") & Time == "T5") %>%
  rownames()

# Extract normalized counts
vsd_cross <- vst(dds_cross)
mat_T5 <- assay(vsd_cross)[top_genes_T5, samples_to_include]

# Center rows
mat_T5 <- mat_T5 - rowMeans(mat_T5)

# Set row names to final labels
rownames(mat_T5) <- labels_T5$final_label

# Create sample annotation
ann_T5 <- as.data.frame(colData(dds_cross)[samples_to_include, c("Media_category", "Time")])

# Plot T5 heatmap
pheatmap(mat_T5,
  annotation_col = ann_T5,
  show_rownames = TRUE,
  cluster_rows = FALSE, # Keep in log2FC order
  cluster_cols = TRUE,
  main = "Top 30 DE genes: Crossfeed vs MMM (T5)",
  fontsize_row = 8
)

# -------------------------
# Repeat for T9
# -------------------------

# Get top 30 genes by absolute log2FoldChange for T9
top30_T9 <- res_cross_T9_df_annotated %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(30)

# Create labels
labels_T9 <- create_gene_labels(top30_T9, gtf_map, protein2pathway, module_counts_T9)

# Get the gene_ids in order
top_genes_T9 <- labels_T9$gene_id

# Filter samples to only Crossfeed and MMM at T9
samples_to_include_T9 <- sample_info %>%
  filter(Media_category %in% c("Crossfeed", "MMM") & Time == "T9") %>%
  rownames()

# Extract normalized counts
mat_T9 <- assay(vsd_cross)[top_genes_T9, samples_to_include_T9]

# Center rows
mat_T9 <- mat_T9 - rowMeans(mat_T9)

# Set row names to final labels
rownames(mat_T9) <- labels_T9$final_label

# Create sample annotation
ann_T9 <- as.data.frame(colData(dds_cross)[samples_to_include_T9, c("Media_category", "Time")])

# Plot T9 heatmap
pheatmap(mat_T9,
  annotation_col = ann_T9,
  show_rownames = TRUE,
  cluster_rows = FALSE, # Keep in log2FC order
  cluster_cols = TRUE,
  main = "Top 30 DE genes: Crossfeed vs MMM (T9)",
  fontsize_row = 8
)

# -------------------------
# Optional: Save the plots
# -------------------------

# Save T5
# pdf("~/Library/CloudStorage/Box-Box/2024_Dome2.0/Transcriptomics/Analysis/250626run/Figures/250922/heatmap_crossfeed_T5.pdf",
#     width = 10, height = 12)
plot_heatmap_mat_T5 <- pheatmap(mat_T5,
  annotation_col = ann_T5,
  show_rownames = TRUE,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  main = "Top 30 DE genes: Crossfeed vs MMM (T5)",
  fontsize_row = 8
)
plot_file <- paste0("~/Library/CloudStorage/Box-Box/2024_Dome2.0/Transcriptomics/Analysis/250626run/Figures/250922/", Sys.Date(), "heatmap_crossfeed_T5.pdf")

# Save T9

plot_heatmap_mat_T9 <- pheatmap(mat_T9,
  annotation_col = ann_T9,
  show_rownames = TRUE,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  main = "Top 30 DE genes: Crossfeed vs MMM (T9)",
  fontsize_row = 8
)
plot_file <- paste0("~/Library/CloudStorage/Box-Box/2024_Dome2.0/Transcriptomics/Analysis/250626run/Figures/250922/", Sys.Date(), "heatmap_crossfeed_T9.pdf")

# -------------------------
# Function: map sig genes → KEGG pathways with directionality
# -------------------------
map_to_pathways_directional <- function(res_df, gtf_map, protein2pathway, padj_cutoff = 0.05) {
  # Convert to data frame if needed and ensure gene_id column exists
  if (!is.data.frame(res_df)) {
    res_df <- as.data.frame(res_df)
    res_df$gene_id <- rownames(res_df)
  }

  # Filter significant genes and add direction
  sig_df <- res_df %>%
    filter(padj < padj_cutoff) %>%
    mutate(direction = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated")) %>%
    select(gene_id, direction)

  # Map to proteins
  sig_proteins_dir <- gtf_map %>%
    inner_join(sig_df, by = "gene_id") %>%
    select(protein_id, direction) %>%
    distinct()

  # Map to pathways
  df <- protein2pathway %>%
    ungroup() %>% # Ensure not grouped
    inner_join(sig_proteins_dir, by = "protein_id") %>%
    select(description, protein_id, direction) %>%
    distinct() %>%
    ungroup() # Ungroup again after operations

  if (nrow(df) == 0) {
    warning("No pathway mapping found for this set of genes")
    return(tibble(description = character(), direction = character(), protein_count = integer()))
  }

  # Count proteins per pathway and direction
  df_counts <- df %>%
    group_by(description, direction) %>%
    summarise(protein_count = n(), .groups = "drop") %>%
    arrange(desc(protein_count))

  return(df_counts)
}

# -------------------------
# Function: map sig genes → KEGG modules with directionality
# -------------------------
map_to_modules_directional <- function(res_df, gtf_map, ko_map, module_ko_clean, module_desc_clean, padj_cutoff = 0.05) {
  # Convert to data frame if needed and ensure gene_id column exists
  if (!is.data.frame(res_df)) {
    res_df <- as.data.frame(res_df)
    res_df$gene_id <- rownames(res_df)
  }

  # Filter significant genes and add direction
  sig_df <- res_df %>%
    filter(padj < padj_cutoff) %>%
    mutate(direction = ifelse(log2FoldChange > 0, "Upregulated", "Downregulated")) %>%
    select(gene_id, direction)

  # Map to proteins
  sig_proteins_dir <- gtf_map %>%
    inner_join(sig_df, by = "gene_id") %>%
    select(protein_id, direction) %>%
    distinct()

  # Proteins → KOs
  protein_ko_dir <- ko_map %>%
    inner_join(sig_proteins_dir, by = "protein_id") %>%
    select(protein_id, ko, direction) %>%
    distinct()

  if (nrow(protein_ko_dir) == 0) {
    warning("No KOs found for these proteins")
    return(tibble(module = character(), description = character(), direction = character(), protein_count = integer()))
  }

  # KOs → modules
  df <- protein_ko_dir %>%
    inner_join(module_ko_clean, by = "ko") %>%
    inner_join(module_desc_clean, by = "module") %>%
    select(module, description, protein_id, direction) %>%
    distinct() %>%
    ungroup() # Ensure not grouped

  if (nrow(df) == 0) {
    warning("No modules found for these proteins")
    return(tibble(module = character(), description = character(), direction = character(), protein_count = integer()))
  }

  # Count proteins per module and direction
  df_counts <- df %>%
    group_by(module, description, direction) %>%
    summarise(protein_count = n(), .groups = "drop") %>%
    arrange(desc(protein_count))

  return(df_counts)
}

# -------------------------
# Generate directional pathway counts for T5 and T9
# -------------------------
pathway_counts_T5_dir <- map_to_pathways_directional(res_cross_T5_df_annotated, gtf_map, protein2pathway)
pathway_counts_T9_dir <- map_to_pathways_directional(res_cross_T9_df_annotated, gtf_map, protein2pathway)

# -------------------------
# Generate directional module counts for T5 and T9
# -------------------------
module_counts_T5_dir <- map_to_modules_directional(res_cross_T5_df_annotated, gtf_map, ko_map, module_ko_clean, module_desc_clean)
module_counts_T9_dir <- map_to_modules_directional(res_cross_T9_df_annotated, gtf_map, ko_map, module_ko_clean, module_desc_clean)

# -------------------------
# Get top 15 pathways/modules by total protein count
# -------------------------
top15_pathways_T5 <- pathway_counts_T5_dir %>%
  group_by(description) %>%
  summarise(total = sum(protein_count)) %>%
  slice_max(total, n = 15) %>%
  pull(description)

top15_pathways_T9 <- pathway_counts_T9_dir %>%
  group_by(description) %>%
  summarise(total = sum(protein_count)) %>%
  slice_max(total, n = 15) %>%
  pull(description)

top15_modules_T5 <- module_counts_T5_dir %>%
  group_by(description) %>%
  summarise(total = sum(protein_count)) %>%
  slice_max(total, n = 15) %>%
  pull(description)

top15_modules_T9 <- module_counts_T9_dir %>%
  group_by(description) %>%
  summarise(total = sum(protein_count)) %>%
  slice_max(total, n = 15) %>%
  pull(description)

# -------------------------
# Plot T5 Pathways with directionality
# -------------------------
p3_dir <- pathway_counts_T5_dir %>%
  filter(description %in% top15_pathways_T5) %>%
  mutate(description = str_wrap(description, width = 30)) %>%
  # Calculate total for ordering
  group_by(description) %>%
  mutate(total = sum(protein_count)) %>%
  ungroup() %>%
  # Make downregulated counts negative for diverging bar chart
  mutate(protein_count_plot = ifelse(direction == "Downregulated", -protein_count, protein_count)) %>%
  ggplot(aes(x = reorder(description, total), y = protein_count_plot, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Upregulated" = "#9F00BB", "Downregulated" = "#00BBB4")) +
  labs(
    title = "KEGG Pathways (T5): Up vs Down Regulation",
    x = "Pathway",
    y = "Protein Count",
    fill = "Direction"
  ) +
  theme_minimal(base_size = 14) +
  theme_pub() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5)

plot_file <- paste0("~/Library/CloudStorage/Box-Box/2024_Dome2.0/Transcriptomics/Analysis/250626run/Figures/250922/", Sys.Date(), "pathway_directional_crossfeedvsMMMT5.pdf")
# ggsave2(p3_dir, filename = plot_file, device=cairo_pdf, height=10, width = 10)

# -------------------------
# Plot T9 Pathways with directionality
# -------------------------
p4_dir <- pathway_counts_T9_dir %>%
  filter(description %in% top15_pathways_T9) %>%
  mutate(description = str_wrap(description, width = 30)) %>%
  group_by(description) %>%
  mutate(total = sum(protein_count)) %>%
  ungroup() %>%
  mutate(protein_count_plot = ifelse(direction == "Downregulated", -protein_count, protein_count)) %>%
  ggplot(aes(x = reorder(description, total), y = protein_count_plot, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Upregulated" = "#9F00BB", "Downregulated" = "#00BBB4")) +
  labs(
    title = "KEGG Pathways (T9): Up vs Down Regulation",
    x = "Pathway",
    y = "Protein Count",
    fill = "Direction"
  ) +
  theme_minimal(base_size = 14) +
  theme_pub() +
  theme(
    axis.text.y = element_text(size = 12, face = "bold"),
    axis.text.x = element_text(size = 12, face = "bold"),
    legend.position = "bottom"
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5)

plot_file <- paste0("~/Library/CloudStorage/Box-Box/2024_Dome2.0/Transcriptomics/Analysis/250626run/Figures/250922/", Sys.Date(), "pathway_directional_crossfeedvsMMMT9.pdf")
# ggsave2(p4_dir, filename = plot_file, device=cairo_pdf, height=10, width = 10)

# -------------------------
# Plot T5 Modules with directionality
# -------------------------
p5_dir <- module_counts_T5_dir %>%
  filter(description %in% top15_modules_T5) %>%
  mutate(description = str_wrap(description, width = 30)) %>%
  group_by(description) %>%
  mutate(total = sum(protein_count)) %>%
  ungroup() %>%
  mutate(protein_count_plot = ifelse(direction == "Downregulated", -protein_count, protein_count)) %>%
  ggplot(aes(x = reorder(description, total), y = protein_count_plot, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Upregulated" = "#9F00BB", "Downregulated" = "#00BBB4")) +
  labs(
    title = "KEGG Modules (T5): Up vs Down Regulation",
    x = "Module",
    y = "Protein Count",
    fill = "Direction"
  ) +
  theme_minimal(base_size = 14) +
  theme_pub() +
  theme(
    axis.text.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold"),
    legend.position = "bottom"
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5)

plot_file <- paste0("~/Library/CloudStorage/Box-Box/2024_Dome2.0/Transcriptomics/Analysis/250626run/Figures/250922/", Sys.Date(), "module_directional_crossfeedvsMMMT5.pdf")
# ggsave2(p5_dir, filename = plot_file, device=cairo_pdf, height=16, width = 10)

# -------------------------
# Plot T9 Modules with directionality
# -------------------------
p6_dir <- module_counts_T9_dir %>%
  filter(description %in% top15_modules_T9) %>%
  mutate(description = str_wrap(description, width = 30)) %>%
  group_by(description) %>%
  mutate(total = sum(protein_count)) %>%
  ungroup() %>%
  mutate(protein_count_plot = ifelse(direction == "Downregulated", -protein_count, protein_count)) %>%
  ggplot(aes(x = reorder(description, total), y = protein_count_plot, fill = direction)) +
  geom_col() +
  coord_flip() +
  scale_fill_manual(values = c("Upregulated" = "#9F00BB", "Downregulated" = "#00BBB4")) +
  labs(
    title = "KEGG Modules (T9): Up vs Down Regulation",
    x = "Module",
    y = "Protein Count",
    fill = "Direction"
  ) +
  theme_minimal(base_size = 14) +
  theme_pub() +
  theme(
    axis.text.y = element_text(size = 15, face = "bold"),
    axis.text.x = element_text(size = 15, face = "bold"),
    legend.position = "bottom"
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5)

plot_file <- paste0("Transcriptomics/Analysis/250626run/Figures/250922/", Sys.Date(), "module_directional_crossfeedvsMMMT9.pdf")

# -------------------------
# Statistical Testing for PCA Groups
# -------------------------

# Assuming you already have 'vsd' object from your code
# Extract PCA data with multiple PCs
pca_data <- plotPCA(vsd, intgroup = c("condition", "Time", "Media"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Get full PCA results for more components
pca_full <- prcomp(t(assay(vsd)))
pca_scores <- as.data.frame(pca_full$x[, 1:5]) # First 5 PCs
pca_scores$sample <- rownames(pca_scores)

# Merge with sample metadata
pca_scores <- pca_scores %>%
  left_join(
    as.data.frame(colData(vsd)) %>%
      rownames_to_column("sample"),
    by = "sample"
  )

# -------------------------
# 1. PERMANOVA (most robust for multivariate data)
# -------------------------
library(vegan)

# Create distance matrix from first 5 PCs
pc_dist <- dist(pca_scores[, 1:5])

# Test effect of each variable
permanova_condition <- adonis2(pc_dist ~ condition,
  data = pca_scores,
  permutations = 999
)

permanova_time <- adonis2(pc_dist ~ Time,
  data = pca_scores,
  permutations = 999
)

permanova_media <- adonis2(pc_dist ~ Media,
  data = pca_scores,
  permutations = 999
)

# Test all variables together
permanova_full <- adonis2(pc_dist ~ Media + Time + Media:Time,
  data = pca_scores,
  permutations = 999
)
