# --- Load libraries ---
library(data.table)
library(dplyr)
library(tidyr)
library(tibble)
library(umap)
library(ggplot2)
library(irlba)

# --- Set working directory ---
setwd("C:/Users/skyun/Desktop/Research/25.Kronos_analysis/03.Capture_analysis/02.GATK_dimension_reduction")

# --- Load genotype data ---
df <- fread("common_EMS_True_per0.05_20837_sites_69.65_mutants.tsv") #ems exom
df <- fread("common_EMS_True_per0.15_56926_sites_208.95_mutants.tsv") #ems promoter
df <- fread("common_EMS_False_per0.05_42518_sites_69.65_mutants.tsv") #nonems exom
df <- fread("common_EMS_False_per0.15_93098_sites_208.95_mutants.tsv") #nonems promoter

# --- Convert genotypes to numeric codes ---
df <- df %>%
  mutate(gt_code = case_when(
    genotype %in% c("0/0", "0|0") ~ 0,
    genotype %in% c("0/1", "1/0", "0|1", "1|0") ~ 1,
    genotype %in% c("1/1", "1|1") ~ 2,
    TRUE ~ NA_real_
  )) %>%
  filter(!is.na(gt_code))

# --- Collapse duplicate sample-site entries by max zygosity ---
df_collapsed <- df %>%
  group_by(sample, site) %>%
  summarise(gt_code = max(gt_code), .groups = "drop")

# --- Pivot to wide matrix ---
geno_matrix <- df_collapsed %>%
  pivot_wider(names_from = site, values_from = gt_code, values_fill = 0)

# --- Convert to numeric matrix with sample as rownames ---
geno_matrix_num <- geno_matrix %>%
  column_to_rownames("sample") %>%
  as.matrix()

storage.mode(geno_matrix_num) <- "numeric"

# --- PCA + UMAP ---
pca_out <- prcomp_irlba(geno_matrix_num, n = 50, center = TRUE)
umap_out <- umap(pca_out$x[, 1:50])

# --- Prepare UMAP dataframe ---
umap_df <- as.data.frame(umap_out$layout)
colnames(umap_df) <- c("UMAP1", "UMAP2")
umap_df$Sample <- rownames(geno_matrix_num)


# --- Visualize without labels ---
ggplot(umap_df, aes(x = UMAP1, y = UMAP2)) +
  geom_point(alpha = 0.8) +
  theme_minimal()

#
# --- Run this part for EC datasets to get labels ---
#
library(dbscan)
library(ggsci)
# Cluster using DBSCAN on non-EMS exome UMAP data
non_ems_exome_df <- umap_df  # Replace with correct filtered dataset if needed
dbscan_result <- dbscan(non_ems_exome_df[, c("UMAP1", "UMAP2")], eps = 1.0, minPts = 5)

non_ems_exome_df$Cluster <- as.factor(dbscan_result$cluster)

# Assign colors to clusters
clusters <- sort(unique(non_ems_exome_df$Cluster))
n_clusters <- length(clusters)

cluster_colors <- setNames(pal_d3("category20")(n_clusters), as.character(clusters))

cluster_lookup <- non_ems_exome_df %>%
  select(Sample, Cluster) %>%
  distinct()

# Create a dataframe with sample ID, cluster, and assigned color
cluster_annotation <- non_ems_exome_df %>%
  dplyr::select(Sample, Cluster) %>%
  mutate(Color = cluster_colors[as.character(Cluster)])

# Write to file
write.table(cluster_annotation,
            file = "mutant_cluster_colors.tsv",
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)

####### Or load the colors
cluster <- read.csv(file = "mutant_cluster_colors.tsv", sep = "\t")

cluster_colors <- cluster %>% 
  distinct(Cluster, Color) %>%   # drop duplicates
  deframe()     

cluster_lookup <- cluster %>%         
  select(Sample, Cluster,) %>%  
  distinct()       

####### Merge color using the commends below
umap_df <- left_join(umap_df, cluster_lookup, by = "Sample")

ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = factor(Cluster))) +
  geom_point(alpha = 0.8) +
  scale_color_manual(values = cluster_colors) +
  theme_minimal()
