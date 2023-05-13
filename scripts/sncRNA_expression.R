# Script maps mt-sncRNAs to F and M mitochondrial genomes of Potamilus streckersoni.
# Then, a heatmap is constructed to display differences in F and M mt-sncRNA
# expression across female and male gonadal samples.



# 0. Load libraries ------------------------------------------------------------

library(BiocManager)   # Bioconductor, manage packages
library(Rsubread)      # Bioconductor, featureCounts read summarization
library(tidyverse)     # CRAN, data manipulation
library(pheatmap)      # CRAN, clustered heat maps
library(ggplotify)     # CRAN, convert pheatmap to ggplot object
library(cowplot)       # CRAN, grid of plots
library(locfit)        # CRAN, used to make density plot for heat map legend
library(grDevices)     # CRAN, make gradient for heatmap



# 1. GTF Annotation Files for sncRNA centroids ----------------------------------
# Commented out b/c manual changes were made to finalize GTFs.
# Final GTF files located at ./sncRNA/mapped_reads_for_expression/ subdirectories


### Prepare data for female centroid GTF ###
# # Geneious text view provides list of centroids and their strand
# get_female_centroid_direction = read.delim('./sncRNA/before_filtering/female_before_filtering/RNA_hybrid_pass_F_with_strand_direction.tsv', header = F)
#
# # Grab centroids that are on reverse strand
# female_centroids_reverse_not_clean = grep('Cluster.+(reverse)', get_female_centroid_direction$V1, value = TRUE)
#
# # Clean names to merge with data set with centroid meta data
# female_centroids_reverse = gsub(' .+', '', female_centroids_reverse_not_clean)
# female_centroids_reverse = as.data.frame(female_centroids_reverse)
# colnames(female_centroids_reverse) = 'centroid'
#
# # Only 13 F centroids that passed RNAHybrid are reverse strand
# female_centroids_reverse = female_centroids_reverse %>% distinct(centroid)
#
# # Import files with filtered centroid meta data
# male_centroids = read.csv('./sncRNA/final_sncRNAs/M_mtsncRNA_filtered_final.csv')
# female_centroids = read.csv('./sncRNA/final_sncRNAs/F_mtsncRNA_filtered_final.csv')
#
# # Check if any of the filtered centroids are reverse (all forward)
# which(female_centroids$centroid %in% female_centroids_reverse)
#
# #Add direction column to female meta data
# female_centroids$direction = 'forward'
#
# ## Construct GTFs ###
# male_centroids_gtf = NULL
# male_scaffold = rep('PohiBra085_Mmt_4Raquel', nrow(male_centroids))
# male_centroids_gtf$scaffold = male_scaffold
# male_centroids_gtf = as.data.frame(male_centroids_gtf)
# male_centroids_gtf$annotation_source = 'ON881148'
# male_centroids_gtf$feature_type = 'gene' # featureCounts looks for gene in GTF
# male_centroids_gtf$start = male_centroids$temp_start
# male_centroids_gtf$stop = male_centroids$temp_stop
# male_centroids_gtf$score = '.'
# male_centroids_gtf$strand = male_centroids$direction
# male_centroids_gtf$strand = gsub('forward', '+', male_centroids_gtf$strand)
# male_centroids_gtf$strand = gsub('reverse', '-', male_centroids_gtf$strand)
# male_centroids_gtf$frame = '0'
# male_gene_id = paste0("'", male_centroids$Number, "';")
# male_centroids_gtf$attribute = paste('gene_id', male_centroids$Number, sep = ' ')
# write_delim(male_centroids_gtf, './sncRNA/mapped_reads_for_expression/M_mtsncRNAs/male_centroids_gtf.GTF', delim = '\t', col_names = F)
# # Note: went into VSCode after and added quotes around sncRNA number in attribute col and gff header
#
# female_centroids_gtf = NULL
# female_scaffold = rep('MW413895', nrow(female_centroids))
# female_centroids_gtf$scaffold = female_scaffold
# female_centroids_gtf = as.data.frame(female_centroids_gtf)
# female_centroids_gtf$annotation_source = 'MW413895'
# female_centroids_gtf$feature_type = 'gene' # featureCounts prob looks for genes?
# female_centroids_gtf$start = female_centroids$start
# female_centroids_gtf$stop = female_centroids$stop
# female_centroids_gtf$score = '.'
# female_centroids_gtf$strand = female_centroids$direction
# female_centroids_gtf$strand = gsub('forward', '+', female_centroids_gtf$strand)
# female_centroids_gtf$strand = gsub('reverse', '-', female_centroids_gtf$strand)
# female_centroids_gtf$frame = '0'
# female_centroids_gtf$attribute = paste('gene_id', female_centroids$centroid, sep = ' ')
# write_delim(female_centroids_gtf, './sncRNA/mapped_reads_for_expression/F_mtsncRNAs/female_centroids_gtf.GTF', delim = '\t', col_names = F)


# For reference, we manually removed rows with duplicate cluster names
# and rows that had same start/stop locations, and added quotes around
# Cluster in the attribute col and gff header



# 2. Map sncRNA reads -----------------------------------------------------------

F_sncRNAs = c(
  './sncRNA/mapped_reads_for_expression/F_mtsncRNAs/PohiBra084-Male_mapped.bam',
  './sncRNA/mapped_reads_for_expression/F_mtsncRNAs/PohiBra088-Female_mapped.bam',
  './sncRNA/mapped_reads_for_expression/F_mtsncRNAs/PohiBra089-Female_S1_L002_mapped.bam',
  './sncRNA/mapped_reads_for_expression/F_mtsncRNAs/PohiBra090-Male_S2_L002_mapped.bam',
  './sncRNA/mapped_reads_for_expression/F_mtsncRNAs/PohiBra091-Male_S3_L002_mapped.bam',
  './sncRNA/mapped_reads_for_expression/F_mtsncRNAs/PohiBra093-Female_S4_L002_mapped.bam'
  )

F_sncRNAs_counts = featureCounts(
  file = F_sncRNAs,
  annot.ext = './sncRNA/mapped_reads_for_expression/F_mtsncRNAs/female_centroids_gtf.gtf',
  isGTFAnnotationFile = T,
  GTF.featureType = 'gene',
  isPairedEnd = F,
  useMetaFeatures = T,
  verbose = T,
  nthreads = 7
  )

M_sncRNAs = c(
  './sncRNA/mapped_reads_for_expression/M_mtsncRNAs/PohiBra084-Male_mapped.bam',
  './sncRNA/mapped_reads_for_expression/M_mtsncRNAs/PohiBra088-Female_mapped.bam',
  './sncRNA/mapped_reads_for_expression/M_mtsncRNAs/PohiBra089-Female_S1_L002_mapped.bam',
  './sncRNA/mapped_reads_for_expression/M_mtsncRNAs/PohiBra090-Male_S2_L002_mapped.bam',
  './sncRNA/mapped_reads_for_expression/M_mtsncRNAs/PohiBra091-Male_S3_L002_mapped.bam',
  './sncRNA/mapped_reads_for_expression/M_mtsncRNAs/PohiBra093-Female_S4_L002_mapped.bam'
  )

M_sncRNAs_counts = featureCounts(
  file = M_sncRNAs,
  annot.ext = './sncRNA/mapped_reads_for_expression/M_mtsncRNAs/male_centroids_gtf.gtf',
  isGTFAnnotationFile = T,
  GTF.featureType = 'gene',
  isPairedEnd = F,
  useMetaFeatures = T,
  verbose = T,
  nthreads = 7
  )




# 3. Merge centroid names with final manuscript IDs ----------------------------

manuscript_sncRNA_IDs = read.csv(
  './sncRNA/final_sncRNAs/manuscript_sncRNA_IDs.csv',
  stringsAsFactors = T
  )

female_centroids = read.csv(
  './sncRNA/final_sncRNAs/F_mtsncRNA_filtered_final.csv',
  stringsAsFactors = T
  )

female_correct_IDs = female_centroids %>%
  select(centroid, sequence) %>%
  inner_join(
    manuscript_sncRNA_IDs,
    by = c('sequence' = 'Sequence')
  ) %>%
  distinct()

male_centroids = read.csv(
  './sncRNA/final_sncRNAs/M_mtsncRNA_filtered_final.csv',
  stringsAsFactors = T
  )

male_centroids$Centroid_Sequence = toupper(male_centroids$Centroid_Sequence)
male_correct_IDs = male_centroids %>%
  select(Number, Centroid_Sequence) %>%
  inner_join(
    manuscript_sncRNA_IDs,
    by = c('Centroid_Sequence' = 'Sequence')
  ) %>%
  distinct() %>%
  rename(
    'centroid' = 'Number',
    'sequence' = 'Centroid_Sequence'
  )

sncRNA_correct_IDs = rbind(female_correct_IDs, male_correct_IDs)


# Combined counts matrix
df.M_sncRNAs_counts = as.data.frame(M_sncRNAs_counts$counts)
df.M_sncRNAs_counts = rownames_to_column(df.M_sncRNAs_counts)
df.F_sncRNAs_counts = as.data.frame(F_sncRNAs_counts$counts)
df.F_sncRNAs_counts = rownames_to_column(df.F_sncRNAs_counts)
df.F_sncRNAs_counts = df.F_sncRNAs_counts %>%
  filter(rowname != 'Cluster1986_size_1015') # low counts across all samples

comb_sncRNAs_counts = rbind(df.M_sncRNAs_counts, df.F_sncRNAs_counts)
comb_sncRNAs_counts = inner_join(
  sncRNA_correct_IDs,
  comb_sncRNAs_counts,
  by = c('centroid' = 'rowname')
  ) %>%
  select(-centroid, -sequence, -Location)

comb_sncRNAs_counts = column_to_rownames(comb_sncRNAs_counts, var = 'Name')




# 4. Counts per Million --------------------------------------------------------

# Counts per million, output is a df for heatmap
CPM = function(count_matrix){
  norm_count_matrix = count_matrix
  for (i in 1:ncol(count_matrix)) {
    norm_count_matrix[,i] = (1000000 * count_matrix[,i]) / sum(count_matrix[,i])
  }
  return(norm_count_matrix)
}

# Counts per million, output is a vector for density
CPM_vector = function(count_matrix){
  norm_vector = c()
  for (i in 1:ncol(count_matrix)) {
    cpm_value = (1000000 * count_matrix[,i]) / sum(count_matrix[,i])
    norm_vector = c(norm_vector, cpm_value)
  }
  return(norm_vector)
}

# Apply CPM functions
cpm_comb_sncRNAs_vec = CPM_vector(comb_sncRNAs_counts)
cpm_comb_sncRNAs = CPM(comb_sncRNAs_counts)

# Log10(CPM + 1)
log_cpm_comb_sncRNAs_vec = log10(cpm_comb_sncRNAs_vec + 1)
log_cpm_comb_sncRNAs = log10(cpm_comb_sncRNAs + 1)



# 5. Summarize raw & CPM counts across female and male samples -----------------

# Summarize raw and CPM in female samples
females_summary_stats_cpm_sncRNAs = cpm_comb_sncRNAs %>%
  select(2, 3, 6) %>%
  pivot_longer(c(1, 2, 3), names_to = 'sample') %>%
  summarize(
    counts = 'CPM',
    sex = 'female',
    mean = mean(value),
    median = median(value),
    IQR = IQR(value)
  )

females_summary_stats_raw_sncRNAs = comb_sncRNAs_counts %>%
  select(2, 3, 6) %>%
  pivot_longer(c(1, 2, 3), names_to = 'sample') %>%
  summarize(
    counts = 'raw',
    sex = 'female',
    mean = mean(value),
    median = median(value),
    IQR = IQR(value)
  )

# Summarize raw and CPM in male samples
males_summary_stats_raw_sncRNAs = comb_sncRNAs_counts %>%
  select(1, 4, 5) %>%
  pivot_longer(c(1, 2, 3), names_to = 'sample') %>%
  summarize(
    counts = 'raw',
    sex = 'male',
    mean = mean(value),
    median = median(value),
    IQR = IQR(value)
  )

males_summary_stats_cpm_sncRNAs = cpm_comb_sncRNAs %>%
  select(1, 4, 5) %>%
  pivot_longer(c(1, 2, 3), names_to = 'sample') %>%
  summarize(
    counts = 'CPM',
    sex = 'male',
    mean = mean(value),
    median = median(value),
    IQR = IQR(value)
  )

# Combine into table
summary_stats_sncRNA_expression_by_sex = rbind(
  females_summary_stats_cpm_sncRNAs,
  females_summary_stats_raw_sncRNAs,
  males_summary_stats_cpm_sncRNAs,
  males_summary_stats_raw_sncRNAs
  )



# 6. sncRNA Heatmap Figure ------------------------------------------------------

### Legend for heat map ###
# Mean-center vector of log10 CPM to match `scale` in `pheatmap`
centered_log_cpm = scale(log_cpm_comb_sncRNAs_vec, scale = TRUE, center = TRUE)

# Convert to continuous density values
# (needed in order to fill same colors as heatmap)
centered_log_density = density(centered_log_cpm, n = 2^10)

# Heatmap colors
cb_palette = palette.colors(palette = "Okabe-Ito")
cb_heatmap_colors = c(cb_palette[7], 'white', cb_palette[4])
cb_heatmap_func = colorRampPalette(cb_heatmap_colors)

# Density plot for heat map legend
heatmap1_legend_plot = data.frame(
  x = centered_log_density$x,
  y = centered_log_density$y
  ) %>%
  ggplot(aes(x = x, y = y)) +
  geom_line(col = 'white') +
  geom_segment(
    aes(
      xend = x,
      yend = 0,
      color = x
      )
    ) +
  scale_color_gradient2(
    low = cb_heatmap_colors[1],
    mid = cb_heatmap_colors[2],
    high = cb_heatmap_colors[3]
    ) +
  scale_x_continuous(
    breaks = c(-1.5, 0, 1.5),
    position = 'top'
    ) +
  scale_y_reverse() +
  theme_minimal() +
  coord_flip() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = 'none',
    axis.text.x = element_blank(),
    axis.title.y = element_blank(),
    axis.title.x = element_blank()
    )

### Heat map ###
# Input needs to be a matrix
heatmap1_input = as.matrix(log_cpm_comb_sncRNAs)

# Change sample column names
column_labels = c(
  'PohiBra084',
  'PohiBra088',
  'PohiBra089',
  'PohiBra090',
  'PohiBra091',
  'PstrBra093'
  )

# Generate heat map with mean-centered values
heatmap1_plot = pheatmap(
  heatmap1_input,
  treeheight_row = 0,
  treeheight_col = 0,
  cutree_cols = 2,
  clustering_distance_rows = 'correlation',
  clustering_distance_cols = 'correlation',
  main = '',
  scale = 'row', # Scale so mean of each sncRNA is zero
  color = cb_heatmap_func(100),
  legend = F,
  labels_col = column_labels
  )

# Convert pheatmap to ggplot object
heatmap1_ggplot = as.ggplot(heatmap1_plot)

# White space
blank = ggplot() +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

# Combine plots as one figure
sncRNA_heatmap = ggdraw() +
  draw_plot(
    blank,
    x = 0,
    y = 0,
    width = 0.05
    ) +
  draw_plot(
    heatmap1_ggplot,
    x = 0.2,
    y = 0,
    width = 0.8
    ) +
  draw_plot(
    heatmap1_legend_plot,
    x = 0.05,
    y = 0,
    width = 0.15,
    height = 1
    ) +
  draw_plot_label(
    'Mean-Centered Log10 Counts Per Million',
    x = 0.05,
    y = 0.25,
    hjust = 0,
    vjust = 0,
    size = 10,
    angle = 90
    )

# Export figure
# ggsave('./figures/sncRNA_heatmap.pdf', width = 8, height = 9, units = 'in')
ggsave(plot = sncRNA_heatmap, './figures/sncRNA_heatmap.pdf', width = 5, height = 6, units = 'in')
