# Script filters putative mt-sncRNAs by base coverage.
# (i.e., valid mt-sncRNAs should have a sharp increase in 5' and 3' coverage)
# Additionally, F mt-sncRNAs are filtered by ΔΔG scores < −9 kJ from RNAup.



# 0. Load libraries and data ---------------------------------------------------
library(tidyverse)

# Import coverage stats
f_coverage = read.delim("./sncRNA/filter_by_coverage/coverage_data/F_base_cov.tsv", header = F)
f_coverage = f_coverage[, -1]
names(f_coverage) = c("position", "depth")

m_cov = read.delim("./sncRNA/filter_by_coverage/coverage_data/M_base_cov.tsv")
m_cov = m_cov[, 2:3]
names(m_cov) = c("position", "depth")



# 1. Clean the read coverage dfs -----------------------------------------------

# F base coverage did not have missing base positions / zero depth

# M read coverage, add in missing base positions that had zero depth
n_bases = 16484
M_mtDNA_size = data.frame(position = 1:n_bases)

# Assign zero depth to previously missing bases
m_coverage = left_join(M_mtDNA_size, m_cov)
missing_depth = which(is.na(m_coverage$depth))
m_coverage[missing_depth, 2] = 0

# Transform all depth + 1 to keep loop running
m_coverage$depth = m_coverage$depth + 1



# 2. Define functions to filter mt-sncRNAs -------------------------------------

# Identify sharp changes in forward and reverse coverage
# df_coverage : dataframe with col of bases (int), col of coverage at that base (int)
# threshold   : threshold value used to determine the depth cutoff (float)
# n_bps       : define the number of previous bps to calculate previous average depth (int)

get_sharp_coverage = function(df_coverage, threshold, n_bps) {

  ### Forward coverage ###
  df_coverage$fwd = rep(NA, each = nrow(df_coverage)) # Initialize with sharp increase in forward coverage = NA

  for (i in 6:nrow(df_coverage)) {  # Loop from the 6th bp to the last bp

    current_bp = df_coverage[i, 2]  # Track of the coverage of the current bp
    previous_bps_fwd = c()          # Track the coverage of the previous n_bps

    for (j in 1:n_bps) {
      # Append the previous n_bps
      previous_bps_fwd = c(previous_bps_fwd, df_coverage[i - j, 2])
    }

    # Is the coverage at the current bp at least a <threshold> increase from the average depth of prev n_bp?
    avg_depth = mean(previous_bps_fwd)  # Calculate average depth of previous n_bps
    depth_cutoff = avg_depth*threshold  # Calculate depth cutoff given threshold value

    if (current_bp >= depth_cutoff) {
      df_coverage$fwd[i] = TRUE
    } else {
      df_coverage$fwd[i] = FALSE
    }

  }


  ### Reverse coverage ###
  df_coverage$rev = rep(NA, each = nrow(df_coverage)) # Initialize with sharp increase in reverse coverage = NA
  reverse_loop_start = nrow(df_coverage) - 6          # Start of the for loop below

  for (i in reverse_loop_start:1) {

    current_bp  = df_coverage[i, 2]  # Track the coverage of the current bp
    previous_bps_rev = c()           # Track the coverage of the previous n_bps

    for (j in 1:n_bps) {
      # Append the previous n_bps
      previous_bps_rev = c(previous_bps_rev, df_coverage[i + j, 2])
    }

    # Is the coverage at the current bp at least a <threshold> increase from the average depth of prev n_bp?
    avg_depth = mean(previous_bps_rev)  # Calculate average depth of previous n_bps
    depth_cutoff = avg_depth*threshold  # Calculate depth cutoff given threshold value

    if (current_bp >= depth_cutoff) {
      df_coverage$rev[i] = TRUE
    } else {
      df_coverage$rev[i] = FALSE
    }

  }

  # Filter by bases with sharp coverage increases in forward/reverse direction
  sharp_coverage = df_coverage %>% filter(fwd == TRUE | rev == TRUE)

  # Return the data frame filtered above.
  return(sharp_coverage)

}



# 3. Filter F mt-sncRNAs -------------------------------------------------------
F_sharp_coverage = get_sharp_coverage(f_coverage, 1.5, 5)

# Filter down sncRNA targets by sharp 5' and 3' changes in coverage,
# where the start is leftmost base, stop is rightmost base
F_targets = read.csv("./sncRNA/filter_by_coverage/female_RNA_hybrid_pass_startstop.csv", header = T)

# Index sncRNA target df if start and stop positions have sharp coverage change
F_starts_targets = which(F_targets$start %in% F_sharp_coverage$position)
F_ends_targets = which(F_targets$stop %in% F_sharp_coverage$position)

F_filtered_targets_start = F_targets[F_starts_targets, ]
F_filtered_targets_ends = F_targets[F_ends_targets, ]

# Only keep sncRNAs that have sharp coverage change in start AND stop positions
F_mtsncRNA_filter_by_coverage = inner_join(
  F_filtered_targets_start,
  F_filtered_targets_ends
  )

# Export final female sncRNA targets
write.csv(
  F_mtsncRNA_filter_by_coverage,
  "./sncRNA/filter_by_coverage/female_centroids_1.5_coverage.csv",
  row.names = F
  )



# 4. Add gene meta data to the filtered F mt-sncRNAs --------------------------
# File contains <F_centroid_name>_<gene_target_ID>_<utr_ID>
female_centroid_info = read.delim("./sncRNA/filter_by_coverage/female_centroids_FUNID_utrs.tsv", header = F)

# Separate F centroid name, gene target ID, and utr ID into separate columns,
# then combine into a df.
centroid = gsub("__FUN.*", "", female_centroid_info$V1)
FUN_ID.utr = gsub("Cluster.*__", "", female_centroid_info$V1)
FUN_ID = gsub(".utr.*", "", FUN_ID.utr)
utr = gsub("FUN_[0-9]+.", "", FUN_ID.utr)
female_centroid_info_df = data.frame(
  centroid = centroid,
  FUN_ID = FUN_ID,
  utr = utr
  )

# Join F mt-sncRNA gene targets with F mt-sncRNA filtered by 1.5x coverage increase
targeted_by_F_mtsncRNA = inner_join(
  female_centroid_info_df,
  F_mtsncRNA_filter_by_coverage,
  by = "centroid"
  )

# Add column with the 3'UTR sequence of the genes targeted by the filtered F mt-sncRNAs
all_3UTRs = read.csv("./sncRNA/3UTRs_w_DNA_sequences.csv", header = TRUE)
all_3UTRs_FUNID = gsub("-T.*", "", all_3UTRs$ID)
all_3UTRs_utrs = gsub("FUN_.*-T[0-9].", "", all_3UTRs$ID)
all_3UTRs_split_IDs = data.frame(
  UTR_sequence = all_3UTRs$Sequence,
  FUN_ID = all_3UTRs_FUNID,
  utr = all_3UTRs_utrs
  )

female_centroids_1.5_meta_info = inner_join(
  targeted_by_F_mtsncRNA,
  all_3UTRs_split_IDs
  )

# Load gene annotations and convert to tidy df
annotations = read.delim("./sncRNA/filter_by_coverage/Pstreck_annotation.tsv", header = F)
genes_and_gene_info = as.data.frame(annotations[annotations$V3 == "gene", c(1, 9)])

# Split gff gene info column by semicolon to find genes with names
genes_IDs_names = cbind(
  genes_and_gene_info,
  str_split_fixed(
    genes_and_gene_info$V9,
    ";", 3
    )
  )
notclean_FUNID = genes_IDs_names$`1`
notclean_genename = genes_IDs_names$`2`
FUN_ID = gsub("ID=", "", notclean_FUNID)
Target_annotation = gsub("Name=", "", notclean_genename)
gene_annotations = data.frame(
  FUN_ID = FUN_ID,
  Target_annotation = Target_annotation,
  Target_Contig = genes_IDs_names$V1
  )

# Join annotation data to F mt-sncRNA and gene target df made above
female_centroids_1.5_meta_info = inner_join(
  female_centroids_1.5_meta_info,
  gene_annotations
  )

# Keep distinct rows only
female_centroids_1.5_meta_info = distinct(female_centroids_1.5_meta_info)

# If target is unannotated, annotate as 'Hypothetical'
female_centroids_1.5_meta_info[female_centroids_1.5_meta_info$Target_annotation == "", "Target_annotation"] = "Hypothetical Protein"

# Export
write.csv(female_centroids_1.5_meta_info, "./sncRNA/filter_by_coverage/female_centroids_1.5_meta_info.csv", row.names = F)



# 5. Add RNAup thermodynamics to female_centroids_1.5_meta_info ----------------

# Import thermodynamic results from RNAUp
F_RNAup = read.csv('./sncRNA/filter_by_coverage/female_centroids_1.5_meta_info_RNAup.csv', header = T)

# Keep interactions with free binding energy < -9 kJ
# Keep interactions that do not have a centroid with unknown nucleotide (N)
F_mtsncRNA_filtered_final = F_RNAup %>%
  filter(FreeBindingEnergy_RNAup_37C < -9 & Pass)

# Export final candidate F mt-sncRNAs
write.csv(
  F_mtsncRNA_filtered_final,
  './sncRNA/final_sncRNAs/F_mtsncRNA_filtered_final.csv',
  row.names = F
  )



# 6. Filter M mt-sncRNAs -------------------------------------------------------
M_sharp_coverage = get_sharp_coverage(m_coverage, 1.5, 5)

# Filter down sncRNA targets by sharp 5' and 3' changes in coverage
M_targets = read.csv("./sncRNA/filter_by_coverage/male_centroids_before_coverage_filter.csv", header = T)

# Make three new columns for loops below
M_targets$direction = M_targets$Start
M_targets$temp_start = M_targets$Start
M_targets$temp_stop = M_targets$Start

# `temp` columns have start and stop locations regardless of strand orientation
for (i in 1:nrow(M_targets)) {
  if (M_targets$Start[i] > M_targets$Stop[i]) {

    M_targets$direction[i] = "reverse"          # Notate reverse strand direction
    M_targets$temp_start[i] = M_targets$Stop[i] # temp_start column has start bp
    M_targets$temp_stop[i] = M_targets$Start[i] # temp_stop column has stop bp
  } else {

    M_targets$direction[i] = "forward"          # Notate forward strand direction

    # Start and stop are already fine, move to new column:
    M_targets$temp_start[i] = M_targets$Start[i]
    M_targets$temp_stop[i] = M_targets$Stop[i]
  }
}

# Index sncRNA target df if start and stop positions have sharp coverage change
M_starts_targets = which(M_targets$Start %in% M_sharp_coverage$position)
M_ends_targets = which(M_targets$Stop %in% M_sharp_coverage$position)

M_filtered_targets_start = M_targets[M_starts_targets, ]
M_filtered_targets_ends = M_targets[M_ends_targets, ]

# Only keep sncRNAs that have sharp coverage change in start AND stop positions
M_mtsncRNA_filtered_final = inner_join(
  M_filtered_targets_start,
  M_filtered_targets_ends
  )

# Export final male sncRNA targets
write.csv(
  M_mtsncRNA_filtered_final,
  "./sncRNA/final_sncRNAs/M_mtsncRNA_filtered_final.csv",
  row.names = F
  )

# Note: M mt-sncRNAs already had gene meta data and were already filtered by RNAup <-9 kJ.
