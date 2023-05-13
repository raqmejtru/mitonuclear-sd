# Script creates a table of genes that are targeted by the
# finalized list of F and M mt-sncRNAs.



# 0. Load libraries and data ---------------------------------------------------
library(tidyverse)

# Final set of F mt-sncRNAs
final_F_sncRNAs = read.csv(
  './sncRNA/final_sncRNAs/F_mtsncRNA_filtered_final.csv',
  header = T
  )

# Final set of M mt-sncRNAs
final_M_sncRNAs = read.csv(
  './sncRNA/final_sncRNAs/M_mtsncRNA_filtered_final.csv',
  header = T
  )



# 1. Table of genes targeted by final F and M mt-sncRNAs -----------------------

# F mt-sncRNA targets
targets_F_sncRNAs = final_F_sncRNAs %>%
  # Filter out these two centroids due to low counts across all samples
  filter(centroid != 'Cluster1986_size_1015') %>%
  select(Target, centroid) %>%
  mutate(sncRNA_type = 'F-sncRNA')


# M mt-sncRNA targets
targets_M_sncRNAs = final_M_sncRNAs %>%
  # Replicates of M-sncRNA-2
  filter(Number != '3' & Number != '7') %>%
  select(Target, Number) %>%
  mutate(sncRNA_type = 'M-sncRNA')

# Both F and M mt-sncRNA targets
all_sncRNA_targets = full_join(
  targets_F_sncRNAs,
  targets_M_sncRNAs,
  by = c(
    'Target' = 'Target',
    'sncRNA_type' = 'sncRNA_type',
    'centroid' = 'Number'
    )
  )

# Export to csv
write.csv(
  all_sncRNA_targets,
  './sncRNA/final_sncRNAs/all_sncRNA_targets.csv',
  row.names = F
  )
