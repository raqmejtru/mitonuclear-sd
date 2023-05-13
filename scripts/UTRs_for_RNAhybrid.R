# Script creates data frames containing 3'UTRs hit by F and M mt-sncRNAs,
# respectively, along with the DNA sequences of the 3'UTRs.

# 0. Load libraries and data ---------------------------------------------------
library(tidyverse)

F_centroids_BLAST = read.csv('./sncRNA/before_filtering/female_before_filtering/F_centroids_BLAST.csv')
M_centroids_BLAST = read.csv('./sncRNA/before_filtering/male_before_filtering/M_centroids_BLAST.csv')
UTRsequences = read.csv('./sncRNA/3UTRs_w_DNA_sequences.csv')


# 1. Clean 3'UTRs  -------------------------------------------------------------

# Split into 3 columns (Sequence, FUN_ID, utr)
UTRsequences.splitID = as.data.frame(
  separate(
    UTRsequences,
    col = 2,
    into = c('FUN', 'utr'),
    sep = '-T[0-9]'
    )
  )

# Merge FUN_ID and utr into a col (removed all variations of mRNA track `-T#`)
UTRsequences$ID.clean = paste0(
  UTRsequences.splitID$FUN,
  UTRsequences.splitID$utr
  )

# 12,762 unique 3'UTRs extracted from annotation
UTRsequences.unique = UTRsequences %>%
  select(1, 3) %>%
  distinct(ID.clean, .keep_all = T)



# 2. Make DF of 3'UTRs hit by F centroids, add DNA sequences -------------------

F_clusters_with_hits = F_centroids_BLAST %>%
  select(
    Document.Name,
    Query.start,
    Query.end,
    ID
  ) %>%
  filter(Query.start < 4) %>%
  distinct(Document.Name)

F_UTR_ID_clean = F_centroids_BLAST %>%
  select(
    Document.Name,
    Type,
    ID
  ) %>%
  filter(Type == "3'UTR") %>%
  separate(
    col = 3,
    into = c('FUN', 'utr'),
    sep = '-T[0-9]'
  ) %>%
  unite(
    'ID.clean',
    3:4,
    sep = ''
  ) %>%
  select(Document.Name, ID.clean) %>%
  distinct(ID.clean, .keep_all = T)

#Female clusters hit to 3,380 unique 3' UTRs with correct query starts
F_hits_and_unique_UTR_IDs = inner_join(
  F_clusters_with_hits,
  F_UTR_ID_clean,
  by = 'Document.Name',
  keep = TRUE
  ) %>%
  select(Document.Name.x, ID.clean) %>%
  rename(F_cluster_ID = Document.Name.x)

F_hits_and_UTR_sequences = inner_join(
  F_hits_and_unique_UTR_IDs,
  UTRsequences.unique,
  by = 'ID.clean'
  )

write.csv(
  F_hits_and_UTR_sequences,
  './sncRNA/before_filtering/female_before_filtering/F_hits_and_UTR_sequences.csv',
  row.names = F
  )


# 3. Make DF of 3'UTRs hit by M centroids, add DNA sequences -------------------
M_clusters_with_hits = M_centroids_BLAST %>%
  select(
    Document.Name,
    Query.start,
    Query.end,
    ID
  ) %>%
  filter(Query.start < 4) %>%
  distinct(Document.Name)


M_UTR_ID_clean = M_centroids_BLAST %>%
  select(
    Document.Name,
    Type,
    ID
  ) %>%
  filter(Type == "3'UTR'") %>%
  separate(
    col = 3,
    into = c('FUN', 'utr'),
    sep = '-T[0-9]'
  ) %>%
  unite('ID.clean', 3:4, sep = '') %>%
  select(Document.Name, ID.clean) %>%
  distinct(ID.clean, .keep_all = T)

#Male clusters hit to 5,233 unique 3' UTRs with correct query starts
M_hits_and_unique_UTR_IDs = inner_join(
  M_clusters_with_hits,
  M_UTR_ID_clean,
  by = 'Document.Name',
  keep = TRUE
  ) %>%
  select(Document.Name.x, ID.clean) %>%
  rename(M_cluster_ID = Document.Name.x)

M_hits_and_UTR_sequences = inner_join(
  M_hits_and_unique_UTR_IDs,
  UTRsequences.unique,
  by = 'ID.clean'
  )

write.csv(
  M_hits_and_UTR_sequences,
  './sncRNA/before_filtering/male_before_filtering/M_hits_and_UTR_sequences.csv',
  row.names = F
  )
