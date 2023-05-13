# Script finds the start and stop base pair positions of all potential
# F mt-sncRNAs that were mapped to the F mitogenome.
# Outputs a file needed in filter_by_coverage.R



# 0. Load libraries and data ---------------------------------------------------
library(tidyverse)

# Mapped reads -> sam file -> converted to tsv
F_mapped_sam = read.delim(
  "./sncRNA/before_filtering/female_before_filtering/RNA_hybrid_pass_F_mapped.tsv",
  header = F
  )

# 1. Get start/stop base pair positions ----------------------------------------

# Each element is potential F mt-sncRNA
# Each element has a list of nucleotides corresponding to sncRNA DNA sequence
nucleotides = strsplit(F_mapped_sam[,10], split = "")

# Get the lengths of each potential F mt-sncRNA DNA sequence
lengths = c()
for (i in 1:length(nucleotides)) {
  lengths[i] = length(nucleotides[[i]])
}

# Get vector of the left-most position of each F mt-sncRNA in the F mitogenome
start = F_mapped_sam[,4]

# Get vector of the right-most position of each F mt-sncRNA in the F mitogenome
stop = F_mapped_sam[,4] + lengths - 1

# Construct and export file used in filter_by_coverage.R
female_RNA_hybrid_pass_startstop = data.frame(
  start = start,
  stop = stop,
  centroid = F_mapped_sam$V1,
  sequence = F_mapped_sam$V10
  )

write.csv(
  female_RNA_hybrid_pass_startstop,
  "./sncRNA/filter_by_coverage/female_RNA_hybrid_pass_startstop.csv",
  row.names = F
  )
