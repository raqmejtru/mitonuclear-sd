# Script contains the analysis workflow to investigate
# mitonuclear sex determination in freshwater mussels.



# ------------------------------------------------------------------------------
setwd("./mitonuclear-sd")


# Script creates data frames containing 3'UTRs hit by F and M mt-sncRNAs,
# respectively, along with the DNA sequences of the 3'UTRs.
source('./scripts/UTRs_for_RNAhybrid.R')


# Script finds the start and stop base pair positions of all potential
# F mt-sncRNAs that were mapped to the F mitogenome.
# Outputs a file needed in filter_by_coverage.R
source('./scripts/get_start_stop_of_F_mtsncRNAs_passed_RNAhybrid.R')


# Script filters putative mt-sncRNAs by base coverage.
# (i.e., valid mt-sncRNAs should have a sharp increase in 5' and 3' coverage)
# Additionally, F mt-sncRNAs are filtered by ΔΔG scores < −9 kJ from RNAup.
source('./scripts/filter_by_coverage.R')


# Script maps mt-sncRNAs to F and M mitochondrial genomes of Potamilus streckersoni.
# Then, a heatmap is constructed to display differences in F and M mt-sncRNA
# expression across female and male gonadal samples.
source('./scripts/sncRNA_expression.R')


# Script creates a table of genes that are targeted by the
# finalized list of F and M mt-sncRNAs.
source('./scripts/get_all_sncRNA_targets.R')


# Script performs differential gene expression analysis across
# female and male gonadal samples.
source('./scripts/DE.R')


# Script performs a weighted gene co-expression network analysis across
# female and male gonadal samples.
source('./scripts/WGCNA.R')


# Script plots the correlation between eigengenes and being female or male.
# Additionally, the intramodular connectivity of each gene is computed to
# identify hub genes.
# WGCNA results are exported, to be used in DE_downstream.R
source('./scripts/WGCNA_downstream.R')


# Script makes volcano plots of gene experession,
# adds mt-sncRNA target info to DESeq2 results,
# performs a gene set enrichment analysis,
# combines DESeq2, WGCNA, and mt-sncRNA results,
# and plots the expression of the GCNT1 protein network predicted by STRING.
source('./scripts/DE_downstream.R')
