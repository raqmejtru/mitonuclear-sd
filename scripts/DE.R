# Script performs differential gene expression analysis across
# female and male gonadal samples.



# 0. Load libraries ------------------------------------------------------------

library(BiocManager)   # Bioconductor, manage packages
library(Rsubread)      # Bioconductor, featureCounts read summarization
library(DESeq2)        # Bioconductor, DGE and normalization
library(tidyverse)     # CRAN, data manipulation
library(factoextra)    # CRAN, PCA
library(pheatmap)      # CRAN, clustered heat maps



# 1. Map transcripts -----------------------------------------------------------

# Make a list of sorted bam file names to copy over to bam_files list below
# for file in *.bam; do echo "'$file'," >> bam_files.txt; done

# Note: bamFiles order should match sample order in sample_info.csv
bam_files = c(
  './RNASeq/PohiBra080_sorted.bam',
  './RNASeq/PohiBra081_sorted.bam',
  './RNASeq/PohiBra083_sorted.bam',
  './RNASeq/PohiBra084_sorted.bam',
  './RNASeq/PohiBra085_sorted.bam',
  './RNASeq/PohiBra086_sorted.bam',
  # 'PohiBra090_sorted.bam', # outlier, removed from DESeq analysis
  './RNASeq/PstrBra093_sorted.bam',
  './RNASeq/PstrBra094_sorted.bam',
  './RNASeq/PstrBra098_sorted.bam'
  )

# Summarize gene counts
featureCounts_output = featureCounts(
  file = bam_files,
  annot.ext = "./RNASeq/Pstreck_annotation.gtf",
  isGTFAnnotationFile = T,
  GTF.attrType = "gene_id",
  isPairedEnd = T,
  useMetaFeatures = T,
  verbose = T,
  nthreads = 64
  )



# 2. Make DESeq dataset --------------------------------------------------------

# Import list of protein coding genes
keep_pcg = as.list(
  read.table(
    "./RNASeq/keep_protein_coding_genes",
    header = T
    )
  )

# 43,193 genes total in raw counts matrix
raw_counts = featureCounts_output$counts
raw_counts %>% nrow()

# Confirm subset of 34,937 protein coding genes
keep_pcg$ID %in% raw_counts$X %>% length()

# Subset raw counts to keep protein coding genes only
pcg_counts = raw_counts %>% filter(X %in% keep_pcg$ID)

# Export pcg counts matrix as csv
write.csv(
  pcg_counts,
  './RNASeq/pcg_counts.csv',
  row.names = F
  )

pcg_counts = read.csv(
  "./RNASeq/pcg_counts.csv",
  header = T
  ) %>%
  column_to_rownames(var = "X") %>%
  as.matrix()


# Note: manually removed sample 90 from sample_info.csv
# Sample condition as factors
sample_info = read.csv(
  "./RNASeq/sample_info.csv",
  header = T
  )

sample_info$condition = factor(sample_info$condition)
# Levels: 0 female, 1 male

# Make DESeq dataset
deseq_input =  DESeqDataSetFromMatrix(
  countData = pcg_counts,
  colData = sample_info,
  design = ~ condition
  )



# 3. Normalize counts ----------------------------------------------------------
# For outlier detection only. Raw counts input for DESeq2.
vst_norm_counts = vst(deseq_input, blind = FALSE)



# 4. Identify sample outliers --------------------------------------------------

# Outliers based on sample distances:
# Distance object, lower tri matrix
lower_tri_distances = dist(t(assay(vst_norm_counts)))

# Distance matrix (rows are samples, cols are genes)
sample_distances = as.matrix(dist(t(assay(vst_norm_counts))))

# Auto exports
pheatmap(
  sample_distances,
  clustering_distance_rows = lower_tri_distances,
  clustering_distance_cols = lower_tri_distances,
  filename = "./RNASeq/DE/DE_outputs/outlier_heatmap.pdf"
  )

# Outliers based on PCA:
pca_results = prcomp(t(assay(vst_norm_counts)))

outlier_pca = fviz_pca_ind(
  pca_results,
  col.ind = "cos2", # Color by quality of representation
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE
  )

# Export PCA plot
ggsave(plot = outlier_pca, "./RNASeq/DE/DE_outputs/outlier_pca.pdf", outlier_pca)



# 5. Differential Gene Expression ----------------------------------------------

# Run DESeq analysis
deseq_output = DESeq(deseq_input)
deseq_results = results(deseq_output, tidy = T)

# Prepare annotation data.frame
annotations         = read.delim("./RNASeq/Pstreck_annotation.tsv", header = F)
genes_and_gene_info = as.data.frame(annotations[annotations$V3 == "gene", c(1, 9)])
genes_IDs_names     = cbind(genes_and_gene_info, str_split_fixed(genes_and_gene_info$V9, ";", 3))
notclean_FUNID      = genes_IDs_names$`1`
notclean_genename   = genes_IDs_names$`2`
FUN_ID              = gsub("ID=", "", notclean_FUNID)
Target_annotation   = gsub("Name=", "", notclean_genename)
gene_annotations    = data.frame(FUN_ID = FUN_ID,
                                 Target_annotation = Target_annotation,
                                 Target_Contig = genes_IDs_names$V1)

# Append column to DESeq results with available gene annotations
deseq_results  = inner_join(deseq_results, gene_annotations, by = c("row" = "FUN_ID"))

# Export DESeq results
write.csv(
  deseq_results,
  file = "./RNASeq/DE/DE_outputs/nuclear_deseq_results.csv",
  row.names = F
  )
