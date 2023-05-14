# Script performs a weighted gene co-expression network analysis across
# female and male gonadal samples.



# 0. Load libraries ------------------------------------------------------------
library(BiocManager)  # Bioconductor, manages packages
library(Rsubread)     # Bioconductor, featureCounts read summarization
library(DESeq2)       # Bioconductor, normalize counts
library(NOISeq)       # Bioconductor, filter low counts
library(WGCNA)        # CRAN, weighted gene co-expression network analysis
library(tidyverse)    # CRAN, data manipulation



# 1. Summarize counts ----------------------------------------------------------

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



# 2. Normalize counts ----------------------------------------------------------

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

# Import sample info
sample_info = read.csv("./RNASeq/sample_info.csv", header = T)

# Filter out genes that have an average expression per condition
# less than 5 CPM and a coefficient of variation per condition higher
# than cv.cutoff (200%) in all conditions.
filtered_counts = filtered.data(
  pcg_counts,
  factor = sample_info$condition,
  norm = FALSE,
  method = 1,
  cv.cutoff = 200,
  cpm = 5
  )

# Normalize counts
countsNorm = vst(filtered_counts)



# 3. Outlier clustering --------------------------------------------------------

countsNorm.t = t(countsNorm)

# Clean data to only include "good" genes & samples
gsg = goodSamplesGenes(countsNorm.t, verbose = 3);
pctGoodSamples = sum(gsg$goodSamples)/length(gsg$goodSamples)
print(paste0('Proportion of good samples: ', pctGoodSamples*100, "%"))
pctGoodGenes = sum(gsg$goodGenes)/length(gsg$goodGenes)
print(paste0('Proportion of good genes: ', pctGoodGenes*100, "%"))

# Filter by good genes only
countsGoodGenes = countsNorm.t[, gsg$goodGenes]

# Export sample clustering to detect outliers
sampleTree = hclust(
  dist(countsGoodGenes),
  method = "average"
  )

pdf("./RNASeq/WGCNA/WGCNA_outputs/3.sample_clustering_outliers.pdf")
plotSampleTree = plot(
  sampleTree,
  main = "Sample clustering to detect outliers",
  sub = "", xlab = "",
  cex.lab = 1.5,
  cex.axis = 1.5,
  cex.main = 2
  )
dev.off()



# 4. Soft threshold ------------------------------------------------------------

# Define vector with sequence of powers
powers = c(c(1:10), seq(from = 12, to = 30, by = 2))

# Choose soft threshold
sft = pickSoftThreshold(
  countsGoodGenes,
  powerVector = powers,
  networkType = "signed",
  verbose = 5
  )


# Scale-free topology fit index as a function of the soft-thresholding power
pdf("./RNASeq/WGCNA/WGCNA_outputs/4.scale_free_topology_model_fit.pdf")
plot(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  type = "n",
  main = paste("Scale independence"))
text(
  sft$fitIndices[,1],
  -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
  labels = powers,
  cex = 0.9,
  col = "red"
  )
abline(h = 0.90, col = "red") # R^2 cut-off
dev.off()


# Mean connectivity as a function of the soft-thresholding power
pdf("./RNASeq/WGCNA/WGCNA_outputs/4.mean_connectivity.pdf")
plot(
  sft$fitIndices[,1],
  sft$fitIndices[,5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = paste("Mean connectivity"))
text(
  sft$fitIndices[,1],
  sft$fitIndices[,5],
  labels = powers,
  cex = 0.9,
  col = "red"
  )
dev.off()



# 5. Step-by-step network construction -----------------------------------------

# 5a. Distances for clustering -------------------------------------------------
# Network adjacency matrix (required for `intramodularConnectivity`)
softPower = 20
adjacency = adjacency(
  countsGoodGenes,
  type = "signed",
  power = softPower
  )

# Topological Overlap Matrix (TOM)
# To minimize effects of noise and spurious associations,
# transform the adjacency into Topological Overlap Matrix,
# and calculate the corresponding dissimilarity.
TOM = TOMsimilarityFromExpr(
  datExpr = countsGoodGenes,
  weights = NULL,
  corType = "pearson",
  networkType = "signed",
  power = softPower,
  TOMType = "signed",
  TOMDenom = "mean",
  maxPOutliers = 1,
  quickCor = 0,
  pearsonFallback = "individual",
  cosineCorrelation = FALSE,
  replaceMissingAdjacencies = FALSE,
  suppressTOMForZeroAdjacencies = FALSE,
  suppressNegativeTOM = FALSE,
  useInternalMatrixAlgebra = FALSE,
  nThreads = 64,
  verbose = 1
  )

# Dissimilarity from TOM
dissTOM = 1 - TOM

# Hierarchical Gene Clustering using dissTOM
geneTree = hclust(
  as.dist(dissTOM),
  method = "average"
  )

# Export dendrogram
pdf("./RNASeq/WGCNA/WGCNA_outputs/5a.dendrogram.pdf")
plot(
  geneTree,
  xlab = "",
  sub = "",
  main = "Gene clustering on TOM-based dissimilarity",
  labels = FALSE,
  hang = 0.04)
dev.off()



# 5b. Module Identification (Cut Branches) -------------------------------------
#  - Table returns modules listed largest to smallest.
#  - Label 0 is reserved for unassigned genes.
minModuleSize = 30
dynamicMods = cutreeDynamic(
  dendro = geneTree,
  distM = dissTOM,
  pamRespectsDendro = FALSE,
  minClusterSize = minModuleSize
  )

dynamicModsTable = table(dynamicMods)

# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
dynamicColorsTable = table(dynamicColors)



# 5c. Module dendrogram --------------------------------------------------------
pdf("./RNASeq/WGCNA/WGCNA_outputs/5c.dynamic_cut_dendrogram.pdf")
plotDendroAndColors(
  geneTree,
  dynamicColors,
  "Dynamic Tree Cut",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene dendrogram and module colors"
  )
dev.off()



# 5d. Merge modules with similar expression profiles ---------------------------

# Calculate eigengenes & dissimilarity of module eigengenes
MEList = moduleEigengenes(countsGoodGenes, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1 - cor(MEs)


# Cluster module eigengenes
METree = hclust(
  as.dist(MEDiss),
  method = "average"
  )

# Export eigengene dendrogram
pdf("./RNASeq/WGCNA/WGCNA_outputs/5d.eigengene_dendrogram.pdf")
plot(
  METree,
  main = "Clustering of Module Eigengenes",
  xlab = "",
  sub = ""
  )
dev.off()



# 5e. Merged dendrogram --------------------------------------------------------

# Height cut is 1 - correlation
MEDissThres = 0.3

merge = mergeCloseModules(
  countsGoodGenes,
  dynamicColors,
  cutHeight = MEDissThres,
  verbose = 3
  )

mergedColors = merge$colors
mergedMEs = merge$newMEs

# Export merged dendrogram
pdf("./RNASeq/WGCNA/WGCNA_outputs/5e.merged_dendrogram_stepByStep.pdf")
plotDendroAndColors(
  geneTree,
  cbind(dynamicColors, mergedColors),
  c("Dynamic Tree Cut", "Merged dynamic"),
  main = "Merged Cluster Dendrogram",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05
  )
dev.off()


# Save merged variables for future analysis
moduleColors = mergedColors

# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder) - 1
MEs = mergedMEs

# Save module colors and labels for use in WGCNA downstream analyses
save(
  adjacency,
  MEs,
  moduleLabels,
  moduleColors,
  geneTree,
  countsGoodGenes,
  dissTOM,
  file = "./RNASeq/WGCNA/WGCNA_outputs/networkConstruction-stepByStep.RData"
  )



# 6. Network heatmap -----------------------------------------------------------
#  - Each row and column of the heatmap correspond to a single gene.
#  - The heatmap can depict adjacency or topological overlap,
#    with light colors denoting low adjacency and darker colors higher adjacency.
#  - In addition, the gene dendrograms and module colors are plotted along the top
#    and left side of the heatmap.

# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^12

# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA

# Export heatmap showing connections between all genes
pdf("./RNASeq/WGCNA/WGCNA_outputs/6.allGene_heatmap.pdf")
TOMplot(
  plotTOM,
  geneTree,
  moduleColors,
  main = "Network heatmap plot, all genes"
  )
dev.off()



# 7. Trait + eigengenes --------------------------------------------------------
#  - Use eigengenes as representative profiles and quantify module similarity
#    by eigengene correlation.
#  - 'plotEigengeneNetworks' generates a summary plot of the eigengene network.
#  - It is usually informative to add a trait (or multiple traits) to the eigengenes
#    to see how the traits fit into the eigengene network

sample_info = read.csv("./RNASeq/sample_info.csv", header = T)

# Import M/F trait (yes = 1, no = 0)
male = as.data.frame(sample_info$male)
names(male) = "male"
female = as.data.frame(sample_info$female)
names(female) = "female"

# Add M/F to existing module eigengenes
MET = orderMEs(cbind(MEs, male, female))

# Plot the relationships among the eigengenes and the trait
pdf("./RNASeq/WGCNA/WGCNA_outputs/7.eigengeneTrait.pdf")
plotEigengeneNetworks(
  MET,
  "",
  marDendro = c(0,4,1,2),
  marHeatmap = c(3,4,1,2),
  xLabelsAngle = 90
  )
dev.off()

# Plot the dendrogram
pdf("./RNASeq/WGCNA/WGCNA_outputs/7.eigengeneTraitDendrogramOnly.pdf")
plotEigengeneNetworks(
  MET, "Eigengene dendrogram",
  marDendro = c(0,4,2,0),
  plotHeatmaps = FALSE)
dev.off()

# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
pdf("./RNASeq/WGCNA/WGCNA_outputs/7.heatmap_eigengenes_sex")
plotEigengeneNetworks(
  MET,
  "Eigengene adjacency heatmap",
  marHeatmap = c(3,4,2,2),
  plotDendrograms = FALSE,
  xLabelsAngle = 90
  )
dev.off()



# 8. Export gene + module information ------------------------------------------
genesWithModuleColors = data.frame(
  FUN_ID = colnames(countsGoodGenes),
  module = moduleColors
  )

write.csv(
  genesWithModuleColors,
  "./RNASeq/WGCNA/WGCNA_outputs/genes_module_colors.csv",
  row.names = F
  )
