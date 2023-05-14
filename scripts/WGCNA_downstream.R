# Script plots the correlation between eigengenes and being female or male.
# Additionally, the intramodular connectivity of each gene is computed to
# identify hub genes.
# WGCNA results are exported, to be used in DE_downstream.R



# 0. Load libraries ------------------------------------------------------------
library(tidyverse)     # CRAN, data manipulation
library(WGCNA)         # CRAN, weighted gene co-expression network analysis
library(grDevices)     # CRAN, make gradient for eigengene heatmap



# 1. Load WGCNA results --------------------------------------------------------
wgcna_results = read.csv(
  './RNASeq/WGCNA/WGCNA_outputs/genes_module_colors.csv',
  stringsAsFactors = F
  )

# 18 modules
wgcna_results$module = factor(wgcna_results$module)
levels(wgcna_results$module)

# Load network results
# (adjacency, MEs, moduleLabels, moduleColors, geneTree, countsGoodGenes, dissTOM)
load("./RNASeq/WGCNA/WGCNA_outputs/networkConstruction-stepByStep.RData")



# 2. Trait and eigengene network -----------------------------------------------

# Add M/F to eigengene matrix
sample_info = read.csv(
  "./RNASeq/sample_info.csv",
  header = T
  )

# Import M/F trait (yes = 1, no = 0)
male = as.data.frame(sample_info$male)
names(male) = "male"
female = as.data.frame(sample_info$female)
names(female) = "female"

# Add M/F to existing module eigengenes
MET = orderMEs(cbind(MEs, male, female))

# Plot dendrogram + heatmap of eigengene and trait network
pdf("./figures/plot_eigengene_sex_network.pdf", width = 8.5, height = 11)
cb_palette = palette.colors(palette = "Okabe-Ito")
cb_heatmap_colors = c(cb_palette[7], 'white', cb_palette[4])
cb_heatmap_func = colorRampPalette(cb_heatmap_colors)
colors = cb_heatmap_func(100)
plotEigengeneNetworks(
  multiME = MET,
  setLabels = "",
  plotDendrograms = T,
  plotHeatmaps = TRUE,
  setMargins = TRUE,
  marDendro = c(1, 5.5, 2, 3.75),
  marHeatmap = c(2, 6, 4, 2),
  colorLabels = TRUE,
  xLabelsAngle = 0,
  xLabelsAdj = 0.5,
  signed = TRUE,
  heatmapColors = colors,
  plotAdjacency = FALSE,         # plots correlation
  printAdjacency = FALSE,
  plotPreservation = "standard",
  zlimPreservation = c(0, 1),
  printPreservation = TRUE,
  showRows = c(2, 15),           # only show male and female
  xLabelsPosition = "top",       # move colors to top
  )
dev.off()



# 3. Merge gene annotations and module info ------------------------------------

# Prepare annotation data.frame
annotations         = read.delim("./RNASeq/Pstreck_annotation.tsv", header = F)
genes_and_gene_info = as.data.frame(annotations[annotations$V3 == "gene", c(1, 9)])
genes_IDs_names     = cbind(genes_and_gene_info, str_split_fixed(genes_and_gene_info$V9, ";", 3))
notclean_FUNID      = genes_IDs_names$`1`
notclean_genename   = genes_IDs_names$`2`
FUN_ID              = gsub("ID=", "", notclean_FUNID)
target_annotation   = gsub("Name=", "", notclean_genename)
gene_annotations    = data.frame(FUN_ID = FUN_ID,
                                 Target_annotation = target_annotation,
                                 Target_Contig = genes_IDs_names$V1)

# Append column to module results with available gene annotations
wgcna_results  = inner_join(
  wgcna_results,
  gene_annotations,
  by = c("FUN_ID" = "FUN_ID")
  )



# 4. Identify hub genes --------------------------------------------------------
# Compute intramodular connectivity
intramodular_connectivity = intramodularConnectivity(
  adjMat = adjacency,
  colors = moduleColors,
  scaleByMax = FALSE
  )

# Rownames to column
intramodular_connectivity = rownames_to_column(
  intramodular_connectivity,
  var = "row"
  )

# Add intramodular connectivity results to wgcna_results
wgcna_results = inner_join(
  intramodular_connectivity,
  wgcna_results,
  by = c('row' = 'FUN_ID')
  )

# Top 5% of kWithin are hub genes
five_pct = nrow(wgcna_results)*0.05
top_five_pct_kWithin = wgcna_results %>%
  arrange(-kWithin) %>%
  dplyr::slice(0:five_pct)

# Export
write.csv(
  top_five_pct_kWithin,
  "./RNASeq/WGCNA/WGCNA_downstream/top_five_pct_kWithin.csv",
  row.names = F
  )


# Identify hub gene from each module
hubs = chooseTopHubInEachModule(
  datExpr = countsGoodGenes,
  colorh = moduleColors,
  omitColors = "grey",
  power = 20,
  type = "signed"
  )

hubs = as.data.frame(hubs)
hubs = rownames_to_column(hubs, var = "module")

# Add WGCNA results to hub genes
module_hub_gene = left_join(
  hubs, wgcna_results,
  by = c('hubs' = 'row', 'module' = 'module')
  )
module_hub_gene = rename(module_hub_gene, 'FUN_ID' = 'hubs')

# Export hub genes data
write.csv(
  module_hub_gene,
  "./RNASeq/WGCNA/WGCNA_downstream/module_hub_gene.csv",
  row.names = F
  )



# 5. Export WGCNA results ------------------------------------------------------
write.csv(
  wgcna_results,
  "./RNASeq/WGCNA/WGCNA_downstream/wgcna_results.csv",
  row.names = F
  )
