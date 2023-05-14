# Script makes volcano plots of gene experession,
# adds mt-sncRNA target info to DESeq2 results,
# performs a gene set enrichment analysis,
# combines DESeq2, WGCNA, and mt-sncRNA results,
# and plots the expression of the GCNT1 protein network predicted by STRING.



# 0. Load libraries ------------------------------------------------------------
library(tidyverse)     # CRAN, data manipulation
library(fgsea)         # Bioconductor, Fast Gene Set Enrichment Analysis
library(igraph)        # CRAN, used by tidygraph
library(ggraph)        # CRAN, ggplot implementation of networks
library(tidygraph)     # CRAN, add tidy node and edge data to network



# 1. Add mt-sncRNA target info to DESeq2 results -------------------------------

# Import DESeq2 results
nuclear_results = read.csv(
  './RNASeq/DE/DE_outputs/nuclear_deseq_results.csv',
  stringsAsFactors = F
  )

# Import all mt-sncRNA target genes
all_sncRNA_targets = read.csv(
  './sncRNA/final_sncRNAs/all_sncRNA_targets.csv',
  header = T
  )

# Index rows that are F mt-sncRNA targets
female_targets = all_sncRNA_targets %>% filter(sncRNA_type == 'F-sncRNA')
index_sncRNA_F_targets = which(nuclear_results$row %in% female_targets$Target)

# Index rows that are M mt-sncRNA targets
male_targets = all_sncRNA_targets %>% filter(sncRNA_type == 'M-sncRNA')
index_sncRNA_M_targets = which(nuclear_results$row %in% male_targets$Target)

# New data.frame with only sncRNA target genes
sncRNA_targets = nuclear_results %>% add_column(sncRNA_type = 'not a target')
sncRNA_targets[index_sncRNA_F_targets, 10] = 'F-sncRNA'
sncRNA_targets[index_sncRNA_M_targets, 10] = 'M-sncRNA'
sncRNA_targets_deseq_results = sncRNA_targets %>% filter(sncRNA_type != 'not a target')

# Add columns with sncRNA ID and type to targets
target_results = inner_join(
  sncRNA_targets_deseq_results,
  all_sncRNA_targets,
  by = c(
    'row' = 'Target',
    'sncRNA_type' = 'sncRNA_type'
    )
  )

# Remove duplicates
target_results = target_results %>% distinct()



# 2. Volcano plots -------------------------------------------------------------

# Expression of all genes
volcano_nuclear = nuclear_results %>%
  mutate(DE = case_when(
    log2FoldChange < -1 & padj < 0.05 ~ 'Female-biased',
    log2FoldChange > 1 & padj < 0.05 ~ 'Male-biased',
    )
  ) %>%
  ggplot() +
  geom_point(
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      fill = DE
      ),
    alpha = 0.5,
    shape = 21,
    size = 1.5,
    col = 'black'
    ) +
  scale_x_continuous(
    breaks = seq(-25, 25, 5),
    limits = c(-25, 25)
    ) +
  scale_fill_manual(
    values = c('#D55E00', '#0072B2', 'light grey'),
    labels = c('Female-biased', 'Male-biased', 'NA')
    ) +
  geom_vline(
    xintercept = c(-1, 1),
    col = 'black',
    linetype = 'dashed'
    ) +
  geom_hline(
    yintercept = -log10(0.05),
    col = 'black',
    linetype = 'dashed'
    ) +
  labs(
    x = expression('Log'[2]*' Fold-change'),
    y = expression('-Log'[10]*'  Adjusted P-value'),
    fill = 'Biased expression'
    ) +
  theme_minimal() +
  theme(
    legend.position = c(0.8, 0.8),
    legend.background = element_rect(fill = 'white', color = 'black'),
    axis.text = element_text(color = 'black', size = 12),
    axis.title = element_text(color = 'black', size = 12),
    )

# Export
ggsave(
  plot = volcano_nuclear,
  './figures/volcano_nuclear.pdf',
  width = 7.5,
  height = 6,
  unit = 'in'
  )


# Expression of sncRNA targets
volcano_mitonuclear = target_results %>%
  mutate(DE = case_when(
    log2FoldChange < -1 & padj < 0.05 ~ 'Female-biased',
    log2FoldChange > 1 & padj < 0.05 ~ 'Male-biased',
    )
  ) %>%
  ggplot() +
  geom_point(
    aes(
      x = log2FoldChange,
      y = -log10(padj),
      fill = DE
      ),
    alpha = 0.8,
    shape = 21,
    size = 1.5,
    col = 'black'
    ) +
  scale_fill_manual(
    values = c('#D55E00', '#0072B2', '#E1E1E1'),
    labels = c('Female-biased', 'Male-biased', 'NA')
    ) +
  facet_grid(~sncRNA_type) +
  geom_vline(
    xintercept = c(-1, 1),
    col = 'black',
    linetype = 'dashed'
    ) +
  geom_hline(
    yintercept = -log10(0.05),
    col = 'black',
    linetype = 'dashed'
    ) +
  labs(
    x = expression('Log'[2]*' Fold-change'),
    y = expression('-Log'[10]*'  Adjusted P-value'),
    fill = 'Biased expression',
    ) +
  scale_x_continuous(
    breaks = seq(-7.5, 2.5, 2.5),
    ) +
  theme(
    legend.position = c(0.65, 0.8),
    axis.text = element_text(color = 'black', size = 10),
    axis.title = element_text(color = 'black', size = 12),
    panel.grid = element_line(color =  '#EBEBEB'),
    panel.background = element_rect(fill = 'white', color = '#EBEBEB'),
    strip.background = element_rect(fill = '#333333'),
    strip.text = element_text(color = 'white')
  )

# Export
ggsave(
  plot = volcano_mitonuclear,
  './figures/volcano_mt-sncRNA_targets.pdf',
  width = 8,
  height = 5,
  units = 'in'
  )



# 3. GSEA ----------------------------------------------------------------------

# Prepare hallmark pathways (human)
hallmark = gmtPathways('./RNASeq/DE/DE_downstream/h.all.v2022.1.Hs.symbols.gmt')
GOBP = gmtPathways('./RNASeq/DE/DE_downstream/c5.go.bp.v2022.1.Hs.symbols.gmt')

# Prepare list of targets that have human gene symbols
nuclear_results$human_ID = gsub('_[0-9]', '', nuclear_results$Target_annotation)
nuclear_results$human_ID = gsub('-', '', nuclear_results$human_ID)
nuclear_results$human_ID = toupper(nuclear_results$human_ID)
fgsea_input_df = nuclear_results %>%
   select(human_ID, stat) %>%
   na.omit() %>%
   group_by(human_ID) %>%
   summarize(stat = mean(stat))

# Convert to named list
fgsea_input = deframe(fgsea_input_df)


# Run FGSEA
fgsea_output_hallmark = fgseaMultilevel(
  pathways = hallmark,
  stats = fgsea_input,
  minSize = 1,
  maxSize = 50
  )

fgsea_output_GOBP = fgseaMultilevel(
  pathways = GOBP,
  stats = fgsea_input,
  minSize = 1,
  maxSize = 50
  )


# Hallmark pathways
plot_hallmark = fgsea_output_hallmark %>%
  mutate(bias = case_when(
    NES < 0 ~ 'Female-biased',
    NES > 0 ~ 'Male-biased'
    )
  ) %>%
  as_tibble() %>%
  filter(padj < 0.05) %>%
  arrange(desc(NES)) %>%
  mutate(
    pathway = gsub('HALLMARK_', '', pathway)
    ) %>%
  mutate(
    pathway = gsub('_', ' ', pathway)
    ) %>%
  ggplot(
    aes(
      x = NES,
      y = reorder(pathway, NES)
      )
    ) +
  geom_point(
    aes(color = bias),
    size = 3
  ) +
  geom_segment(
    aes(
      x = 0,
      xend = NES,
      y = reorder(pathway, NES),
      yend = reorder(pathway, NES),
      color = bias
      ),
    size = 1.2
    ) +
  labs(
    x = 'Normalized Enrichment Score',
    y = '',
    title = 'Gene Set Enrichment of Hallmark Pathways',
    color = 'Biased enrichment'
    ) +
  scale_color_manual(
    values = c('#D55E00', '#0072B2')
    ) +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(
  plot = plot_hallmark,
  './figures/gsea_plot_hallmark.pdf',
  width = 7.5,
  height = 6,
  unit = 'in'
  )

# Number of significant hallmark pathways:
fgsea_output_hallmark %>% filter(padj < 0.05) %>% nrow()

# GO: Biological Pathways (too many for legible plot)
# plot_GOBP = fgsea_output_GOBP %>%
#   as_tibble() %>%
#   arrange(desc(NES)) %>%
#   mutate(pathway = gsub('GOBP_', '', pathway)) %>%
#   filter(padj < 0.05) %>%
#   filter(size > 10) %>%
#   ggplot(aes(reorder(pathway, NES), NES)) +
#   geom_point() +
#   coord_flip() +
#   labs(x = '',
#        y = 'Normalized Enrichment Score (NES)',
#        title = 'Gene Set Enrichment of GO: Biological Pathways',
#        subtitle = 'Female-biased: NES < 0\nMale-biased: NES > 0') +
#   theme_minimal() +
#   theme(plot.subtitle = element_text(hjust = 0),
#         plot.title = element_text(hjust = 0),
#         legend.position = 'none')
# png('gsea_plot_GOBP.png',
#     width = 1000,
#     height = 2000)


# Combine hallmark and GOBP fgsea results
fgsea_results_complete = full_join(
  fgsea_output_hallmark,
  fgsea_output_GOBP
  ) %>%
  as.data.frame()

# Convert gene sets to character to be able to export as csv
fgsea_results_complete$leadingEdge = as.character(fgsea_results_complete$leadingEdge)

# Export hallmark and GOBP fgsea results
write.csv(
  fgsea_results_complete,
  './RNASeq/DE/DE_downstream/fgsea_results_hallmark_and_GOBP.csv',
  row.names = F
)

# Filter by significantly enriched pathways (p < 0.05)
fgsea_results_p0.05 = fgsea_results_complete %>%
  filter(padj < 0.05) %>%
  as.data.frame()

# Export significantly enriched pathways (p < 0.05)
write.csv(
  fgsea_results_p0.05,
  './RNASeq/DE/DE_downstream/fgsea_results_hallmark_and_GOBP_p0.05.csv',
  row.names = F
  )



# 4. Combine DESeq2, WGCNA, and mt-sncRNA results ------------------------------

# 4a. Join sncRNA modules with WGCNA modules ------------------------------------
wgcna_results = read.csv(
  './RNASeq/WGCNA/WGCNA_downstream/wgcna_results.csv',
  header = T,
  stringsAsFactors = F
  )

wgcna_results = wgcna_results[,1:6]

# Some targets were not assigned as part of modules in WGCNA
complete_target_results = left_join(target_results, wgcna_results)

# Which M-sncRNAs in yellow (female) module?
complete_target_results %>%
  filter(sncRNA_type == 'M-sncRNA') %>%
  filter(module == 'yellow')

# Which M-sncRNAs in greenyellow (male) module?
complete_target_results %>%
  filter(sncRNA_type == 'M-sncRNA') %>%
  filter(module == 'greenyellow')

# Which F-sncRNAs in greenyellow (male) module?
complete_target_results %>%
  filter(sncRNA_type == 'F-sncRNA') %>%
  filter(module == 'greenyellow')

# Which F-sncRNAs in yellow (female) module?
complete_target_results %>%
  filter(sncRNA_type == 'F-sncRNA') %>%
  filter(module == 'yellow')


# 4b. Merge sncRNA names with final manuscript IDs -----------------------------

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


# 4c. sncRNA metadata, target expression, WGCNA kWithin & module ---------------

sncRNA_targets_meta_info_l2fc_WGCNA = inner_join(
  sncRNA_correct_IDs,
  complete_target_results,
  by = c('centroid' = 'centroid')
  ) %>%
  select(
    Name,
    sequence,
    Location,
    row,
    Target_annotation,
    log2FoldChange,
    padj,
    kWithin,
    module
  ) %>%
  rename(
    'mt-sncRNA' = 'Name',
    'mt-sncRNA_sequence' = 'sequence',
    'mt-sncRNA_location' = 'Location',
    'target_ID' = 'row',
    'target_annotation' = 'Target_annotation',
  )

# Export table
write.csv(
  sncRNA_targets_meta_info_l2fc_WGCNA,
  './figures/sncRNA_targets_meta_info_l2fc_wgcna_MISSING_KEGG.csv',
  row.names = F
  )

# Then, manually add in KEGG annotations for un-annotated targets.



# 5. Expression of GCNT1 protein network --------------------------------------

# Import GCNT1 network generated by STRING protein-protein interaction prediction
interactions = read.delim(
  './RNASeq/DE/DE_downstream/pstreck_gcnt1_edges.tsv',
  sep = '\t',
  header = T
  )

# Clean column names
interactions = interactions %>%
  select(1, 2) %>%
  rename(
    from = 'X.node1',
    to = 'node2'
    )

interactions$from = gsub('-.+', '', interactions$from)
interactions$to = gsub('-.+', '', interactions$to)


# Join interactions with DESeq2 results
interactions_expression = inner_join(
  interactions,
  nuclear_results,
  by = c('from' = 'row')
  )

# Boolean vector of p-adj, by < 0.05 or > 0.05
significance = interactions_expression %>%
  distinct(from, .keep_all = T) %>%
  select(padj) %>%
  mutate(significance = case_when(
    padj <= 0.05 ~ TRUE,
    padj > 0.05 ~ FALSE
    )
  )
significance = significance$significance

# Vector of l2fc, should be in increasing order by FUN_ID
l2fc = interactions_expression %>%
  distinct(from, .keep_all = T) %>%
  select(log2FoldChange)
l2fc = l2fc$log2FoldChange

# Vector of annotated names, should be in increasing order by FUN_ID
annotations = interactions_expression %>%
  distinct(from, .keep_all = T) %>%
  mutate(Target_annotation = case_when(
    Target_annotation == '' ~ from,                    # Copy FUN_ID if not annotated
    Target_annotation != '' ~ Target_annotation        # Else use annotation
    )
  ) %>%
  mutate(Target_annotation = case_when(
   Target_annotation == 'FUN_010176' ~ 'MUC13-like',  # KEGG
   Target_annotation == 'FUN_016811' ~ 'GALNT9-like', # KEGG
   Target_annotation == 'FUN_036891' ~ 'Hypothetical',
   Target_annotation != '' ~ Target_annotation,       # Else use annotation,
   )
  ) %>%
  mutate(Target_annotation = case_when(
   Target_annotation == 'GCNT1' ~ 'GCNT1*',
   Target_annotation != '' ~ Target_annotation
   )
  ) %>%
  select(Target_annotation)
annotations = annotations$Target_annotation

# Add expression and annotation data to graph
graph = as_tbl_graph(interactions) %>%
  mutate(l2fc = l2fc) %>%
  mutate(annotations = annotations)

# Specify undirected graph
graph = as.undirected(graph, mode = 'collapse')

# Plot network, colored by expression level.
ggraph(graph, layout = 'gem') +
  geom_edge_diagonal(
    edge_width = 0.1,
    strength = 0.8,
    edge_color = 'black',
    alpha = 0.3
    ) +
  geom_node_label(
    aes(
      label = annotations,
      fill = l2fc
      ),
    size = 4.5
    ) +
  scale_fill_gradient2(
    low = '#D55E00',
    mid = 'white',
    high = '#0072B2',
    midpoint = 0,
    limits = c(-8, 8),
    na.value = "grey50",
    ) +
  labs(
    title = 'Expression of GCNT1 Predicted Protein Network',
    fill = 'log2fold change',
    subtitle = '*p-adj < 0.05'
    ) +
  theme(
    panel.background = element_blank(),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.327)
  )

# Export graph
ggsave(
  './figures/gcnt1_network_graph.pdf',
  width = 10,
  height = 5,
  units = 'in'
  )
