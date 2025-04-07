library(dplyr)
library(tidyr)
Microglia <- c("TMEM119", "CX3CR1", "P2RY12")
Oligodendrocytes <- c("SOX10", "OLIG2","MBP","PLP1","CNP")
OPCs <- c("SOX10", "OLIG1", "OLIG2", "PDGFRA")
Astrocytes <- c("GFAP", "SOX9", "ALDH1L1", "SLC2A1","S100B")
Fibroblasts <- c("COL1A1", "COL6A1", "DCN")
Pericytes <- c("KCNJ8", "HIGD1B", "ABCC9", "RGS5")
Ependymal_cells <- c("CFAP126", "TMEM212", "PIFO", "TEKT1")
Endothelial_cells <- c("PECAM1", "KDR", "FLT1")
Neutrophils <- c("MMP9", "S100A9", "CD177", "LTF")
Lymphocytes <- c("CD3E", "CD19", "NKG7")
Neurons <- c("RBFOX3", "MAP2", "TUBB3", "NEUROD1", "SYP", "CAMK2A")
Early_neurons <-c("TUBB3","STMN2","NEFL","NEFH","NCAM1","POU3F2","MAP2","MAPT","ASCL1","TUBB51","TUBA1A","PAX6")
Dividing_NPCs <- c("SOX2", "NES", "DCX", "MALAT1","ANXA2","VCP","MSX1", "PTN")
Schwann_cells <- c("S100B", "PLP1", "NF1")
Progenitors <- c("SOX1", "SOX2","NES","MKI67","POU5F1","SOX11","LIN28A")
Motor_neurons <- c("ISL1","MNX1","CHAT","FOXP1","FOXP2","CHODL","TENM2")
Sensory_neurons <- c("NTRK1","NTRK2","TRPM3","TRPC1","RUNX1","TRPV1")
gene_sets <- list(
  Microglia = Microglia,
  Monocytes = Monocytes,
  Macrophages = Macrophages,
  Neutrophils = Neutrophils,
  Lymphocytes = Lymphocytes,
  OPCs = OPCs,
  Endothelial_cells = Endothelial_cells,
  Fibroblasts = Fibroblasts,
  Pericytes = Pericytes,
  Ependymal_cells = Ependymal_cells,
  Astrocytes = Astrocytes,
  Oligodendrocytes = Oligodendrocytes,
  Neurons = Neurons,
  Dividing_NPCs = Dividing_NPCs,
  Schwann_cells = Schwann_cells,
  Progenitors = Progenitors,
  Motor_neurons = Motor_neurons,
  Sensory_neurons = Sensory_neurons,
  Pan_neurons = Pan_neurons
)

# Convert the list to a data.frame with equal number of rows
# The number of rows will be determined by the maximum length of any gene list
max_length <- max(sapply(gene_sets, length))

# Make all lists the same length, pad with NA if necessary
gene_sets_padded <- lapply(gene_sets, function(x) {
  length(x) <- max_length
  return(x)
})

# Create the data.frame
gene_df <- data.frame(gene_sets_padded)

# View the resulting data.frame
head(gene_df)

# Check which markers are missing in the Seurat object
all_markers <- unique(unlist(gene_sets))

missing_markers <- setdiff(all_markers, rownames(combined_seurat_umap))

if(length(missing_markers) > 0) {
  print("Missing markers:")
  print(missing_markers)
} else {
  print("All markers are present in the Seurat object.")
}# Create a list of gene sets

#FindMarkersGeneofEachCluster
markers <- FindAllMarkers(combined_seurat_umap, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.5)
head(markers)

assign_oligo_subtype <- function(cluster_markers) {
  if (any(cluster_markers %in% Microglia)) return("Microglia")
  if (any(cluster_markers %in% Oligodendrocytes)) return("Oligodendrocytes")
  if (any(cluster_markers %in% OPCs)) return("OPCs")
  if (any(cluster_markers %in% Astrocytes)) return("Astrocytes")
  if (any(cluster_markers %in% Ependymal_cells)) return("Ependymal_cells")
  if (any(cluster_markers %in% Neurons)) return("Neurons")
  if (any(cluster_markers %in% Dividing_NPCs)) return("Dividing_NPCs")
  if (any(cluster_markers %in% Schwann_cells)) return("Schwann_cells")
  if (any(cluster_markers %in% Progenitors)) return("Progenitors")
  if (any(cluster_markers %in% Motor_neurons)) return("Motor_neurons")
  if (any(cluster_markers %in% Sensory_neurons)) return("Sensory_neurons")
  return("Unknown")
}

#Assign the celltype using the pre-written function
combined_seurat_umap$celltype <- sapply(combined_seurat_umap$seurat_clusters, function(cluster) {
  cluster_genes <- markers %>% filter(cluster == !!cluster) %>% pull(gene)
  assign_oligo_subtype(cluster_genes)
})
unique(combined_seurat_umap$celltype)
#Check the annotation
unique(combined_seurat_umap$seurat_clusters[combined_seurat_umap$celltype == "Unknown"])
write.csv(markers, file = "screeningpurpose.csv")

#Extract the metadata
metadata <- cbind(seurat_obj$condition,seurat_obj$celltype)
metadata <- as.data.frame(metadata)
metadata <- metadata %>% rename(condition = V1, cell_type = V2)
# Count the occurrences of each cell type within each condition
cell_type_condition_counts <- metadata %>%
  group_by(condition, cell_type) %>%
  summarise(count = n()) %>%
  ungroup()

# Calculate the percentage of each cell type within each condition
cell_type_condition_counts <- cell_type_condition_counts %>%
  group_by(condition) %>%
  mutate(percentage = (count / sum(count)) * 100)

# Generate a custom palette with 20 colors by combining Set1, Set2, and Set3
library(RColorBrewer)
custom_palette <- c(
  brewer.pal(8, "Set2"),        # 8 colors from Set1        
  brewer.pal(8, "Set3")         # 8 colors from Set3
)

#Rearrange the order of the plot order
cell_type_condition_counts$condition <- factor(cell_type_condition_counts$condition, levels = c("Unsorted", "Small","Large"))
cell_type_condition_counts$cell_type <- factor(cell_type_condition_counts$cell_type, levels = c("Progenitors", "Dividing_NPCs","Astrocytes","Neurons"))
# Create the stacked bar chart with annotations
y <- ggplot(cell_type_condition_counts, aes(x = condition, y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.5) +
  labs(title = "Cell Type Distribution Across Conditions",
       x = "Condition", y = "Percentage (%)") +
  theme_minimal() +
  scale_fill_manual(values = custom_palette) +  # Choose a color palette for the bars
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size =14),
        legend.title = element_blank(),
        legend.position = "top")   # Ensure full visibility of the stacked bars
  
ggsave(y, file = "Composition_condition.jpg", width = 6, height = 6, dpi = 600)  # Add text labels to each section

#Plot the raw UMAP with the cell type annotation
u <- DimPlot(combined_seurat_umap, reduction = "umap", group.by = "celltype") + # You can change the colors here as needed
  labs(title = "UMAP - Cell Types") +
  theme_minimal()
ggsave(u, file="rawcluster_celltype.jpg", width =6, height =6)