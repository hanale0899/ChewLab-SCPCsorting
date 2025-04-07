##This code is to prepare merged RDS files representing for each biological group for shinyApp## #nolint
##Version 09.12.2024 by Hannah Le##
##Have fun making science works!##

BiocManager::install('hdf5r')

#Load libraries##nolint
library(hdf5r)
library(Seurat)
library(SeuratObject)
library(Matrix)
library(ggplot2)
library(harmony)


#Unsorted samples
data_dir1 <- "C://Users/Chew Lab Leica/Documents/R Data Analysis/RNA sequencing analysis/Jerome/Unsorted_counts_matrix" # nolint
data_dir2 <- "C://Users/Chew Lab Leica/Documents/R Data Analysis/RNA sequencing analysis/Jerome/Unsorted_2_counts_matrix"# nolint

#Small samples
data_dir3 <- "C://Users/Chew Lab Leica/Documents/R Data Analysis/RNA sequencing analysis/Jerome/Small_counts_matrix" #nolint
data_dir4 <- "C://Users/Chew Lab Leica/Documents/R Data Analysis/RNA sequencing analysis/Jerome/Small_2_counts_matrix" # nolint

#Large
data_dir5 <- "C://Users/Chew Lab Leica/Documents/R Data Analysis/RNA sequencing analysis/Jerome/Large_counts_matrix" # nolint
data_dir6 <- "C://Users/Chew Lab Leica/Documents/R Data Analysis/RNA sequencing analysis/Jerome/Large_2_counts_matrix" # nolint

list.files(data_dir1) # Should show barcodes.tsv, genes.tsv, and matrix.mtx

#Read files using Read10X from Seurat# #nolint
data1 <- Read10X(data.dir = data_dir1)
data2 <- Read10X(data.dir = data_dir2)

data3 <- Read10X(data.dir = data_dir3)
data4 <- Read10X(data.dir = data_dir4)

data5 <- Read10X(data.dir = data_dir5)
data6 <- Read10X(data.dir = data_dir6)

#Creat Seurat Object for the analysis# #nolint
seurat_object1 = CreateSeuratObject(counts = data1)
seurat_object2 = CreateSeuratObject(counts = data2)

seurat_object3 = CreateSeuratObject(counts = data3)
seurat_object4 = CreateSeuratObject(counts = data4)

seurat_object5 = CreateSeuratObject(counts = data5)
seurat_object6 = CreateSeuratObject(counts = data6)

#Append conditions/sample groups
seurat_object1$condition = "Unsorted"
seurat_object2$condition = "Unsorted"
seurat_object3$condition = "Small"
seurat_object4$condition = "Small"
seurat_object5$condition = "Large"
seurat_object6$condition = "Large"

#Append batch information
seurat_object1$batch = "1"
seurat_object2$batch = "2"
seurat_object3$batch = "1"
seurat_object4$batch = "2"
seurat_object5$batch = "1"
seurat_object6$batch = "2"

##QC for percentage of mitochondrial genes
seurat_object1[["percent.mt"]] <- PercentageFeatureSet(seurat_object1, pattern = "^MT-") # nolint

# Filter out low-quality cells
seurat_obj_filtered1 <- subset(seurat_object1, 
                     subset = nFeature_RNA > 200 &  # At least 500 genes detected # nolint
                       nFeature_RNA < 100000 &          # Less than 100000 genes (to remove doublets) # nolint
                       percent.mt < 10)                # Less than 10% mitochondrial genes for human genome # nolint


seurat_object2[["percent.mt"]] <- PercentageFeatureSet(seurat_object2, pattern = "^MT-") # nolint

# Filter out low-quality cells
seurat_obj_filtered2 <- subset(seurat_object2, 
                     subset = nFeature_RNA > 200 &  # At least 500 genes detected # nolint
                       nFeature_RNA < 100000 &          # Less than 100000 genes (to remove doublets) # nolint
                       percent.mt < 10)                # Less than 10% mitochondrial genes for human genome # nolint

seurat_object3[["percent.mt"]] <- PercentageFeatureSet(seurat_object3, pattern = "^MT-") # nolint: line_length_linter.

# Filter out low-quality cells
seurat_obj_filtered3 <- subset(seurat_object3, 
                     subset = nFeature_RNA > 200 &  # At least 500 genes detected # nolint: line_length_linter.
                       nFeature_RNA < 100000 &          # Less than 100000 genes (to remove doublets) # nolint: line_length_linter.
                       percent.mt < 5)                # Less than 10% mitochondrial genes for human genome # nolint: line_length_linter.

seurat_object4[["percent.mt"]] <- PercentageFeatureSet(seurat_object4, pattern = "^MT-") # nolint: line_length_linter.

# Filter out low-quality cells
seurat_obj_filtered4 <- subset(seurat_object4, 
                     subset = nFeature_RNA > 200 &  # At least 500 genes detected # nolint: line_length_linter.
                       nFeature_RNA < 100000 &          # Less than 100000 genes (to remove doublets) # nolint: line_length_linter.
                       percent.mt < 5)                # Less than 10% mitochondrial genes for human genome # nolint: line_length_linter.

seurat_object5[["percent.mt"]] <- PercentageFeatureSet(seurat_object5, pattern = "^MT-") # nolint: line_length_linter.

# Filter out low-quality cells
seurat_obj_filtered5 <- subset(seurat_object5, 
                     subset = nFeature_RNA > 200 &  # At least 500 genes detected # nolint: line_length_linter.
                       nFeature_RNA < 100000 &          # Less than 100000 genes (to remove doublets) # nolint: line_length_linter.
                       percent.mt < 5)                # Less than 10% mitochondrial genes for human genome # nolint: line_length_linter.

seurat_object6[["percent.mt"]] <- PercentageFeatureSet(seurat_object6, pattern = "^MT-") # nolint: line_length_linter.
# Filter out low-quality cells
seurat_obj_filtered6 <- subset(seurat_object6, 
                     subset = nFeature_RNA > 200 &  # At least 500 genes detected # nolint: line_length_linter.
                       nFeature_RNA < 100000 &          # Less than 100000 genes (to remove doublets) # nolint: line_length_linter.
                       percent.mt < 5)                # Less than 10% mitochondrial genes for human genome # nolint: line_length_linter.


# Normalize the data
seurat_obj_norm1 <- NormalizeData(seurat_obj_filtered1, normalization.method = "LogNormalize", scale.factor = 10000) # nolint: line_length_linter.
seurat_obj_norm2 <- NormalizeData(seurat_obj_filtered2, normalization.method = "LogNormalize", scale.factor = 10000) # nolint: line_length_linter.
seurat_obj_norm3 <- NormalizeData(seurat_obj_filtered3, normalization.method = "LogNormalize", scale.factor = 10000) # nolint: line_length_linter.
seurat_obj_norm4 <- NormalizeData(seurat_obj_filtered4, normalization.method = "LogNormalize", scale.factor = 10000) # nolint: line_length_linter.
seurat_obj_norm5 <- NormalizeData(seurat_obj_filtered5, normalization.method = "LogNormalize", scale.factor = 10000) # nolint: line_length_linter.
seurat_obj_norm6 <- NormalizeData(seurat_obj_filtered6, normalization.method = "LogNormalize", scale.factor = 10000) # nolint: line_length_linter.

#export raw counts to csv
# Assuming seurat_obj is your Seurat object
# Extract the count matrix from the Seurat object
count_matrix1 <- GetAssayData(seurat_obj_norm1, layer = "counts")
count_matrix2 <- GetAssayData(seurat_obj_norm2, layer = "counts")
count_matrix3 <- GetAssayData(seurat_obj_norm3, layer = "counts")
count_matrix4 <- GetAssayData(seurat_obj_norm4, layer = "counts")
count_matrix5 <- GetAssayData(seurat_obj_norm5, layer = "counts")
count_matrix6 <- GetAssayData(seurat_obj_norm6, layer = "counts")

# Write to rds (with gene names as row names)
saveRDS(count_matrix1, file = "scRNA_count_matrix_unsorted.rds")
saveRDS(count_matrix2, file = "scRNA_count_matrix_unsorted_2.rds")
saveRDS(count_matrix3, file = "scRNA_count_matrix_small.rds")
saveRDS(count_matrix4, file = "scRNA_count_matrix_small_2.rds")
saveRDS(count_matrix5, file = "scRNA_count_matrix_large.rds")
saveRDS(count_matrix6, file = "scRNA_count_matrix_large_2.rds")

# Merge the multiple Seurat objects into one
combined_seurat <- merge(seurat_obj_norm1, y = c(seurat_obj_norm2, seurat_obj_norm3, seurat_obj_norm4, seurat_obj_norm5, seurat_obj_norm6), 
                         add.cell.ids = c("Unsorted", "Unsorted2","Small", "Small2", "Large", "Large2"))


combined_seurat <- JoinLayers(combined_seurat)

# Normalize the combined object (optional)
combined_seurat_norm <- NormalizeData(combined_seurat)

# Scale and run PCA (optional)
combined_seurat_norm <- FindVariableFeatures(combined_seurat_norm)
combined_seurat_norm <- ScaleData(combined_seurat_norm, features = VariableFeatures(combined_seurat_norm)) #nolint
combined_seurat_norm <- RunPCA(combined_seurat_norm, features = VariableFeatures(combined_seurat_norm)) #nolint
combined_seurat_norm <- RunHarmony(combined_seurat_norm, group.by.vars = "batch")#Remove batch effects
# Perform clustering (optional)
combined_seurat_norm <- FindNeighbors(combined_seurat_norm, dims = 1:20)
combined_seurat_norm <- FindClusters(combined_seurat_norm, resolution = 1.2)

combined_seurat_umap <- RunUMAP(combined_seurat_norm, reduction = "harmony", dims = 1:20)
a <- DimPlot(combined_seurat_umap , reduction = "umap")
ggsave(a, file="rawcluster.jpg", width=6, height=6)

# Save the merged Seurat object
saveRDS(combined_seurat_umap, file = "combined.rds")

# Now let's assume that 'condition' is the metadata column you want to highlight by
# Create a new column in metadata if you want to highlight certain conditions specifically
highlight_condition <- 'Unsorted'  # Replace this with the condition you want to highlight

# Create a new factor that highlights the condition of interest
combined_seurat_umap$highlight <- ifelse(combined_seurat_umap$condition == highlight_condition, highlight_condition, "Other")

# UMAP plot with highlighting specific condition
p <- DimPlot(combined_seurat_umap, reduction = "umap", group.by = "highlight", 
             cols = c("gray", "purple")) +  # You can change the colors here as needed
  labs(title = paste("UMAP -", highlight_condition)) +
  theme_minimal()

ggsave(p, file="rawcluster_unsorted.jpg", width =6, height =6)


t <- DimPlot(combined_seurat_umap, reduction = "umap", group.by = "condition", 
             cols = c("red", "darkblue","darkgreen")) +  # You can change the colors here as needed
  labs(title = "UMAP - Groups") +
  theme_minimal()
ggsave(t, file="rawcluster_allgroups.jpg", width =6, height =6)


