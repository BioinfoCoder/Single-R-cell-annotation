# script to annotate cell types from 20k Human PBMCs from a healthy female donor
# setwd("~/Desktop/demo/singleCell_singleR_part2/scripts")

library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
library(pheatmap)

# 10X CellRanger .HDF5 format ---------
hdf5_obj <- Read10X_h5(filename = "C:/Users/doaad/Downloads/20k_PBMC_3p_HT_nextgem_Chromium_X_raw_feature_bc_matrix.h5",
                       use.names = TRUE,
                       unique.features = TRUE)
pbmc.seurat <- CreateSeuratObject(counts = hdf5_obj)
pbmc.seurat

# QC and Filtering ----------
# explore QC
pbmc.seurat$mitoPercent <- PercentageFeatureSet(pbmc.seurat, pattern = '^MT-')
VlnPlot(pbmc.seurat, features = c("nFeature_RNA", "nCount_RNA", "mitoPercent"), ncol = 3)
pbmc.seurat.filtered <- subset(pbmc.seurat, subset = nCount_RNA > 800 &
                                 nFeature_RNA > 500 &
                                 mitoPercent < 10)

# It is a good practice to filter out cells with non-sufficient genes identified and genes with non-sufficient expression across cells.


# pre-process standard workflow ---------------
pbmc.seurat.filtered <- NormalizeData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindVariableFeatures(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- ScaleData(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunPCA(object = pbmc.seurat.filtered)
ElbowPlot(pbmc.seurat.filtered)
pbmc.seurat.filtered <- FindNeighbors(object = pbmc.seurat.filtered, dims = 1:20)
pbmc.seurat.filtered <- FindClusters(object = pbmc.seurat.filtered)
pbmc.seurat.filtered <- RunUMAP(object = pbmc.seurat.filtered, dims = 1:20)

# running steps above to get clusters
DimPlot(pbmc.seurat.filtered, reduction = "umap")
View(pbmc.seurat.filtered@meta.data)

# run SingleR with multiple reference datasets (default mode) ---------

# for pbmc data, we will use two datasets
hpca <- celldex::HumanPrimaryCellAtlasData()#blood subpopulation
dice <- celldex::DatabaseImmuneCellExpressionData()#bulk immune cell 

# ...1. Strategy 1: Using reference-specific labels ----------
hpca$label.main
dice$label.main#the problem that there is common labels between referances

# adding ref info to labels 
#add label to each refer to distinguish where it is come from 
hpca$label.main <- paste0('HPCA.', hpca$label.main)
dice$label.main <- paste0('DICE.', dice$label.main)

# create a combined ref based on shared genes THAT CONTAIN THE SHARED LABELS BY THE INTERSECTION 
shared <- intersect(rownames(hpca), rownames(dice))
combined <- cbind(hpca[shared,], dice[shared,])#combined the refrences based on the shared genes
combined #and it is distinguish by the string we paste earlier
combined$label.main

# run singleR using combined ref
# savings counts into a separate object
pbmc_counts <- GetAssayData(pbmc.seurat.filtered, slot = 'counts')#get the count from seurat object to run single R 

com.res1 <- SingleR(test = pbmc_counts, ref = combined, labels = combined$label.main)
table(com.res1$labels)
#add the labels to the metadata of the seurat object 
pbmc.seurat.filtered$com.res1.labels <- com.res1[match(rownames(pbmc.seurat.filtered@meta.data), rownames(com.res1)), 'labels']
View(pbmc.seurat.filtered@meta.data)

DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = 'com.res1.labels', label = TRUE)

# ...2. Strategy 2: Comparing scores across references ----------

hpca$label.main
dice$label.main
#If you want to compare cell types across both references (hpca and dice), it's easier when the labels are identical — without the dataset-specific prefixes
hpca$label.main <- gsub('HPCA\\.','', hpca$label.main)
dice$label.main <- gsub('DICE\\.','', dice$label.main)

com.res2 <- SingleR(test = pbmc_counts, 
                    ref = list(HPCA = hpca, DICE = dice),
                    labels = list(hpca$label.main, dice$label.main))

# Check the final label from the combined assignment.
table(com.res2$labels)

# which reference scored best for which label?
grouping <- paste0(com.res2$labels,'.', com.res2$reference)
best_ref <- as.data.frame(split(com.res2, grouping))

# get de. genes from each individual references
metadata(com.res2$orig.results$HPCA)$de.genes
metadata(com.res2$orig.results$DICE)$de.genes
# Combined diagnostics
plotScoreHeatmap(com.res2)#it first visualise the score in rhe combined results and for each ref 
pbmc.seurat.filtered$com.res2.labels <- com.res2[match(rownames(pbmc.seurat.filtered@meta.data), rownames(com.res2)), 'labels']

DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = 'com.res2.labels', label = TRUE)

# ...3. Strategy 3: Using Harmonized Labels ---------- Harmonized Labels is a standard vocublary we will use cell ontology it will use the same term across ref

hpca.ont <- celldex::HumanPrimaryCellAtlasData(cell.ont = 'nonna')#only samples with nonany ontology will retrieve 
dice.ont <- celldex::DatabaseImmuneCellExpressionData(cell.ont = 'nonna')

# Using the same sets of genes:
shared <- intersect(rownames(hpca.ont), rownames(dice.ont))#retain genes shared by both ref 
hpca.ont <- hpca.ont[shared,]
dice.ont <- dice.ont[shared,]

# Showing the top 10 most frequent terms:
tail(sort(table(hpca.ont$label.ont)),10)#explore the ten most frequent terms
tail(sort(table(dice.ont$label.ont)), 10)

# using label.ont instead on label.main while running SingleR

com.res3 <- SingleR(test = pbmc_counts,
                    ref = list(HPCA = hpca.ont, DICE = dice.ont),
                    labels = list(hpca.ont$label.ont, dice.ont$label.ont))#labels are ontology terms 


table(com.res3$labels)#the chalenge is to get what these gene ontology means and what cell types correspond to it 



# How to map cell ontology terms? ----------------

colData(hpca.ont)#by get the col data for each ref as these ontology terms map to label name or fine label 
colData(dice.ont)

hpca.fle <- system.file("mapping","hpca.tsv", package = "celldex")#files provide mapping between ontology terms and cell names #provide a path for where this mapping files present 

hpca.mapping <- read.delim(hpca.fle, header = F)#read this mapping file using read.delim

pbmc.seurat.filtered$com.res3.labels <- com.res3[match(rownames(pbmc.seurat.filtered@meta.data), rownames(com.res3)), 'labels']


DimPlot(pbmc.seurat.filtered, reduction = 'umap', group.by = 'com.res3.labels', label = TRUE)


# Rename for clarity
colnames(hpca.mapping) <- c("main.label", "ontology")

# Match ontology labels to main labels
mapped.labels <- hpca.mapping$main.label[match(com.res3$labels, hpca.mapping$ontology)]

# Add to Seurat metadata
pbmc.seurat.filtered$main.label <- mapped.labels

DimPlot(pbmc.seurat.filtered, reduction = "umap", group.by = "main.label", label = TRUE)

