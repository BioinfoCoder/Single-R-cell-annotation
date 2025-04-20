BiocManager::install("SingleR")
BiocManager::install("Seurat")
BiocManager::install("tidyverse")
BiocManager::install("celldex")
##################################
# ---------------------- #
# SingleR_tutorial.R
# ---------------------- #

# Setting up environment ===================================================
# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation

# Set seed
set.seed(42)
# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(celldex)
library(SingleR)
library(Seurat)
BiocManager::install("scRNAseq")
library(scRNAseq)

# We'll use a PBMC dataset from the R package scRNAseq
#PBMC are immune cella in the blood
sce <- scRNAseq::KotliarovPBMCData(mode = c('rna'))
seu <- CreateSeuratObject(counts = counts(sce), meta.data = as.data.frame(colData(sce)))
rm(sce)
seu <- NormalizeData(object = seu)#adjusting for sequencing depth
# Additional Seurat preprocessing steps - 
# This is not necessary for this tutorial but I will use it for visualisation later
seu <- FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000)
gc() #free up some memory
seu <- ScaleData(seu, features = VariableFeatures(seu))#centring the data mean =0 SD=1
seu <- RunPCA(seu, features = VariableFeatures(object = seu))
seu <- FindNeighbors(seu, dims = 1:10)
seu <- FindClusters(seu, resolution = 0.5)
seu <- RunUMAP(seu, dims = 1:10) 
raw_counts <- LayerData(seu, assay = "RNA", layer = 'counts') #extract raw count from seurat object using the function layerdata
raw_counts[c('VIM', 'BCL2', 'TP53', 'CD4'),1:5]
norm_counts <- LayerData(seu, assay = "RNA", layer = 'data') #
gh<-norm_counts[c('VIM', 'BCL2', 'TP53', 'CD4'),1:5]
scaled_genes <- GetAssayData(seu, slot = "scale.data")

# 2. Get reference dataset(very general)
#this pakage for download single cell annotated dataset
ref <- celldex::HumanPrimaryCellAtlasData()
unique(ref$label.main) #label main major
unique(ref$label.fine)#label fine specific
# Subset to include only relevant cell types (CAREFUL!)
ref <- ref[,grepl('DC|B_cell|Neutrophils|T_cells|Monocyte|Erythroblast|
                 Macrophage|NK_cell|Platelets|Myelocyte', ref$label.main)]
unique(ref$label.main)
BiocManager::install("scran")
library(scran)
# 3. Run SingleR
ct_ann <- SingleR(test = norm_counts, # we could also use sce or raw_counts
                  ref = ref, 
                  labels = ref$label.main,
                  de.method = 'wilcox')#wilcoxon ranked test for differential expression recommended for single cell data , for bulk keep the default
#this is the output of singleR (ct_ann) in form of dataframe with predicted labels 
#if the label in the test dataset not found in ref singleR still able to asign label for it(for thid reasons it assign incorest labels)
#to avoid this we remove low quality prediction with low scores
#so yhey assigned delta values which is a gap or diff. between assigned label and the median scores accross labels 
#if the delta is small this indicate that the cell matches all labels with the same confidance so the assigned label not very meanigful 
#so single R can discard cells with low delta values caused by ambigeous assignment with closlely related referance labels and incorect assignment that matches very poorly to all referances 
#so in prune label colonm we find cleaner or more reliable labels 
#we can compare prune labels with the intial labels
#there is a slight change between the numbers of pruned labels and the number of labels 
# Inspect quality of the predictions
#cells with too small delta will clean up 
plotScoreHeatmap(ct_ann) #to inspect the quality of the prediction 
#on y accesss we have the referance labels which is the rows and each colunm is a cell 
#example cells with high score of T cells mostly clasified as t cells 
plotDeltaDistribution(ct_ann, ncol = 4, dots.on.top = FALSE)#plot the distripution of the gap between the score of the assignrd labels and the score of all the remaining labels accross cells assigned to each referance label
#we want high deltas as it gives more confident to the assigned labels 
# Add to seurat object
rownames(ct_ann)[1:5] # make sure you have cell IDs as the row names that the function will use to match the prediction 
#the rownames(ct_ann) must be equal to rownames(seu) they schould have the same cell ids 
seu <- AddMetaData(seu, ct_ann$pruned.labels, col.name = 'SingleR_HCA')#to add this cell annotation to our seurat object 
#here we will save the pruned labels and give it a colunm name SingleR_HCA(human cell atlas)
# Visualise them on the UMAP
seu <- SetIdent(seu, value = "SingleR_HCA")#color the umap based on SingleR predictions
DimPlot(seu, label = T , repel = T, label.size = 3) + NoLegend()
####################################################################################################
#tips and tricks for singleR#
#make sure that the gene annotation in your referance is the same at your data set(make sure the gene names are the same between referance and test dataset
#choose ref data set curfully it well determine wich cells will be annotated ion your dataset(search for publication that use similar methods and tissue with similar origin)
#recommendation use it with different reference dataset and find a consensus between them this will give more reliable cell type annotation 
#cell types depend on gene expression level and the combination of genes