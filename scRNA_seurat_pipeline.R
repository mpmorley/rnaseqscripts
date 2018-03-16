#install.packages("Seurat")
library(Seurat)
library(dplyr)
library(Matrix)
library(openxlsx)

# Load the dataset
pbmc.data <- Read10X(data.dir = "DF_E13.5_Nkx2.1/outs/filtered_gene_bc_matrices/mm10/")

# Examine the memory savings between regular and sparse matrices
dense.size <- object.size(x = as.matrix(x = pbmc.data))
dense.size
sparse.size <- object.size(x = pbmc.data)
sparse.size
dense.size/sparse.size

# Initialize the Seurat object with the raw (non-normalized data).  
pbmc <- CreateSeuratObject(raw.data = pbmc.data, min.cells = 3, min.genes = 200,project = "10X_PBMC")

# calculate the percent.mito values.
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pbmc@data), value = TRUE)
percent.mito <- Matrix::colSums(pbmc@raw.data[mito.genes, ])/Matrix::colSums(pbmc@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
pbmc <- AddMetaData(object = pbmc, metadata = percent.mito, col.name = "percent.mito")
VlnPlot(object = pbmc, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)

# Plot GenePlot
par(mfrow = c(1, 2))
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = pbmc, gene1 = "nUMI", gene2 = "nGene")


# Filter out cells that have unique gene counts less than 500 
pbmc <- FilterCells(object = pbmc, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.05))

#################################
#Define cell types in meta data based on a certain criteria (if present. Else proceed without this section)
my.data=FetchData(pbmc,c("ident","nGene","Hopx","Sftpc"))
my.data$celltype=ifelse(my.data$Hopx>=1 & my.data$Sftpc <10, "AT1", ifelse(my.data$Hopx<1 & my.data$Sftpc >=10,"AT2",ifelse(my.data$Hopx>=1 & my.data$Sftpc >=10,"AT1/AT2","NULL")))
my.data=my.data %>% select(celltype)
my.data2=as.character(my.data$celltype)
names(my.data2)=rownames(my.data)
pbmc <- AddMetaData(object = pbmc, metadata = my.data2, col.name = "celltype")
###################################

#normalize data
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize",scale.factor = 10000)

#detection of variable genes
#calculates the average expression and dispersion for each gene, places these genes into bins, and then calculates a z-score for dispersion within each bin
pbmc <- FindVariableGenes(object = pbmc, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = pbmc@var.genes)

# Also read in a list of cell cycle markers, from Tirosh et al, 2015
cc.genes <- readLines(con = "data/regev_lab_cell_cycle_genes.txt")
cc.genes= tolower(cc.genes)
cc.genes=paste0(toupper(substr(cc.genes, 1, 1)), substr(cc.genes, 2, nchar(cc.genes)))
# Segregate this list into markers of G2/M phase and markers of Sphase
s.genes <- cc.genes[1:43]
g2m.genes <- cc.genes[44:97]

#Assign scores in the CellCycleScoring function.Stores S and G2/M scores in object@meta.data, along with the predicted classification of each cell in either G2M, S or G1 phase
pbmc <- CellCycleScoring(object = pbmc, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)

# view cell cycle scores and phase assignments
head(x = pbmc@meta.data)

#Scaling the data and removing unwanted sources of variation
pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito","S.Score", "G2M.Score"))
#pbmc <- ScaleData(object = pbmc, vars.to.regress = c("nUMI", "percent.mito"))

#Perform linear dimensional reduction (Note: performed on the variable genes)
pbmc <- RunPCA(object = pbmc, pc.genes = pbmc@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)
VizPCA(object = pbmc, pcs.use = 1:2)
PCAPlot(object = pbmc, dim.1 = 1, dim.2 = 2)
PCHeatmap(object = pbmc, pc.use = 1:15, cells.use = 500, do.balanced = TRUE, 
          label.columns = FALSE, use.full = FALSE)


#Determine statistically significant principal components by randomly permuting a subset of the data (1% by default) and rerunning PCA
pbmc <- JackStraw(object = pbmc, num.replicate = 100, do.print = FALSE)

JackStrawPlot(object = pbmc, PCs = 1:14)
PCElbowPlot(object = pbmc)

#Clustering
pbmc <- FindClusters(object = pbmc, reduction.type = "pca", dims.use = 1:14, resolution = 0.6, print.output = 0, save.SNN = TRUE)

PrintFindClustersParams(object = pbmc)

#Run tsne
pbmc <- RunTSNE(object = pbmc, dims.use = 1:14, do.fast = TRUE)
TSNEPlot(object = pbmc,pt.size=2,group.by = "ident")
#TSNEPlot(object = pbmc,pt.size=2,group.by = "celltype",colors.use = c("coral2","lightgray","turquoise3","mediumorchid"))

#Create Feature Plot of genes of choice
FeaturePlot(object = pbmc, features.plot = c("Hopx","Sftpc","Sox2","Sox9"), cols.use = c("grey", "blue"), 
            reduction.use = "tsne")

#find cluster markers
cluster1_2.markers <- FindMarkers(object = pbmc, ident.1 = 1,ident.2 =2, min.pct = -Inf)
print(x = head(x = cluster1_2.markers, n = 5))
cluster1_2.markers=cluster1_2.markers[order(-cluster1_2.markers$avg_logFC),]

cluster1_all.markers <- FindMarkers(object = pbmc, ident.1 = 1,ident.2 = c(2,3), min.pct = 0.25)
print(x = head(x = cluster1_all.markers, n = 5))

#write output
result=list(AT1_vs_AT2=cluster1_2.markers,AT1_vs_all=cluster1_all.markers)

write.xlsx(result,file="results.xlsx",row.names=T)
