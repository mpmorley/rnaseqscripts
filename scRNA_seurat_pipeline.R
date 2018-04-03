#install.packages("Seurat")
library(Seurat)
library(dplyr)
library(Matrix)
library(openxlsx)


################################### Set input parameters ##################################
#
dir='DF_E15.5_Nkx2.1'
project="DF_E15.5_Nkx2.1"
type="single" #choose either single or aggregate
regresscellcycle="no" #choose either yes or no
#
##########################################################################################################

# Load the dataset
if(type=="single"){
inputdata <- Read10X(data.dir = paste(dir,"/outs/filtered_gene_bc_matrices/mm10/",sep=""))
}else{
  inputdata <- Read10X(data.dir = paste(dir,"/outs/filtered_gene_bc_matrices_mex/mm10/",sep=""))
}


# Initialize the Seurat object with the raw (non-normalized data).  
scrna <- CreateSeuratObject(raw.data = inputdata, min.cells = 3, min.genes = 200,project = project)

# calculate the percent.mito values.
mito.genes <- grep(pattern = "^mt-", x = rownames(x = scrna@data), value = TRUE)
percent.mito <- Matrix::colSums(scrna@raw.data[mito.genes, ])/Matrix::colSums(scrna@raw.data)

# AddMetaData adds columns to object@meta.data, and is a great place to stash QC stats
scrna <- AddMetaData(object = scrna, metadata = percent.mito, col.name = "percent.mito")

pdf(file=paste(dir,"/plots/",project,"_violinplot.pdf",sep=""),height = 7,width = 11)
VlnPlot(object = scrna, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3)
dev.off()

# Plot GenePlot
pdf(file=paste(dir,"/plots/",project,"_geneplot.pdf",sep=""),height = 7,width = 11)
par(mfrow = c(1, 2))
GenePlot(object = scrna, gene1 = "nUMI", gene2 = "percent.mito")
GenePlot(object = scrna, gene1 = "nUMI", gene2 = "nGene")
dev.off()

# Filter out cells that have unique gene counts less than 500 
scrna <- FilterCells(object = scrna, subset.names = c("nGene", "percent.mito"), 
                    low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.05))

#################################
#Define cell types in meta data based on a certain criteria (if present. Else proceed without this section)
# my.data=FetchData(scrna,c("ident","nGene","Hopx","Sftpc"))
# my.data$celltype=ifelse(my.data$Hopx>=1 & my.data$Sftpc <10, "AT1", ifelse(my.data$Hopx<1 & my.data$Sftpc >=10,"AT2",ifelse(my.data$Hopx>=1 & my.data$Sftpc >=10,"AT1/AT2","NULL")))
# my.data=my.data %>% select(celltype)
# my.data2=as.character(my.data$celltype)
# names(my.data2)=rownames(my.data)
# scrna <- AddMetaData(object = scrna, metadata = my.data2, col.name = "celltype")
###################################

#normalize data
scrna <- NormalizeData(object = scrna, normalization.method = "LogNormalize",scale.factor = 10000)

#detection of variable genes
#calculates the average expression and dispersion for each gene, places these genes into bins, and then calculates a z-score for dispersion within each bin
scrna <- FindVariableGenes(object = scrna, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = scrna@var.genes)

if(regresscellcycle=="yes"){
  # Read in a list of cell cycle markers, from Tirosh et al, 2015
    cc.genes <- readLines(con = "/Users/bapoorva/Downloads/cell_cycle_vignette_files/regev_lab_cell_cycle_genes.txt")
    cc.genes= tolower(cc.genes)
    cc.genes=paste0(toupper(substr(cc.genes, 1, 1)), substr(cc.genes, 2, nchar(cc.genes)))
  
  # Segregate this list into markers of G2/M phase and markers of Sphase
    s.genes <- cc.genes[1:43]
    g2m.genes <- cc.genes[44:97]
  
  #Assign scores in the CellCycleScoring function.Stores S and G2/M scores in object@meta.data, along with the predicted classification of each cell in either G2M, S or G1 phase
    scrna <- CellCycleScoring(object = scrna, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)
  
  # view cell cycle scores and phase assignments
    head(x = scrna@meta.data)
  
  #Scaling the data and removing unwanted sources of variation
    scrna <- ScaleData(object = scrna, vars.to.regress = c("nUMI", "percent.mito","S.Score", "G2M.Score"))
}else{
  scrna <- ScaleData(object = scrna, vars.to.regress = c("nUMI", "percent.mito"))
}


#Perform linear dimensional reduction (Note: performed on the variable genes)
scrna <- RunPCA(object = scrna, pc.genes = scrna@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5)

  pdf(file=paste(dir,"/plots/",project,"_vizplot.pdf",sep=""),height = 7,width = 11)
  VizPCA(object = scrna, pcs.use = 1:2)
  dev.off()

  pdf(file=paste(dir,"/plots/",project,"_pcaplot.pdf",sep=""),height = 7,width = 11)
  PCAPlot(object = scrna, dim.1 = 1, dim.2 = 2)
  dev.off()

  pdf(file=paste(dir,"/plots/",project,"_pcheatmap.pdf",sep=""),height = 7,width = 11)
  PCHeatmap(object = scrna, pc.use = 1:15, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
  dev.off()

#Determine statistically significant principal components by randomly permuting a subset of the data (1% by default) and rerunning PCA
scrna <- JackStraw(object = scrna, num.replicate = 100, do.print = FALSE)

  pdf(file=paste(dir,"/plots/",project,"_jackstrawplot.pdf",sep=""),height = 7,width = 11)
  JackStrawPlot(object = scrna, PCs = 1:14)
  dev.off()

  pdf(file=paste(dir,"/plots/",project,"_pcelbow.pdf",sep=""),height = 7,width = 11)
  PCElbowPlot(object = scrna)
  dev.off()

#Clustering
scrna <- FindClusters(object = scrna, reduction.type = "pca", dims.use = 1:15, resolution = 0.6, print.output = 0, save.SNN = TRUE)

PrintFindClustersParams(object = scrna)

#Run tsne
scrna <- RunTSNE(object = scrna, dims.use = 1:15, do.fast = TRUE)

# scrna_eset <- newCellDataSet(scrna@scale.data,
#                        phenoData =scrna@meta.data, featureData = fd)


  pdf(file=paste(dir,"/plots/",project,"_TSNEplot.pdf",sep=""),height = 9,width = 9)
  TSNEPlot(object = scrna,pt.size=2,group.by = "ident")
  dev.off()
  
save(scrna,file=paste(dir,"/",project,".RData",sep=""))
#TSNEPlot(object = scrna,pt.size=2,group.by = "celltype",colors.use = c("coral2","lightgray","turquoise3","mediumorchid"))

#Create Feature Plot of genes of choice
# FeaturePlot(object = scrna, features.plot = c("Hopx","Sftpc","Sox2","Sox9"), cols.use = c("grey", "blue"), 
#             reduction.use = "tsne")

#find cluster markers
# cluster1_2.markers <- FindMarkers(object = scrna, ident.1 = 1,ident.2 =2, min.pct = -Inf)
# print(x = head(x = cluster1_2.markers, n = 5))
# cluster1_2.markers=cluster1_2.markers[order(-cluster1_2.markers$avg_logFC),]
# 
# cluster1_all.markers <- FindMarkers(object = scrna, ident.1 = 1,ident.2 = c(2,3), min.pct = 0.25)
# print(x = head(x = cluster1_all.markers, n = 5))
# 
# #write output
# result=list(AT1_vs_AT2=cluster1_2.markers,AT1_vs_all=cluster1_all.markers)
# 
# write.xlsx(result,file="results.xlsx",row.names=T)
# 


