#install.packages("Seurat")
library(Seurat)
library(dplyr)
library(Matrix)
library(openxlsx)

#Specify input
dir='DF_E15.5_Nkx2.1' #specify directory
project='DF_E15.5_Nkx2.1' # specify project name. 
maxdim = 40 # Maximum number of dimensions (principle components) to compute
mindim = 4 # Minimum number of dimensions (principle components) to compute
maxres = 1.6 # Maximum number of resolutions to use to find clusters
minres = 0.4 # Minimum number of resolutions to use to find clusters
markergenes=NA #If you have a list of markergenes, enter here as a character vector or leave it as NA. Keep in mind the organism and use the right genenames
#markergenes <- c("NKX2-5","SOX2","SOX17")

#Note: Occasionally the program might throw an error if one or multiple marker genes entered is not present in the dataset (Probably filtered out due to low expression).
# In that case, check genename for typos (also case as the program is case-sensitive) or else remove it from the list. 


##################### Create funtion to process the data ##################################
processExper <- function(dir,name,ccscale=F,type='single',org='human'){
#create directory to save the plots  
dir.create(paste0(dir,"/plots",sep=""),recursive = T)
  
# Load the dataset
if(type=="single"){
inputdata <- Read10X(data.dir = paste(dir,"/outs/filtered_gene_bc_matrices/mm10/",sep=""))
}else{
  inputdata <- Read10X(data.dir = paste(dir,"/outs/filtered_gene_bc_matrices_mex/mm10/",sep=""))
}

# Initialize the Seurat object with the raw (non-normalized data).  
scrna <- CreateSeuratObject(raw.data = inputdata, min.cells = 3, min.genes = 200,project = name)

# calculate the percent.mito values.
if(org='mouse'){
mito.genes <- grep(pattern = "^mt-", x = rownames(x = scrna@data), value = TRUE)
}else{
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = scrna@data), value = TRUE)
}
percent.mito <- Matrix::colSums(scrna@raw.data[mito.genes, ])/Matrix::colSums(scrna@raw.data)

# AddMetaData adds columns to object@meta.data. metadata is a great place to stash QC stats
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

#normalize data
scrna <- NormalizeData(object = scrna, normalization.method = "LogNormalize",scale.factor = 10000)

#detection of variable genes
#calculates the average expression and dispersion for each gene, places these genes into bins, and then calculates a z-score for dispersion within each bin
scrna <- FindVariableGenes(object = scrna, mean.function = ExpMean, dispersion.function = LogVMR,x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
length(x = scrna@var.genes)

if(ccscale==T){
  if(org=='human'){
  #Assign scores in the CellCycleScoring function.Stores S and G2/M scores in object@meta.data, along with the predicted classification of each cell in either G2M, S or G1 phase
  scrna <- CellCycleScoring(object = scrna, s.genes = cc.genes$s.genes, g2m.genes = cc.genes$g2m.genes)
  }else{
    m2h <- read_csv('~/NGSshare/homologs/mouse_human.csv')
    cc.genes$s.genes <- m2h %>% filter(human_name %in% cc.genes$s.genes) %>% pull(mouse_name)
    cc.genes$g2m.genes <- m2h %>% filter(human_name %in% cc.genes$g2m.genes) %>% pull(mouse_name)
    scrna <- CellCycleScoring(object = scrna, s.genes = cc.genes$s.genes, g2m.genes = cc.genes$g2m.genes)
  }
  #Scaling the data and removing unwanted sources of variation
    scrna <- ScaleData(object = scrna, vars.to.regress = c("nUMI", "percent.mito","S.Score", "G2M.Score"))
}else{
  scrna <- ScaleData(object = scrna, vars.to.regress = c("nUMI", "percent.mito"))
}
}

#Setup seurat objects for each of the single cell experiments
################################################################################################################################################
##  dir - directory where you have your input files                                                                                           ##
##  name -project name                                                                                                                        ## 
##  type - choose either single or aggregate. Default is single                                                                               ##
##  ccscale - choose either TRUE (T) ot FALSE(F) based on whether or not you want to regress out cell cycle. Default is FALSE                 ## 
##  org - specify if mouse/human. default is human                                                                                            ##
##  mergedata - choose either TRUE (T) ot FALSE(F) based on whether or not you want to merge datasets. If true, specify                       ##
################################################################################################################################################
scrna <- processExper(dir='DF_E15.5_Nkx2.1',name='DF_E15.5_Nkx2.1',ccscale=T,type="single",org="mouse")


#Perform linear dimensional reduction (Note: performed on the variable genes)
scrna <- RunPCA(object = scrna, pc.genes = scrna@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5,pcs.compute=maxdim)

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
scrna <- JackStraw(object = scrna, num.replicate = 100, do.print = FALSE,num.pc=maxdim)

  pdf(file=paste(dir,"/plots/",project,"_jackstrawplot.pdf",sep=""),height = 7,width = 11)
  JackStrawPlot(object = scrna, PCs = 1:14)
  dev.off()

  pdf(file=paste(dir,"/plots/",project,"_pcelbow.pdf",sep=""),height = 7,width = 11)
  PCElbowPlot(object = scrna)
  dev.off()

#Run tSNE,uMAP and Clustering in a loop and save the results in seperate folders within the project directory
  
  for(dim in seq(mindim,maxdim,2)){
    #Create dir for each dimension
    dir.create(paste0(dir,"/",dim,'/'))
    
    #Run tsne and umap
    scrna <- RunTSNE(object = scrna, dims.use = 1:dim, do.fast = TRUE)
    scrna <- RunUMAP(object = scrna, dims.use = 1:dim)
    
    
    for(res in seq(minres,maxres,0.2)){
      scrna <- FindClusters(scrna, reduction.type = "pca",dims.use = 1:dim, save.SNN = T, resolution = res)
      
      # Visualization make 2 plots.. TSNE and UMAP 
      plot_grid(
        TSNEPlot(scrna, do.label = T,do.return=T),
        DimPlot(scrna, reduction.use = "umap", do.label = T,do.return=T)
      )
      
      ggsave(paste0(dir,"/",dim,'/tSNE_umap_dim',dim,'_res',res,'.png'), width=15,height=9)
    }
    
    if(is.na(markergenes)==FALSE){
      fp <- FeaturePlot(object = scrna, features.plot = markergenes, cols.use = c("lightgrey", "blue"),do.return = T)
      cowplot::ggsave(paste0(dir,"/",dim,'/TSNE_dim',dim,'_markergenes.png'),plot_grid(plotlist = fp),width=16,height=16)
      
      fp_umap <- FeaturePlot(object = scrna, features.plot = markergenes, cols.use = c("lightgrey", "blue"),do.return = T,reduction.use = "umap")
      cowplot::ggsave(paste0(dir,"/",dim,'/UMAP_dim',dim,'_markergenes.png'),plot_grid(plotlist = fp),width=16,height=16)
    }
    save(scrna,file=paste0(dir,"/",dim,"/",project,"_dim",dim,'.RData'))
  }

              ################################################## STOP HERE ##############################################
  ################################################################################################################################################
  ################################################################################################################################################
  ##  Look at the tsne and umap plots. The resolutions vary between minimum and maximum resolutions specified incremented by 0.2.               ##
  ##  There will be one plot for each resolution                                                                                                ## 
  ##  The dimensions vary between minimum and maximum dimension specified incremented by 2. The above code will generate a dir for each         ##
  ##  dimension and within each directory will be n number of plots of varying resolution and a seurat object. Go through the plots and select  ## 
  ##  the dimension and resolution that makes the most sense and then proceed.                                                                  ##
  ##  If you have entered markergenes, the Featureplots will show the gene expresssion using both tSNE and umap. If you know the cell type,     ##
  ##  you can add the cell type as a variable as well and run FindMarkers fuction.Defaults to cluster.                                          ##
  ################################################################################################################################################
  ################################################################################################################################################

  #specify dim and res
  dim=22
  res=1.2
  addcelltype="yes" #choose between yes or no

  #Load the right data
  load(paste0(dir,"/",dim,"/",project,"_dim",dim,'.RData'))
  
  #Rename the right resolution colum in the meta data to var_cluster
  scrna@meta.data=scrna@meta.data %>% rename("var_cluster"=paste0("res.",res))

#Add celltypes if you have the information 
  #Create and excel sheet names celltypes.csv with the first column having the cluster id and second column having the celltype
  #Note: Every cluster should have a celltype associated with it. If there isnt any, use "novel_celltype_<num>"
if(addcelltype=="yes"){
  celltype=as.data.frame(read.csv("celltype.csv"))
  colnames(celltype)=c("clust","var_celltype")
  scrna@meta.data=left_join(scrna@meta.data,celltype,by=c("var_cluster"="clust"))
  scrna <- SetAllIdent(object = scrna, id = "var_celltype")
}else{
  #Else set the cluster of selected resolution as the identity/main group of comparison
  scrna <- SetAllIdent(object = scrna, id = var_cluster)
  }

  #For each group, run find markers in a loop and save it in the seurat object
  scrna@misc=NA
  scrna@misc <-  vector(mode="list", length=length(levels(scrna@ident)))
  names(scrna@misc)=levels(scrna@ident)
  for(c in levels(scrna@ident)){
    scrna@misc[[c]] <- FindMarkers(scrna,ident.1 = c) %>% tibble::rownames_to_column('gene_name')
    rownames(scrna@misc[[c]])=scrna@misc[[c]]$gene_name
  } 

#Save Seurat object  
save(scrna,file=paste(dir,"/",project,".RData",sep=""))
