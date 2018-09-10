library(Seurat)
library(dplyr)
library(readr)
cpallette=c("#64B2CE", "#DA5724", "#74D944", "#CE50CA", "#C0717C", "#CBD588", "#5F7FC7", 
            "#673770", "#D3D93E", "#8569D5", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
            "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
            "#8A7C64", "#599861")
#Specify input
dir='sampleD' #specify directory
name='MB_human_lung_SampleD' # specify project name. 
files=c("sampleD/rep1/filtered_gene_bc_matrices/GRCh38",
          "sampleD/rep2/filtered_gene_bc_matrices/GRCh38")
org='human'
maxdim = 75 # Maximum number of dimensions (principle components) to compute
mindim = 10 # Minimum number of dimensions (principle components) to compute
maxres = 1.2 # Maximum number of resolutions to use to find clusters
minres = 0.4 # Minimum number of resolutions to use to find clusters
markergenes=NA #If you have a list of markergenes, enter here as a character vector or leave it as NA. Keep in mind the organism and use the right genenames
#markergenes <- c("PDGFRA","DCN","ACTA2","CD34","COL3A1","MYH11","ITGA8","PDGFRB","DES")
#mouseorthologfile <- '~/dsdata/NGSshare/homologs/mouse_human.csv'
mouseorthologfile <- '~/NGSshare/homologs/mouse_human.csv'
paramsweep=F

#Note: Occasionally the program might throw an error if one or multiple marker genes entered is not present in the dataset (Probably filtered out due to low expression).
# In that case, check genename for typos (also case as the program is case-sensitive) or else remove it from the list. 
#create directory to save the plots  

dir.create(paste0(dir,"/seurat/plots",sep=""),recursive = T)


##################### Create funtion to process the data ##################################
processExper <- function(dir,name,ccscale=T,org,files){
  try(if(length(files)==0) stop("No files"))
  
  if(length(files)==1){
    # Load the dataset
    inputdata <- Read10X(data.dir =files[1])
    # Initialize the Seurat object with the raw (non-normalized data).  
    scrna <- CreateSeuratObject(raw.data = inputdata, min.cells = 10, min.genes = 200,project = name)
  }else{

    #Initialize the first object with the raw (non-normalized data) and add rest of the data 
    inputdata <- Read10X(data.dir =files[1])
    scrna <- CreateSeuratObject(raw.data = inputdata, min.cells = 10, min.genes = 200, project = 'Rep1')
    cat('Rep1', length(scrna@cell.names), "\n")
    for(i in 2:length(files)){
      tmp.data <- Read10X(data.dir =files[i])
      tmp.scrna <- CreateSeuratObject(raw.data = tmp.data, min.cells = 10, min.genes = 200, project = paste0('Rep',i))
      cat('Rep', i, ": ", length(tmp.scrna@cell.names), "\n", sep="")
      scrna <- MergeSeurat(scrna, tmp.scrna, do.normalize = FALSE, min.cells = 0, min.genes = 0, add.cell.id2 = paste0('Rep',i))
    }
    cat("merged: ", length(scrna@cell.names), "\n", sep="")
  }
  
  # calculate the percent.mito values.
  if(org=='mouse'){
    mito.genes <- grep(pattern = "^mt-", x = rownames(x = scrna@data), value = TRUE)
  }else{
    mito.genes <- grep(pattern = "^MT-", x = rownames(x = scrna@data), value = TRUE)
  }
  percent.mito <- Matrix::colSums(scrna@raw.data[mito.genes, ])/Matrix::colSums(scrna@raw.data)
  
  # AddMetaData adds columns to object@meta.data. metadata is a great place to stash QC stats
  scrna <- AddMetaData(object = scrna, metadata = percent.mito, col.name = "percent.mito")
  
  pdf(file=paste0(dir,"/seurat/plots/",name,"_violinplot.pdf"),height = 7,width = 11)
  VlnPlot(object = scrna, features.plot = c("nGene", "nUMI", "percent.mito"),do.return = T)
  dev.off()
  
  # Plot GenePlot
  pdf(file=paste(dir,"/seurat/plots/",name,"_geneplot.pdf",sep=""),height = 7,width = 11)
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
      m2h <- read_csv(mouseorthologfile)
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



scrna <- processExper(dir=dir,name=name,ccscale=T,org=org,files=files)


#Perform linear dimensional reduction (Note: performed on the variable genes)
scrna <- RunPCA(object = scrna, pc.genes = scrna@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5,pcs.compute=maxdim)

pdf(file=paste(dir,"/seurat/plots/",name,"_vizplot.pdf",sep=""),height = 7,width = 11)
VizPCA(object = scrna, pcs.use = 1:2)
dev.off()

pdf(file=paste(dir,"/seurat/plots/",name,"_pcaplot.pdf",sep=""),height = 7,width = 11)
PCAPlot(object = scrna, dim.1 = 1, dim.2 = 2)
dev.off()

pdf(file=paste(dir,"/seurat/plots/",name,"_pcheatmap_1_15.pdf",sep=""),height = 16,width = 11)
PCHeatmap(object = scrna, pc.use = 1:25, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
dev.off()

pdf(file=paste(dir,"/seurat/plots/",name,"_pcheatmap_26_50.pdf",sep=""),height = 16,width = 11)
PCHeatmap(object = scrna, pc.use = 26:50, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
dev.off()
pdf(file=paste(dir,"/seurat/plots/",name,"_pcheatmap_51_75.pdf",sep=""),height = 16,width = 11)
PCHeatmap(object = scrna, pc.use = 51:75, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
dev.off()


#Determine statistically significant principal components by randomly permuting a subset of the data (1% by default) and rerunning PCA
scrna <- JackStraw(object = scrna, num.replicate = 100,num.pc=maxdim)

pdf(file=paste(dir,"/seurat/plots/",name,"_jackstrawplot.pdf",sep=""),height = 20,width = 11)
scrna <- JackStraw(object = scrna, num.replicate = 100,num.pc=maxdim,do.par=T)
dev.off()

pdf(file=paste(dir,"/seurat/plots/",name,"_pcelbow.pdf",sep=""),height = 7,width = 11)
PCElbowPlot(object = scrna,num.pc = maxdim)
dev.off()

#Save the pre clsutered and dim reduced object, we will load this after we decided on which params 

#save(scrna,file=paste0(dir,"/seurat/",name,'step1_.RData'))

#Run tSNE,uMAP and Clustering in a loop and save the results in seperate folders within the name directory
if (paramsweep==T){
  for(dim in seq(mindim,maxdim,2)){
    #Create dir for each dimension
    dir.create(paste0(dir,"/seurat/",dim,'/'))
    
    #Run tsne and umap
    scrna <- RunTSNE(object = scrna, dims.use = 1:dim, do.fast = TRUE)
    scrna <- RunUMAP(object = scrna, dims.use = 1:dim, min_dist=0.5,n_neighbors = 15,metric = 'correlation')
    
    
    for(res in seq(minres,maxres,0.2)){
      scrna <- FindClusters(scrna, reduction.type = "pca",dims.use = 1:dim, resolution = res)
      p1 <- plot_grid(
        TSNEPlot(scrna, do.label = T,do.return=T,colors.use = cpallette),
        DimPlot(scrna, reduction.use = "umap", do.label = T,do.return=T,cols.use = cpallette)
      )
      #dev.off()
      save_plot(paste0(dir,"/seurat/",dim,'/tSNE_umap_dim',dim,'_res',res,'.pdf'), p1,base_height=8,base_width=12)
      
      
    }
    
    if(is.na(markergenes)==FALSE){
      
      fp_tsne <- FeaturePlot(object = scrna, features.plot = markergenes, cols.use = c("lightgrey", "blue"),do.return = T)
      save_plot(paste0(dir,"/seurat/",dim,'/TSNE_dim',dim,'_markergenes.pdf'),plot_grid(plotlist = fp_tsne),base_width=16,base_height=16)
      fp_umap <- FeaturePlot(object = scrna, features.plot = markergenes, cols.use = c("lightgrey", "blue"),do.return = T,reduction.use = "umap")
      save_plot(paste0(dir,"/seurat/",dim,'/UMAP_dim',dim,'_markergenes.pdf'),plot_grid(plotlist = fp_umap),base_width=16,base_height=16)
      
    }
    
  }
  quit()
}

stop()
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
dim=39
res=.6
addcelltype="no" #choose between yes or no

#Load the right data
load(file=paste0(dir,"/seurat/",name,'step1_.RData'))

scrna <- RunTSNE(object = scrna, dims.use = 1:dim, do.fast = TRUE)
scrna <- RunUMAP(object = scrna, dims.use = 1:dim, min_dist=0.5,n_neighbors = 15,metric = 'correlation')
scrna <- RunDiffusion(object = scrna, dims.use = 1:dim)

scrna <- FindClusters(scrna, reduction.type = "pca",dims.use = 1:dim, resolution = res)


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
  scrna <- AddMetaData(scrna, scrna@ident,col.name= "var_cluster")
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
save(scrna,file=paste(dir,"/",name,".RData",sep=""))
