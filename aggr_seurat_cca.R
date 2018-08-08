library(Seurat)
library(dplyr)
library(Matrix)
library(cowplot)
library(readr)

#Specify input
dir='CS_Endoflu' #specify directory
projectdir=paste(dir,'/seurat_cca/')
project='CS_LungEndo' # specify project name. 
maxdim = 40 # Maximum number of dimensions (principle components) to compute
mindim = 4 # Minimum number of dimensions (principle components) to compute
maxres = 1.6 # Maximum number of resolutions to use to find clusters
minres = 0.4 # Minimum number of resolutions to use to find clusters
#If you have a list of markergenes, enter here as a character vector or leave it as NA. Keep in mind the organism and use the right genenames
markergenes <- c('Pecam1','Nkx2-1','Acta2','Wt1','Ptprc','Hba-a1','Epcam','Wnt2','Cd86','Gypa','Sftpc','Pdgfrb','Hopx','Pdgfra')			

#Note: Occasionally the program might throw an error if one or multiple marker genes entered is not present in the dataset (Probably filtered out due to low expression).
# In that case, check genename for typos (also case as the program is case-sensitive) or else remove it from the list. 

################################################################################################
#     Map the cc.genes from human to mouse. 
################################################################################################
m2h <- read_csv('~/NGSshare/homologs/mouse_human.csv')
cc.genes$s.genes <- m2h %>% filter(human_name %in% cc.genes$s.genes) %>% pull(mouse_name)
cc.genes$g2m.genes <- m2h %>% filter(human_name %in% cc.genes$g2m.genes) %>% pull(mouse_name)
################################################################################################

 #Create project dir
dir.create(projectdir,recursive = T)
#create directory to save the plots  
dir.create(paste0(dir,"/plots",sep=""),recursive = T)

#create function to set up seurat objects
processExper <- function(dir,name,ccscale=F,org='human'){
  #Read in input data
  a.raw <- Read10X(data.dir = dir)
  dim(a.raw)
  colnames(a.raw) <- paste0(colnames(a.raw), name)
  
  #Create seurat object
  a <- CreateSeuratObject(raw.data = a.raw, project = name)
  a@meta.data$var_sample <- name
  
  #Calculate mito-gene percentage
  if(org='mouse'){
    mito.genes <- grep(pattern = "^mt-", x = rownames(x = scrna@data), value = TRUE)
  }else{
    mito.genes <- grep(pattern = "^MT-", x = rownames(x = scrna@data), value = TRUE)
  }
  percent.mito <- Matrix::colSums(a@raw.data[mito.genes, ])/Matrix::colSums(a@raw.data)
  a <- AddMetaData(object = a, metadata = percent.mito, col.name = "percent.mito")
  a@meta.data %>% filter(percent.mito > 0.05) %>% summarise(n=n())
  a@meta.data %>% filter(nGene < 1000) %>% summarise(n=n())
  
  #Save QC plots
  png(paste0(dir,"/plots/",name,'_QC_scatterplot.png'))
  par(mfrow = c(1, 2))
  GenePlot(object = a, gene1 = "nUMI", gene2 = "percent.mito")
  GenePlot(object = a, gene1 = "nUMI", gene2 = "nGene")
  dev.off()
  
  vln=VlnPlot(object = a, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3,do.return=TRUE, return.plotlist = T)
  cowplot::ggsave(paste0(dir,"/plots/",name,'_QC_vlnplot.png'),plot_grid(plotlist = vln),width=16,height=16)
  
  #Filter cells
  a <- FilterCells(a, subset.names = c("nGene","percent.mito"), low.thresholds = c(200,-Inf), 
                   high.thresholds = c(Inf,0.05))
  #normalize data
  a <- NormalizeData(a)
  
  #Find the most variable genes
  a <- FindVariableGenes(a, do.plot = F)
  
  #Scaling the data and removing unwanted sources of variation
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
    
    scrna <- ScaleData(object = scrna, vars.to.regress = c("nUMI", "percent.mito","S.Score", "G2M.Score"))
  }else{
    scrna <- ScaleData(object = scrna, vars.to.regress = c("nUMI", "percent.mito"))
  }
  
  return(a)
  
}

#Setup seurat objects for each of the single cell experiments
################################################################################################################################################
##  dir - directory where you have your input files                                                                                           ##
##  name -project name                                                                                                                        ## 
##  ccscale - choose either TRUE (T) ot FALSE(F) based on whether or not you want to regress out cell cycle. Default is FALSE                 ## 
##  org - specify if mouse/human. default is human                                                                                            ##
################################################################################################################################################

ctrl <- processExper(dir="data/control/mm10/",name='ctrl',ccscale=F,org='mouse')
flu2 <- processExper(dir="data/Flu2/mm10/",name='Flu2',ccscale=F,org='mouse')
flu5 <- processExper(dir="data/Flu5/mm10/",name='Flu5',ccscale=F,org='mouse')

#create a list of the seurat objects for each expt
ob.list <- list(ctrl,flu2,flu5)

#Gene selection for input to cca
#select top 1000 variable genes from each expt
genes.use <- c()
for (i in 1:length(ob.list)) {
  ob.list[[i]] = FindVariableGenes(ob.list[[i]], do.plot = F)
  genes.use <- c(genes.use, head(rownames(ob.list[[i]]@hvg.info), 1000))
}
#retain genes that appear more than once and are present in the scaled data object
genes.use <- names(which(table(genes.use) > 1))
for (i in 1:length(ob.list)) {
  genes.use <- genes.use[genes.use %in% rownames(ob.list[[i]]@scale.data)]
}

#perform CCA to find common sources of variation between the datasets
combined <- RunMultiCCA(ob.list, genes.use = genes.use, num.ccs = maxdim)

p1 <- DimPlot(object = combined, reduction.use = "cca", group.by = "var_sample", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = combined, features.plot = "CC1", group.by = "var_sample", 
              do.return = TRUE,pt.size = "NA")
pg=plot_grid(p1, p2)
ggsave(pg,filename=paste0(dir,'/plots/',"dimplot_vln.pdf"),height=9,width=12,units = "in")

#create plot to measure of correlation strength for each CC
p <- MetageneBicorPlot(combined, grouping.var = "var_sample", dims.eval = 1:50)
cowplot::ggsave(paste0(dir,'/plots/','MetageneBiCorr.png'), width=16,height=8)


for(dim in seq(mindim,maxdim,2)){

  # Alignment
  dir.create(paste0(projectdir,dim,'/'))
  
  combined <- CalcVarExpRatio(object = combined, reduction.type = "pca",grouping.var = "var_sample", dims.use = 2:dim)
  
  scrna <- SubsetData(combined, subset.name = "var.ratio.pca",accept.low = 0.5)
  
  scrna <- AlignSubspace(scrna,reduction.type = "cca",grouping.var = "var_sample",dims.align = 1:dim)
  
  # t-SNE and Clustering
  scrna <- RunTSNE(scrna,reduction.use = "cca.aligned",dims.use = 1:dim,check_duplicates = FALSE,nthreads = 10,max_iter = 2000)
  
  ### run UMAP
  scrna <- RunUMAP(scrna, reduction.use = "cca.aligned", dims.use = 1:dim)
  
  ### Create a new loop to loop over min to max resolutions incremented by 0.2
  for(res in seq(minres,maxres,0.2)){
    scrna <- FindClusters(scrna, reduction.type = "cca.aligned",dims.use = 1:dim, save.SNN = T, resolution = res)
    
    # Visualization make 4 plots.. TSNE with cluster, TSNE with var_sample, UMAP cluster nd var_sample
    
    plot_grid(
      TSNEPlot(scrna, do.label = T,do.return=T),
      DimPlot(scrna,reduction.use = "tsne",group.by='var_sample',do.return = T),
      DimPlot(scrna, reduction.use = "umap", do.label = T,do.return=T),
      DimPlot(scrna, reduction.use = "umap",group.by='var_sample', do.label = T,do.return=T)
    )
    
    
    ggsave(paste0(projectdir,"/",dim,'/TSNE_dim',dim,'_res',res,'.png'), width=16,height=13)
  }
  
  #IF marker genes are specified, plot the gene expression of the marker genes into a feature plot
  if(is.na(markergenes)==FALSE){
    fp <- FeaturePlot(object = scrna, features.plot = markergenes, cols.use = c("lightgrey", "blue"),do.return = T)
    cowplot::ggsave(paste0(projectdir,"/",dim,'/TSNE_dim',dim,'_markergenes.png'),plot_grid(plotlist = fp),width=16,height=16)
    
    fp_umap <- FeaturePlot(object = scrna, features.plot = markergenes, cols.use = c("lightgrey", "blue"),do.return = T,reduction.use = "umap")
    cowplot::ggsave(paste0(projectdir,"/",dim,'/UMAP_dim',dim,'_markergenes.png'),plot_grid(plotlist = fp),width=16,height=16)
  }
  
  save(scrna,file=paste0(projectdir,"/",dim,"/",project,"_dim",dim,'.RData'))
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
load(paste0(projectdir,"/",dim,"/",project,"_dim",dim,'.RData'))

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
