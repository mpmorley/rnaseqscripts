library(Seurat)
library(dplyr)
library(readr)
#Pallette for nicer looking plots
cpallette=c("#64B2CE", "#DA5724", "#74D944", "#CE50CA", "#C0717C", "#CBD588", "#5F7FC7", 
            "#673770", "#D3D93E", "#8569D5", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", 
            "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", 
            "#8A7C64", "#599861")
### Specify the name of the outpur directory. 
outdir='P15'
## Spefiicy the name of the analysis 
projectname='JZ_lung_P12' # specify project name. 

## Currently this script is coded to use 10x genomics output. 
files=c("P15/mm10/")

# Currently only mouse and human 
org='mouse'

#Set the min and max of PCA dimensions to test. 
maxdim = 80 
mindim = 10 
#If you have a list of markergenes, enter here as a character vector or leave it as NA. Keep in mind the organism and use the right genenames
# example is markergenes <- c("PDGFRA","DCN","ACTA2","CD34","COL3A1","MYH11","ITGA8","PDGFRB","DES")
markergenes=NA 


## In order for the cellc cycle calculation to work with mouse, we need an file with mouse/human orthologs
mouseorthologfile <- '~/dsdata/NGSshare/homologs/mouse_human.csv'


#Note: Occasionally the program might throw an error if one or multiple marker genes entered is not present in the dataset (Probably filtered out due to low expression).
# In that case, check genename for typos (also case as the program is case-sensitive) or else remove it from the list. 
#create directory to save the plots  
plotdir=paste0(outdir,"/plots",sep="")
dir.create(plotdir,recursive = T)

################################################################################################################################################
##  dir - directory where you have your input files                                                                                           ##
##  name -project name                                                                                                                        ## 
##  type - choose either single or aggregate. Default is single                                                                               ##
##  ccscale - choose either TRUE (T) ot FALSE(F) based on whether or not you want to regress out cell cycle. Default is FALSE                 ## 
##  org - specify if mouse/human. default is human                                                                                            ##
##  mergedata - choose either TRUE (T) ot FALSE(F) based on whether or not you want to merge datasets. If true, specify                       ##
################################################################################################################################################
processExper <- function(dir,name,ccscale=T,org,files,filtergenes=NULL){
  try(if(length(files)==0) stop("No files"))
  
  if(length(files)==1){
    # Load the dataset
    inputdata <- Read10X(data.dir =files[1])
    colnames( inputdata) <- paste0(colnames(inputdata), '_',name)
    
    # Initialize the Seurat object with the raw (non-normalized data).  
    scrna <- CreateSeuratObject(raw.data = inputdata, min.cells = 10, min.genes = 200,project = name)
  }else{
    
    #Initialize the first object with the raw (non-normalized data) and add rest of the data 
    inputdata <- Read10X(data.dir =files[1])
    scrna <- CreateSeuratObject(raw.data = inputdata, min.cells = 10, min.genes = 200, project = name)
    cat('Rep1', length(scrna@cell.names), "\n")
    for(i in 2:length(files)){
      tmp.data <- Read10X(data.dir =files[i])
      tmp.scrna <- CreateSeuratObject(raw.data = tmp.data, min.cells = 10, min.genes = 200, project = name)
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
  
  pdf(file=paste0(plotdir,"/violinplot.pdf"),height = 7,width = 11)
  VlnPlot(object = scrna, features.plot = c("nGene", "nUMI", "percent.mito"),do.return = T)
  dev.off()
  
  # Plot GenePlot
  pdf(file=paste(plotdir,"/geneplot.pdf",sep=""),height = 7,width = 11)
  par(mfrow = c(1, 2))
  GenePlot(object = scrna, gene1 = "nUMI", gene2 = "percent.mito")
  GenePlot(object = scrna, gene1 = "nUMI", gene2 = "nGene")
  dev.off()
  
  # Filter out cells that have unique gene counts less than 500 
  scrna <- FilterCells(object = scrna, subset.names = c("nGene", "percent.mito"), 
                       low.thresholds = c(500, -Inf), high.thresholds = c(Inf, 0.05))
  
  
  if(length(filtergenes)>0){
    for(gene in filtergenes){
      try(
      scrna <- FilterCells(object = scrna, subset.names = gene, 
                         low.thresholds = c(-Inf), high.thresholds = c(1))
      )
    }
  }
  
  
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



## Preprocess the data and create a seurat object. If using he view website the object has to be named scrna. 

scrna <- processExper(dir=outdir,name=projectname,ccscale=T,org=org,files=files)


#Perform linear dimensional reduction (Note: performed on the variable genes)
scrna <- RunPCA(object = scrna, pc.genes = scrna@var.genes, do.print = TRUE, pcs.print = 1:5, genes.print = 5,pcs.compute=maxdim)

pdf(file=paste(plotdir,"/vizplot.pdf",sep=""),height = 7,width = 11)
VizPCA(object = scrna, pcs.use = 1:2)
dev.off()

pdf(file=paste(plotdir,"/pcaplot.pdf",sep=""),height = 7,width = 11)
PCAPlot(object = scrna, dim.1 = 1, dim.2 = 2)
dev.off()

#can set number of cores to speed it up
scrna <- JackStraw(object = scrna, num.replicate = 100,num.pc=maxdim,num.cores = 3)


for(l in seq(0,(maxdim %/%20-1)*20,20)+1){
  h=ifelse(l+19>=maxdim,maxdim,l+19)
  pdf(file=paste(plotdir,"/pcheatmap_",l,"_",h,".pdf",sep=""),height = 18,width = 11)
  PCHeatmap(object = scrna, pc.use = l:h, cells.use = 500, do.balanced = TRUE, label.columns = FALSE, use.full = FALSE)
  dev.off()
}

pdf(file=paste(plotdir,"/jackstrawplot.pdf",sep=""),height = 20,width = 11)
JackStrawPlot(object = scrna, PCs = 1:maxdim)
dev.off()


pdf(file=paste(plotdir,"/pcelbow.pdf",sep=""),height = 7,width = 11)
PCElbowPlot(object = scrna,num.pc = maxdim)
dev.off()



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

try(if(exists("dim")) stop("Dim isn't define, please set the number of dims"))

dim=
addcelltype="no" #choose between yes or no



### Loop over a few resolution choices, 0.6 is last, it's the default and typically works well. 
for(res in c(0.4,0.8,1.0,1.2,0.6)){
  scrna <- FindClusters(scrna, reduction.type = "pca",dims.use = 1:dim, resolution = res)
  scrna <- AddMetaData(scrna, scrna@ident,col.name= paste0("var_cluster_altres_",res))
}

scrna@meta.data <-select(scrna@meta.data, -starts_with("res"))
scrna <- AddMetaData(scrna, scrna@ident,col.name= "var_cluster")


### Run various Dimension reduction techinques, Diffusion Maps takes while so I typically do not run it. 
scrna <- RunTSNE(object = scrna, dims.use = 1:dim, do.fast = TRUE)
scrna <- RunUMAP(object = scrna, dims.use = 1:dim, min_dist=0.5,n_neighbors = 30,metric = 'correlation')
#scrna <- RunDiffusion(object = scrna, dims.use = 1:dim)

 
###################################################################
# Look over the UMAP and make sure params used make a decent plot, 
# you can change the n_neihbor and min_dist param 
# 
###################################################  
DimPlot(scrna,reduction.use = 'umap',cols.use = cpallette,group.by = 'var_cluster_altres_0.4')

  
####################################################################  
  
  
#Add celltypes if you have the information 
#Create and excel sheet names celltypes.csv with the first column having the cluster id and second column having the celltype
#Note: Every cluster should have a celltype associated with it. If there isnt any, use "novel_celltype_<num>"
if(addcelltype=="yes"){
  celltype=as.data.frame(read.csv("celltype.csv"))
  colnames(celltype)=c("clust","var_celltype")
  scrna@meta.data=left_join(scrna@meta.data,celltype,by=c("var_cluster"="clust"))
  scrna <- SetAllIdent(object = scrna, id = "var_celltype")
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
save(scrna,file=paste0(outdir,"/",projectname,".RData"))



    
