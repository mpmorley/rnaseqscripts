library(Seurat)
library(dplyr)
library(Matrix)
library(cowplot)
library(readr)

setwd("/Users/bapoorva/Desktop/Server/Morrisey/Colin/scRNA/EndoFlu")
#################################################################################
#     Need to map the cc.genes from human to mouse. 
#
#################################################################################
m2h <- read_csv('~/NGSshare/homologs/mouse_human.csv')
cc.genes$s.genes <- m2h %>% filter(human_name %in% cc.genes$s.genes) %>% pull(mouse_name)
cc.genes$g2m.genes <- m2h %>% filter(human_name %in% cc.genes$g2m.genes) %>% pull(mouse_name)
#################################################################################

 #Create project dir
projectdir='aggr/seurat_cca/'
dir.create(projectdir,recursive = T)
markergenes <- c('Pecam1','Nkx2-1','Acta2','Wt1','Ptprc','Hba-a1','Epcam','Wnt2','Cd86','Gypa','Sftpc','Pdgfrb','Hopx','Pdgfra')			

#create function to set up seurat objects
processExper <- function(dir,name,ccscale=F){
  a.raw <- Read10X(data.dir = dir)
  dim(a.raw)
  colnames(a.raw) <- paste0(colnames(a.raw), name)
  
  a <- CreateSeuratObject(raw.data = a.raw, project = name)
  a@meta.data$var_sample <- name
  
  mito.genes <- grep(pattern = "^mt-", x = rownames(x = a@data), value = TRUE)
  percent.mito <- Matrix::colSums(a@raw.data[mito.genes, ])/Matrix::colSums(a@raw.data)
  a <- AddMetaData(object = a, metadata = percent.mito, col.name = "percent.mito")
  a@meta.data %>% filter(percent.mito > 0.05) %>% summarise(n=n())
  a@meta.data %>% filter(nGene < 1000) %>% summarise(n=n())
  
  png(paste0(projectdir,name,'_QC_scatterplot.png'))
  par(mfrow = c(1, 2))
  GenePlot(object = a, gene1 = "nUMI", gene2 = "percent.mito")
  GenePlot(object = a, gene1 = "nUMI", gene2 = "nGene")
  dev.off()
  
  vln=VlnPlot(object = a, features.plot = c("nGene", "nUMI", "percent.mito"), nCol = 3,do.return=TRUE, return.plotlist = T)
  cowplot::ggsave(paste0(projectdir,name,'_QC_vlnplot.png'),plot_grid(plotlist = vln),width=16,height=16)
  
  #this is a comment
 
  #We assign scores in the CellCycleScoring function, which stores S and G2/M scores in object@meta.data, along with the predicted classification of each cell in either G2M, S or G1 phase
  a <- CellCycleScoring(object = a, s.genes = cc.genes$s.genes, g2m.genes = cc.genes$g2m.genes)
  
  a <- FilterCells(a, subset.names = c("nGene","percent.mito"), low.thresholds = c(200,-Inf), 
                   high.thresholds = c(Inf,0.05))
  
  a <- NormalizeData(a)
  a <- FindVariableGenes(a, do.plot = F)
  if(ccscale==T){
  a <- ScaleData(object =a, vars.to.regress = c("nUMI", "percent.mito","S.Score", "G2M.Score"))}
  else{a <- ScaleData(a, display.progress = F, vars.to.regress=c("nUMI", "percent.mito"))}
  return(a)
  
}

#Setup seurat objects for each of the sc experiments
ctrl <- processExper("data/control/mm10/",'ctrl',ccscale=F)
flu2 <- processExper("data/Flu2//mm10/",'Flu2',ccscale=F)
#flu4 <- processExper("data/Flu4//mm10/",'Flu5')
flu5 <- processExper("data/Flu5/mm10/",'Flu5',ccscale=F)

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
combined <- RunMultiCCA(ob.list, genes.use = genes.use, num.ccs = 40)

p1 <- DimPlot(object = combined, reduction.use = "cca", group.by = "var_sample", 
              pt.size = 0.5, do.return = TRUE)
p2 <- VlnPlot(object = combined, features.plot = "CC1", group.by = "var_sample", 
              do.return = TRUE,pt.size = "NA")
pg=plot_grid(p1, p2)
ggsave(pg,filename="dimplot_vln.pdf",height=9,width=12,units = "in")
#DimPlot(object = combined, reduction.use = "cca", group.by = "var_sample", pt.size = 0.5)


p <- MetageneBicorPlot(combined, grouping.var = "var_sample", dims.eval = 1:50)
cowplot::ggsave(paste0(projectdir,'MetageneBiCorr.png'), width=16,height=8)


#alnmetric <- data.frame(dim=1:50,metric=0)
for(dim in seq(4,40,2)){

  # Alignment
  
  dir.create(paste0(projectdir,dim,'/'))
  
  combined <- CalcVarExpRatio(object = combined, reduction.type = "pca",
                              grouping.var = "var_sample", dims.use = 2:dim)
  
  scrna <- SubsetData(combined, 
                      subset.name = "var.ratio.pca",
                      accept.low = 0.5)
  
  
  
  scrna <- AlignSubspace(scrna,
                         reduction.type = "cca",
                         grouping.var = "var_sample",
                         dims.align = 1:dim
  )
  
  
  #alnmetric$metric[dim] <- CalcAlignmentMetric(scrna, reduction.use = "cca.aligned", dims.use=1:dim,
  #grouping.var='var_sample')
  
  
  
  #Create TSNE
  
  # t-SNE and Clustering
  scrna <- RunTSNE(scrna,
                   reduction.use = "cca.aligned",
                   dims.use = 1:dim,
                   check_duplicates = FALSE,
                   nthreads = 10,
                   max_iter = 2000)
  
  ### run UMAP.. 
  scrna <- RunUMAP(scrna, reduction.use = "cca.aligned", dims.use = 1:dim)
  
  ### Create a new loop to loop over res .4 to 1.6 by .2 
  for(res in seq(0.4,1.6,0.2)){
    scrna <- FindClusters(scrna, reduction.type = "cca.aligned",
                          dims.use = 1:dim, save.SNN = T, resolution = res)
    
    
    
    # Visualization make 4 plots.. TSNE with cluster, TSNE with var_sample, UMAP cluster nd var_sample
    
    plot_grid(
      TSNEPlot(scrna, do.label = T,do.return=T),
      DimPlot(scrna,reduction.use = "tsne",group.by='var_sample',do.return = T),
      DimPlot(scrna, reduction.use = "umap", do.label = T,do.return=T),
      DimPlot(scrna, reduction.use = "umap",group.by='var_sample', do.label = T,do.return=T)
    )
    
    
    ggsave(paste0(projectdir,dim,'/TSNE_dim',dim,'_res',res,'.png'), width=16,height=13)
  }
  fp <- FeaturePlot(object = scrna, features.plot = markergenes, cols.use = c("lightgrey", "blue"),do.return = T)
  cowplot::ggsave(paste0(projectdir,dim,'/TSNE_dim',dim,'_markergenes.png'),plot_grid(plotlist = fp),width=16,height=16)
  
  save(scrna,file=paste0(projectdir,dim,"/Lung_time_series",dim,'.RData'))
}

