library(OligoTools)
library(mogene20sttranscriptcluster.db)
require(sva)
require(dplyr)
library(limma)
require(EnsDb.Mmusculus.v75)
require(devtools)
install_github('mpmorley/ExpressExtras')
library(ExpressExtras)
library(plyr)
library(dplyr)
library(tidyr)
library(limma)
library(NMF)
library(RColorBrewer)
library(edgeR)
library(ggplot2)
require(org.Mm.eg.db)
require(EnsDb.Mmusculus.v75)
require(SPIA)

projectname='Falcor_Naphthalene'

#read phenodata
pData<- read.csv('data/PhenoData.csv')
rownames(pData)=pData$sample_name
phenoData <- new("AnnotatedDataFrame", data=pData)

#read cel files
celFiles <- list.celfiles("./CEL", full.names = TRUE)
affyGeneFS <- read.celfiles(celFiles,phenoData=phenoData,verbose=TRUE)

#normalize cel data
eset = rma(affyGeneFS,target='core')

#annotate microarray data, get gene_biotype, seq_name and create new featuresdata 
eset.bk <- EnsemblAnnotate(eset,mogene20sttranscriptcluster.db)
eset.bk <-BackgrdFilter(eset.bk)
eset.bk <-GetMainProbes(eset.bk)

voom=exprs(eset.bk)
pData=pData(eset.bk)
genenames=pData(featureData(eset.bk))

genes <- GeneAnnotate(as.character(genenames$ENSEMBL))

final_res <- left_join(genenames,genes,by=c('ENSEMBL'='ENSEMBL')) %>% select(ID,SYMBOL.x,GENENAME,ENTREZID.x,ENSEMBL,GenestoProbe,biotype,geneloc) %>% dplyr::rename(SYMBOL=SYMBOL.x,ENTREZID=ENTREZID.x)
rownames(final_res)=make.names(final_res$ID,unique=TRUE)
p=strsplit(rownames(final_res),"X")
m=sapply(p,"[",2)
rownames(final_res)=m
genenames=final_res


#create design matrix
(design <-model.matrix(~0+maineffect,data=pData))
colnames(design)=c("FalcorKO","WildType")
############### Create the eset #######################
## Create phenoData
rownames(pData) <- (pData$sample_name)
phenoData <- new("AnnotatedDataFrame",
                 data=pData)
fData <- new("AnnotatedDataFrame",
             data=genenames)
all(rownames(genenames)==rownames(voom))

eset<- ExpressionSet(assayData=as.matrix(voom),phenoData=phenoData,featureData=fData,annotation="mm9")

############# Create contrast matrix and fit models ################

# Use the combn functio to make all possible contrasts 
f <-as.vector(unlist(combn(colnames(design),2,function(x)paste(x,collapse="-"))))
(contrast.matrix <- makeContrasts(contrasts = f,levels=design))
fit <- lmFit(eset,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

######## load and prepare all the MSigDB sets for camera ######

load('/Users/bapoorva/Desktop/ANALYSIS/msigdb/mouse_H_v5p1.rdata')
h.indices <- ids2indices(Mm.H,genenames$ENTREZID)
load('/Users/bapoorva/Desktop/ANALYSIS/msigdb/mouse_c2_v5p1.rdata')
c2.indices <- ids2indices(Mm.c2,genenames$ENTREZID)
load('/Users/bapoorva/Desktop/ANALYSIS/msigdb/mouse_c3_v5p2.rdata')
c3.indices <- ids2indices(Mm.c3,genenames$ENTREZID)
load('/Users/bapoorva/Desktop/ANALYSIS/msigdb/mouse_c4_v5p2.rdata')
c4.indices <- ids2indices(Mm.c4,genenames$ENTREZID)
load('/Users/bapoorva/Desktop/ANALYSIS/msigdb/mouse_GO.rdata')
GO.indices <- ids2indices(Mm.GO,genenames$ENTREZID)



##################################################################

#Remove the '-' from the constrat name, it will cause issues down stream
(contrastnames <-gsub('-','_vs_',colnames(contrast.matrix)))
#Create list to hold the results for limma,togo and camera for all contrasts
limma <-  vector(mode="list", length=length(contrastnames))
names(limma) <- contrastnames
camera <- vector(mode="list", length=length(contrastnames))
names(camera) <- contrastnames
topgo <-vector(mode="list", length=length(contrastnames))
names(topgo) <- contrastnames
spia<-  vector(mode="list", length=length(contrastnames))
names(spia) <- contrastnames
# eplot<-  vector(mode="list", length=length(contrastnames))
# names(eplot) <- contrastnames
# gsea<-  vector(mode="list", length=length(contrastnames))
# names(gsea) <- contrastnames


i=1

for(i in 1:length(contrastnames)){
  print(contrastnames[i])
  limma[[contrastnames[i]]] <- Cleanup(topTable(fit2,coef=i,n=Inf,p.value=1)) 
  topgo[[contrastnames[i]]] <- runTopGO(limma[[contrastnames[i]]])
  
#   #for each limma data (corresponding to the contrast), run ReactomePA
   k=limma[[contrastnames[i]]]
#   genes = k$fc
#   names(genes) = k$ENTREZID 
#   genes = genes[complete.cases(names(genes))]
#   genes = genes[unique(names(genes))] 
#   ss=genes[order(-genes)]
#   y <- gsePathway(ss,organism="mouse")
#   gsea[[contrastnames[i]]]=y
#   if(nrow(summary(y))>0){
#     eplot[[contrastnames[i]]] <-summary(y)}
#   else{
#     eplot[[contrastnames[i]]] <- data.frame()
#   }
  
  #for each limma data (corresponding to the contrast), run SPIA
  limma_sel <- k[which(abs(k$fc) > 2 & k$adj.P.Val < 0.05),]
  if(nrow(limma_sel)>0){
    all_genes = as.numeric(k$ENTREZID)
    sig_genes = limma_sel$fc
    names(sig_genes) = limma_sel$ENTREZID 
    sig_genes = sig_genes[complete.cases(names(sig_genes))]
    sig_genes = sig_genes[unique(names(sig_genes))] 
    spia[[contrastnames[i]]] <- spia(de=sig_genes, all=all_genes, organism="mmu")}
  else{
    spia[[contrastnames[i]]] <- data.frame()
  }
  
  #run camera
  res.h <- camera(voom, h.indices, design,contrast.matrix[,i],inter.gene.cor=0.01)
  res.c2 <- camera(voom, c2.indices, design,contrast.matrix[,i],inter.gene.cor=0.01)
  res.GO <- camera(voom, GO.indices, design,contrast.matrix[,i],inter.gene.cor=0.01)
  #res.c3 <- camera(v, c3.indices, design,i,inter.gene.cor=0.01)
  #res.c4 <- camera(v, c4.indices, design,i,inter.gene.cor=0.01)
  camera[[contrastnames[i]]] <- list(Hallmark=list(camera_result=res.h,indices=h.indices),Curated=list(camera_result=res.c2,indices=c2.indices),GO=list(camera_result=res.GO,indices=GO.indices))
}


############ Save list of results ##############
  results <- list(eset=eset,limma=limma,camera=camera, topgo=topgo,spia=spia)
save(results,file=paste(projectname, ".RData",sep=''))


