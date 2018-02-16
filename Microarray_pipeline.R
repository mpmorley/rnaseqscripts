library(OligoTools)
require(sva)
require(dplyr)
library(limma)
require(devtools)
install_github('mpmorley/ExpressExtras')
library(ExpressExtras)
library(plyr)
library(dplyr)
library(tidyr)
library(limma)
library(NMF)
library(RColorBrewer)
library(ggplot2)
require(SPIA)
library(oligo)
library(readr)
library(topGO)


################################### Please Set the following Paramtrers ##################################
#
projectname='CCAM_array'
constrastmaker='auto' #set to either file or auto
#
##########################################################################################################

#read phenodata
(pData<- read.csv('data/phenodata.csv'))
rownames(pData)=pData$sample_name
phenoData <- new("AnnotatedDataFrame", data=pData)

hum=c("human","Human","Hs","Homo sapiens","Homo_sapiens")
mouse=c("mouse","Mouse","Mm","Mus musculus","Mus_musculus")
if(unique(pData$organism) %in% hum){
  library(hugene20sttranscriptcluster.db)
  annotationdb=hugene20sttranscriptcluster.db
}else if (unique(pData$organism) %in% mouse){
  library(mogene20sttranscriptcluster.db)
  annotationdb=mogene20sttranscriptcluster.db
}else{
  stop("No proper orgasnism found, please edit the pData file")
}


#read cel files
celFiles <- list.celfiles("./CEL", full.names = TRUE)
affyGeneFS <- read.celfiles(celFiles,phenoData=phenoData,verbose=TRUE)

#normalize cel data
eset = rma(affyGeneFS,target='core')

#annotate microarray data, get gene_biotype, seq_name and create new featuresdata 
eset.bk <- EnsemblAnnotate(eset,annotationdb)
eset.bk <-BackgrdFilter(eset.bk)
eset.bk <-GetMainProbes(eset.bk)

#remove all genes without annotation
eset.bk <- eset.bk[!is.na(eset.bk@featureData@data$ENTREZID),]

#Update annotation to include biotype and geneloc
genenames=pData(featureData(eset.bk))
genes <- GeneAnnotate(as.character(genenames$ENSEMBL),organism = pData$organism)
final_res <- left_join(genenames,genes,by=c('ENSEMBL'='ENSEMBL')) %>% select(ID,SYMBOL.x,GENENAME,ENTREZID.x,ENSEMBL,GenestoProbe,biotype,geneloc) %>% dplyr::rename(SYMBOL=SYMBOL.x,ENTREZID=ENTREZID.x)
rownames(final_res)=make.names(final_res$ID,unique=TRUE)
p=strsplit(rownames(final_res),"X")
m=sapply(p,"[",2)
rownames(final_res)=m
genenames=final_res

eset=eset.bk
#substitute new annotation data as featuresdata in the eset
fData <- new("AnnotatedDataFrame",data=genenames)
eset@featureData=fData

#create design matrix
(design <-model.matrix(~0+maineffect,data=pData))
colnames(design)=gsub('maineffect','',colnames(design))


# Use the combn functio to make all possible contrasts 
#f <-as.vector(unlist(combn(colnames(design),2,function(x)paste(x,collapse="-"))))
if(constrastmaker=='auto'){
  f<-as.vector(unlist(combn(colnames(design),2,function(x)paste(x,collapse="-"))))
}else{
  f<-read.csv('data/contrastlist.csv') %>% 
    mutate(c=paste(treatment,control,sep="-")) %>%
    .$c
}
(contrast.matrix <- makeContrasts(contrasts = f,levels=design))
fit <- lmFit(eset,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

######## load and prepare all the MSigDB sets for camera ######
hum=c("human","Human","Hs","Homo sapiens")
mouse=c("mouse","Mouse","Mm","Mus musculus")
rat=c("rat","Rat","Rn","Rattus norvegicus")
if(unique(pData$organism) %in% hum){
  data('human_H_v5',package="ExpressExtras")
  h.indices <- ids2indices(Hs.H,genenames$ENTREZID)
  data('human_c2_v5',package="ExpressExtras")
  c2.indices <- ids2indices(Hs.c2,genenames$ENTREZID)
  data('human_c3_v5',package="ExpressExtras")
  c3.indices <- ids2indices(Hs.c3,genenames$ENTREZID)
  data('human_c5_v5',package="ExpressExtras")
  GO.indices <- ids2indices(Hs.c4,genenames$ENTREZID)
}else if(unique(pData$organism) %in% mouse){
  data('mouse_H_v5',package="ExpressExtras")
  h.indices <- ids2indices(Mm.H,genenames$ENTREZID)
  data('mouse_c2_v5',package="ExpressExtras")
  c2.indices <- ids2indices(Mm.c2,genenames$ENTREZID)
  data('mouse_c3_v5',package="ExpressExtras")
  c3.indices <- ids2indices(Mm.c3,genenames$ENTREZID)
  data('mouse_c4_v5',package="ExpressExtras")
  c4.indices <- ids2indices(Mm.c4,genenames$ENTREZID)
  data('mouse_GO',package="ExpressExtras")
  GO.indices <- ids2indices(Mm.GO,genenames$ENTREZID)
}else if(unique(pData$organism) %in% rat){
  data('Rat_Hallmark',package="ExpressExtras")
  h.indices <- ids2indices(Rn.H,genenames$ENTREZID)
  data('Rat_C2_v4p2',package="ExpressExtras")
  c2.indices <- ids2indices(Rn.c2,genenames$ENTREZID)
  data('Rat_GO',package="ExpressExtras")
  GO.indices <- ids2indices(Rn.GO,genenames$ENTREZID)
}else {
  print("Incorrect organism name")
}



##################################################################

#Remove the '-' from the constrat name, it will cause issues down stream
(contrastnames <-gsub('-','_vs_',colnames(contrast.matrix)))
(contrastnames <-gsub(' ','',contrastnames))
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



for(i in 1:length(contrastnames)){
  print(contrastnames[i])
  limma[[contrastnames[i]]] <- Cleanup(topTable(fit2,coef=i,n=Inf,p.value=1)) 
  topgo[[contrastnames[i]]] <- runTopGO(limma[[contrastnames[i]]],organism =unique(pData$organism))
  
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
    spia[[contrastnames[i]]] <- spia(de=sig_genes, all=all_genes, organism=ifelse(unique(pData$organism)=='mouse',"mmu",'hsa'))
  }else{
    spia[[contrastnames[i]]] <- data.frame()
  }
  res.h <- camera(v, h.indices, design,contrast.matrix[,i],inter.gene.cor=0.01)
  res.c2 <- camera(v, c2.indices, design,contrast.matrix[,i],inter.gene.cor=0.01)
  res.GO <- camera(v, GO.indices, design,contrast.matrix[,i],inter.gene.cor=0.01)
  #res.c3 <- camera(v, c3.indices, design,i,inter.gene.cor=0.01)
  #res.c4 <- camera(v, c4.indices, design,i,inter.gene.cor=0.01)
  camera[[contrastnames[i]]] <- list(Hallmark=list(camera_result=res.h,indices=h.indices),Curated=list(camera_result=res.c2,indices=c2.indices),GO=list(camera_result=res.GO,indices=GO.indices))
}

############ Save list of results ##############
results <- list(eset=eset,limma=limma,camera=camera, topgo=topgo,spia=spia)
save(results,file=paste(projectname, ".RData",sep=''))


