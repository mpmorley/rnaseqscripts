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

dir='mm9/STAR'
projectname='NANCI_Effects_v2'

#Get a list of all of teh Star count files
(file_list <- list.files(dir,pattern="ReadsPerGene.out.tab$",full.names=F))

#This function will return a data frame of with gene and reads column. I need to modify this function to choose the proper col of counts, right now it assume
# coulmn 4 which is Strand paired reads. 
parse <- function(x, dir){
  d<-read.table(paste(dir,'/',x,sep=''),header=F, sep="\t")
  
  return(data.frame(sample=sub('ReadsPerGene.out.tab','',x),gene=d$V1,
                    signal=d$V4))
}

#the ldply loops over the filelists and appends data in a long format data.frame, this is piped to a filter to keep only 
#Ensembl ids and then formated into wide format. 
dataset <- ldply(file_list, parse, dir) %>% 
  filter(grepl('E',gene)) %>% 
  spread(sample,signal) %>% 
  select(-WT2,-Homozygous1,-Nkx4,-Double2)

#Create a matrix by dropiiing hte genes column
cts <- as.matrix(dataset[,2:dim(dataset)[2]])
rownames(cts)<-dataset$gene

#Get a data.frame of gene annotations, this may be not be the same size of input IDs. 
genenames <- GeneAnnotate(as.character(dataset$gene))
rownames(genenames)=genenames$ENSEMBL

############### Remove some Samples 


(pData <-read.csv('data/phenodata.csv'))
#given the genenames object might have less genes, limit the orginal counts abject to just that are annoated. 
dge <- DGEList(counts=cts[rownames(cts) %in% genenames$ENSEMBL,], genes=genenames)
cpms = cpm(dge)
#Filter low expressed genes, must have 25% with a cpm>1
keep = rowSums(cpms>1)>=.25*dim(dge)[2]
dge <- dge[keep,]
genenames = genenames[keep,]
dge <- calcNormFactors(dge)
plotMDS(dge)
boxplot(log2(cts))


(design <-model.matrix(~0+maineffect,data=pData))
#Clean up the colnames of the design matrix, don't need the colname in. 
colnames(design)<-gsub('maineffect','',colnames(design))
v <- voom(dge,design,plot=TRUE,normalize="quantile")
plotMDS(v)
v$E <- v$E + abs(min(v$E))+1

dge2 <- estimateDisp(dge, design, robust=TRUE)
plotBCV(dge2)



############### Create the eset #######################
## Create phenoData
rownames(pData) <- (pData$sample_name)
if(all(rownames(pData)!=colnames(v$E))){
  stop('Sample names not equal')
}
  
phenoData <- new("AnnotatedDataFrame",
data=pData)
fData <- new("AnnotatedDataFrame",
data=genenames)
all(rownames(genenames)==rownames(v$E))

eset<- ExpressionSet(assayData=v$E,phenoData=phenoData,featureData=fData,annotation="mm9")


############# Create contrast matrix and fit models ################

# Use the combn functio to make all possible contrasts 
#(f <-as.vector(unlist(combn(colnames(design),2,function(x)paste(x,collapse="-")))))
#read contrastlist.csv to create all possible contrasts
(contrastlist <-read.csv('data/contrastlist.csv'))
contrastlist$x_vs_y=paste(contrastlist$x,contrastlist$y,sep="-")
f=as.vector(contrastlist$x_vs_y)

(contrast.matrix <- makeContrasts(contrasts = f,levels=design))
fit <- lmFit(v,design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

######## load and prepare all the MSigDB sets for camera ######

if(pData$organism=="human"){
  load( '/fujfs/d1/projects/data_public/MSigDB/human_H_v5.rdata')
  h.indices <- ids2indices(Mm.H,genenames$ENTREZID)
  load(' /fujfs/d1/projects/data_public/MSigDB/human_c2_v5.rdata')
  c2.indices <- ids2indices(Mm.c2,genenames$ENTREZID)
  load(' /fujfs/d1/projects/data_public/MSigDB/human_c3_v5.rdata')
  c3.indices <- ids2indices(Mm.c3,genenames$ENTREZID)
  load(' /fujfs/d1/projects/data_public/MSigDB/human_c4_v5.rdata')
  c4.indices <- ids2indices(Mm.c4,genenames$ENTREZID)
}else
{
  load( '/fujfs/d1/projects/data_public/MSigDB/mouse_H_v5.rdata')
  h.indices <- ids2indices(Mm.H,genenames$ENTREZID)
  load(' /fujfs/d1/projects/data_public/MSigDB/mouse_c2_v5.rdata')
  c2.indices <- ids2indices(Mm.c2,genenames$ENTREZID)
  load(' /fujfs/d1/projects/data_public/MSigDB/mouse_c3_v5.rdata')
  c3.indices <- ids2indices(Mm.c3,genenames$ENTREZID)
  load(' /fujfs/d1/projects/data_public/MSigDB/mouse_c4_v5.rdata')
  c4.indices <- ids2indices(Mm.c4,genenames$ENTREZID)
  load('/Users/bapoorva/Desktop/ANALYSIS/msigdb/mouse_GO.rdata')
  GO.indices <- ids2indices(Mm.GO,genenames$ENTREZID)
}
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


i=1

for(i in 1:length(contrastnames)){
  print(contrastnames[i])
  limma[[contrastnames[i]]] <- Cleanup(topTable(fit2,coef=i,n=Inf,p.value=1)) 
  topgo[[contrastnames[i]]] <- runTopGO(limma[[contrastnames[i]]])
  
  k=limma[[contrastnames[i]]]

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


