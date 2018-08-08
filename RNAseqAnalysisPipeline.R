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


################################### Please Set the following Paramtrers ##################################
#
dir='RNAseq_1n2nCardio/STAR'
projectname='RNAseq_1n2nCardio'
#MSigDB_path='~/dsdata/projects/data_public/MSigDB/'
librarytype='unstranded'
constrastmaker='auto' #set to either file or auto
#
##########################################################################################################

#read phenodata
pData<- read.csv('data/phenodata.csv') %>% arrange(sample_name)
rownames(pData)=pData$sample_name
phenoData <- new("AnnotatedDataFrame", data=pData)


if(unique(pData$organism)=='human'){
  library(hugene20sttranscriptcluster.db)
  annotationdb=hugene20sttranscriptcluster.db
}else if (unique(pData$organism)=='mouse'){
  library(mogene20sttranscriptcluster.db)
  annotationdb=mogene20sttranscriptcluster.db
}else{
  stop("No proper orgasnism found, please edit the pData file")
}



#Get a list of all of the Star count files
(file_list <- list.files(dir,pattern="ReadsPerGene.out.tab$",full.names=F))


#This function will return a data frame of with gene and reads column. I need to modify this function to choose the proper col of counts, right now it assume
# coulmn 4 which is Strand paired reads. 

parse <- function(x, dir,librarytype){
  d<-read.table(paste(dir,'/',x,sep=''),header=F, sep="\t")
  switch(librarytype,
         unstranded=return(data.frame(sample=sub('ReadsPerGene.out.tab','',x),gene=d$V1,signal=d$V4)),
         fwdstrand=return(data.frame(sample=sub('ReadsPerGene.out.tab','',x),gene=d$V1,signal=d$V3)),
         revstrand=return(data.frame(sample=sub('ReadsPerGene.out.tab','',x),gene=d$V1,signal=d$V2))
  )
}

#the ldply loops over the filelists and appends data in a long format data.frame, this is piped to a filter to keep only 
#Ensembl ids and then formated into wide format. 
dataset <- ldply(file_list, parse, dir,librarytype) %>% 
  filter(grepl('E',gene)) %>% 
  spread(sample,signal)

######### Create a matrix of count data #######################
cts <- as.matrix(dataset[,2:dim(dataset)[2]])
rownames(cts)<-dataset$gene

############# Get a data.frame of gene annotations #######################
genenames <- GeneAnnotate(as.character(dataset$gene),organism = unique(pData$organism))
rownames(genenames)=genenames$ENSEMBL

######### Add protein type to the fdata #######################
#ptntype=read.csv("~/NGSshare/mm9_data/ptntype.csv")
data("ptntype",package="ExpressExtras")
genenames=left_join(genenames,ptntype,by=c("ENSEMBL"="ENSEMBL")) %>% dplyr::select(-gene)
rownames(genenames)=genenames$ENSEMBL
#################### Create a count matrix and filter CPM ###########

dge <- DGEList(counts=cts[rownames(cts) %in% genenames$ENSEMBL,], genes=genenames)
cpms = cpm(dge)
#Filter low expressed genes, must have 25% with a cpm>1
keep = rowSums(cpms>1)>=.25*dim(dge)[2]
dge <- dge[keep,]
genenames = genenames[keep,]
dge <- calcNormFactors(dge)
plotMDS(dge)
boxplot(log2(cts))


################ VOOM transform the data ############################

(design <-model.matrix(~0+maineffect,data=pData))
#Clean up the colnames of the design matrix, don't need the colname in. 
colnames(design)<-gsub('maineffect','',colnames(design))
v <- voom(dge,design,plot=TRUE)
plotMDS(v)
pData$minexpr=abs(min(v$E))+1
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

eset<- ExpressionSet(assayData=v$E,phenoData=phenoData,featureData=fData,annotation=unique(as.character(pData$organism)))


############# Create contrast matrix and fit models ################
if(constrastmaker=='auto'){
  f<-as.vector(unlist(combn(colnames(design),2,function(x)paste(x,collapse="-"))))
}else{
  f<-read.csv('data/contrastlist.csv') %>% 
    mutate(c=paste(treatment,control,sep="-")) %>%
    .$c
}


(contrast.matrix <- makeContrasts(contrasts = f,levels=design))
fit <- lmFit(v,design)
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
  print("Incorrect organism name. ")
}

##################################################################
#Loop over all contrasts and run limma, camera, topgo and spia for each one and save each result as Rdata object containing list of lists
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

for(i in 1:length(contrastnames)){
  print(contrastnames[i])
  limma[[contrastnames[i]]] <- Cleanup2(topTable(fit2,coef=i,n=Inf,p.value=1)) 
  topgo[[contrastnames[i]]] <- runTopGO(limma[[contrastnames[i]]],organism =unique(pData$organism))
  
  k=limma[[contrastnames[i]]]
  
  #for each limma data (corresponding to the contrast), run SPIA
  limma_sel <- k[which(abs(k$fc) > 2 & k$adj.P.Val < 0.05),]
  mm=c("mouse","Mouse","Mm","Mus musculus","Mus_musculus")
  hs=c("human","Human","Hs","Homo sapiens","Homo_sapiens")
  org=ifelse(unique(pData$organism) %in% mm,"mmu",ifelse(unique(pData$organism) %in% hs,'hsa',"NA"))
  if(nrow(limma_sel)>0){
    all_genes = as.numeric(k$ENTREZID)
    sig_genes = limma_sel$fc
    names(sig_genes) = limma_sel$ENTREZID 
    sig_genes = sig_genes[complete.cases(names(sig_genes))]
    sig_genes = sig_genes[unique(names(sig_genes))] 
    spia[[contrastnames[i]]] <- spia(de=sig_genes, all=all_genes, organism=org)
  }else{
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


