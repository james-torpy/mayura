
library("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome)
library(reshape2)
library(ggplot2)
library(edgeR)
library(rtracklayer)
library(RColorBrewer)
library(org.Mm.eg.db)
library(KEGGREST)
library(pathview)


timeStamp <- format(Sys.time(), "%Y_%m_%d")
species <- "mouse"
ensVer <- 84



homedir="/share/ClusterShare/biodata/contrib/nenbar"
homedir="../../../../"
inPath=paste0(homedir,"/projects/claudia/project_results/erv.rsem/")


######## directory structure #######
projectDir=paste0(homedir,"/projects/claudia")
resultsDir=paste0(projectDir,"/project_results")
imageDir=paste0(resultsDir,"/figures")
annotationDir=paste0(projectDir,"/annotation/API/")
projectname="erv"
robjectsDir = paste(resultsDir,"/erv.Robjects/",sep="")
scriptsPath=paste(projectDir,"/scripts/repeats")
logDir=paste0(scriptsPath,"/logs")
sizeDir=paste(resultsDir,"/erv.libSize/",sep="")

outPath=paste0(resultsDir,"/erv.repeatOverlap/")
system(paste("mkdir",outPath))
system(paste("mkdir",sizeDir))

chrs=seqlengths(Mmusculus)[!grepl("_",names(seqlengths(Mmusculus)))]
#names(chrs)=gsub("chr","",names(chrs))
#names(chrs)[25]="MT"
gr<-GRanges(seqnames=names(chrs),IRanges(start=1,end=chrs))

#Load both type and class lists
#load(paste0(annotationDir,"class.Rdata"))
#grClass<-class
#load(paste0(annotationDir,"type.Rdata"))
#grType<-type



#########################################
########## 1. Load in the files
#########################################

geneFiles<-list.files(inPath,pattern="genes",full.names=T)

geneResults<-list()

for(file in geneFiles){
  sampleName<-basename(file)
  sampleName<-gsub(".transcriptome.*","",sampleName)
  cat(sampleName)
  cat("\n")
  data<-read.table(file,header=T)  
  geneResults[[sampleName]]<-as.integer(data$expected_count)
}
temp<-geneResults

df<-as.data.frame(geneResults)
row.names(df)<-data$gene_id


#########################################
########## 1. Load in the files
#########################################



cols<-brewer.pal(6,"Paired")
pdf("../../project_results/figures/pca.pdf",width=12,height=8)
pca<-princomp(df)
plot(pca$loading,pch=19, cex=2,col=cols)
text(pca$loading, names(geneResults),pos = 1)
dev.off()

pdf("../../project_results/figures/pca_components.pdf",width=12,height=8)
plot(pca)
dev.off()

#Annotate with symbols, aggregate

#Annotate with entrez, aggregate
egENSEMBL <- toTable(org.Mm.egENSEMBL)
row.names(df)=gsub("\\..*","",row.names(df))
m <- match(row.names(df), egENSEMBL$ensembl_id)
df$EntrezGene<-egENSEMBL$gene_id[m]

egSYMBOL <- toTable(org.Mm.egSYMBOL)
m <- match(df$EntrezGene, egSYMBOL$gene_id)
df$symbol<-egSYMBOL$symbol[m]

#eliminate duplicated symbols
o <- order(rowSums(df[,1:6]), decreasing=TRUE)
df=df[o,]
d<-duplicated(df$symbol)
df<-df[!d,]

#add the gene length
#data$gene_id=gsub("\\..*","",data$gene_id)
#df_merged<-merge(df,data[,c(1,4)],by.x=0,by.y="gene_id")
#df<-df_merged
#row.names(df)<-df[,1]
#df<-df[,-1]
#eliminate lowly expressed
include<-apply(df[,1:6],1,function(x){sum(x>=50)>=2})
df<-df[include,]
df=df[!is.na(df$symbol),]
row.names(df)<-df$symbol

#change the id to entrez
o <- order(rowSums(df[,1:6]), decreasing=TRUE)
df=df[o,]
d<-duplicated(df$EntrezGene)
df<-df[!d,]
row.names(df)<-df$EntrezGene





################## DESEQ of the individual samples

grouping <- colnames(df)[1:6]
#grouping <- factor(type)

countTable=df[,1:6]
condition=grouping
cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )
cds = estimateDispersions( cds, method='blind',sharingMode="fit-only" )

results<-list()
for(i in c(2,4,6)){
    res = nbinomTest( cds, colnames(df)[i-1], colnames(df)[i] )
    sampleName=gsub("_.*","",colnames(df)[i])

    res[res$baseMeanA==0,"log2FoldChange"]<-log2(res[res$baseMeanA==0,"baseMeanB"])
    res[(res$baseMeanA==0&res$baseMeanB==0),"log2FoldChange"]<-0
    res[res$baseMeanB==0,"log2FoldChange"]<-log2(res[res$baseMeanB==0,"baseMeanB"])
    res[(res$baseMeanA==0&res$baseMeanB==0),"log2FoldChange"]<-0
    results[[sampleName]]<-res

#Fix this bit for the p-values of the pathways for the individual comparisons
#    go <- goana(lrt, species="Mm")
#    topGO(go)
#    write.table(topGO(go,number=Inf),"GO_categories.txt",row.names=T,quote=F,sep="\t")
#
#    keg <- kegga(lrt, species="Mm")
#    topKEGG(keg)
#    write.table(topKEGG(keg,number=Inf),"KEGG_categories.txt",row.names=T,quote=F,sep="\t")



}


diffExprs<-do.call(cbind,results)
colnames(diffExprs)<-paste0(rep(names(results),each=8),"_",colnames(res))
colnames(diffExprs)[1]<-"id"
diffExprs<-diffExprs[,-grep("_id",colnames(diffExprs))]
row.names(diffExprs)<-diffExprs$id

types=c("B6","Nod","SJL")
pathways=c("04060","05168","05322","04062","05164","04612","04620","05416","04940","04623","04621","04622")
for(type in types){
    for(pathway in pathways){
        small=diffExprs[,paste0(type,"_log2FoldChange")]
        names(small)=row.names(diffExprs)

        pv.out <- pathview(gene.data = small, pathway.id = pathway,
            species = "mmu", out.suffix = paste0(type), kegg.native = T, same.layer = F)




    }
}




treatment <- rep(c("control","inflammed" ), 3)
type <- c(rep(c("B6", "Nod", "SJL"), each=2))

design <- model.matrix(~type + treatment)
df=df[,c(2,1,4,3,6,5,7,8)]
#one record has a length of 0
#df=df[row.names(df)!="ENSMUSG00000064945",]
expr <- DGEList(counts=df[1:6])
expr <- calcNormFactors(expr)
expr <- estimateDisp(expr,design)
fit <- glmFit(expr,design)
lrt <- glmLRT(fit)
topTags(lrt)


d=topTags(lrt,n=10000,p.value=0.01)
d=as.data.frame(d)


types=c("common")
pathways=c("04060","05168","05322","04062","05164","04612","04620","05416","04940","04623","04621","04622")
for(type in types){
    for(pathway in pathways){
        small=d[,"logFC"]
        small=-small
        names(small)=row.names(d)

        pv.out <- pathview(gene.data = small, pathway.id = pathway,
            species = "mmu", out.suffix = paste0(type), kegg.native = T, same.layer = F)
    }
}




