
library(GenomicRanges)
library(ShortRead)
#library(R.utils)
library("BSgenome.Mmusculus.UCSC.mm10")
library(BSgenome)
library(reshape2)
library(ggplot2)
library(edgeR)
library(rtracklayer)
library(RColorBrewer)
library(RUVSeq)
library(org.Mm.eg.db)
library(DESeq)
library(pheatmap)
library(GSEABase)
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

outPath=paste0(resultsDir,"/GSEA/")
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
include<-apply(df[,1:6],1,function(x){sum(x>=5)>=2})
df<-df[include,]
df=df[!is.na(df$symbol),]
row.names(df)<-df$symbol




treatment <- rep(c("inflammed", "control"), 3)
type <- c(rep(c("B6", "Nod", "SJL"), each=2))

design <- model.matrix(~type + treatment)

#one record has a length of 0
#df=df[row.names(df)!="ENSMUSG00000064945",]
expr <- DGEList(counts=df[1:6])
expr <- calcNormFactors(expr)
expr <- estimateDisp(expr,design)
fit <- glmFit(expr,design)
lrt <- glmLRT(fit)


logFCs<-as.data.frame(topTags(lrt,n=20000))
logFCs$symbol<-rownames(logFCs)
forGSEA<-logFCs[,c(6,1,5)]
write.table(forGSEA,"../../project_results/tables/forGSEA.rnk",quote=F,sep="\t",row.names=F)



plus<-logFCs[logFCs$logFC>0,]
minus<-logFCs[logFCs$logFC<0,]
significantIDs<-c(row.names(plus[1:50,]))
fittedValues<-fitted.values(lrt)
fittedValues50<-fittedValues[row.names(fittedValues) %in% significantIDs,]
colnames(fittedValues50)<-gsub("_R1","",colnames(fittedValues50))
fittedValues50=fittedValues50[,c(2,1,4,3,6,5)]
cols=rev(colorRampPalette(rev(brewer.pal(n = 7, name ="Blues")))(100))
pdf("../../project_results/figures/logFC_genes_50plus.pdf",width=12,height=8)
pheatmap(log10(fittedValues50),cluster_cols=F,cluster_rows=T,main="log10 counts")
dev.off()

#expr_norm <- rpkm(expr, log=TRUE, gene.length=df$effective_length)
#save(expr_norm,file="normalized_expression.Rdata")
#pdf("../../project_results/figures/pca_fpkm.pdf",width=12,height=8)
#pca<-princomp(expr_norm)
#plot(pca$loading,pch=19, cex=2,col=cols)
#text(pca$loading, names(geneResults),pos = 1)
#dev.off()

#pdf("../../project_results/figures/pca_components_fpkm.pdf",width=12,height=8)
#plot(pca)
#dev.off()

go <- goana(lrt, species="Mm")
topGO(go)
write.table(topGO(go,number=Inf),"GO_categories.txt",row.names=T,quote=F,sep="\t")

keg <- kegga(lrt, species="Mm")
topKEGG(keg)
write.table(topKEGG(keg,number=Inf),"KEGG_categories.txt",row.names=T,quote=F,sep="\t")
break()

row.names(df)<-df$symbol
expr <- DGEList(counts=df[1:6])
expr <- calcNormFactors(expr)
expr <- estimateDisp(expr,design)
fit <- glmFit(expr,design)
lrt <- glmLRT(fit)
topTags(lrt)

d=topTags(lrt,n=10000,p.value=0.01)
write.table(d,"all_inf.txt",quote=F,sep="\t",row.names=T)




################## DESEQ of the individual samples

grouping <- colnames(df)[1:6]
#grouping <- factor(type)

countTable=df[,1:6]
condition=grouping
cds = newCountDataSet( countTable, condition )
cds = estimateSizeFactors( cds )
cds = estimateDispersions( cds, method='blind',sharingMode="fit-only" )

results<-list()
for(i in c(1,3,5)){
    res = nbinomTest( cds, colnames(df)[i+1],colnames(df)[i])
    sampleName=gsub("_.*","",colnames(df)[i])
    res[res$baseMeanA==0,"log2FoldChange"]<-log2(res[res$baseMeanA==0,"baseMeanB"])
    res[(res$baseMeanA==0&res$baseMeanB==0),"log2FoldChange"]<-0
    res[res$baseMeanB==0,"log2FoldChange"]<-log2(res[res$baseMeanB==0,"baseMeanB"])
    res[(res$baseMeanA==0&res$baseMeanB==0),"log2FoldChange"]<-0
    results[[sampleName]]<-res

}

diffExprs<-do.call(cbind,results)
colnames(diffExprs)<-paste0(rep(names(results),each=8),"_",colnames(res))
colnames(diffExprs)[1]<-"id"
diffExprs<-diffExprs[,-grep("_id",colnames(diffExprs))]

write.table(diffExprs,"pairwise_diff_expression.txt",row.names=F,quote=F,sep="\t")



normalized.count <- counts(cds, normalized = TRUE)
normalized.count<-normalizeQuantiles(normalized.count)
#save(normalized.count,file="../project_results//endo_1.tables/normalized.counts.Rdata")

############ heatmap

geneSelection<-normalized.count[grep("Ifnb1|Ifna2|Aim2|Nod2|Mb21d1",row.names(normalized.count)),]
dataM<-melt(geneSelection)
dataM$sample<-gsub("_.*","",dataM$Var2)
dataM$inflammed<-rep(rep(c("inflammed","noninflammed"),each=5),3)
colnames(dataM)[1]="gene"
colnames(dataM)[3]="normalized_count"
dataM$sampleNames<-paste(dataM$sample,dataM$inflammed,sep="_")
cols<-rev(brewer.pal(6,"Paired"))
names(cols)<-sampelNames
dataM$cols<-rep(cols,each=5)

pdf("../../project_results/figures/normalized_counts_genes.pdf",width=12,height=8)
p<-ggplot(dataM,aes(gene,normalized_count,group=dataM$sampleNames))
p<-p+geom_line(aes(group=sampleNames),colour=dataM$cols)
p<-p+scale_color_manual(values=cols)
p
dev.off()

annotation<-read.table("../project_results//endo_1.tables/annotation.txt",header=T)
#combine the names and put into degust
normalized.count<-as.data.frame(normalized.count)
normalized.count$id<-row.names(normalized.count)
merged<-merge(normalized.count,annotation,by.x="id",by.y="gene_id")
write.table(merged,"../project_results/endo_1.tables/normalized.count.txt",row.names=T,quote=F,sep="\t")











