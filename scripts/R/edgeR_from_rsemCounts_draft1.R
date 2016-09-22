##### edgeR_from_rsemCounts #####

# This script takes a transcriptome counts file outputted by RSEM and performs a differential analysis:

# load packages:
library(RColorBrewer)
library(org.Mm.eg.db)
library(edgeR)
library(stringr)
library(dplyr)

library(GenomicRanges)
library(ShortRead)
library(reshape2)
library("BSgenome.Mmusculus.UCSC.mm10")
library(R.utils)
library(ggplot2)
library(rtracklayer)
library(RUVSeq)
library(DESeq)

# set up directory structure:
projectname <- "mayura"

# if fetching all files from project (not just one sample type), leave 'samplename' blank (i.e. "")
samplename <- "MW_310816"
inType <- ".rsem"


# use following homeDir for R shell, files on cluster:
#homeDir <- "/home/jamtor"
# use following homeDir for R Studio, files on cluster:
#homeDir <- "/Users/jamestorpy/clusterHome"
# use following homeDir for files on local drive:
homeDir <- "/Users/jamestorpy/Documents/Garvan/phd"

projectDir <- paste0(homeDir, "/projects/", projectname)
resultsDir <- paste0(projectDir, "/results/")
inDir <- paste0(resultsDir, samplename, "/", samplename, inType)
paste0("The inDir is: ", inDir)
setwd(inDir)

plotDir <- paste0(resultsDir, samplename, "/plots")
tableDir <- paste0(resultsDir, samplename, "/tables")

system(paste("mkdir", plotDir))
system(paste("mkdir", tableDir))

##### 1. Load in the files #####

# fetch inFiles:
inFiles = list.files(inDir, pattern = ".transcriptome.genes.results", recursive = T, full.names=T)
inFiles = grep("subset", inFiles, inv=T, value=T)
paste0("The inFiles are: ")
inFiles

# create empty list for geneResults:
geneResults <- list()

for (file in inFiles) {
  # define samplenames:
  sampleName <- basename(file)
  sampleName <- sub(".transcriptome.*","",sampleName)
  cat(sampleName)
  cat("\n")
  # read in samples and add to geneResults list:
  data <- read.table(file, header = T)
  geneResults[[sampleName]]<-as.integer(data$expected_count)
}
temp <- geneResults

# define df as data frame of geneResults with gene ID row names:
df <- as.data.frame(geneResults)
row.names(df)<-data$gene_id

# define plot colour scheme:
colrs <- brewer.pal(6, "Paired")
# perform pca and visualise on plot, save as pdf:
pdf(paste0(plotDir, "/pca.pdf"), width=12, height=8)
pca <- princomp(df)
plot(pca$loading, pch=19, cex=2, col=colrs)
text(pca$loading, names(geneResults), pos=1)
dev.off()

# visualise pca components on plot, save as pdf:
pdf(paste0(plotDir, "/pca_components.pdf"), width=12, height=8)
plot(pca)
dev.off()

# annotate df with entrez IDs:
egENSEMBL <- toTable(org.Mm.egENSEMBL)
row.names(df) <- gsub("\\..*", "", row.names(df))
m <- match(row.names(df), egENSEMBL$ensembl_id)
df$EntrezGene <- egENSEMBL$gene_id[m]

# annotate df with gene symbols:
egSYMBOL <- toTable(org.Mm.egSYMBOL)
m <- match(df$EntrezGene, egSYMBOL$gene_id)
df$symbol<-egSYMBOL$symbol[m]

# order rows from highest to lowest row sums:
o <- order(rowSums(df[,1:8]), decreasing=TRUE)
df <- df[o,]

# eliminate duplicated symbols:
d <- duplicated(df$symbol)
df <- df[!d,]
# eliminate lowly expressed genes (rows where there are less than 2 counts >=50):
include <- apply(df[1,1:8], 1, function(x){sum(x>=50)>=2})
df <- df[include,]
# eliminate genes with symbol 'NA':
df <- df[!is.na(df$symbol),]

# make row names gene symbols:
row.names(df) <- df$symbol

# order rows from highest to lowest row sums:
### why do this again? ###
o <- order(rowSums(df[,1:8]), decreasing=TRUE)
df <- df[o,]

# eliminate duplicated entrez gene IDs:
d <- duplicated(df$EntrezGene)
df <- df[!d,]

# change row names to entrez gene IDs:
### why not do this to start instead of symbols? ###
row.names(df) <- df$EntrezGene

# add the EC (enzyme) numbers:
egENZYME <- toTable(org.Mm.egENZYME)
m <- match(row.names(df), egENZYME$gene_id)
df[,"EC number"]<-egENZYME$ec_number[m]

# sanity check - check how many enzyme entries:
enzyme_no <- length(df[complete.cases(df["EC number"]), "EC number"])
print(paste0("There are ", enzyme_no, " enzyme entries"))

# define plot colour scheme:
colrs <- brewer.pal(6, "Paired")
# visualise clean pca on plot, save as pdf:
pdf(paste0(plotDir, "/pca_clean.pdf"), width=12, height=8)
#pca <- princomp(df[1:8])
pca <- princomp(df[1:8])
plot(pca$loading, pch=19, cex=2, col=colrs)
text(pca$loading, names(geneResults), pos=1)
dev.off()

# check pca without MW1:
df_1 <- cbind(df[,3:4], df[5:8])
geneResults_1 <- c(geneResults[3:4], geneResults[5:8])
# visualise clean pca on plot, save as pdf:
pdf(paste0(plotDir, "/pca_clean_no_MW1.pdf"), width=12, height=8)
pca <- princomp(df_1[1:6])
plot(pca$loading, pch=19, cex=2, col=colrs)
text(pca$loading, names(geneResults_1), pos=1)
dev.off()

# check pca without MW2:
df_2 <- cbind(df[,1:2], df[5:8])
geneResults_2 <- c(geneResults[1:2], geneResults[5:8])
# visualise clean pca on plot, save as pdf:
pdf(paste0(plotDir, "/pca_clean_no_MW2.pdf"), width=12, height=8)
pca <- princomp(df_2[1:6])
plot(pca$loading, pch=19, cex=2, col=colrs)
text(pca$loading, names(geneResults_2), pos=1)
dev.off()

# check pca without MW3:
df_3 <- cbind(df[,1:2], df[7:8])
geneResults_3 <- c(geneResults[1:2], geneResults[7:8])
# visualise clean pca on plot, save as pdf:
pdf(paste0(plotDir, "/pca_clean_no_MW3.pdf"), width=12, height=8)
pca <- princomp(df_1[1:6])
plot(pca$loading, pch=19, cex=2, col=colrs)
text(pca$loading, names(geneResults_1), pos=1)
dev.off()

# check pca without MW4:
df_4 <- df[,1:6]
geneResults_4 <- c(geneResults[1:6])
# visualise clean pca on plot, save as pdf:
pdf(paste0(plotDir, "/pca_clean_no_MW4.pdf"), width=12, height=8)
pca <- princomp(df_4[1:6])
plot(pca$loading, pch=19, cex=2, col=colrs)
text(pca$loading, names(geneResults_4), pos=1)
dev.off()

# visualise clean pca components on plot, save as pdf:
pdf(paste0(plotDir, "/pca_components_clean.pdf"),width=12,height=8)
plot(pca)
dev.off()

# make alternate arrangement of df for following edgeR tut, with all ctls first followed
#by all treated:
df_alt <- df
# store all treated columns in additional columns:
df_alt[,12] <- df_alt[,1]
colnames(df_alt)[12] <- colnames(df_alt)[1]
df_alt[,13] <- df_alt[,3]
colnames(df_alt)[13] <- colnames(df_alt)[3]
df_alt[,14] <- df_alt[,5]
colnames(df_alt)[14] <- colnames(df_alt)[5]
df_alt[,15] <- df_alt[,7]
colnames(df_alt)[15] <- colnames(df_alt)[7]

# arrange all columns into new order:
df_alt[,1] <- df_alt[,2]
colnames(df_alt)[1] <- colnames(df_alt)[2]
df_alt[,2] <- df_alt[,4]
colnames(df_alt)[2] <- colnames(df_alt)[4]
df_alt[,3] <- df_alt[,6]
colnames(df_alt)[3] <- colnames(df_alt)[6]
df_alt[,4] <- df_alt[,8]
colnames(df_alt)[4] <- colnames(df_alt)[8]
df_alt[,5] <- df_alt[,12]
colnames(df_alt)[5] <- colnames(df_alt)[12]
df_alt[,6] <- df_alt[,13]
colnames(df_alt)[6] <- colnames(df_alt)[13]
df_alt[,7] <- df_alt[,14]
colnames(df_alt)[7] <- colnames(df_alt)[14]
df_alt[,8] <- df_alt[,15]
colnames(df_alt)[8] <- colnames(df_alt)[15]

# remove extra columns:
df_alt <- df_alt[,1:8]

# change names of treated back to original names:
new_names <- c()
i=1
for (n in colnames(df_alt)) {
  new_names[i] <- sub("2.1", "2", n)
  i=i+1
}
colnames(df_alt) <- new_names

##### edgeR tutorial: #####
# from manual section 4.2 #
# note - df_alt has Entrez gene IDs as rownames

# specify each group type as factors:
Type1 <- str_split(colnames(df_alt)[1:4], "_")[[1]][2]
Type2 <- str_split(colnames(df_alt)[5:8], "_")[[2]][2]
Type <- factor(rep(c(Type1, Type2), each=4))
# make type1 the first level of 'type'
Type <- relevel(Type, ref=Type1)
# specify each animal number:
Animal <- factor(substring(colnames(df_alt), 3, 3))

# put df data in DGEList:
expr_alt <- DGEList(counts=df_alt[1:8], group=type)

# perform a Trimmed Mean of M-values (TMM) normalisation to account for RNA
# compositional biases:
expr_alt <- calcNormFactors(expr_alt)

# sanity check - print norm. factors for each sample library:
print("The norm. factors of each sample are:")
expr_alt$samples

# create an MDS plot to show relative similarities of the samples and save to plotDir:
pdf(paste0(plotDir, "/MDS.pdf"),width=16,height=12)
plotMDS(expr_alt)
dev.off()

# compute predictive log2-fold changes to assess consistency of Egr2 replicates
# and save to tableDir:
design <- model.matrix(~Animal+Animal:Type)
logFC <- predFC(expr_alt, design, prior.count=1, dispersion=0.06)
logFC_cor <- cor(logFC[,5:8])
write.table(logFC_cor, file = paste0(tableDir, "/logFC_correlations.txt", sep="\t"))

# define design matrix for an additive linear model:
design <- model.matrix(~Animal+Type)
rownames(design) <- colnames(expr_alt)
print("The matrix design to be used is:")
design

save.image(paste0(projectDir, "/Robjects/", samplename, "_edgeR_from_rsemCounts_temp.RData"))


###
# load each counts file as a data frame:
#for (file in inFiles) {
#  split <- str_split(file, ".tr")
#  uID <- split[[1]][1]
#  print(paste0("The unique ID of the file used is: ", uID))
#}
###