##### edgeR_from_rsemCounts #####

# This script takes a transcriptome counts file outputted by RSEM and performs a differential analysis:

# load packages:
library(RColorBrewer)
library(org.Mm.eg.db)
library(edgeR)
library(stringr)
library(dplyr)
library(biomaRt)

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

# eliminate genes with symbol 'NA':
df <-dplyr::filter(df, !is.na(EntrezGene) )

# order rows, grouping rows with the same gene ids, and take the mean
# of rows with the same gene ids:
df <- dplyr::group_by(df, EntrezGene) %>% summarise_each(funs(mean))

# eliminate lowly expressed genes (rows where there are less than 1 counts >=4):
df <- dplyr::filter(df, rowSums(df > 1) >= 4 )

# make the row names entrez gene IDs
df <- as.data.frame(df)
rownames(df) <- df$EntrezGene 
df <- dplyr::select(df, -EntrezGene)

# define PCA plot colour scheme:
colrs <- brewer.pal(6, "Paired")
# create PCA plot for clean data, save as pdf:
pdf(paste0(plotDir, "/pca_clean.pdf"), width=12, height=8)
#pca <- princomp(df[1:8])
pca <- princomp(df[1:8])
plot(pca$loading, pch=19, cex=2, col=colrs)
text(pca$loading, names(geneResults), pos=1)
dev.off()

# define sample types for DGEList:
samples <- colnames(df)
type <- gsub('MW\\d+_' , '', samples )
type <- as.factor(type)

# make a DGEList from dr for differential analysis:
y <- DGEList(counts = df, group = type)

# normalise for library size:
y <- calcNormFactors(y)

# create an MDS plot to show relative similarities of the samples and save to plotDir:
pdf(paste0(plotDir, "/MDS.pdf"),width=16,height=12)
plotMDS(y)
dev.off()

# design matrix reflecting experimental design on, to base
# differential expression on:
design <- model.matrix(~0 + type)

# calculate dispersion of data:
y <- estimateDisp(y, design = design)
print("The common dispersion is:")
y$common.dispersion
par(mar=c(2, 2, 2, 2))
pdf(paste0(plotDir, "/dispersion.pdf"),width=7,height=6)
plotBCV(y)
dev.off()

### fit data to model (which model?) and plot ###
fit <- glmFit(y, design=design, robust=TRUE)

### calculate differential expression (how?) ###
lrt <- glmLRT(fit, contrast = c(-1,1))

# display counts of differentially expressed genes:
summary(dt <- decideTestsDGE(lrt))

# display differentially expressed gene IDs:
isDE <- as.logical(dt)
DEnames <- rownames(y)[isDE]

# plot logFCs against average count size:
plotSmear(lrt, de.tags=DEnames)
abline(h=c(-1,1), col="blue")

# determine down- and up-regulated genes:
down_genes_temp <- row.names(y$counts)[dt == 1]
up_genes_temp <- row.names(y$counts)[dt == -1]

# create table of gene symbols:
egSYMBOL <- toTable(org.Mm.egSYMBOL)

# display the down- and up-regulated genes:
down_genes <- dplyr::filter(egSYMBOL, gene_id %in% down_genes)
print("Downregulated genes are:")
down_genes
up_genes <- dplyr::filter(egSYMBOL, gene_id %in% down_genes)
print("Upregulated genes are:")
up_genes

# display the top differentially expressed genes:
top_genes <- as.data.frame(topTags(lrt))
top_genes[,5] <- rownames(top_genes)
colnames(top_genes)[5] <- "gene_id"
symbols <- dplyr::filter(egSYMBOL, gene_id %in% top_genes$gene_id)
top_genes <- merge(top_genes, symbols, by.x='gene_id', by.y='gene_id')
rownames(top_genes) <- top_genes$symbol
top_genes <- top_genes[1:5]
print("Top differentially expressed genes are:")
top_genes

# add descriptions to top_genes from bioMart:
mart = useMart("ENSEMBL_MART_ENSEMBL")
mart_dset = useDataset("mmusculus_gene_ensembl", mart = mart)

### fix function!!! ###
add_descr <- function(x, ENSEMBL = FALSE) {
  print("Behold z:")
  print(x)
  mart = useMart("ENSEMBL_MART_ENSEMBL")
  mart_dset = useDataset("mmusculus_gene_ensembl", mart = mart)
  print("Behold mart dataset:")
  mart_dset
  
  # define number of columns in x, and the number of a potential extra column:
  col_no <- ncol(x)
  extra_col <- col_no+1
  extra_col2 <- extra_col+1
  print("Behold column number, extra column number and extra extra column number:")
  print(paste0(col_no, ", ", extra_col, ", ", extra_col2))
  
  i=1
  for (id in x$gene_id) {
    # get attributes for each gene_id:
    bm <- getBM(attributes=c("entrezgene", "ensembl_gene_id", "description"), filters = "entrezgene", values=id, mart = dset)
    print("Behold bm:")
    print(bm)
    
    if (!ENSEMBL) {
      # add gene descriptions to extra column:
      x[i,extra_col] <- bm$description
      print("Behold new z:")
      print(x)
      
      # name extra columns:
      colnames(x)[extra_col] <- "description"
      
      i=i+1
      
    } else {
      # add ENSEMBL IDs to extra column:
      x[i,extra_col] <- bm$ensembl_gene_id
      print("Behold new z:")
      print(x)
      
      # add gene descriptions to extra column:
      x[i,extra_col2] <- bm$description
      print("Behold new z:")
      print(x)
      
      # name extra columns:
      colnames(x)[extra_col] <- "ensembl_gene_id"
      colnames(x)[extra_col2] <- "description"
      print("Behold new z:")
      print(x)
      
      i=i+1
    }
  }
}


### fix function!!! ###
add_descr <- function(x, ENSEMBL = FALSE) {
  # specify mart and dataset to use:
  mart = useMart("ENSEMBL_MART_ENSEMBL")
  mart_dset = useDataset("mmusculus_gene_ensembl", mart = mart)
  
  # define number of columns in x, and the number of a potential extra column:
  col_no <- ncol(x)
  extra_col <- col_no+1
  extra_col2 <- extra_col+1
  
  i=1
  for (id in x$gene_id) {
    # get attributes for each gene_id:
    bm <- getBM(attributes=c("entrezgene", "ensembl_gene_id", "description"), filters = "entrezgene", values=id, mart = dset)

    if (!ENSEMBL) {
      # add ENSEMBL IDs to extra column:
      x[i,extra_col] <- bm$ensembl_gene_id
      # add gene descriptions to extra column:
      x[i,extra_col2] <- bm$description
      # name extra columns:
      colnames(x)[extra_col] <- "ensembl_gene_id"
      colnames(x)[extra_col2] <- "description"
      i=i+1
      
    } else {
      # add gene descriptions to extra column:
      x[i,extra_col] <- bm$description
      # name extra columns:
      colnames(x)[extra_col] <- "description"
      i=i+1
    }
  }
}  
###

i=1
for (id in top_genes$gene_id) {
  print("Info on gene:")
  bm <- getBM(attributes=c("entrezgene", "ensembl_gene_id", "hgnc_symbol", "mgi_symbol", "description"), filters = "entrezgene", values=id, mart = dset)
  print(bm)
  top_genes[i,6] <- bm$description
  i=i+1
}
colnames(top_genes)[6] <- "description"

# display the cpm of top differentially expressed genes in
# individual samples:
top <- rownames(topTags(lrt))
top_cpm <- cpm(y)[top,]
rownames(top_cpm) <- symbols$symbol
top_cpm

for (id in top_genes$gene_id) {
  descr <- getGene()
}

save.image(paste0(projectDir, "/Robjects/", samplename, "_edgeR_from_rsemCounts_temp.RData"))
