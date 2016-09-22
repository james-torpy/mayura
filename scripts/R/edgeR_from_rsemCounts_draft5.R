##### edgeR_from_rsemCounts #####

# This script takes a transcriptome counts file outputted by RSEM and performs a differential analysis:

# load packages:
library(RColorBrewer)
library(org.Mm.eg.db)
library(edgeR)
library(stringr)
library(dplyr)
library(biomaRt)
library(DT)
library(pheatmap)

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
type <- gsub('MW\\d+_' , '', samples)
type <- as.factor(type)
type <- relevel(type, "Wt")

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

### fit data to matrix model and plot:
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
down_genes_temp <- row.names(y$counts)[dt == -1]
up_genes_temp <- row.names(y$counts)[dt == 1]

# create table of gene symbols:
egSYMBOL <- toTable(org.Mm.egSYMBOL)

# display the down- and up-regulated genes:
down_genes <- dplyr::filter(egSYMBOL, gene_id %in% down_genes_temp)
print("Downregulated genes are:")
down_genes
up_genes <- dplyr::filter(egSYMBOL, gene_id %in% up_genes_temp)
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

# display the cpm of top differentially expressed genes in
# individual samples:
top <- rownames(topTags(lrt))
top_cpm <- as.data.frame(cpm(y)[top,])
col_no <- ncol(top_cpm)
top_cpm[,col_no+1] <- rownames(top_cpm)
top_cpm <- merge(top_cpm, symbols, by.x='V9', by.y='gene_id')
rownames(top_cpm) <- top_cpm$symbol
top_cpm <- top_cpm[,2:col_no+1]

# write top_cpm to table:
write.table(top_cpm, file = paste0(tableDir, "/top_cpm.txt"), sep = "\t")
top_cpm_table <- datatable(top_cpm)
saveWidget(top_cpm_table, paste0(tableDir, "/top_cpm.html"))

# add function for adding descriptions and/or ENSEMBL IDs from bioMart to gene dataframes:
add_descr <- function(x, ENSEMBL = FALSE) {
  mart = useMart("ENSEMBL_MART_ENSEMBL")
  mart_dset = useDataset("mmusculus_gene_ensembl", mart = mart)
 
  # define number of columns in x, and the number of a potential extra column:
  col_no <- ncol(x)
  extra_col <- col_no+1
  extra_col2 <- extra_col+1
  
  i=1
  for (id in x$gene_id) {
    # get attributes for each gene_id:
    bm <- getBM(attributes=c("entrezgene", "ensembl_gene_id", "description"), filters = "entrezgene", values=id, mart = mart_dset)
    if (nrow(bm)>0 & nrow(bm)<2) {
      print(i)
      print(bm)
      if (!ENSEMBL) {
        # add gene descriptions to extra column:
        x[i,extra_col] <- bm$description
        # name extra columns:
        colnames(x)[extra_col] <- "description"
      
        i=i+1
      
      } else {
        # add ENSEMBL IDs to extra column:
        x[i,extra_col] <- bm$ensembl_gene_id
        # add gene descriptions to extra column:
        x[i,extra_col2] <- bm$description
        # name extra columns:
        colnames(x)[extra_col] <- "ensembl_gene_id"
        colnames(x)[extra_col2] <- "description"
      
        i=i+1
      }
    }
  }
  return(x)
}


# add descriptions and ENSEMBL IDs to top_genes from bioMart:
top_genes <- add_descr(top_genes, ENSEMBL = TRUE)

# add a straight FC column:
top_genes[,8] <- 2^top_genes$logFC
colnames(top_genes)[8] <- "fold count"
top_genes <- top_genes[,c(1, 8, 2:7)]
print("Top differentially expressed genes are:")
top_genes

# create function to count the number of CPM values equal to 0 for each row of gene dataframes, and make separate lists of genes
# with one, two, three, four counts of CPM values equal to 0
count_0s <- function(x, type = "Wt") {
  if (type == "Wt") {
    Wt_0_count=0
    for (value in x[1:4]) {
      if (value == 0) {
        Wt_0_count = Wt_0_count + 1
      }
    }
    return(Wt_0_count)
  }
  
  if (type == "Egr2") {
    Egr2_0_count=0
    for (value in x[5:8]) {
      if (value == 0) {
        Egr2_0_count <- Egr2_0_count + 1
      }
    }
    return(Egr2_0_count)
  }
}

# remove log values derived from genes with zero CPMs in up_genes:
# create a dataframe of all cpm values and group Wt and Egr2 columns together:
all_cpms <- cpm(y)
all_cpms <- all_cpms[,c(2, 4, 6, 8, 1, 3, 5, 7)]

# make vector of CPM=0 counts for each gene for Wt and Egr2 samples:
Wt_0_counts <- c()
i=1
for (row in rownames(all_cpms)) {
  Wt_0_counts[i] <- count_0s(all_cpms[row,], type = "Wt")
  i=i+1
}

Egr2_0_counts <- c()
i=1
for (row in rownames(all_cpms)) {
  Egr2_0_counts[i] <- count_0s(all_cpms[row,], type = "Egr2")
  i=i+1
}

# add CPM=0 counts to dataframes with a gene_id column:
Wt_0_counts <- data.frame(rownames(all_cpms), Wt_0_counts)
colnames(Wt_0_counts)[1] <- "gene_id"

Egr2_0_counts <- data.frame(rownames(all_cpms), Egr2_0_counts)
colnames(Egr2_0_counts)[1] <- "gene_id"

# add column with CPM=0 counts to top_genes:
top_genes <- merge(top_genes, Wt_0_counts, by.x='gene_id', by.y='gene_id')
top_genes <- merge(top_genes, Egr2_0_counts, by.x='gene_id', by.y='gene_id')

# make function to remove log values from gene dataframes of rows with CPM=0 counts of >=3 for Wt, Egr2, or both:
remove_bad_log <- function(x) {
  for (row in rownames(x)) {
    if (x[row, "Wt_0_counts"] >= 3 | x[row, "Egr2_0_counts"] >=3) {
      x[row, "fold count"] <- NA
      x[row, "logFC"] <- NA
      x[row, "logCPM"] <- NA
    }
  }
  return(x)
}

# remove log values from top_genes of rows with CPM=0 counts of >=3 for Wt, Egr2, or both:
top_genes <- remove_bad_log(top_genes)

# write top_genes to table:
toptable <- datatable(top_genes)
write.table(top_genes, file = paste0(tableDir, "/top_genes.txt"), sep = "\t")
saveWidget(toptable, paste0(tableDir, "/top_genes.html"))


# add descriptions and ENSEMBL IDs to up_genes from bioMart:
up_genes <- add_descr(up_genes, ENSEMBL = TRUE)

# add logFC and other table info to up_genes:
alltable <- as.data.frame(lrt$table)
alltable[,5]  <- rownames(alltable)
colnames(alltable)[5] <- "gene_id"
up_table <- dplyr::filter(alltable, gene_id %in% up_genes$gene_id)
up_genes <- merge(up_table, up_genes, by.x='gene_id', by.y='gene_id')
rownames(up_genes) <- up_genes$symbol
up_genes <- up_genes[,c(1:5, 7:8)]

# add a straight FC column:
up_genes[,8] <- 2^up_genes$logFC
colnames(up_genes)[8] <- "fold count"
up_genes <- up_genes[,c(1, 8, 2:7)]
print("Top differentially expressed genes are:")
up_genes

# remove log values derived from genes with zero CPMs in up_genes:
# add column with CPM=0 counts to up_genes:
up_genes <- merge(up_genes, Wt_0_counts, by.x='gene_id', by.y='gene_id')
up_genes <- merge(up_genes, Egr2_0_counts, by.x='gene_id', by.y='gene_id')

# sort by logFC:
up_genes_sorted <- up_genes[with(up_genes, order(-logFC)), ]

# remove log values from up_genes of rows with CPM=0 counts of >=3 for Wt, Egr2, or both:
up_genes <- remove_bad_log(up_genes)

# write up_genes to table:
write.table(up_genes_sorted, file = paste0(tableDir, "/up_genes.txt"), sep = "\t")
uptable <- datatable(up_genes_sorted)
saveWidget(uptable, paste0(tableDir, "/up_genes.html"))


# add descriptions and ENSEMBL IDs to up_genes from bioMart:
down_genes <- add_descr(down_genes, ENSEMBL = TRUE)

# add logFC and other table info to down_genes:
alltable <- as.data.frame(lrt$table)
alltable[,5]  <- rownames(alltable)
colnames(alltable)[5] <- "gene_id"
down_table <- dplyr::filter(alltable, gene_id %in% down_genes$gene_id)
down_genes <- merge(down_table, down_genes, by.x='gene_id', by.y='gene_id')
rownames(down_genes) <- down_genes$symbol
down_genes <- down_genes[,c(1:5, 7:8)]

# add a straight FC column:
down_genes[,8] <- 2^down_genes$logFC
colnames(down_genes)[8] <- "fold count"
down_genes <- down_genes[,c(1, 8, 2:7)]
print("Top differentially expressed genes are:")
down_genes

# remove log values derived from genes with zero CPMs in down_genes:
# add column with CPM=0 counts to down_genes:
down_genes <- merge(down_genes, Wt_0_counts, by.x='gene_id', by.y='gene_id')
down_genes <- merge(down_genes, Egr2_0_counts, by.x='gene_id', by.y='gene_id')

# sort by logFC:
down_genes_sorted <- down_genes[with(down_genes, order(logFC)), ]

# remove log values from down_genes of rows with CPM=0 counts of >=3 for Wt, Egr2, or both:
down_genes <- remove_bad_log(down_genes)

# write down_genes to table:
write.table(down_genes_sorted, file = paste0(tableDir, "/down_genes.txt"), sep = "\t")
downtable <- datatable(down_genes_sorted)
saveWidget(downtable, paste0(tableDir, "/down_genes.html"))


# create dataframe of CPM values for those genes without log values due to CPM=0:
both_0_counts <- merge(Wt_0_counts, Egr2_0_counts, by.x = "gene_id", by.y = "gene_id")
rownames(both_0_counts) <- as.vector(both_0_counts[,1])

test_cpms_for_no_logs <- data.frame(column1 = numeric(0), column2 = numeric(0), column3 = numeric(0), column4 = numeric(0),
                               column5 = numeric(0), column6 = numeric(0), column7 = numeric(0), column8 = numeric(0))
colnames(test_cpms_for_no_logs) <- colnames(test_all_cmps)
i=1
for (row in rownames(test_both_0_counts)) {
  if (test_both_0_counts[row, "Wt_0_counts"] >= 3 | test_both_0_counts[row, "Egr2_0_counts"] >=3) {
    test_cpms_for_no_logs[i, 9] <- row
    test_cpms_for_no_logs[i, 1:8] <- all_cpms[row,]
  }
  i=i+1
}
rownames(test_cpms_for_no_logs) <- test_cpms_for_no_logs[,9]
colnames(test_cpms_for_no_logs)[9] <- "gene_id"

up_cpms_for_no_logs <- merge(test_cpms_for_no_logs, up_genes, by.x = "gene_id", by.y="gene_id")
down_cpms_for_no_logs <- merge(test_cpms_for_no_logs, down_genes, by.x = "gene_id", by.y="gene_id")

# create pheatmap of down_genes cpm:
down <- down_genes$gene_id
down_cpm <- as.matrix(cpm(y)[down,])
rownames(down_cpm) <- rownames(down_genes)
down_cpm <- down_cpm[,c(2, 4, 6, 8, 1, 3, 5, 7)]

pheatmap(down_cpm)

save.image(paste0(projectDir, "/Robjects/", samplename, "_edgeR_from_rsemCounts_temp.RData"))

### tests ###
test <- c(2, 3, 0, 3, 0, 0, 0, 5)
test_count_0s <- function(x) {
  i=0
  for (value in x) {
    if (value == 0)
      i=i+1
    count <- i
  }
  return(count)
}
test_count_0s(test)

###
test2 <- c(2, 3, 0, 3, 0, 0, 0, 5)

test2_count_0s <- function(x, type = "Wt") {
  if (type == "Wt") {
    Wt_0_count=0
    for (value in x[1:4]) {
      if (value == 0) {
        Wt_0_count = Wt_0_count + 1
      }
    }
    return(Wt_0_count)
  }
  
  if (type == "Egr2") {
    Egr2_0_count=0
    for (value in x[5:8]) {
      if (value == 0) {
        Egr2_0_count <- Egr2_0_count + 1
      }
    }
    return(Egr2_0_count)
  }
}
test2_count_0s(test2, "Wt")
test2_count_0s(test2, "Egr2")

###
test3 <- data.frame(c(1, 2, 3), c(0, 2, 0), c(3, 0 ,0), c(0, 0, 0), c(1, 2, 3), c(0, 6, 3), c(0, 0, 2), c(3, 4, 5))
colnames(test3) <- c("Wt1", "Wt2", "Wt3", "Wt4", "Egr2_1", "Egr2_2", "Egr2_3", "Egr2_4")

test3_count_0s <- function(x, type = "Wt") {
  if (type == "Wt") {
    Wt_0_count=0
    for (value in x[1:4]) {
      if (value == 0) {
        Wt_0_count = Wt_0_count + 1
      }
    }
    return(Wt_0_count)
  }
  
  if (type == "Egr2") {
    Egr2_0_count=0
    for (value in x[5:8]) {
      if (value == 0) {
        Egr2_0_count <- Egr2_0_count + 1
      }
    }
    return(Egr2_0_count)
  }
}

for (row in rownames(test3)) {
  Wt_count <- test3_count_0s(test3[row,], type = "Wt")
  print(Wt_count)
}

for (row in rownames(test3)) {
  Egr2_count <- test3_count_0s(test3[row,], type = "Egr2")
  print(Egr2_count)
}

###
