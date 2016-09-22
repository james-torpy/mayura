### featureCounts_for_edgeR ###
### incomplete ###

# This script takes a .bam file from a STAR alignment and counts overlaps with
# exons from an annotation file

# load packages:
library(Rsubread)
library(stringr)

# set up directory structure:
projectname = "mayura"
genome = "GRCm38_p4/"
annotation = "gencode.vM4.annotation"

# if fetching all files from project (not just one sample type), leave 'samplename' blank (i.e. "")
samplename = "MW_310816"
inType = ".star"


# use following homeDir for R shell, files on cluster:
#homeDir = "/home/jamtor"
# use following homeDir for R Studio, files on cluster:
#homeDir = "/Users/jamestorpy/clusterHome"
# use following homeDir for files on local drive:
homeDir = "/Users/jamestorpy/Documents/Garvan/phd"
genomeDir = paste0(homeDir, "/genomes/", genome)
annotationFile = paste0(genomeDir, annotation, ".gtf")

projectDir = paste0(homeDir, "/projects/", projectname)
resultsDir = paste0(projectDir, "/results/")
inDir = paste0(resultsDir, samplename, inType)
paste("The annotationFile is:", annotationFile)
paste0("The inDir is: ", inDir)
setwd(inDir)

# fetch inFiles:
inFiles = list.files(inDir, pattern = ".transcriptome.bam", recursive = T)
inFiles = grep("subset", inFiles, inv=T, value=T)
paste0("The inFiles are: ")
inFiles

# count overlaps between inFile reads and annotationFile:
for (file in inFiles) {
  split <- str_split(file, "_")
  uID <- split[[1]][1]
  print(paste0("The unique ID of the file used is: ", uID))
  
  counts <- featureCounts(files=paste0(inDir, file),annot.ext=annotationFile, isGTFAnnotationFile=TRUE
  ,GTF.featureType="exon",GTF.attrType="gene_id")
  
  row_no <- nrow(counts)
  
  assign(paste0(uID, "_exon_counts"))
}

save.image(paste0(projectDir, "/Robjects/", samplename, "_featureCounts_temp.RData"))

###subset test:###
#counts <- featureCounts(files="/Users/jamestorpy/Documents/Garvan/phd/projects/mayura/results/subset.star/MW1A_subset.transcriptome.bam"
#          ,annot.ext="/Users/jamestorpy/Documents/Garvan/phd/genomes/GRCm38_p4/gencode.vM4.annotation.gtf", isGTFAnnotationFile=TRUE
#          ,GTF.featureType="exon",GTF.attrType="gene_id")

# amalgamate all counts into a data frame with columns = each data set, rows = gene IDs:
# create empty data frame to put counts into:
all_counts <- data.frame(column1 = numeric(row_no))

i = 1
for (file in inFiles) {
  # fetch unique IDs of sample data frames:
  split <- str_split(file, "_")
  uID <- split[[1]][1]
  print(paste0("The unique ID of the data frame used is: ", uID))
  # fetch each data frame:
  df <- get(paste0(uID, "_exon_counts")$counts)
  
  # name rows of all_counts data frame:
  rownames(all_counts) <- rownames(df)
  # name column of data frame for counts for current sample:
  colnames(all_counts)[i] <- uID
  
  # merge counts into column on all_counts:
  all_counts[,i] <- df
}
  


all_counts <- data.frame()

i = 1
for (file in testFiles) {
  baseFile <- basename(file)
  # fetch unique IDs of sample data frames:
  split <- str_split(baseFile, "_")
  uID <- split[[1]][1]
  print(paste0("The unique ID of the data frame used is: ", uID))
  # fetch each data frame:
  df <- counts$counts
  
  # merge counts into column on all_counts:
  all_counts[,1] <- df[,1]
  
  # name rows of all_counts data frame:
  rownames(all_counts) <- rownames(df)
  # name column of data frame for counts for current sample:
  colnames(all_counts)[i] <- uID
  
  i = i+1
}
