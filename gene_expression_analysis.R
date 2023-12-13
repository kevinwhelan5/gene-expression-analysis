# ----------------------------------------
# Breast Cancer Gene Expression Analysis 
# Author: Kevin Whelan
# Date: 01/12/23
# ----------------------------------------

# load libraries
library(ggplot2)
library(DESeq2)
library(clusterProfiler)
library(genefilter)
library(BiocParallel)
library(EnhancedVolcano)
register(MulticoreParam(4)) # multicore processing for DESeq2


# --------------
# 2. Untar folder and extract files
# --------------

#folder_name = "brca_tcga_pan_can_atlas_2018.tar"
#untar(folder_name)

new_dir = paste(getwd(), "brca_tcga_pan_can_atlas_2018", sep ="/" )
# set working directory to brca data folder
# setwd(new_dir)



# --------------
# 3.Read the RNASeq file: data_mrna_seq_v2_rsem.txt - provides counts for each gene, per patient
# --------------
data_Rnaseq = read.delim("data_mrna_seq_v2_rsem.txt")

# Delete genes with duplicated Hugo symbol
keep = !duplicated(data_Rnaseq[,1])
data_Rnaseq = data_Rnaseq[keep,]

# set rownames of data_Rnaseq to Hugo symbols
rownames(data_Rnaseq)  = data_Rnaseq[,1]



# --------------
# 4.Read the Patient Data file: data_clinical_patient.txt 
# --------------
data_patient = read.delim("data_clinical_patient.txt")



# --------------
# 5.Read the Copy Number Aberrations Data: data_cna.txt
# --------------

# we are interested in values > 0 which indicate amplification of gene expression
# VALUES: -2 = homozygous deletion; -1 = hemizygous deletion; 0 = neutral / no change; 1 = gain; 2 = high level amplification.
data_cna = read.delim("data_cna.txt")

# find ERBB2 in data_cna - this row will be used to determine overexpression
erbb2_indx = which(data_cna[,1] == 'ERBB2')


# plot histogram to see distribution of cna values
cna_values <- data.frame(value=as.numeric(data_cna[erbb2_indx,-c(1,2)]))

ggplot(cna_values, aes(x=value)) +
  geom_histogram(binwidth=0.5, fill = "aquamarine3") +
  labs(x = "CNA Value", y = "Count", title = "Histogram of CNA values for ERBB2 Gene") +
  theme_light(base_size = 14)



# --------------
# 6.Match the RNASeq patient ids with the CNA ids and the Patient Data ids.
# --------------

# match patients in data_Rnaseq to patients in data_cna using patient ID
rna_cna_id = which(is.element(colnames(data_Rnaseq[,-c(1,2)]), colnames(data_cna[,-c(1,2)])))

# select only the rna cases which have cna data.
rna_cna_sub = data_Rnaseq[,2+rna_cna_id]

# check all patients in rna_can_sub are in cna
no_pats_in_rna_cna_sub_and_cna = sum(is.element(colnames(data_Rnaseq[,2+rna_cna_id]), colnames(data_cna[,-c(1,2)]))) 
# sanity check.This will print an error if the result is not the same.
sanity_check = no_pats_in_rna_cna_sub_and_cna == dim(rna_cna_sub)[2]



# --------------
# 7.Create metadata using the CNA level of ERBB2+ (greater than 0 means amplified).
# --------------

# Create matrix for meta data
meta_erbb2 = matrix(0,length(rna_cna_id),1)

# build meta data based on amplification of ERBB2 gene
# meta_erbb2 = 1 if cna value > 0; meta_erbb2 = 0 otherwise
for (i in 1:length(rna_cna_id)){
  # access the colnames of i
  col_i = colnames(rna_cna_sub)[i]
  # get the index in cna for the same patient
  col_cna = which(colnames(data_cna)==col_i)
  # store if they're amplified.
  meta_erbb2[i,] = 1*(data_cna[erbb2_indx,col_cna]>0)
}


# simple checks to make sure. 
col_i = colnames(rna_cna_sub)[1]

col_cna = which(colnames(data_cna)==col_i)

# sanity check

(data_cna[erbb2_indx,col_cna]>0) == meta_erbb2[1,1]

# see now if a positive meta_erbb2 is amplified.
pos_example = which(meta_erbb2==1)[1]

col_i = colnames(rna_cna_sub)[pos_example]
col_cna = which(colnames(data_cna)==col_i)

# sanity check
(data_cna[erbb2_indx,col_cna]>0) == meta_erbb2[pos_example,1]

# botch checks should print true.

# We will add a title to the metadata.
colnames(meta_erbb2) = 'ERBB2Amp'

# change labels to more user-friendly names
# DESeq2 will choose a reference level for factors based on alphabetical order
# labelling below means 0_not_amplified will be chosen as reference group
meta_erbb2[meta_erbb2 == 1] <- "1_amplified"
meta_erbb2[meta_erbb2 == 0] <- "0_not_amplified"


# --------------
# 8.Normalize data using DESeq2.
# --------------

# transform into integers
rna_cna_sub = round(rna_cna_sub)

# convert to matrix for input to DESeq
assay = as.matrix(rna_cna_sub)
assay[is.na(assay)] = 0  # Impute with zeros the NA
assay[assay<0] = 0 # Set negative values to zero

# Construct DESeqDataSet object from matrix of counts and meta data
# design formula tells which columns in colData specify the experimental design
# for our design we use ERBB2Amp
dds <- DESeqDataSetFromMatrix(countData = assay,
                              colData = meta_erbb2,
                              design = ~ ERBB2Amp)

# pre-filtering to keep only rows that have a count of at least 10 for a minimal number of samples.
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Normalize
dds <- DESeq(dds, parallel = TRUE)



# --------------
# 9. Obtain Differentially Expressed Genes.
# --------------

# Get Results - set significance level to 0.05
res <- results(dds, alpha=0.05)

# results sorted according to p-value
res_padj_ordered <- res[order(res$padj, decreasing = FALSE),]
# show sample of first 10
res_padj_ordered[1:10, ]

# Summary
summary(res)
rownames(res) = data_Rnaseq[keep,1]

# Significantly Differentially Expressed - padj > 0.05 threshold
signif = which(res$padj<0.05)
deg = res[signif,]

# Order DE genes by log2 FC  
lfc_ranked <- deg[order(abs(deg$log2FoldChange), decreasing = TRUE),]
lfc_ranked[1:10, ]

# Separate them - into up regulated and down regulated
dup = deg[deg[,2]>0.,]
# order by decreasing log2FoldChange
dupOrdered <- dup[order(dup$log2FoldChange, decreasing = TRUE),]
dupOrdered[1:5,]

ddown = deg[deg[,2]<0.,]
# order by decreasing log2FoldChange
ddownOrdered <- ddown[order(ddown$log2FoldChange, decreasing = FALSE),]
ddownOrdered[1:5,]

# Volcano plot to show top DE genes
EnhancedVolcano(res,
               lab = rownames(res),
               x = 'log2FoldChange',
               y = 'pvalue',
               pCutoff = 0.05,
               FCcutoff = 1,
               xlim = c(-5, 5),
               ylim = c(0, 80))


# plot the counts for the most significantly expressed gene
plotCounts(dds, gene=which.min(res$padj), intgroup="ERBB2Amp")


# --------------
# 10. Perform a Pathway Enrichment Analysis
# --------------

# For Pathway Enrichment we need Entrez IDs
entrez_all = data_Rnaseq[signif,2]

# selecting rows which are significant, getting geneID for those rows
entrez_up = data_Rnaseq[signif[deg[,2]>0.],2]
entrez_down = data_Rnaseq[signif[deg[,2]<0.],2]

# Do a KEGG pathway over-representation analysis
all_paths =   enrichKEGG(gene = entrez_all, organism = 'hsa', pvalueCutoff = 0.1)
head(all_paths,10)

# dot plot of top 10 over-represented pathways
dotplot(all_paths, showCategory=10)

# Optionally you can divide between up and down.
up_paths = enrichKEGG(gene = entrez_up, organism = 'hsa', pvalueCutoff = 0.2)
down_paths = enrichKEGG(gene = entrez_down, organism = 'hsa', pvalueCutoff = 0.1)

# --------------
# 11. Get the variance stabilised transformed expression values.
# --------------

# Transform the data to visualize
vsd <- vst(dds, blind=FALSE)

# --------------
# 12. With the vst values obtain a PCA plot.
# --------------

# DESeq2 PCA plot
plotPCA(vsd, intgroup=c("ERBB2Amp"), ntop = 500)

# --------------
# 13. Cluster the data and show in PCA
# --------------


# --------------
# 14. With the vst values of the DE genes generate an overall survival model.
# --------------



# --------------
# 15. Use lasso cross validation to find the set of genes which predict survival the best.
# --------------



# Instructions ------------------------------------------------------------

# Download the dataset on: https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018
# Untar the folder and extract the files.
# Read the RNASeq file: data_mrna_seq_v2_rsem.txt
# Read the Patient Data file: data_clinical_patient.txt
# Read the Copy Number Aberrations Data: data_cna.txt
# Match the RNASeq patient ids with the CNA ids and the Patient Data ids.
# Create metadata using the CNA level of ERBB2+ (greater than 0 means amplified).
# Normalize data using DESeq2.
# Obtain Differentially Expressed Genes.
# Perform a Pathway Enrichment Analysis
# Get the variance stabilised transformed expression values.
# With the vst values obtain a PCA plot.


# Optional for Additional Marks
# Cluster the data and show in PCA
# With the vst values of the DE genes generate an overall survival model.
# Use lasso cross validation to find the set of genes which predict survival the best.
