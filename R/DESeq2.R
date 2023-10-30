# # -----------------------------------------------------------------------
# Differential gene expression analysis with DESeq2 package
# Authour: Ahkyeong Choe
# # -----------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

library(DESeq2)
library(corrplot)
library(stringr)
#library(rhdf5)
#library(vctrs)
#library(dplyr)
#library(devtools)

#getwd()
#setwd("C:/Users/USER/OneDrive - Wageningen University & Research/WUR Lecture/Advanced Bioinformatics/group/DESeq2/araport11_wo_outlier")

db <- 'araport11'

# Read expression file ----------------------------------------------------
expression_data <- read.csv2("gene_count_matrix.csv", sep =',', row.names = 1)


# re-setting raw sample name to Leaf number -------------------------------

sample_info <- read.table("sample_info.txt", header = TRUE)
new_sample_names <- sprintf('%s_%s', sample_info$sample_name, sample_info$sample_rep)
colnames(expression_data) <- new_sample_names[order(sample_info$file_name)]
expression_data <- expression_data[,order(colnames(expression_data))]

write.table(expression_data, sprintf('gene_count_matrix_%s.csv', 
                                      db), quote = FALSE, sep = ',')

# Check total reads of each samples ---------------------------------------

apply(expression_data, 2, sum)

# remove outlier ----------------------------------------------------------

out <- c('Axenic_R5',
         'Fr1_R3',
         'Leaf15_R4',
         'Leaf48_R4',
         'Leaf51_R1',
         'Leaf53_R3',
         'Leaf70_R5',
         'Leaf117_R4',
         'Leaf130_R5',
         'Leaf131_R4',
         'Leaf434_R5')

expression_data <- expression_data[, !(colnames(expression_data) %in% out)]

# Remove genes (under 10 expression) --------------------------------------

mx = apply( expression_data, 1, max)
expression_data = expression_data[ mx > 10, ] 

write.table(expression_data, sprintf('gene_count_matrix_%s_trim10.csv', db),
            quote = FALSE, sep = ',')

# correlation -------------------------------------------------------------

cor_matrix <- cor(expression_data)

pdf(sprintf('correlation.pdf'), width = 23, height = 23)
corrplot(cor_matrix, method='number',pch.cex = 0.2)
dev.off()

# Set the sample names ----------------------------------------------------

# control
cont <- 'Axenic'
idx <- str_detect(colnames(expression_data), cont)
matrix_Ncont <- expression_data[,idx]

#cont <- 'Fr1'
#idx <- str_detect(colnames(expression_data), cont)
#matrix_Pcont <- expression_data[,idx]
unique(sample_info$sample_name)
for (treat in unique(sample_info$sample_name)){
  #treat <- unique(sample_info$sample_name)[3]
  if (treat == 'Axenic'){
    next
  } 
  if (treat == 'Fr1'){
    next
  }
  
  idx <- str_detect(colnames(expression_data), treat)
  matrix_treat <- expression_data[,idx]
  subexp_data <- cbind(matrix_Ncont, matrix_treat)
  condition  <-  factor(c(rep(cont,ncol(matrix_Ncont)),
                          rep(treat,ncol(matrix_treat))),
                        c(cont,treat))
  col_data <- data.frame(condition)
  dds <- DESeqDataSetFromMatrix(subexp_data, col_data, ~condition)
  
  # Estimate size factores --------------------------------------------------
  
  dds = estimateSizeFactors(dds) 
  sizeFactors(dds)
  
  rld = rlog(dds)
  
  plot(density(assay(dds)[,1]), main="counts")
  plot(density(assay(rld)[,1]), main="log counts")
  
  # Distance between samples ------------------------------------------------
  
  dists = dist(t(assay(rld)))
  
  pdf(sprintf('dendrogram/dendrogram_%s_%s-%s.pdf', db, cont, treat)
      , width = 4, height = 4)
  plot(hclust(dists))  
  dev.off()
  
  # DEseq result ------------------------------------------------------------
  dds = estimateDispersions(dds)
  
  pdf(sprintf('dispersion_estimates/dispersion_estimates_%s_%s-%s.pdf', 
              db, cont, treat)
      , width = 10, height = 10)
  plotDispEsts(dds) 
  dev.off()
  
  dds = nbinomWaldTest(dds)
  res = results(dds)
  
  # Replace na ---------------------------------------------------------------
  res$padj = ifelse(is.na(res$padj), 1, res$padj)
  write.table(res, col.names=NA, row.names=T, 
              file =sprintf("expressions/expressions_%s_%s-%s.tsv", 
                            db, cont, treat), 
              sep ="\t", quote = FALSE)
  
  # MA plot -----------------------------------------------------------------
  
  pdf(sprintf('MAplot/MAplot_%s_%s-%s.pdf', db, cont, treat)
      , width = 5, height = 5)
  plotMA(res, main="MA plot",ylim=c(-8,8),alpha=0.01)
  dev.off()
  
  # DESeq2- deg -------------------------------------------------------------
  
  deseq2_significant = data.frame(res[res$padj<=0.05,])
  
  ## ------------------------------------------------------------------------
}
