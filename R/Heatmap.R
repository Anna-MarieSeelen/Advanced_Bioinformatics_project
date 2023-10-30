# # -----------------------------------------------------------------------
# Heatmap of log2(FC)
# Author: Ahkyeong Choe
# # -----------------------------------------------------------------------

library(RColorBrewer)
library(pheatmap)

#setwd("C:/Users/USER/OneDrive - Wageningen University & Research/WUR Lecture/Advanced Bioinformatics/group/DESeq2/araport11_wo_outlier")
#getwd()

# read expression data ----------------------------------------------------

data15  <- read.table("expressions/expressions_araport11_Axenic-Leaf15.tsv", sep ='\t',header = TRUE, row.names = 1)
data48  <- read.table("expressions/expressions_araport11_Axenic-Leaf48.tsv", sep ='\t',header = TRUE, row.names = 1)
data51  <- read.table("expressions/expressions_araport11_Axenic-Leaf51.tsv", sep ='\t',header = TRUE, row.names = 1)
data53  <- read.table("expressions/expressions_araport11_Axenic-Leaf53.tsv", sep ='\t',header = TRUE, row.names = 1)
data70  <- read.table("expressions/expressions_araport11_Axenic-Leaf70.tsv", sep ='\t',header = TRUE, row.names = 1)
data130  <- read.table("expressions/expressions_araport11_Axenic-Leaf130.tsv", sep ='\t',header = TRUE, row.names = 1)
data131  <- read.table("expressions/expressions_araport11_Axenic-Leaf131.tsv", sep ='\t',header = TRUE, row.names = 1)
data434  <- read.table("expressions/expressions_araport11_Axenic-Leaf434.tsv", sep ='\t',header = TRUE, row.names = 1)

# merged data -------------------------------------------------------------

data_padj <- data.frame(data15$padj, 
                        data48$padj, 
                        data51$padj, 
                        data53$padj, 
                        data70$padj, 
                        data130$padj, 
                        data131$padj, 
                        data434$padj)
rownames(data_padj) <- rownames(data51)

data_padj_wo70 <- data.frame(data15$padj, 
                           data48$padj, 
                           data51$padj, 
                           data53$padj, 
                           #data70$padj, 
                           data130$padj, 
                           data131$padj, 
                           data434$padj)

data_logfc <- data.frame(Leaf_15=data15$log2FoldChange,
                         Leaf_48=data48$log2FoldChange,
                         Leaf_51=data51$log2FoldChange,
                         Leaf_53=data53$log2FoldChange,
                         Leaf_70=data70$log2FoldChange,
                         Leaf_130=data130$log2FoldChange,
                         Leaf_131=data131$log2FoldChange,
                         Leaf_434=data434$log2FoldChange)

data_logfc_wo70 <- data.frame(data15$log2FoldChange,
                            data48$log2FoldChange,
                            data51$log2FoldChange,
                            data53$log2FoldChange,
                            #data70$log2FoldChange,
                            data130$log2FoldChange,
                            data131$log2FoldChange,
                            data434$log2FoldChange)
#rownames(data_logfc_wo70) <- rownames(data51)

# DEGs only ---------------------------------------------------------------

cutoff_padj <- apply(data_padj_wo70, 1, max)<0.01
cutoff_logfc <- apply(abs(data_logfc_wo70), 1, min)>1

data <- data_logfc[cutoff_padj & cutoff_logfc, ]

pheatmap(data,
         color=colorRampPalette(c('white','white','yellow','red', 'red'))(50))
