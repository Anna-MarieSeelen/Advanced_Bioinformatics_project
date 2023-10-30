# # -----------------------------------------------------------------------
# Generate full output file from the output of DESeq2 with DEGs cutoff
# Author: Ahkyeong Choe
# # -----------------------------------------------------------------------

library(dplyr)
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


# padj matrix -------------------------------------------------------------

data_padj <- data.frame(data15$padj, 
                        data48$padj, 
                        data51$padj, 
                        data53$padj, 
                        data70$padj, 
                        data130$padj, 
                        data131$padj, 
                        data434$padj)

data_padj_wo70 <- data.frame(data15$padj, 
                           data48$padj, 
                           data51$padj, 
                           data53$padj, 
                           #data70$padj, 
                           data130$padj, 
                           data131$padj, 
                           data434$padj)

# logfc matrix ------------------------------------------------------------

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



# cutoff ------------------------------------------------------------------

ID <- rownames(data15)

data_padj_05 <- data.frame(ID, padj_05=apply(data_padj, 1, max)<0.05)
data_padj_01 <- data.frame(ID, padj_01=apply(data_padj, 1, max)<0.01)
data_logfc_1 <- data.frame(ID, logfc_1=apply(abs(data_logfc), 1, min)>1)
data_logfc_2 <- data.frame(ID, logfc_2=apply(abs(data_logfc), 1, min)>2)

data_padj_05_wo70 <- data.frame(ID, padj_05_wo70=apply(data_padj_wo70, 1, max)<0.05)
data_padj_01_wo70 <- data.frame(ID, padj_01_wo70=apply(data_padj_wo70, 1, max)<0.01)
data_logfc_1_wo70 <- data.frame(ID, logfc_1_wo70=apply(abs(data_logfc_wo70), 1, min)>1)
data_logfc_2_wo70 <- data.frame(ID, logfc_2_wo70=apply(abs(data_logfc_wo70), 1, min)>2)

# get 24 genes ------------------------------------------------------------

gene24_id <- c("AT2G19190","AT5G64120","AT1G65500","AT3G46280","AT2G39200",
               "AT2G43620","AT1G26380","AT1G26410","AT1G02920","AT1G02930",
               "AT1G21110","AT1G26420","AT1G21120","AT4G12500","AT4G12490",
               "AT2G30750","AT5G24110","AT2G25470","AT4G23140","AT1G35230",
               "AT4G23220","AT1G76930","AT1G64170","AT4G28420")

data_gnsr <- data.frame(ID=gene24_id, GNSR='GNSR')

# merge -------------------------------------------------------------------

data_padj_id <- data.frame(ID, data_padj)
data_logfc_id <- data.frame(ID, data_logfc)

data <- full_join(data_padj_id, data_logfc_id)
data <- full_join(data, data_padj_05)
data <- full_join(data, data_padj_01)
data <- full_join(data, data_logfc_1)
data <- full_join(data, data_logfc_2)
data <- full_join(data, data_gnsr)
data <- full_join(data, data_padj_05_wo70)
data <- full_join(data, data_padj_01_wo70)
data <- full_join(data, data_logfc_1_wo70)
data <- full_join(data, data_logfc_2_wo70)

write.table(data, 'DESeq2_fulldata.csv', quote = FALSE, sep = ',', row.names = FALSE)
