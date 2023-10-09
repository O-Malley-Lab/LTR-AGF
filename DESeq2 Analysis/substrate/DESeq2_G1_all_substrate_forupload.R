# script to perform differential gene expression analysis using DESeq2 package
# setwd("~/Desktop/demo/DESeq2_tutorial/data")

# load libraries
library(DESeq2)
library(tidyverse)
library(airway)
library(EnhancedVolcano) #not needed if you just are looking for the DESeq outputs without Volcanoplot

packageVersion("DESeq2")
#change this line to YOUR working directory
setwd("/Users/tejasn/Documents/OMalley_lab/Code/G1_all_substrates")

# Step 1: preparing count data ----------------
counts_data_GvsAll <- read.csv('G1_substrate_raw_expectedcounts.csv', row.names = NULL)
row.names(counts_data_GvsAll) <- counts_data_GvsAll[, 1]
counts_data_GvsAll <- counts_data_GvsAll[, -1]

counts_data_GvsAll<-round(counts_data_GvsAll)

# read in sample info
colData_GvsAll <- read.csv('all_substrate_columns.csv')
row.names(colData_GvsAll) <- colData_GvsAll[, 1]
colData_GvsAll <- colData_GvsAll[, -1]

# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data_GvsAll) %in% rownames(colData_GvsAll))

# are they in the same order?
all(colnames(counts_data_GvsAll) == rownames(colData_GvsAll))


# Step 2: construct a DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = counts_data_GvsAll,
                       colData = colData_GvsAll,
                       design = ~ substrate)

dds

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds


# set the factor level
dds$substrate <- relevel(dds$substrate, ref = "G")

# NOTE: collapse technical replicates

# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)

resA = results(dds, contrast = c("substrate","A","G"))
resAS = results(dds, contrast = c("substrate","AS","G"))
resCB = results(dds, contrast = c("substrate","CB","G"))
resCS = results(dds, contrast = c("substrate","CS","G"))
resM = results(dds, contrast = c("substrate","M","G"))
resRCG = results(dds, contrast = c("substrate","RCG","G"))
resSG = results(dds, contrast = c("substrate","SG","G"))

write.csv(resA, file = "AvsGresults_fromall.csv")
write.csv(resAS, file = "ASvsGresults_fromall.csv")
write.csv(resCB, file = "CBvsGresults_fromall.csv")
write.csv(resCS, file = "CSvsGresults_fromall.csv")
write.csv(resM, file = "MvsGresults_fromall.csv")
write.csv(resRCG, file = "RCGvsGresults_fromall.csv")
write.csv(resSG, file = "SGvsGresults_fromall.csv")

#stop here if you don't need merged csv files (with raw data and DESeq analysis together)

data_A = data.frame(resA)
  colnames(data_A)[colnames(data_A) == "log2FoldChange"] ="log2FoldChange_A"
  colnames(data_A)[colnames(data_A) == "padj"] ="padj_A"
data_AS = data.frame(resAS)
  colnames(data_AS)[colnames(data_AS) == "log2FoldChange"] ="log2FoldChange_AS"
  colnames(data_AS)[colnames(data_AS) == "padj"] ="padj_AS"
data_CB = data.frame(resCB)
    colnames(data_CB)[colnames(data_CB) == "log2FoldChange"] ="log2FoldChange_CB"
    colnames(data_CB)[colnames(data_CB) == "padj"] ="padj_CB"
data_CS = data.frame(resCS)
    colnames(data_CS)[colnames(data_CS) == "log2FoldChange"] ="log2FoldChange_CS"
    colnames(data_CS)[colnames(data_CS) == "padj"] ="padj_CS"
data_M = data.frame(resM)
    colnames(data_M)[colnames(data_M) == "log2FoldChange"] ="log2FoldChange_M"
    colnames(data_M)[colnames(data_M) == "padj"] ="padj_M"
data_RCG = data.frame(resRCG)
    colnames(data_RCG)[colnames(data_RCG) == "log2FoldChange"] ="log2FoldChange_RCG"
    colnames(data_RCG)[colnames(data_RCG) == "padj"] ="padj_RCG"
data_SG = data.frame(resSG)
    colnames(data_SG)[colnames(data_SG) == "log2FoldChange"] ="log2FoldChange_SG"
    colnames(data_SG)[colnames(data_SG) == "padj"] ="padj_SG"


        
subset_A = subset(data_A,select = -c(baseMean,lfcSE,stat,pvalue))    
subset_AS = subset(data_AS,select = -c(baseMean,lfcSE,stat,pvalue)) 
subset_CB = subset(data_CB,select = -c(baseMean,lfcSE,stat,pvalue))
subset_CS = subset(data_CS,select = -c(baseMean,lfcSE,stat,pvalue))
subset_M = subset(data_M,select = -c(baseMean,lfcSE,stat,pvalue))    
subset_RCG = subset(data_RCG,select = -c(baseMean,lfcSE,stat,pvalue))
subset_SG = subset(data_SG,select = -c(baseMean,lfcSE,stat,pvalue))

subset_A$locustag = rownames(subset_A)
subset_AS$locustag = rownames(subset_A)
subset_CB$locustag = rownames(subset_A)
subset_CS$locustag = rownames(subset_A)
subset_M$locustag = rownames(subset_A)
subset_RCG$locustag = rownames(subset_A)
subset_SG$locustag = rownames(subset_A)

#A,AS  CB, CS.  M.  RCG
#A+AS+SG CB + CS + M + RCG
#all combination

subset_A_AS = merge(subset_A,subset_AS,by = "locustag",all=T)
subset_CS_M = merge(subset_CS,subset_M,by = "locustag",all=T)
subset_RCG_SG = merge(subset_RCG,subset_SG,by = "locustag", all=T)

subset_A_AS_CB = merge(subset_A_AS,subset_CB,by = "locustag",all=T)
subset_CS_M_RCG_SG = merge(subset_CS_M,subset_RCG_SG,by = "locustag",all=T)
subset_ALL = merge(subset_A_AS_CB,subset_CS_M_RCG_SG,by = "locustag",all=T)



merged_data_A<-merge(data_A,counts_data_GvsAll,by = "row.names")
write.csv(merged_data_A, file="DGE_analysis_G1_substrate_A.csv")

merged_data_AS<-merge(data_AS,counts_data_GvsAll,by = "row.names")
write.csv(merged_data_AS, file="DGE_analysis_G1_substrate_AS.csv")

merged_data_CB<-merge(data_CB,counts_data_GvsAll,by = "row.names")
write.csv(merged_data_CB, file="DGE_analysis_G1_substrate_CB.csv")

merged_data_CS<-merge(data_CS,counts_data_GvsAll,by = "row.names")
write.csv(merged_data_CS, file="DGE_analysis_G1_substrate_CS.csv")

merged_data_M<-merge(data_M,counts_data_GvsAll,by = "row.names")
write.csv(merged_data_M, file="DGE_analysis_G1_substrate_M.csv")

merged_data_RCG<-merge(data_RCG,counts_data_GvsAll,by = "row.names")
write.csv(merged_data_RCG, file="DGE_analysis_G1_substrate_RCG.csv")

merged_data_SG<-merge(data_SG,counts_data_GvsAll,by = "row.names")
write.csv(merged_data_SG, file="DGE_analysis_G1_substrate_SG.csv")

