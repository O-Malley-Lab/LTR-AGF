# script to perform differential gene expression analysis using DESeq2 package
# setwd("~/Desktop/demo/DESeq2_tutorial/data")

# load libraries
library(DESeq2)
library(tidyverse)
library(airway)

# Step 1: preparing count data ----------------
counts_data <- read.csv('t15.csv', row.names = NULL)
row.names(counts_data) <- counts_data[, 1]
counts_data <- counts_data[, -1]

counts_data<-round(counts_data)

# read in sample info
colData <- read.csv('t15col.csv')
row.names(colData) <- colData[, 1]
colData <- colData[, -1]

# making sure the row names in colData matches to column names in counts_data
all(colnames(counts_data) %in% rownames(colData))

# are they in the same order?
all(colnames(counts_data) == rownames(colData))


# Step 2: construct a DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = counts_data,
                       colData = colData,
                       design = ~ heatshock)

dds

# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds

# set the factor level
dds$heatshock <- relevel(dds$heatshock, ref = "untreated")

# NOTE: collapse technical replicates

# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)
res <- results(dds)

res



# Explore Results ----------------

summary(res)
write.csv(res, file = "t15results.csv")
res0.01 <- results(dds, alpha = 0.01)
summary(res0.01)

# contrasts
resultsNames(dds)

# e.g.: treated_4hrs, treated_8hrs, untreated


# MA plot
plotMA(res)



