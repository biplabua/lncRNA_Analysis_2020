library(DESeq2)
library(tidyverse)
library(biomaRt)
library(psych)
library(ComplexHeatmap)
library(circlize)
library(corrplot)
library(RColorBrewer)

#Make a table that we can use to convert gene ID to gene symbol.
lncRNA_id_conversion <- read.table("./input_files/lncRNA_ensemble2symble.tsv")
lncRNA_id_conversion$biotype <- "lncRNA"

mRNA_id_converstion <- read.table("./input_files/mRNA_2_symble.tsv")
mRNA_id_converstion$biotype <- "mRNA"

id_conversion <- rbind(lncRNA_id_conversion,mRNA_id_converstion)
colnames(id_conversion) <- c("Geneid", "GeneSymbol", "biotype")

id_conversion <- unique(id_conversion)
#id_conversion <- apply(id_conversion,2,as.character)

#Combined featureCount output file for mRNA and lncRNA
lncRNA_count = readLines("./input_files/count_file.csv")[-1]
lncRNA_count <- read.table(textConnection(lncRNA_count), sep = "\t", row.names = 1, header = TRUE)
lncRNA_count <- lncRNA_count[,6:ncol(lncRNA_count)]
colnames(lncRNA_count) <- c('O3_1', 'C_1', 'O3_2', 'C_2', 'OH_1', 'O3_3', 
                  'OH_2', 'OH_3', 'C_3')


mRNA_counts = readLines("./input_files/mRNA_gene_counts.csv")[-1]
mRNA_counts <- read.csv(textConnection(mRNA_counts), sep = "\t", row.names = 1)
mRNA_counts <- mRNA_counts[,6:ncol(mRNA_counts)]
colnames(mRNA_counts) <- c('O3_1', 'C_1', 'O3_2', 'C_2', 'OH_1', 'O3_3', 
                   'OH_2', 'OH_3', 'C_3')


count <- rbind(mRNA_counts, lncRNA_count)

#Organize metadata for input in colData
meta_data <- data.frame(sample=c('O3_1', 'C_1', 'O3_2', 'C_2', 'OH_1', 'O3_3','OH_2', 'OH_3', 'C_3'),
                        condition=c('O3', 'C', 'O3', 'C', 'OH', 'O3', 'OH', 'OH', 'C'))

#Load count matrix to the DESeqDataSetFromMatrix
dds <- DESeqDataSetFromMatrix(countData=count, 
                              colData=meta_data, 
                              design=~condition)
#Filter raw count data 
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#Run DESeq
dds <- DESeq(dds, test="Wald")

#Compare between condition
CvsO3 <- results(dds, contrast = c("condition", "O3", "C"))
CvsO3 <- CvsO3 %>% as.data.frame() %>% rownames_to_column("Geneid") 
CvsO3 <- CvsO3[,c(1,3,7)]
colnames(CvsO3) <- c("Geneid", "log2FC_CvsO3", "padj_CvsO3")
CvsO3 <- merge(CvsO3, id_conversion, by="Geneid", all=F)
CvsO3$Geneid <- gsub("\\.[0-9]", replacement = "", CvsO3$Geneid)

CvsOH <- results(dds, contrast = c("condition", "OH", "C"))
CvsOH <- CvsOH %>% as.data.frame() %>% rownames_to_column("Geneid")
CvsOH <- CvsOH[,c(1,3,7)]
colnames(CvsOH) <- c("Geneid", "log2FC_CvsOH", "padj_CvsOH")
CvsOH <- merge(CvsOH, id_conversion, by="Geneid", all=F)
CvsOH$Geneid <- gsub("\\.[0-9]", replacement = "", CvsOH$Geneid)

