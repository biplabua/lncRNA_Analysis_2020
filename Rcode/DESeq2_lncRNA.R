library(DESeq2)
library(tidyverse)
library(biomaRt)
library(psych)
library(ComplexHeatmap)
library(circlize)
library(corrplot)
library(RColorBrewer)


setwd("~/Desktop/Des/lncRNA_Analysis_07_16_2020/")
#Combined id and gene symbol of lncRNA and mRNA
id_conversion <- read.table("./input_files/lncRNA_ensemble2symble.tsv")
cds_id_converstion <- read.table("./input_files/mRNA_2_symble.tsv")

colnames(id_conversion) <- c("Geneid", "GeneSymble")
colnames(cds_id_converstion) <- c("Geneid", "GeneSymble")

id_conversion$biotype <- 'lncRNA'
cds_id_converstion$biotype <- "mRNA"
id_symbol <- rbind(id_conversion,cds_id_converstion)
id_symbol <- unique(id_symbol)
id_symbol <- apply(id_symbol,2,as.character)

#Combined featureCount output file for mRNA and lncRNA
all_content = readLines("./input_files/count_file.csv")[-1]
df <- read.table(textConnection(all_content), sep = "\t", row.names = 1, header = TRUE)#[,6:ncol(df)]

df <- df[,5:ncol(df)]
colnames(df) <- c('length', 'O3_1', 'C_1', 'O3_2', 'C_2', 'OH_1', 'O3_3', 
                  'OH_2', 'OH_3', 'C_3')

all_content1 = readLines("./input_files/mRNA_gene_counts.csv")[-1]

df1 <- read.csv(textConnection(all_content1), sep = "\t", row.names = 1)#[,6:ncol(df)]
df1 <- df1[,5:ncol(df1)]
colnames(df1) <- c('length', 'O3_1', 'C_1', 'O3_2', 'C_2', 'OH_1', 'O3_3', 
                   'OH_2', 'OH_3', 'C_3')

#Count_with_legth has geneid and their length information which will be used for coExpression matrix
count_with_legth <- rbind(df1,df)

#Remove the length information to load into DESeq
count <- count_with_legth[,2:ncol(count_with_legth)]

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
CvsOH <- results(dds, contrast = c("condition", "OH", "C"))
#-------------------------------------------------------------------------------
#Merge with id to symbol conversion table 

CvsO3 <- as.data.frame(CvsO3, stringsAsFactors=FALSE) %>% rownames_to_column("Geneid") #%>% filter(padj < 0.01)
CvsOH <- CvsOH %>% as.data.frame() %>% rownames_to_column("Geneid")  #%>% filter(padj < 0.01)
CvsO3 <- merge(CvsO3, id_symbol, by="Geneid", all=F)
lncRNA_O3 <- CvsO3 %>% filter(biotype=="lncRNA")
CvsOH <- merge(CvsOH, id_symbol, by="Geneid", all=F)
lncRNA_OH <- CvsOH %>% filter(biotype=="lncRNA")
aal <- merge(lncRNA_O3, lncRNA_OH, by="Geneid", all=T)
plot(aal$log2FoldChange.x, aal$log2FoldChange.y, xlab = "O3", ylab = "OH")
##Volcano plot with "significant" genes labeled
par(mar = c(5, 5, 5, 5))
with(lncRNA_O3, plot(log2FoldChange, -log10(padj), pch=20, xlim=c(-5,4), cex = 1.0, cex.axis=2.0, cex.lab=2.0))
with(subset(lncRNA_O3, padj<.05 & log2FoldChange > 1), points(log2FoldChange, -log10(padj), pch=20, col="red", cex = 1.0))
with(subset(lncRNA_O3, padj<.05 & log2FoldChange < -1), points(log2FoldChange, -log10(padj), pch=20, col="green", cex = 1.0))

with(lncRNA_OH, plot(log2FoldChange, -log10(padj), pch=20,  xlim=c(-5,4), cex = 1.0, cex.axis=2.0, cex.lab=2.0))
with(subset(lncRNA_OH, padj<.05 & log2FoldChange > 1), points(log2FoldChange, -log10(padj), pch=20, col="red", cex = 1.0))
with(subset(lncRNA_OH, padj<.05 & log2FoldChange < -1), points(log2FoldChange, -log10(padj), pch=20, col="green", cex = 1.0))

#---------------------------------------------------------------------------------
#Vendiagram
library(VennDiagram)
lncRNA_O3_V <- lncRNA_O3 %>% filter(abs(log2FoldChange) > 1 & padj<.05)
lncRNA_OH_V <- lncRNA_OH %>% filter(abs(log2FoldChange) > 1 & padj<.05)
lncRNA_O3_V_up <- lncRNA_O3 %>% filter(log2FoldChange > 1 & padj<.05)
lncRNA_O3_V_down <- lncRNA_O3 %>% filter(log2FoldChange < -1 & padj<.05)
lncRNA_OH_V_up <- lncRNA_OH %>% filter(log2FoldChange > 1 & padj<.05)
lncRNA_OH_V_down <- lncRNA_OH %>% filter(log2FoldChange < -1 & padj<.05)

lncRNA_O3_V <- lncRNA_O3 %>% filter(abs(log2FoldChange) > 0 & padj<.05)
lncRNA_OH_V <- lncRNA_OH %>% filter(abs(log2FoldChange) > 0 & padj<.05)
lncRNA_O3_V_up <- lncRNA_O3 %>% filter(log2FoldChange > 0 & padj<.05)
write.csv(lncRNA_O3_V_up, file = "./lncRNA_O3_up.csv")
lncRNA_O3_V_down <- lncRNA_O3 %>% filter(log2FoldChange < 0 & padj<.05)
write.csv(lncRNA_O3_V_down, file = "./lncRNA_O3_down.csv")

lncRNA_OH_V_up <- lncRNA_OH %>% filter(log2FoldChange > 0 & padj<.05)
write.csv(lncRNA_OH_V_up, file = "./lncRNA_OH_up.csv")

lncRNA_OH_V_down <- lncRNA_OH %>% filter(log2FoldChange < 0 & padj<.05)
write.csv(lncRNA_OH_V_down, file = "./lncRNA_OH_down.csv")


colors <- c("#6b7fff", "#c3db0f")
venn.diagram(x=list(lncRNA_O3_V$Geneid, lncRNA_OH_V$Geneid), category.names = c("O3", "OH"), 
             filename = "fig1_venn.png", col="black", fill =  colors, cat.col = colors,
             cex = 3, cat.cex = 3)
colors <- c("#6b7fff", "#c3db0f")
venn.diagram(x=list(lncRNA_O3_V_up$Geneid, lncRNA_OH_V_up$Geneid), category.names = c("O3", "OH"), 
             filename = "fig2_venn_up.png", col="black", fill =  colors, cat.col = colors,
             cex = 3, cat.cex = 3, sub.col = c("black", "black"))

colors <- c("#6b7fff", "#c3db0f")
venn.diagram(x=list(lncRNA_O3_V_down$Geneid, lncRNA_OH_V_down$Geneid), category.names = c("O3", "OH"), 
             filename = "fig3_venn_down.png", col="black", fill =  colors, cat.col = colors,
             cex = 3, cat.cex = 3)


#Remove version information from the geneid
CvsO3$Geneid1 <- gsub("\\.[0-9]", replacement = "", CvsO3$Geneid)
CvsOH$Geneid1 <- gsub("\\.[0-9]", replacement = "", CvsOH$Geneid)

#Saving log2FC file as csv file
write.csv(CvsO3, file = "./log2FC/CvsO3.csv")
write.csv(CvsOH, file = "./log2FC/CvsOH.csv")
#----------------------------------------------------------------------------------------------------------
library(fgsea)
#Gene set enrichment analysis of cancer related genes
negative_lncRNA <- read.table("./input_files/negative.lncRNA.glist.xls.txt")
positiveative_lncRNA <- read.table("./input_files/positive.lncRNA.glist.xls.txt")
all_lncRNA <- rbind(negative_lncRNA, positiveative_lncRNA)
colnames(all_lncRNA) <- c("Geneid1")
genesets <- list(Cancer_Related=all_lncRNA$Geneid)
#gsea_analysis <- function()
gseaDat <- CvsO3 %>% filter(!is.na(Geneid)) %>% filter(log2FoldChange > 0)
ranks <- gseaDat$log2FoldChange
names(ranks) <- gseaDat$Geneid
ranks <- sort(ranks, decreasing = T)
fgseaRes <- fgsea(genesets, ranks, minSize=15, maxSize = 500, nperm=1000)
plotEnrichment(genesets[["Cancer_Related"]], ranks)

gseaDat <- CvsO3 %>% filter(!is.na(Geneid)) %>% filter(log2FoldChange < 0)
ranks <- gseaDat$log2FoldChange
names(ranks) <- gseaDat$Geneid1
ranks <- sort(ranks, decreasing = T)
fgseaRes <- fgsea(genesets, ranks, minSize=15, maxSize = 500, nperm=1000)
plotEnrichment(genesets[["Cancer_Related"]], ranks)

gseaDat <- CvsOH %>% filter(!is.na(Geneid1)) %>% filter(log2FoldChange>0)
ranks <- gseaDat$log2FoldChange
names(ranks) <- gseaDat$Geneid1
ranks <- sort(ranks, decreasing = T)
barplot(sort(ranks, decreasing = T))
fgseaRes <- fgsea(genesets, ranks, minSize=15, maxSize = 500, nperm=1000)
plotEnrichment(genesets[["Cancer_Related"]], ranks)

gseaDat <- CvsOH %>% filter(!is.na(Geneid1)) %>% filter(log2FoldChange<0)
ranks <- gseaDat$log2FoldChange
names(ranks) <- gseaDat$Geneid1
ranks <- sort(ranks, decreasing = T)
barplot(sort(ranks, decreasing = T))
fgseaRes <- fgsea(genesets, ranks, minSize=15, maxSize = 500, nperm=1000)
plotEnrichment(genesets[["Cancer_Related"]], ranks)


#------------------------------------------------------------------------------------------------

#Analysis of genes regulated by lncRNA in cis
cis_mRNA <- read.table("./input_files/genelncRNA.bed")[,c(4,10)]
colnames(cis_mRNA) <- c("Geneid", "lncRNA_gene")

CvsO3_lncRNA <- CvsO3 %>% filter(biotype=="lncRNA") %>% filter(Geneid %in% cis_mRNA$Geneid)
CvsO3_lncRNA <- merge(CvsO3_lncRNA, cis_mRNA, by="Geneid", all = FALSE) 

CvsO3_lncRNA <- CvsO3_lncRNA %>% rename(lncRNA=Geneid, Geneid=Gene)

CvsO3_mRNA <- CvsO3 %>% filter(biotype=="mRNA") %>% filter(Geneid %in% cis_mRNA$Gene) 
colnames(cis_mRNA) <- c("Geneid", "lncRNA")

CvsO3_mRNA <- merge(CvsO3_mRNA, cis_mRNA, by="Geneid", all = FALSE) 
mRNA_lncRNA_lg2FC <- merge(CvsO3_mRNA, CvsO3_lncRNA, by="Geneid", all=F)
mRNA_lncRNA_lg2FC <- mRNA_lncRNA_lg2FC %>% filter(padj.x < 0.01 & padj.y < 0.01)%>% distinct()
write.csv(mRNA_lncRNA_lg2FC, file = "mRNA_lncRNA_lg2FC_cis.csv")
mRNA_lncRNA_lg2FC1 <- mRNA_lncRNA_lg2FC[,c(1, 3, 14)] %>% distinct()
colnames(mRNA_lncRNA_lg2FC1) <- c("Geneid", "lncRNA_FC_O3", "mRNA_FC_O3") 
write.csv(mRNA_lncRNA_lg2FC1, file = "./mRNA_lncRNA_lg2FC_O3.csv")


cis_mRNA <- read.table("./input_files/genelncRNA.bed")[,c(4,10)]
colnames(cis_mRNA) <- c("Gene", "Geneid")

CvsOH_lncRNA <- CvsOH %>% filter(biotype=="lncRNA") %>% filter(Geneid %in% cis_mRNA$Geneid)
CvsOH_lncRNA <- merge(CvsOH_lncRNA, cis_mRNA, by="Geneid", all = FALSE) 
CvsOH_lncRNA <- CvsOH_lncRNA %>% rename(lncRNA=Geneid, Geneid=Gene)
CvsOH_mRNA <- CvsOH %>% filter(biotype=="mRNA") %>% filter(Geneid %in% cis_mRNA$Gene) 
colnames(cis_mRNA) <- c("Geneid", "lncRNA")

CvsOH_mRNA <- merge(CvsOH_mRNA, cis_mRNA, by="Geneid", all = FALSE) 

mRNA_lncRNA_lg2FC_OH <- merge(CvsOH_mRNA, CvsOH_lncRNA, by="Geneid", all=F)
mRNA_lncRNA_lg2FC_OH <- mRNA_lncRNA_lg2FC_OH %>% filter(padj.x < 0.01 & padj.y < 0.01)
write.csv(mRNA_lncRNA_lg2FC_OH, file = "mRNA_lncRNA_lg2FC_OH_cis.csv")

mRNA_lncRNA_lg2FC1_OH <- mRNA_lncRNA_lg2FC_OH[,c(1, 3, 14)] %>% distinct()
colnames(mRNA_lncRNA_lg2FC1_OH) <- c("Geneid", "lncRNA_FC_OH", "mRNA_FC_OH") 
write.csv(mRNA_lncRNA_lg2FC1_OH, file = "./mRNA_lncRNA_lg2FC_OH.csv")


colors <- c("#6b7fff", "#c3db0f")
venn.diagram(x=list(mRNA_lncRNA_lg2FC1$Geneid, mRNA_lncRNA_lg2FC1_OH$Geneid), category.names = c("O3", "OH"), 
             filename = "fig3_venn.png", col="black", fill =  colors, cat.col = colors, cex = 3, cat.cex = 3)

df_fig3 <- merge(mRNA_lncRNA_lg2FC1, mRNA_lncRNA_lg2FC1_OH, by="Geneid", all = FALSE)
colnames(df_fig3) <- c("Geneid", "lncRNA O3", "Nearest Gene O3", "lncRNA OH", "Nearest Gene OH")
df_fig3 <- df_fig3 %>% distinct(Geneid, .keep_all = T)
mat = as.matrix(df_fig3[,c(2,3,4,5)])
Heatmap(mat, name = "log2FC", cluster_columns = FALSE)
common_trans_gene <- intersect(mRNA_lncRNA_lg2FC1$Geneid, mRNA_lncRNA_lg2FC1_OH$Geneid)
library(biomaRt)
library(fgsea)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
genes <- getBM(
  filters="ensembl_gene_id",
  attributes=c("ensembl_gene_id", "entrezgene_id"),
  values=common_trans_gene,
  mart=mart)
colnames(genes) <- c("Geneid1", "Entrez")

library(goseq)
Go_CvsO3 <- CvsO3 %>% filter(biotype=="mRNA")
Go_CvsO3 <- Go_CvsO3 %>% distinct(Geneid, .keep_all = T)

supportedOrganisms() %>% filter(str_detect(Genome, "hg19"))
isSigGene <- unique(Go_CvsO3$Geneid) %in% unique(common_trans_gene )
names(isSigGene) <- unique(Go_CvsO3$Geneid)
pwf <- nullp(isSigGene, "hg19", "ensGene")
goResults <- goseq(pwf, "hg19","ensGene")
select_go <- c("GO:0036462", "GO:1902253", "GO:0097190", "GO:0043067", "GO:0006977", "GO:2001233", "GO:0008219", 
               "GO:0012501", "GO:0006950", "GO:0043408", "GO:1902175", "GO:0043410")
select_go1 <- subset(goResults, goResults$category %in% select_go) %>% mutate(hitsPerc=numDEInCat*100/numInCat)
select_go1 <- select_go1 %>% rename(number_of_gene=numDEInCat, p_value=over_represented_pvalue)
ggplot(select_go1, aes(hitsPerc, category, colour=p_value,size=number_of_gene)) + 
  geom_point() + xlab("Percent of genes hits") + ylab("") + scale_size_continuous(range = c(3, 12))+
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), text = element_text(size = 18))

enriched.GO=goResults$category[p.adjust(goResults$over_represented_pvalue)<.05]

library(goseq)
Go_CvsO3 <- CvsO3 %>% filter(biotype=="mRNA")
Go_CvsO3 <- Go_CvsO3 %>% distinct(Geneid, .keep_all = T)

supportedOrganisms() %>% filter(str_detect(Genome, "hg19"))
isSigGene <- unique(Go_CvsO3$Geneid) %in% unique(common_trans_gene )
names(isSigGene) <- unique(Go_CvsO3$Geneid)
pwf <- nullp(isSigGene, "hg19", "ensGene")
goResults <- goseq(pwf, "hg19","ensGene", test.cats=c("GO:BP"))
select_go <- c("GO:0036462", "GO:1902253", "GO:0097190", "GO:0043067", "GO:0006977", "GO:2001233", "GO:0008219", 
               "GO:0012501", "GO:0006950", "GO:0043408", "GO:1902175", "GO:0043410")
select_go1 <- subset(goResults, goResults$category %in% select_go) %>% mutate(hitsPerc=numDEInCat*100/numInCat)
select_go1 <- select_go1 %>% rename(number_of_gene=numDEInCat, p_value=over_represented_pvalue)
ggplot(select_go1, aes(hitsPerc, category, colour=p_value,size=number_of_gene)) + 
  geom_point() + xlab("Percent of genes hits") + ylab("") + scale_size_continuous(range = c(3, 12))+
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), text = element_text(size = 18))

#--------------------------------------------------------------------------------------------------


#Analysis of gene regulated by lncRNA in trans
CvsO3_lncRNA <- CvsO3 %>% filter(biotype=="lncRNA") %>% arrange(padj)
CvsO3_lncRNA_20 <- CvsO3_lncRNA[1:25,]
write.csv(CvsO3_lncRNA_20, "./CvsO3_lncRNA_25.csv", quote = F)
CvsOH_lncRNA <- CvsOH %>% filter(biotype=="lncRNA") %>% arrange(padj)
CvsOH_lncRNA_20 <- CvsOH_lncRNA[1:25,]
write.csv(CvsOH_lncRNA_20, "./CvsOH_lncRNA_25.csv", quote = F)

lncRNA_mRNA <- read.csv("./lncRNA_RNA interaction_final.csv", header = TRUE)
lncRNA_mRNA <- lncRNA_mRNA[complete.cases(lncRNA_mRNA), ] %>% rename(GeneSymble=Name)
CvsO3_interaction <- merge(CvsO3, lncRNA_mRNA, by="GeneSymble", all=F) %>% rename(mRNA=Geneid)
CvsO3_interaction <- CvsO3_interaction %>% rename(Geneid=Name_lncRNA)
CvsO3_interaction <- merge(CvsO3_interaction, id_symbol, by="Geneid", all=F)

CvsO3_interaction1 <- CvsO3_interaction %>% filter(abs(log2FoldChange) > 1 & padj < 0.01)
write.csv(CvsO3_interaction1, "./CvsO3_network.csv", quote = F)
CvsO3_interaction1 <- CvsO3_interaction1[,c(1,7)]

CvsOH_interaction <- merge(CvsOH, lncRNA_mRNA, by="GeneSymble", all=F) %>% rename(mRNA=Geneid)
CvsOH_interaction <- CvsOH_interaction %>% rename(Geneid=Name_lncRNA)
CvsOH_interaction <- merge(CvsOH_interaction, id_symbol, by="Geneid", all=F)
CvsOH_interaction1 <- CvsOH_interaction[,c(2,5,9,10,14,18,20,21)]
CvsOH_interaction1 <- CvsOH_interaction1 %>% filter(abs(log2FoldChange) > 1 & padj < 0.01)
write.csv(CvsOH_interaction1, "./CvsOH_network.csv", quote = F)
CvsOH_interaction1 <- CvsOH_interaction1[,c(1,7)]

library(GGally)
library(ggnet)
library(network)
library(sna)
library(ggplot2)

levs <- unique(unlist(CvsOH_interaction1, use.names = FALSE))
net <- table(lapply(CvsOH_interaction1, factor, levs))
label <- colnames(net)
net = network(net, directed = FALSE)
attribute <- as.data.frame(label)
colnames(attribute) <- c("GeneSymble")
OH <- inner_join(attribute, CvsOH)
OH$status <- factor(ifelse(OH$log2FoldChange >0, "Up", "Down")) 
network.vertex.names(net) = label
net %v% "Status" <- as.character(OH$status)
color = c(rep(c("mRNA"), 203),rep(c("lncRNA"), 20))
net %v% "Type" <- color
ggnet2(net, node.size = 4, label = TRUE, label.size=2.5, 
       color = "Type", shape = "Status", palette = "Set1")


levs <- unique(unlist(CvsOH_interaction1, use.names = FALSE))
net <- table(lapply(CvsOH_interaction1, factor, levs))
label <- colnames(net)
net = network(net, directed = FALSE)
attribute <- as.data.frame(label)
colnames(attribute) <- c("GeneSymble")
OH <- inner_join(attribute, CvsOH)
OH$status <- factor(ifelse(OH$log2FoldChange >0, "Up", "Down")) 
network.vertex.names(net) = label
net %v% "Status" <- as.character(OH$status)
color = c(rep(c("mRNA"), 201),rep(c("lncRNA"), 20))
net %v% "Type" <- color
ggnet2(net, node.size = 4, label = TRUE, label.size=2.5, 
       color = "Type", shape = "Status", palette = "Set1")

CvsO3_lncRNA_20 <- CvsO3_lncRNA_20 %>% rename(GeneSymble.y=GeneSymble)
CvsO3_interaction <-  CvsO3_interaction%>% filter(abs(log2FoldChange) > 1 & padj < 0.01)
SI_Table4 <- merge(CvsO3_interaction, CvsO3_lncRNA_20, by="GeneSymble.y", all=F)
write.csv(SI_Table4, "SI_Table4.csv")
CvsOH_lncRNA_20 <- CvsOH_lncRNA_20 %>% rename(GeneSymble.y=GeneSymble)
CvsOH_interaction <-  CvsOH_interaction%>% filter(abs(log2FoldChange) > 1 & padj < 0.01)
SI_Table4 <- merge(CvsO3_interaction, CvsO3_lncRNA_20, by="GeneSymble.y", all=T)
write.csv(SI_Table4, "SI_Table5.csv")







