#Analysis of gene regulated by lncRNA in trans
CvsO3_lncRNA <- CvsO3 %>% filter(biotype=="lncRNA") %>% arrange(padj_CvsO3)
CvsO3_lncRNA_20 <- CvsO3_lncRNA[1:25,]
#write.csv(CvsO3_lncRNA_20, "./CvsO3_lncRNA_25.csv", quote = F)
CvsOH_lncRNA <- CvsOH %>% filter(biotype=="lncRNA") %>% arrange(padj_CvsOH)
CvsOH_lncRNA_20 <- CvsOH_lncRNA[1:25,]
write.csv(CvsOH_lncRNA_20, "./CvsOH_lncRNA_25.csv", quote = F)

lncRNA_mRNA <- read.csv("./input_files/lncRNA_RNA interaction_final.csv", header = TRUE)
lncRNA_mRNA <- lncRNA_mRNA[complete.cases(lncRNA_mRNA), ]
lncRNA_mRNA <- lncRNA_mRNA %>% dplyr::select(Name, Name_lncRNA)
colnames(lncRNA_mRNA) <- c("GeneSymbol", "lncRNA_id")
lncRNA_mRNA$lncRNA_id <- gsub("\\.[0-9]", replacement = "", lncRNA_mRNA$lncRNA_id)
id_conversion$Geneid <- gsub("\\.[0-9]", replacement = "", id_conversion$Geneid)
id_conversion <- id_conversion %>% filter(biotype=="lncRNA")
colnames(id_conversion) <- c("lncRNA_id", "lncRNA_symbol", "biotype")
CvsO3_interaction <- merge(CvsO3, lncRNA_mRNA, by="GeneSymbol", all=F) 

CvsO3_interaction <- merge(CvsO3_interaction, id_conversion, by="lncRNA_id", all=F)

CvsO3_interaction1 <- CvsO3_interaction[,c(2,4,5,7)]
CvsO3_interaction1 <- CvsO3_interaction1 %>% filter(abs(log2FC_CvsO3) > 2 & padj_CvsO3 < 0.01)
write.csv(CvsO3_interaction1, "./CvsO3_network.csv", quote = F)
CvsO3_interaction1 <- CvsO3_interaction1[,c(1,4)]

CvsOH_interaction <- merge(CvsOH, lncRNA_mRNA, by="GeneSymbol", all=F) 
CvsOH_interaction <- merge(CvsOH_interaction, id_conversion, by="lncRNA_id", all=F)
CvsOH_interaction1 <- CvsOH_interaction[,c(2,4,5,7)]
CvsOH_interaction1 <- CvsOH_interaction1 %>% filter(abs(log2FC_CvsOH) > 2 & padj_CvsOH < 0.01)
write.csv(CvsOH_interaction1, "./CvsOH_network.csv", quote = F)
CvsOH_interaction1 <- CvsOH_interaction1[,c(1,4)]

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
colnames(attribute) <- c("GeneSymbol")
OH <- inner_join(attribute, CvsOH)
OH$status <- factor(ifelse(OH$log2FC_CvsOH > 0, "Up", "Down")) 
network.vertex.names(net) = label
net %v% "Status" <- as.character(OH$status)
color = c(rep(c("mRNA"), 79),rep(c("lncRNA"), 19))
net %v% "Type" <- color

png(filename="figure_4A_log2FC2.png", width=600,height = 600, bg="white")
par(mar=c(5,6,4,1)+.1)
ggnet2(net, node.size = 4, label = TRUE, label.size=2.5, 
       color = "Type", shape = "Status", palette = "Set1")
dev.off()

levs <- unique(unlist(CvsO3_interaction1, use.names = FALSE))
net <- table(lapply(CvsO3_interaction1, factor, levs))
label <- colnames(net)
net = network(net, directed = FALSE)
attribute <- as.data.frame(label)
colnames(attribute) <- c("GeneSymbol")
O3 <- inner_join(attribute, CvsO3)
O3$status <- factor(ifelse(O3$log2FC_CvsO3 > 0, "Up", "Down")) 
network.vertex.names(net) = label
net %v% "Status" <- as.character(O3$status)
color = c(rep(c("mRNA"), 73),rep(c("lncRNA"), 18))
net %v% "Type" <- color
png(filename="figure_4B.png", width=600,height = 600, bg="white")
par(mar=c(5,6,4,1)+.1)
ggnet2(net, node.size = 4, label = TRUE, label.size=2.5, 
       color = "Type", shape = "Status", palette = "Set1")
dev.off()

