library(fgsea)

#Gene set enrichment analysis of cancer related genes
negative_lncRNA <- read.table("./input_files/negative.lncRNA.glist.xls.txt")
positive_lncRNA <- read.table("./input_files/positive.lncRNA.glist.xls.txt")
all_lncRNA <- rbind(negative_lncRNA, positive_lncRNA)
colnames(all_lncRNA) <- c("Geneid")
genesets <- list(Cancer_Related=all_lncRNA$Geneid)

gseaDat <- CvsO3 %>% filter(!is.na(Geneid)) %>% filter(log2FC_CvsO3 > 0)
ranks <- gseaDat$log2FC_CvsO3
names(ranks) <- gseaDat$Geneid
ranks <- sort(ranks, decreasing = T)
fgseaRes <- fgsea(genesets, ranks, minSize=15, maxSize = 500, nperm=1000)
png(filename="figure_2A.png", width=400,height = 200, bg="white")
par(mar=c(5,6,4,1)+.1)
plotEnrichment(genesets[["Cancer_Related"]], ranks, ticksSize = 0.3)
dev.off()


gseaDat <- CvsO3 %>% filter(!is.na(Geneid)) %>% filter(log2FC_CvsO3 < 0)
ranks <- gseaDat$log2FC_CvsO3
names(ranks) <- gseaDat$Geneid
ranks <- sort(ranks, decreasing = T)
fgseaRes <- fgsea(genesets, ranks, minSize=15, maxSize = 500, nperm=1000)

png(filename="figure_2B.png", width=400,height = 200, bg="white")
par(mar=c(5,6,4,1)+.1)
plotEnrichment(genesets[["Cancer_Related"]], ranks)
dev.off()

gseaDat <- CvsOH %>% filter(!is.na(Geneid)) %>% filter(log2FC_CvsOH > 0)
ranks <- gseaDat$log2FC_CvsOH
names(ranks) <- gseaDat$Geneid
ranks <- sort(ranks, decreasing = T)
fgseaRes <- fgsea(genesets, ranks, minSize=15, maxSize = 500, nperm=1000)

png(filename="figure_2C.png", width=400,height = 200, bg="white")
par(mar=c(5,6,4,1)+.1)
plotEnrichment(genesets[["Cancer_Related"]], ranks)
dev.off()

gseaDat <- CvsOH %>% filter(!is.na(Geneid)) %>% filter(log2FC_CvsOH < 0)
ranks <- gseaDat$log2FC_CvsOH
names(ranks) <- gseaDat$Geneid
ranks <- sort(ranks, decreasing = T)
fgseaRes <- fgsea(genesets, ranks, minSize=15, maxSize = 500, nperm=1000)

png(filename="figure_2D.png", width=400,height = 200, bg="white")
par(mar=c(5,6,4,1)+.1)
plotEnrichment(genesets[["Cancer_Related"]], ranks)
dev.off()


