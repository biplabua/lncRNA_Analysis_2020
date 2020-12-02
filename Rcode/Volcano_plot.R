

lncRNA_O3 <- CvsO3 %>% filter(biotype=="lncRNA")
lncRNA_OH <- CvsOH %>% filter(biotype=="lncRNA")

##Volcano plot with "significant" genes labeled
png(filename="figure_1A.png", width=400, bg="white")
par(mar=c(5,6,4,1)+.1)
with(lncRNA_O3, plot(log2FC_CvsO3, -log10(padj_CvsO3), pch=20, xlim=c(-5,4), cex = 1.0, cex.axis=2.0, cex.lab=2.0))
with(subset(lncRNA_O3, padj_CvsO3 < 0.05 & log2FC_CvsO3 > 1), points(log2FC_CvsO3, -log10(padj_CvsO3), pch=20, col="red", cex = 1.0))
with(subset(lncRNA_O3, padj_CvsO3 < 0.05 & log2FC_CvsO3 < -1), points(log2FC_CvsO3, -log10(padj_CvsO3), pch=20, col="green", cex = 1.0))
dev.off()

png(filename="figure_1B.png", width=400, bg="white")
par(mar=c(5,6,4,1)+.1)
with(lncRNA_OH, plot(log2FC_CvsOH, -log10(padj_CvsOH), pch=20,  xlim=c(-5,4), cex = 1.0, cex.axis=2.0, cex.lab=2.0))
with(subset(lncRNA_OH, padj_CvsOH < 0.05 & log2FC_CvsOH > 1), points(log2FC_CvsOH, -log10(padj_CvsOH), pch=20, col="red", cex = 1.0))
with(subset(lncRNA_OH, padj_CvsOH < 0.05 & log2FC_CvsOH < -1), points(log2FC_CvsOH, -log10(padj_CvsOH), pch=20, col="green", cex = 1.0))
dev.off()
