#Vendiagram
library(VennDiagram)
lncRNA_O3_V <- lncRNA_O3 %>% filter(abs(log2FC_CvsO3) > 1 & padj_CvsO3 < 0.05)
lncRNA_O3_V_up <- lncRNA_O3 %>% filter(log2FC_CvsO3 > 1 & padj_CvsO3 < 0.05)
lncRNA_O3_V_down <- lncRNA_O3 %>% filter(log2FC_CvsO3 < -1 & padj_CvsO3 < 0.05)

lncRNA_OH_V <- lncRNA_OH %>% filter(abs(log2FC_CvsOH) > 1 & padj_CvsOH < 0.05)
lncRNA_OH_V_up <- lncRNA_OH %>% filter(log2FC_CvsOH > 1 & padj_CvsOH < 0.05)
lncRNA_OH_V_down <- lncRNA_OH %>% filter(log2FC_CvsOH < -1 & padj_CvsOH < 0.05)


write.csv(lncRNA_O3_V_up, file = "./lncRNA_O3_up.csv")
write.csv(lncRNA_O3_V_down, file = "./lncRNA_O3_down.csv")

write.csv(lncRNA_OH_V_up, file = "./lncRNA_OH_up.csv")
write.csv(lncRNA_OH_V_down, file = "./lncRNA_OH_down.csv")


colors <- c("darkgreen", "coral3")

venn.diagram(x=list(lncRNA_O3_V_up$Geneid, lncRNA_OH_V_up$Geneid), category.names = c("O3", "OH"), 
             filename = "figure_1C.png", col="black", fill =  colors, cat.col = colors,
             cex = 3, cat.cex = 3, sub.col = c("black", "black"))

venn.diagram(x=list(lncRNA_O3_V_down$Geneid, lncRNA_OH_V_down$Geneid), category.names = c("O3", "OH"), 
             filename = "figure_1D.png", col="black", fill =  colors, cat.col = colors,
             cex = 3, cat.cex = 3)
