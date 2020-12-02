#Analysis of genes regulated by lncRNA in cis
cis_mRNA <- read.table("./input_files/genelncRNA.bed")[,c(4,10)]
colnames(cis_mRNA) <- c("mRNA_id", "lncRNA_id")
cis_mRNA$lncRNA_id <- gsub("\\.[0-9]", replacement = "", cis_mRNA$lncRNA_id)

CvsO3_lncRNA <- CvsO3 %>% filter(biotype=="lncRNA") %>% dplyr::select(Geneid, log2FC_CvsO3, padj_CvsO3, 
                                                          GeneSymbol) #filter(Geneid %in% cis_mRNA$lncRNA_id)
colnames(CvsO3_lncRNA) <- c("lncRNA_id", "lncRNA_log2FC_O3", "lncRNA_padj_O3", "lncRNA_GeneSymbol")

CvsO3_lncRNA <- merge(CvsO3_lncRNA, cis_mRNA, by="lncRNA_id", all = FALSE) 


CvsO3_mRNA <- CvsO3 %>% filter(biotype=="mRNA") %>% dplyr::select(Geneid, log2FC_CvsO3, padj_CvsO3, 
                                                                  GeneSymbol) #%>% filter(Geneid %in% cis_mRNA$Gene) 
colnames(CvsO3_mRNA) <- c("mRNA_id", "mRNA_log2FC_O3", "mRNA_padj_O3", "mRNA_GeneSymbol")


CvsO3_mRNA_lncRNA <- merge(CvsO3_mRNA, CvsO3_lncRNA, by="mRNA_id", all = FALSE) 


CvsOH_lncRNA <- CvsOH %>% filter(biotype=="lncRNA") %>% dplyr::select(Geneid, log2FC_CvsOH, padj_CvsOH, 
                                                                      GeneSymbol) #filter(Geneid %in% cis_mRNA$lncRNA_id)
colnames(CvsOH_lncRNA) <- c("lncRNA_id", "lncRNA_log2FC_OH", "lncRNA_padj_OH", "lncRNA_GeneSymbol")

CvsOH_lncRNA <- merge(CvsOH_lncRNA, cis_mRNA, by="lncRNA_id", all = FALSE) 


CvsOH_mRNA <- CvsOH %>% filter(biotype=="mRNA") %>% dplyr::select(Geneid, log2FC_CvsOH, padj_CvsOH, 
                                                                  GeneSymbol) #%>% filter(Geneid %in% cis_mRNA$Gene) 
colnames(CvsOH_mRNA) <- c("mRNA_id", "mRNA_log2FC_OH", "mRNA_padj_OH", "mRNA_GeneSymbol")

CvsOH_mRNA_lncRNA <- merge(CvsOH_mRNA, CvsOH_lncRNA, by="mRNA_id", all = FALSE) 

df_fig3 <- merge(CvsO3_mRNA_lncRNA, CvsOH_mRNA_lncRNA, by="mRNA_id", all = FALSE)
df_fig3 <- df_fig3[,c(1,6,2,13,9)]
colnames(df_fig3) <- c("Geneid", "lncRNA O3", "Nearest Gene O3", "lncRNA OH", "Nearest Gene OH")
df_fig3 <- df_fig3 %>% distinct(Geneid, .keep_all = T)

mat = as.matrix(df_fig3[,c(2,3,4,5)])
png(filename="figure_3A.png", width=200,height = 400, bg="white")
par(mar=c(5,6,4,1)+.1)
Heatmap(mat, name = "log2FC", cluster_columns = FALSE)
dev.off()




CvsO3_mRNA_lncRNA <- CvsO3_mRNA_lncRNA %>% filter(mRNA_padj_O3 < 0.01 & lncRNA_padj_O3 < 0.01) %>% 
  filter(abs(mRNA_log2FC_O3) > 1.0 & abs(lncRNA_log2FC_O3) > 1.0) %>% distinct()
write.csv(CvsO3_mRNA_lncRNA, file = "CvsO3_mRNA_lncRNA_lg2FC_cis.csv")

CvsOH_mRNA_lncRNA <- CvsOH_mRNA_lncRNA %>% filter(mRNA_padj_OH < 0.01 & lncRNA_padj_OH < 0.01) %>% 
  filter(abs(mRNA_log2FC_OH) > 1.0 & abs(lncRNA_log2FC_OH) > 1.0) %>% distinct()
write.csv(CvsOH_mRNA_lncRNA, file = "CvsOH_mRNA_lncRNA_lg2FC_cis.csv")



colors <- c("darkgreen", "coral3")
venn.diagram(x=list(CvsO3_mRNA_lncRNA$mRNA_id, CvsOH_mRNA_lncRNA$mRNA_id), category.names = c("O3", "OH"), 
             filename = "figure_3B.png", col="black", fill =  colors, cat.col = colors, cex = 3, cat.cex = 3)

