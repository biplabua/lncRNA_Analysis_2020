source("./Rcode/get_coexpression_matrix.R")
#Now we have generated a table containing the statistically significant correlations for every pair of gene differentially regulated.
#The following commands can perform the mentioned normalization. 
count <- rbind(df1,df)
count <- as.data.frame(count) %>% rownames_to_column("Geneid")
count <- merge(id_symbol, count, by="Geneid")
count <- count %>% distinct(GeneSymble, .keep_all = TRUE)
gene_length <- count$length
#count <- count[,c(1,2,3,4,5,7,10)]
count <- count[,c(2,5:ncol(count))]
count <- column_to_rownames(count, 'GeneSymble')
library(edgeR)
y <- DGEList(counts=count,genes=data.frame(Length=gene_length))
y <- calcNormFactors(y)
RPKM <- rpkm(y)
NormData.log <- log2(RPKM+1)
#3.12 Network Stats
CvsO3_top_500_down1 <- get_coexpressionMatrix(NormData.log, CvsO3_top_500_up_down$GeneSymble)
CvsO3_top_500_down <- merge(CvsO3_top_500_down1, id_symbol, by="GeneSymble", all=T )
CvsO3_top_500_down <- filter(CvsO3_top_500_down, biotype=="lncRNA")
CvsO3_top_500_down <- rename(CvsO3_top_500_down,c('GeneSymble'='lncRNA'))
write.table(CvsO3_top_500_down, "./CoExpressionMatrix/CvsO3_Cor.table.filter_down.txt", sep="\t", row.names=F, quote=F)

CvsO3_top_500_up <- get_coexpressionMatrix(NormData.log, CvsO3_top_500_up$GeneSymble)
CvsO3_top_500_up <- merge(CvsO3_top_500_up, id_symbol, by="GeneSymble")
write.table(CvsO3_top_500_up, "./CoExpressionMatrix/CvsO3_Cor.table.filter_up.txt", sep="\t", row.names=F, quote=F)

CvsOH_top_500_down <- get_coexpressionMatrix(NormData.log, CvsOH_top_500_down$GeneSymble)
CvsOH_top_500_down <- merge(CvsOH_top_500_down, id_symbol, by="GeneSymble")
write.table(CvsOH_top_500_down, "./CoExpressionMatrix/CvsOH_Cor.table.filter_down.txt", sep="\t", row.names=F, quote=F)

CvsOH_top_500_up <- get_coexpressionMatrix(NormData.log, CvsOH_top_500_up$GeneSymble)
CvsOH_top_500_up <- merge(CvsOH_top_500_up, id_symbol, by="GeneSymble")
write.table(CvsOH_top_500_up, "./CoExpressionMatrix/CvsOH_Cor.table.filter_up.txt", sep="\t", row.names=F, quote=F)
