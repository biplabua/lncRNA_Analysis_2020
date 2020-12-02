get_coexpressionMatrix <- function(NormData.log, GeneList){
  #3.11 Calculating correlation between genes pairs
  
  #Since the goal of this guide is to generate a gene co-expression network, it is necessary to determinate the genes that co-express in our determined conditions. 
  #First, the user should select the regulated genes from the normalized table generated previously by the DESeq2 and extract their normalized counts, using the following command:
  
  Norm.interest <- subset(NormData.log, rownames(NormData.log) %in% GeneList)
  
  Norm.interest.corr <- corr.test( t(Norm.interest), method="pearson", ci=F)
  
  #Among the many results of this function, there are two triangular matrices with the needed data. 
  #One matrix contains the correlation values and the other contain the p-values in the upper part and the adjusted p-values in the lower part. 
  #To generate a table comprising the data organized properly to build the network, the user should type the following commands: 
  
  
  Norm.interest.corr$p[lower.tri( Norm.interest.corr$p,diag=TRUE)]=NA
  
  Pval.adj <- as.data.frame(as.table(Norm.interest.corr$p))
  
  Norm.interest.corr$r [lower.tri( Norm.interest.corr$r,diag=TRUE)]=NA
  
  Correlation <- as.data.frame(as.table(Norm.interest.corr$r))
  
  Cor.table <- na.exclude(cbind( Correlation, Pval.adj))[,c(1,2,3,6)]
  
  colnames(Cor.table) <- c("GeneSymble","gene2","cor","p.adj")
  
  #The generated table, can be filer according absolute correlation (0.9) and adjusted p-value (0.01) thresholds:
  
  Cor.table.filt <- Cor.table [(abs(Cor.table[,3])>0.9 & Cor.table[,4] <0.01 ),]
  
  #To facilitate graphs generation in future network, adjusted p-value can be converted to logarithmic base 10 and to avoid problems with log of 0,
  #we propose to sum the minimal adjusted p-value found to the entire column.  
  
  Log.p.adj <- log10(Cor.table.filt[,4]+min(Cor.table.filt [Cor.table.filt[,4]!=0,4]))
  
  Cor.table.filt<-cbind(Cor.table.filt,Log.p.adj)
  return(Cor.table.filt)
}
