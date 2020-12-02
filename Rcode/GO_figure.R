library(tidyverse)
setwd("~/Desktop/Des/lncRNA_Analysis_07_16_2020/")
df <- read.csv("input_files/GO_from_CPDB.csv")
df <- df[,c(2,6,8,9 )]
colnames(df) <- c("q_value", "GO_term", "Number_of_DE_genes", "Percent_of_Gene")

png(filename="figure_3C.png", width=600,height = 600, bg="white")
par(mar=c(5,6,4,1)+.1)
ggplot(df, aes(Percent_of_Gene, GO_term, colour=q_value,size=Number_of_DE_genes)) + 
  geom_point() + xlab("Percent of genes hits") + ylab("") + scale_size_continuous(range = c(3, 12))+
  theme(axis.text = element_text(size = 18), axis.title = element_text(size = 20), text = element_text(size = 18))
dev.off()
