awk 'FS=";"{print $1}' gencode.v19.long_noncoding_RNAs.gtf | awk '{if($3=="gene") print $1"\t"$4"\t"$5"\t"$10"\t"$6"\t"$7}' | sed 's/"//g' | sed 's/chr//g' > lncRNA_genes.bed
awk 'FS=";"{print $1}' Homo_sapiens.GRCh37.75.gtf | awk '{if($2=="protein_coding" && $3=="gene") print $1"\t"$4"\t"$5"\t"$10"\t"$6"\t"$7}' | sed 's/"//g' | awk '{if($1!="MT")print $0}' | sort -k1,1 -k2,2n > codingGene.bed
bedtools closest -a codingGene.bed -b lncRNA_genes.bed -d | awk '{if($13 > 20000) print $0}' > genelncRNA.bed
