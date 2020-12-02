awk '{if($3=="exon") print $10 "\t" $16}' ../input_files/Homo_sapiens.GRCh37.75.gtf | sed 's/";//g' | sed 's/"//g'| uniq > ../input_files/mRNA_2_symble.tsv
lncRNA_Analysis % awk '{if($3=="exon") print $10 "\t" $16}' ../input_files/gencode.v19.long_noncoding_RNAs.gtf | sed 's/";//g' | sed 's/"//g'| uniq > ./input_files/lncRNA_ensemble2symble.tsv

