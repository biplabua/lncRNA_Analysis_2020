from Bio import SeqIO
top_lncRNA = []
with open("/Users/biplab/Desktop/lncRNA_Analysis/CvsO3_lncRNA_25.csv", "r") as fh:
	for line in fh:
		geneName20 = line.split(",")[1]
		#print(geneName20.strip())
		top_lncRNA.append(geneName20)
top_fasta = []
for seq_record in SeqIO.parse("/Users/biplab/Downloads/gencode.v34.lncRNA_transcripts.fa", "fasta"):
    geneName = seq_record.id.split("|")[1]
    #print(geneName)
    if geneName in top_lncRNA:
    	top_fasta.append(seq_record)

SeqIO.write(top_fasta, "CvsO3_short_seqs.fasta", "fasta")
top_lncRNA = []
with open("/Users/biplab/Desktop/lncRNA_Analysis/CvsOH_lncRNA_25.csv", "r") as fh:
	for line in fh:
		geneName20 = line.split(",")[1]
		#print(geneName20.strip())
		top_lncRNA.append(geneName20)
top_fasta = []
for seq_record in SeqIO.parse("/Users/biplab/Downloads/gencode.v34.lncRNA_transcripts.fa", "fasta"):
    geneName = seq_record.id.split("|")[1]
    #print(geneName)
    if geneName in top_lncRNA:
    	top_fasta.append(seq_record)

SeqIO.write(top_fasta, "CvsOH_short_seqs.fasta", "fasta")
   	
    
