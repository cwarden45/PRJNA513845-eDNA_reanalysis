import os
import sys
import re
from Bio import SeqIO
from Bio.Seq import Seq

inputFA = "FLASH_combined_unique_seqs-min_2_reads--FILTER_OTU_MIN_100_READS-Swarm_OTU_with_counts.fa"
input_table = "FLASH_combined_unique_seqs-min_2_reads--FILTER_OTU_MIN_100_READS-Swarm_OTU_with_counts.megablast"
output_table = "FLASH_combined_unique_seqs-min_2_reads--FILTER_OTU_MIN_100_READS-Swarm_OTU_with_counts.megablast-simple.txt"

#megablast hash

megablastHash = {}

inHandle = open(input_table)
line = inHandle.readline()
		
while line:
	line = re.sub("\n","",line)
	line = re.sub("\r","",line)
			

	lineInfo = line.split("\t")
	query_seq = lineInfo[13]
	query_name = lineInfo[0]
	OTU_name = re.sub("_\d+$","",query_name)
	#print(OTU_name)
	#sys.exit()
	hit_accession = lineInfo[1]
	hit_kingdom = lineInfo[8]
	if hit_kingdom == "N/A":
		hit_kingdom = "No Kingdom Info"
	hit_info = lineInfo[11]
	hit_length = lineInfo[4]

	hash_id = query_seq + "\t" + OTU_name
	
	megablastHash[hash_id]=hit_length + "\t" + hit_accession + "\t" + hit_info + "\t" + hit_kingdom
				
	line = inHandle.readline()
inHandle.close()

#compare to input FASTA and create new sample

missing_hit_count = 0

outHandle = open(output_table, 'w')

text = "OTU.seq\tOTU.name\tMEGABLAST.length\tMEGABLAST.sseqid\tMEGABLAST.sblastnames\tMEGABLAST.sskingdoms\n"
outHandle.write(text)
	
fasta_parser = SeqIO.parse(inputFA, "fasta")
	
for fasta in fasta_parser:
	seqID = fasta.id
	seq = str(fasta.seq)
	
	OTU_name = re.sub("_\d+$","",seqID)
			
	hash_id = seq+ "\t" + OTU_name

	BLAST_info = "NA\tNo MEGABLAST Hit\tNo MEGABLAST Hit\tNo MEGABLAST Hit"
	if hash_id in megablastHash:
		BLAST_info=megablastHash[hash_id]
	else:
		missing_hit_count +=1

	text = hash_id + "\t" + BLAST_info +"\n"
	outHandle.write(text)
		
outHandle.close()

print("Filtered Sequences without MEGABLAST hits: " +str(missing_hit_count))