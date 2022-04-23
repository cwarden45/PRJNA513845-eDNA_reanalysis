input.file = "FLASH_combined_unique_seqs-min_2_reads--FILTER_OTU_MIN_100_READS-Swarm_OTU_with_counts.megablast"
output.plot = "FLASH_combined_unique_seqs-min_2_reads--FILTER_OTU_MIN_100_READS-Swarm_OTU_with_counts.megablast-bar_plots.png"

input.table = read.delim(input.file, head=F, sep="\t")
otu_counts = as.character(input.table$V1)
otu_counts = gsub("OTU\\d+_","",otu_counts)
print(otu_counts[1:5])
otu_counts = as.numeric(otu_counts)
print(otu_counts[1:5])

count_per_kingdom = tapply(otu_counts, input.table$V9, sum)

png(output.plot, height=400, width=800)
par(mfcol=c(1,2))
barplot(table(input.table$V9), main = "Unique, Filtered Sequences with BLAST Hits",
			ylab="Counts")
barplot(count_per_kingdom, main = "Summed Counts per Filtered OTU with BLAST Hit",
			ylab="Counts")
dev.off()