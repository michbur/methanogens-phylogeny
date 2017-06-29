

seq_type <- "aa"


#For nucleotide sequences
sek <- process_aln_files(readAAStringSet(fasta_file, format="fasta"))



full_aln <- msa(sek, chosen_aln_method)

conv_aln <- msaConvert(full_aln, type=c("seqinr::alignment"))


write.fasta(sequences=conv_aln[["seq"]], names=seq_aln$nam[i], file.out="seq.aln", open= "a")
conv_aln()


chosen_file <- "glob.fasta.aln"

# zmienna wewnętrzna, wartości - dna dla DNA i aa dla aminokwasów
seq_type <- "dna"


