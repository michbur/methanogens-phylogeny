
library(msa)
library(ape)
library(phangorn)
library(phytools)
library(ggtree)

#For nucleotide sequences
sek <- process_aln_files(readAAStringSet(fasta_file, format="fasta"))

process_aln_files <- function(x) {
  names(slot(x, "ranges")) <- lapply(strsplit(names(slot(x, "ranges")), "|", fixed = TRUE), 
                                       function(i) paste0(i[1L:4], collapse = "|"))
  names(slot(x, "ranges")) <- gsub("|[|]+", "|", names(slot(x, "ranges")), fixed = TRUE)
  x
}

full_aln <- msa(sek, chosen_aln_method)

conv_aln <- msaConvert(full_aln, type=c("seqinr::alignment"))

dir_name <- tempdir()

write.fasta(sequences=conv_aln[["seq"]], names=seq_aln$nam[i], file.out="seq.aln", open= "a")
conv_aln()

chosen_file <- "homo.fasta.aln"
chosen_file <- "glob.fasta.aln"

# zmienna wewnętrzna, wartości - dna dla DNA i aa dla aminokwasów
seq_type <- "dna"
seq_type <- "aa"

source("./functions/phyl.R")

seq_aln <- try(seq_fun[["read"]](chosen_file), silent = TRUE)
res_tree_cmp <- seq_fun[["phyl"]](seq_aln)
plot_tree(res_tree_cmp)

svg("tmp_name.svg", width = 7, height = 0.45*sum(!is.na(fort_tree[["label"]])), pointsize = 12)
plot_tree(res_tree_cmp)
dev.off()