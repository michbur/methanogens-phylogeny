library(msa)
library(seqinr)
library(ape)
library(phangorn)
library(phytools)
library(ggtree)

source("aln_params.R")
source("phyl_params.R")

# poziomy: aa lub dna
seq_type <- "aa"

source("./functions/phyl.R")
source("./functions/aln.R")

tmp_dir <- tempdir()

raw_seqs <- process_aln_files(read_fun("mrca_prot.fasta"))

aln <- msa(raw_seqs, chosen_aln_method)

# plik do sciagniecia aln.pdf
msaPrettyPrint(aln, 
               output="pdf",
               file = paste0(tmp_dir, "/aln.pdf"),
               showNames="left", 
               showLogo="top", 
               shadingModeArg="chemical", 
               consensusThreshold=50, 
               logoColors="chemical", 
               shadingMode="functional")

conv_aln <- msaConvert(aln, type=c("seqinr::alignment"))

# drugi plik do sciagniecia aln.fasta
write.fasta(sequences = as.list(conv_aln[["seq"]]), names=conv_aln[["seq"]], 
            file.out = paste0(tmp_dir, "/aln.fasta"))

# --------------------------------------------------------

# tutaj filogeneza
aln_file <- paste0(tmp_dir, "/aln.fasta")

seq_aln <- try(seq_fun[["read"]](aln_file), silent = TRUE)
res_tree_cmp <- seq_fun[["phyl"]](seq_aln)

svg(paste0(tmp_dir, "/tmp_name.svg"), width = 7, height = 0.45*sum(!is.na(fort_tree[["label"]])), pointsize = 12)
plot_tree(res_tree_cmp)
dev.off()