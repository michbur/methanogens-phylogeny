library(msa)

source("aln_params.R")
source("./functions/aln.R")

seq_type <- "aa"

raw_seqs <- process_aln_files(read_fun("mrca_prot.fasta"))

aln <- msa(raw_seqs, chosen_aln_method)

msaPrettyPrint(aln, output="pdf", 
               showNames="left", 
               showLogo="top", 
               shadingModeArg="chemical", 
               consensusThreshold=50, 
               logoColors="chemical", 
               shadingMode="functional")

