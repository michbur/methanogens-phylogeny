#Select a nucleotide substitution model
# model:"raw", "N", "TS", "TV", "JC69", "K80", "F81", "K81", "F84", "BH87", "T92", "TN93", "GG95", "logdet", "paralin", "indel", "indelblock".
chosen_model_dna <- "F84"

#Select a amino acid substitution model
# model: "JC69", "F81", "WAG", "JTT", "LG", "Dayhoff", "cpREV", "mtmam", "mtArt", "MtZoa", "mtREV24", "VT","RtREV", "HIVw", "HIVb", "FLU", "Blossum62", "Dayhoff_DCMut", "JTT_DCMut"
chosen_model_aa <- "F81"


# Gamma: okienko do wpisywania wartosci domyslnie nic
chosen_gamma <- FALSE

#SN1987
#Select a method to reconstruct phylogenetic tree
#the neighbor-joining method of Saitou and Nei (1987)

#G1997
#the BIONJ algorithm of Gascuel (1997)

#DG2002
#the the minimum evolution algorithm of Desper and Gascuel (2002)

# tree_type: lista trzech mozliwości (SN1987, G1997, DG2002)
chosen_tree_type = "SN1987"


# Okienko do wyboru opcji
# Dla pairwise.deletion checkbox, nazwa: Delete sites with at least one missing data for all sequences
chosen_deletion <- FALSE


# Opcja do wyboru: Bootstrap z okienkiem do wartosci parametru B (możliwe wartości: 100, 200, 300)
chosen_B <- 100