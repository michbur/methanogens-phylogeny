read_fun <- switch(seq_type, 
                   aa = function(x) readAAStringSet(x),
                   dna = function(x) readDNAStringSet(x))

process_aln_files <- function(x) {
  names(slot(x, "ranges")) <- lapply(strsplit(names(slot(x, "ranges")), "|", fixed = TRUE), 
                                     function(i) paste0(i[1L:2], collapse = "|"))
  names(slot(x, "ranges")) <- gsub("|[|]+", "|", names(slot(x, "ranges")), fixed = TRUE)
  x
}
