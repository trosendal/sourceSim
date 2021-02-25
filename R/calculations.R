sequences <-
    sapply(list.files("data/MLST_alleles", full.names = T), readLines)

combined <- unlist(sequences)

bases <- paste(combined[!grepl("^>", combined)],
               collapse = "")

n_bases <- nchar(bases)

frequencies <- sapply(c("A", "T", "C", "G"),
                      function(x)
                          lengths(regmatches(bases,
                                             gregexpr(x, bases)))) / n_bases

frequencies <- paste0(names(frequencies), ": ", frequencies)

writeLines(frequencies, con = "data/nucleotide_frequencies.txt")
