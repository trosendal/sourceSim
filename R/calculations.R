
## Nucleotide frequencies
sequences <-
    sapply(list.files("data/MLST_alleles", full.names = T),
           function(x) {
               seq <- readLines(x)
               paste(seq[!grepl("^>", seq)], collapse = "")
           })

combined <- unlist(sequences)

bases <- paste(sequences, collapse = "")

n_bases <- nchar(bases)

frequencies <- sapply(c("A", "T", "C", "G"),
                      function(x)
                          lengths(regmatches(bases,
                                             gregexpr(x, bases)))) / n_bases

params <- paste0("PRO", names(frequencies), ":\t", frequencies)

## Mutation rate for C. jejuni
## Estimates from Wilson et al. (2009, Mol Biol Evol)
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2639114/
generation_length <- 2.79e-4    # years
mutation_rate <- 3.23e-2 *  # per kilobase per year
    generation_length * 1000    # per base per generation

params <- c(params, paste0("MUTR:\t", mutation_rate))

## Recombination rate relative to mutations
## Estimate from Yu et al. (2012, J Mol Evol)
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3985069/
rr <- 6.96

params <- c(params, paste0("RECR:\t", rr))

writeLines(params, con = "data/params.txt")
