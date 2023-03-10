library(sourceSim)

## Test with migration
sequences <- internal_sequences("Campylobacter")
bases <- paste(sequences, collapse = "")
frequencies <- base_frequencies(bases)
proa <- frequencies[names(frequencies) == "A"]
prot <- frequencies[names(frequencies) == "T"]
prog <- frequencies[names(frequencies) == "G"]
proc <- frequencies[names(frequencies) == "C"]
nloc <- length(sequences)
opfn <- 1
lole <- 500
nbac <- 50000
npop <- 3
migr <- 0
migp <- 0
generation_length <- 2.79e-4    # years
gen_per_year <- 1 / generation_length
mutation_rate <- 3.23e-2 # per kilobase per year
mutation_rate <- mutation_rate / 1000 # per base per year
mutation_rate <- mutation_rate / gen_per_year # per base per generation
## multiply mutation rate by a large number to speed up the
## simulation:
mutation_rate <- mutation_rate * 10
rr <- 6.96
rl <- 0
iseq <- 1
seqs <- 0.01
## set the generations to 1000 to speed up
seqi <- 100
genr <- 1000
parameters <- list(
    OPFN = opfn,
    LOLE = lole,
    GENR = genr,
    NLOC = nloc,
    NBAC = nbac,
    NPOP = npop,
    MIGR = 0,
    MIGP = 0.5,
    MUTR = mutation_rate,
    RECR = rr,
    RECL = rl,
    PROA = proa,
    PROT = prot,
    PROG = prog,
    PROC = proc,
    SEQS = seqs,
    SEED = 123
)
mig_mat <- matrix(c(0,    0, 0,
                    0.01,  0, 0,
                    0,    0, 0),
                  nrow = 3,
                  byrow = TRUE)
result <- simu(input = parameters, migration = mig_mat)

ex <- paste("The data must contain 'Pop_human'.",
                   "You need to run the 'sample_humans'",
                   "function to create this.")
ob <- tools::assertError(isource(result))[[1]]$message
stopifnot(identical(ob, ex))

frequency <- c(0.2, 0.5, 0.3)
result <- sample_humans(x = result,
                        attribution = frequency,
                        n = 1000)

sa <- isource(result, simplify = FALSE)

## sourceSim:::plot.isource_output(sa, type = "area")
