library(sourceSim)
result <- simu(list(MIGR = 0, MIGP = 0))

## Make sure there is no overlap in the populations
stopifnot(!any(rowSums(result$population[, 1:4] > 0) > 1))


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
    MIGR = migr,
    MIGP = migp,
    MUTR = mutation_rate,
    RECR = rr,
    RECL = rl,
    PROA = proa,
    PROT = prot,
    PROG = prog,
    PROC = proc,
    SEQS = seqs
)
parameters$MIGI <- 1
mig_mat <- matrix(c(0,    0, 0,
                    0.7,  0, 0,
                    0.02, 0, 0),
                  nrow = 3,
                  byrow = TRUE)
result2 <- simu(input = list(MIGI = 1, NPOP = 3), migration = mig_mat)
## confirm that the simulation generates overlapping metapopulations
stopifnot(any(rowSums(result2$population[, 1:3] > 0) > 1))
