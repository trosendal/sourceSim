library(sourceSim)

## Nucleotide frequencies
sequences <- get_sequences()

bases <- paste(sequences, collapse = "")

frequencies <- base_frequencies(bases)

# The number of loci is the number of sequences
nloc <- length(sequences)

# output filename modifier (must be alphanumeric)
opfn <- 1

# delete for performance
rm(sequences)

# length of each locus
lole <- 500

## Number of bacteria in each population
nbac <- 50000

# The number of populations to simulate
npop <- 3

# Migration rate scaler
migr <- 0

# Migration probability
migp <- 0

## Mutation rate for C. jejuni
## Estimates from Wilson et al. (2009, Mol Biol Evol)
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2639114/
generation_length <- 2.79e-4    # years
gen_per_year <- 1/generation_length
mutation_rate <- 3.23e-2 # per kilobase per year
mutation_rate <- mutation_rate/1000 # per base per year
mutation_rate <- mutation_rate/gen_per_year # per base per generation

## Recombination rate relative to mutations
## Estimate from Yu et al. (2012, J Mol Evol)
## https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3985069/
rr <- 6.96

## Recombintation length was determined to be 3000 bases therefore
## longer than the gene length so we will set the parameter to 0 to
## indicate the whole allele.
rl <- 0

# Save initial genome (1 == yes)
iseq <- 1

# sequence saving: sample size as proportion of population
seqs <- 0.01

# sequence saving: generation interval
seqi <- 10000

## Insertion and deletion set to 0 to avaoid alignment at the end fo
## the simulation.

parameters <- list(
  OPFN = opfn,
  LOLE = lole,
  NLOC = nloc,
  NBAC = nbac,
  NPOP = npop,
  MIGR = migr,
  MIGP = migp,
  MUTR = mutation_rate,
  RECR = rr,
  RECL = rl,
  PROA = frequencies[names(frequencies) == "A"],
  PROT = frequencies[names(frequencies) == "T"],
  PROG = frequencies[names(frequencies) == "G"],
  PROC = frequencies[names(frequencies) == "C"],
  ISEQ = iseq,
  SEQS = seqs,
  SEQI = seqi
)

dir.create("simu1")

simu_input <- create_simu.input(params = parameters,
                                out_path = "simu1",
                                suffix = 1)

res_1 <- simu(input = simu_input)

## Now generate a set of sequences that have known migration between
## reservoirs. We have 3 reservoirs (populations) so we will set the
## migration as a matrix rows are from and columns are to.
##
## First set the model parameter to use the matrix:

parameters$MIGI <- 1
parameters$OPFN <- 2
parameters$MIGR <- 0.01
parameters$MIGP <- 0.01

dir.create("simu2")

mat <- matrix(c(0,    0, 0,
                0.7,  0, 0,
                0.02, 0, 0),
              nrow = 3,
              byrow = TRUE)

migration_2 <- create_migration.input(n_populations = 3,
                       rates = mat,
                       out_path = "simu2",
                       suffix = 2)

simu_2 <- create_simu.input(params = parameters, out_path = "simu2", suffix = 2)

res_2 <- simu(input = simu_2, migration = migration_2, out_path = "simu2")
