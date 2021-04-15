## Nucleotide frequencies
sequences <-
    sapply(list.files("../data/MLST_alleles", full.names = T),
           function(x) {
               seq <- readLines(x)
               paste(seq[!grepl("^>", seq)], collapse = "")
           })

bases <- paste(sequences, collapse = "")

n_bases <- nchar(bases)

frequencies <- sapply(c("A", "T", "C", "G"),
                      function(x)
                          lengths(regmatches(bases,
                                             gregexpr(x, bases)))) / n_bases
## Number of bacteria in each population
nbac <- 50000


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

## Insertion and deletion set to 0 to avaoid alignment at the end fo
## the simulation.

parameters <-
    c("OPFN:\t 1", 			#Run id: Output file name modifier, alphanumeric. Former files with same modifier will be overwritten!
      "SEED:\t 123456789",		#Set seed or give zero to have random seed from time.
      "GENR:\t 20000",			#Number of generations to run
      "MEAN:\t 0",				#Compute and save mean distance summary every 1000 generation (1) or not (0)
      "LOLE:\t 500",			#Length of locus
      "NLOC:\t 7",				#Number of loci in bacteria
      paste0("NBAC:\t", nbac),			#Number of bacteria in population
      "NPOP:\t 3",				#Number of populations in metapopulation
      "PRIG:\t 0",		#Debugging output generation, prints sequences and metadata to terminal.
      "PRIP:\t 0",				#Debugging output population, prints sequences and metadata to terminal.
      "MIGI:\t 0",				#Using migration rates from "migration.input"-file (1) or universal value given with MIGR (0)
      "MIGR:\t 0",			#Mean migration rate scaler
      "MIGP:\t 0",			#Migration probability
      "MICA:\t 0.0", 			#Mean microepidemic amount
      "MICS:\t 0.0",			#Mean microepidemic size scaler
      paste0("MUTR:\t ", mutation_rate),			#Mutation rate per nucleotide per population per generation
      paste0("RECR:\t ", rr),			#Recombination rate in relation to mutations
      paste0("RECL:\t ", rl),			#Recombination mean length. Give (0) if whole loci are to be recombined
      "RECA:\t 18",			#Recombination acceptance parameter for similarity test
      "RECS:\t 0",				#Gather recombination site metadata (1) or not (0)
      "INSR:\t 0",			#Insertion rate in relation to mutations
      "INSL:\t 1.7",			#Insertion length parameter for Zipf distribution
      "DELR:\t 0",			#Deletion rate in relation to mutations
      "DELL:\t 1.7",			#Deletion length parameter for Zipf distribution
      "INDM:\t 0.02",			#Maximum indel length as a proportion of the loci length
      paste0("PROA:\t ", frequencies[names(frequencies) == "A"]),			#Proportion of base A
      paste0("PROT:\t ", frequencies[names(frequencies) == "T"]),			#Proportion of base T
      paste0("PROG:\t ", frequencies[names(frequencies) == "G"]),			#Proportion of base G
      paste0("PROC:\t ", frequencies[names(frequencies) == "C"]),			#Proportion of base C
      "ATOT:\t 0.3334",		#P(a->t) i.e. Probability of mutation from base A to T
      "ATOG:\t 0.3333",		#P(a->g)
      "ATOC:\t 0.3333",		#P(a->c)
      "TTOA:\t 0.3334",		#P(t->a)
      "TTOG:\t 0.3333",		#P(t->g)
      "TTOC:\t 0.3333",		#P(t->c)
      "GTOA:\t 0.3333",		#P(g->a)
      "GTOT:\t 0.3333",		#P(g->t)
      "GTOC:\t 0.3334",		#P(g->c)
      "CTOA:\t 0.3333",		#P(c->a)
      "CTOC:\t 0.3333",		#P(c->t)
      "CTOG:\t 0.3334",		#P(c->g)
      "RSEL:\t 1",				#Randomly select bacteria for new generation (1) or keep as is (0)
      "RORD:\t 1",				#Keep order of event types same (0) or randomize each generation (1)
      "SUMI:\t 20000",			#Summary info: Generation interval of computing summaries.
      "GWDS:\t 0.01",			#Loci-wise mismatches and mutation-only loci-wise mismatches: Size as proportion of the population.
      "MGWD:\t 0",				#Do mutation-only loci-wise mismatches along regular (1) or not (0). Activates recombination site metadata.
      "GWDI:\t -20000",			#Loci-wise mismatches and mutation-only mismatches: Generation interval of computing and saving to file.
      "PGWS:\t 0.01",			#Interpopulation: Loci-wise mismatches and mutation-only loci-wise mismatches: Size as proportion of the population.
      "PMGW:\t 0",				#Interpopulation: Do mutation-only loci-wise mismatches along regular (1) or not (0). Activates recombination site metadata.
      "PGWI:\t -20000",			#Interpopulation: Loci-wise mismatches and mutation-only mismatches: Generation interval of computing and saving to file.
      "ISEQ:\t 1",				#Save initial genome (1) or not (0)
      "SEQS:\t 0.1",			#Sequence saving: Sample size as proportion of population
      "SEQI:\t 10000",			#Sequence saving: Generation interval
      "STRA:\t -20000",		#Strain-id composition documenting interval
      "MUTD:\t -20000",			#Do mutation documenting with this interval, give 0 to use recombination documenting interval RECI
      "RECI:\t -20000",			#Do recombination count documenting with this interval
      "RECT:\t 0",			#Do recombination event documenting with this hd threshold or not (0). Outputfile can be big! Use small values cautiously.
      "INDD:\t 0",					#Do indel documenting (1) or not (0)
      "MIGD:\t 0",					#Do migration documenting (1) or not (0)
      "MICD:\t 0",					#Do microepidemic documenting (1) or not (0)
      "MLST:\t 0",					#Do MLST documenting (1) or not (0).
      "STAL:\t 0")					#Save sequence-allele type pairs (1) or not (0). Outputfile can be big!


writeLines(parameters, con = "../data/params.txt")
writeLines(parameters, con = "simu1.input")

## Now run the simulation with the these rates without migration
## between populations. The results will be in the output directory
## with the suffix "1":
system("./simu -p 1")

## Now generate a set of sequences that have known migration between
## reservoirs. We have 3 reservoirs (populations) so we will set the
## migration as a matrix rows are from and columns are to.
##
## First set the model parameter to use the matrix:
parameters[grep("MIGI", parameters)] <-  "MIGI:\t 1"
parameters[grep("OPFN", parameters)] <-  "OPFN:\t 2"
##
## There appears to be a bug in the bactmeta software that requires an
## extra column that is ignored. If I run this with just the 3 columns
## the software stops and incidated the matrix is too small:
mat <- matrix(c(0,    0, 0, "*",
                0.7,  0, 0, "*",
                0.02, 0, 0, "*"), nrow = 3, byrow = TRUE)
## So this should be 70% migration from source 2 --> 1 and 5%
## migration from 3 --> 1 and other migration rates are 0.
write.table(mat, file = "migration.input", row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
writeLines(c(readLines("migration.input"), "", "##Commments here", "", "***"), "migration.input")
writeLines(parameters, con = "simu2.input")

## Now run this simulation. The results will also land in the output directory but will have the suffix "2"
system("./simu -p 2")
