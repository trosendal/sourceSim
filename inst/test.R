library(sourceSim)
## Parameters reasonable to campylobacter
proa <- c(A = 0.322069733492811)
prot <- c(T = 0.314661475155537)
prog <- c(G = 0.213327503265256)
proc <- c(C = 0.149941288086396)
nloc <- 7
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
rr <- 6.96
rl <- 0
iseq <- 1
seqs <- 0.01
seqi <- 10000
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
  PROA = proa,
  PROT = prot,
  PROG = prog,
  PROC = proc,
  SEQS = seqs
)

result <- lapply(1:30, function(i) {
    cat(i, "\n")
    ## Sample the migration rates
    ab <- runif(1)
    ac <- runif(1)
    bc <- runif(1)
    ## Sample the actual attrbution fractions
    attriba <- runif(1)
    attribb <- runif(1)
    attribc <- runif(1)
    frequency <- c(attriba, attribb, attribc)
    ## Assume the migration rate are balanced forth and back
    mig_mat <- matrix(c(0,  ab, ac,
                        ab, 0,  bc,
                        ac, bc, 0),
                      nrow = 3,
                      byrow = TRUE)
    result <- simu(input = parameters, migration = mig_mat)
    result <- sample_humans(x = result,
                            attribution = frequency,
                            n = 1000)
    res <- isource(result)
    list(attribution = res, migration = mig_mat, sampling = frequency)
})

if(!dir.exists("results"))
    dir.create("results")
filename <- tempfile(tmpdir = "results", fileext = ".Rds")
saveRDS(result, file = filename)
