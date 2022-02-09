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
seed <- 0
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
  SEQS = seqs,
  SEED = seed
)

rdirichlet <- function(n = 1, alpha = c(1, 1)) {
    g <- do.call("cbind", lapply(alpha, function(x) {
        rgamma(n, shape = x)
    }))
    g / matrix(rowSums(g), nrow = n, ncol = length(alpha))
}

result <- lapply(1:5, function(i) {
    cat("\n###################\n########",
        i,
        "########\n###################\n")

    ## Sample the actual attribution fractions
    frequency <- rdirichlet(1, c(1, 1, 1))

    ## Assume the migration rate from pop 0 to 1
    mig0 <- runif(1, 0, 0.1)
    mig_mat <- matrix(c(0,    0, 0,
                        mig0, 0, 0,
                        0,    0, 0),
                      nrow = 3,
                      byrow = TRUE)

    result <- simu(input = parameters, migration = mig_mat)

    result <- sample_humans(x = result,
                            attribution = frequency,
                            n = 1000)

    res <- isource(result)

    list(attribution = res,
         result = result,
         migration = mig_mat,
         sampling = frequency)
})

if(!dir.exists("results"))
    dir.create("results")
filename <- tempfile(tmpdir = "results", fileext = ".Rds")
saveRDS(result, file = filename)
