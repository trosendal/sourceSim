library(sourceSim)
## No migration at all

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

frequency <- rdirichlet(1, c(1, 1, 1))

result <- simu(input = parameters)

result_h <- sample_humans(x = result,
                          attribution = frequency,
                          n = 1000)

res <- hald(result_h, iter = 10000,
            burnin = 5000,
            thinning = 10)

## we expect to get values of a that reflect the 'frequency' object
## above which should be the attribution for each source.... But we don't:

a <- res$sims.list$a

par(mfrow = c(1, 3))
apply(a, 2, function(x) {
    plot(density(x))
})

q <- res$sims.list$q

par(mfrow = c(8, 5),
    mar = c(0,0,0,0))
apply(q, 2, function(x) {
    plot(density(x))
})
