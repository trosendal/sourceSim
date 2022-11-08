library(sourceSim)
## check that the rosum is always very close to 1
stopifnot(all(round(rowSums(rdirichlet(100, c(10, 1, 1))), 5) == 1))
stopifnot(all(round(rowSums(rdirichlet(100, c(1, 1, 1, 1, 1))), 5) == 1))
