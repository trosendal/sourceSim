## Experiment 6: Calculate the expected error given a random model output,
## i.e. the "expected" are sampled the same way as the observed attribution
## frequencies. This value is the maximum error that a model needs to
## outperform.

library(sourceSim)

expected <-  rdirichlet(1000, c(1, 1, 1))

errors <- apply(expected, 1, function(x) {
    observed <- rdirichlet(1000, c(1,1,1))
    sapply(seq_len(ncol(observed)), function(i) {
        rmse(observed[, i], x[i])
    })
})

apply(errors, 1, mean)
mean(errors)
