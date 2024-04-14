## Experiment 6: Calculate the expected error given a random model output,
## i.e. the "expected" are sampled the same way as the observed attribution
## frequencies. This value is the maximum error that a model needs to 
## outperform.

library(sourceSim)

errors <- sapply(1:10000, function(i) {
    frequencies <- rdirichlet(2, c(1, 1, 1))
    apply(frequencies, 2, function(x) abs(diff(x)))
})

mean(errors)
