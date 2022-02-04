# Simulation of data for source attribution analysis

The idea of this package is to use the methods described in this
paper: https://doi.org/10.1093/bioinformatics/bty093 and its
accompanying software:
https://bitbucket.org/aleksisipola/bacmeta/src/master/ to simulate a
bacterial population across several artificial partially overlapping
sources. This will be used to test the efficacy of the existing
methods of source attribution for at least the microbial subtyping
methods.

## Install the sourceSim software

To install the software you can install it in R directly from github like
this:

```R
library(remotes)
remotes::install_github(repo = "trosendal/sourceSim", build_vignettes = TRUE)
```

## Using the sourceSim to generate artificial data

An example of running sourceSim is provided in the package as a
vignette. If you run the following you can work through an example:

```R
library('sourceSim')
vignette('run_simulation')
```
