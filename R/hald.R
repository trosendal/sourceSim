##' hald functions
##'
##' @param ... other arguments
##' @param x An object to run the hald method on
##' @export
##' @return A character vector of the proportions
hald <- function(x, ...) UseMethod("hald")

##' run hald model on a \code{sourceSim_result} object
##' @param x the \code{sourceSim_result} object
##' @param iter the number of iterations to run the simulation for
##' @param burnin the number of burnin iterations
##' @param thinning the thinning rate
##' @param n_chains the number of simulation chains
##' @param q_MaxRange the maximum value of the Hald model 'q' parameter
##'        source-specific coefficient)
##' @param a_MaxRange the maximum value of the Hald model 'a' parameter
##'        (MLST-specific coefficient)
##' @param others_cutoff For MLST's that have less than or equal to
##'        \code{others_cutoff} individuals, reclassify them as type "others"
##' @param simplify When TRUE the return is a matrix of
##'        summary stats. When FALSE the return is a list with the same
##'        matrix and the mcmc results.
##' @param ... other arguments
##' @return proportions of the attribution for each population
##' @export
hald.sourceSim_result <- function(x,
                                  iter = 10000,
                                  burnin = iter / 2,
                                  thinning = 10,
                                  n_chains = 3,
                                  q_max = 2000,
                                  a_max = 100,
                                  others_cutoff = 4,
                                  simplify = TRUE,
                                  ...) {

    pops <- x$population

    if (!("Pop_human" %in% names(pops))) {
        stop("The data must contain 'Pop_human'. ",
             "You need to run the 'sample_humans' ",
             "function to create this.")
    }

    pops <- pops[, c("seqID", names(pops)[grepl("^Pop", names(pops))])]

    pop_names <- names(pops)[-1]
    pop_names <- pop_names[-length(pop_names)]

    names(pops) <- gsub("^Pop\\_", "", names(pops))

    pops[pops$human <= others_cutoff, ]$seqID <- "others"

    if ("others" %in% pops$seqID) {
        pops_others <- pops[pops$seqID == "others", ]
        pops <- rbind(
            pops[pops$seqID != "others", ],
            aggregate(pops_others, . ~ seqID, sum)
        )
    }

    human_pops <- pops$human
    names(human_pops) <- pops$seqID

    ## Don't use 'others' as most prevalent genotype. This indexing works
    ## because 'others' is always last.
    max_human <- names(which.max(human_pops[human_pops != "others"]))
    human_pops <- c(human_pops[max_human],
                    human_pops[names(human_pops) != max_human])

    source_pops <- as.data.frame(
        t(pops[, !names(pops) %in% c("seqID", "human")])
    )
    colnames(source_pops) <- pops$seqID

    source_pops <- cbind(data.frame(max_human = source_pops[, max_human]),
                         source_pops[, colnames(source_pops) != max_human])
    colnames(source_pops)[colnames(source_pops) == "max_human"] <- max_human

    source_pops <- as.matrix(source_pops)

    source_pops[source_pops == 0] <- 1e-6
    human_pops[human_pops == 0] <- 1e-6

    n_sources <- nrow(source_pops)
    n_serotypes <- ncol(source_pops)

    sim_data <- list(
        pop_names = pop_names,
        humans = human_pops,
        sources = source_pops,
        FoodSourceCount = n_sources,
        SeroVarCount = n_serotypes,
        q_MaxRange = q_max,
        a_MaxRange = a_max
    )

    hald(sim_data,
         iter,
         burnin,
         thinning,
         n_chains,
         simplify,
         ...)
}

##' run hald model on a \code{list} object using OpenBUGS
##'
##' @param x the \code{list} object
##' @param iter the number of iterations to run the simulation for
##' @param burnin the number of burnin iterations
##' @param thinning the thinning rate
##' @param n_chains the number of simulation chains
##' @param simplify When TRUE the return is a matrix of
##'        sumamry stats. When FALSE the return is a list with the same
##'        matrix and the simulation results.
##' @param ... other arguments
##' @return proportions of the attribution for each population
hald.list <- function(
        x,
        iter = 10000,
        burnin = iter / 2,
        thinning = 10,
        n_chains = 3,
        simplify = TRUE,
        ...) {

    stopifnot(setequal(
        names(x),
        c(
            "pop_names",
            "humans",
            "sources",
            "FoodSourceCount",
            "SeroVarCount",
            "q_MaxRange",
            "a_MaxRange"
        )
    ))

    pop_names <- x$pop_names
    x <- x[names(x) != "pop_names"]

    hald_dir <- tempdir()
    wd <- setwd(hald_dir)
    on.exit(setwd(wd))

    model_output <- runjags::run.jags(system.file("hald/Bugs_model_code.txt",
                                                  package = "sourceSim"),
                                      monitor = c("lambdaji", "lambdaexp", "a", "q"),
                                      data = x,
                                      n.chains = n_chains,
                                      inits = NA,
                                      sample = iter,
                                      burnin = burnin,
                                      thin = thinning
                                      )

    ## 'lambdaji' is the number of samples for every combination of serovar i
    ## and food source j. here we sum along i to get the number of samples
    ## per food source (for every iteration separately)
    sources <- apply(
        model_output$sims.list$lambdaji, 1, function(x) apply(x, 1, sum)
    )

    ## the sum of every iteration is the total number of samples
    sums <- apply(sources, 2, sum)

    ## divide sample counts by total number to obtain attribution fractions
    source_fractions <- sapply(seq_along(sums),
                               function(i) sources[, i] / sums[i])

    ## calculate summary statistics for each source population
    pe <- apply(source_fractions, 1, function(x) {
        c(
            "mean" = mean(x),
            "median" = stats::median(x),
            "sd" = stats::sd(x),
            stats::quantile(x, c(.025, .975))
        )
    })

    colnames(pe) <- pop_names
    class(pe) <- c("hald_result_table", class(pe))

    if (isFALSE(simplify)) {
        x$sources[x$sources < 1] <- 0
        x$humans[x$humans < 1] <- 0

        pe <- list(sim_data = x, model = model_output, pe = pe, dir = hald_dir)

        class(pe) <- c("hald_output", class(pe))
    }

    pe
}
