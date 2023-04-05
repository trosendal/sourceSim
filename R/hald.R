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

    names(pops) <- gsub("^Pop\\_", "", names(pops))

    pops_long <- stats::reshape(
        pops,
        direction = "long",
        varying = list(names(pops)[-1]),
        idvar = "seqID",
        v.names = "n",
        timevar = "population",
        times = names(pops)[-1]
    )

    pops_long[pops_long$n <= others_cutoff, ]$seqID <- "others"

    pops_long_agg <- stats::aggregate(
        pops_long$n,
        by = list(seqID = pops_long$seqID,
                  population = pops_long$population),
        FUN = sum
    )

    pops_wide <- stats::reshape(
        pops_long_agg,
        idvar = "seqID",
        timevar = "population",
        direction = "wide"
    )
    pops_wide[is.na(pops_wide)] <- 0
    names(pops_wide) <- names(pops)

    human_pops <- pops_wide$human
    names(human_pops) <- pops_wide$seqID
    human_pops <- c(human_pops[names(human_pops) != "others"],
                    human_pops["others"])
    max_human <- names(which(human_pops == max(human_pops)))
    human_pops <- c(human_pops[max_human],
                    human_pops[names(human_pops) != max_human])

    source_pops <-
        as.data.frame(t(pops_wide[, !names(pops_wide) %in% c("seqID",
                                                            "human")]))
    colnames(source_pops) <- pops_wide$seqID
    source_pops <-
        cbind(source_pops[, colnames(source_pops) != "others"],
              data.frame(others = source_pops[, "others"]))

    source_pops <- cbind(data.frame(max_human = source_pops[, max_human]),
                         source_pops[, colnames(source_pops) != max_human])
    colnames(source_pops)[colnames(source_pops) == "max_human"] <- max_human

    source_pops <- as.matrix(source_pops)

    n_sources <- nrow(source_pops)
    n_serotypes <- ncol(source_pops)

    sim_data <- list(
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
            "humans",
            "sources",
            "FoodSourceCount",
            "SeroVarCount",
            "q_MaxRange",
            "a_MaxRange"
        )
    ))

    hald_dir <- tempdir()
    wd <- setwd(hald_dir)
    on.exit(setwd(wd))

    model_output <- R2OpenBUGS::bugs(
        x,
        inits = NULL,
        parameters.to.save = "lambdaji",
        model.file = system.file("hald/Bugs_model_code.txt",
                                 package = "sourceSim"),
        n.chains = n_chains,
        n.iter = iter,
        n.burnin = burnin,
        n.thin = thinning,
        working.directory = hald_dir
    )

    ## 'lambdaji' is the number of samples for every combination of serovar i
    ## and food source j. here we sum along i to get the number of samples
    ## per food source (for every iteration separately)
    sources <- apply(
        model_output$sims.list$lambdaji, 1, function(x) apply(x, 1, sum))

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

    class(pe) <- c("hald_result_table", class(pe))

    if (isFALSE(simplify)) {
        pe <- list(pe = pe, model = model_output)

        class(pe) <- c("hald_output", class(pe))
    }

    pe
}
