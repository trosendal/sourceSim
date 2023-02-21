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
##' @param ... other arguments
##' @return a \code{hald_output} object, the result of an OpenBUGS run with
##'         the Hald model.
##' @export
hald.sourceSim_result <- function(x,
                                  iter = 20000,
                                  burnin = 1000,
                                  thinning = 50,
                                  n_chains = 3,
                                  q_MaxRange = 2000,
                                  a_MaxRange = 100,
                                  others_cutoff = 4,
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
        as.data.frame(t(pops_wide[,!names(pops_wide) %in% c("seqID",
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
        q_MaxRange = q_MaxRange,
        a_MaxRange = a_MaxRange
    )

    hald(sim_data,
         iter,
         burnin,
         thinning,
         n_chains,
         ...)
}

##' run hald model on a \code{list} object
##'
##' @param x the \code{list} object
##' @param iter the number of iterations to run the simulation for
##' @param burnin the number of burnin iterations
##' @param thinning the thinning rate
##' @param n_chains the number of simulation chains
##' @param ... other arguments
##' @return a \code{hald_output} object, the result of an OpenBUGS run with
##'         the Hald model.
##' @export
hald.list <- function(
        x,
        iter = 20000,
        burnin = 1000,
        thinning = 50,
        n_chains = 3,
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

    output <- R2OpenBUGS::bugs(
        x,
        inits = NULL,
        parameters.to.save = c("a", "q"),
        model.file = system.file("hald/Bugs_model_code.txt",
                                 package = "sourceSim"),
        n.chains = n_chains,
        n.iter = iter,
        n.burnin = burnin,
        n.thin = thinning,
        working.directory = hald_dir
    )

    class(output) <- c("hald_output", class(output))

    output
}
