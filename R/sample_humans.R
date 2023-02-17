##' functions to sample humans
##'
##' @param ... other arguments
##' @param x An object to sample humans on
##' @return A character vector of the proportions
##' @export
sample_humans <- function(x, ...) UseMethod("sample_humans")

##' sample_humans
##'
##' Sample the human cases from the bacterial sources given an
##' attribution fraction.
##'
##' @param x A result of a \code{simu()} function, a \code{sourceSim_result}
##'          object
##' @param attribution the attribution fraction for each of the source
##'        populations
##' @param n The number of human cases
##' @param overwrite If humans have already been sampled for this simulation,
##'        should a new sampling be run? Will overwrite the old one. Default
##'        is \code{FALSE}, in which case nothing happens.
##' @export
##' @return A sourceSim_result object with added humans
sample_humans.sourceSim_result <- function(x,
                          attribution = 1,
                          n = 100,
                          overwrite = FALSE) {
    df <- x$population

    if ("Pop_human" %in% names(df) && isFALSE(overwrite)) {
        cat("Humans have already been sampled and 'overwrite' is FALSE. ",
            "Nothing will be done.")
        return(x)
    }

    if (length(attribution) == 1)
        attribution <- rep(attribution, x$parameters$NPOP)
    stopifnot(length(attribution) == x$parameters$NPOP)

    pops <- seq_len(x$parameters$NPOP)

    ## Given the attribution, how many cases should come from each population?
    n_pop <- table(factor(sample(x = pops,
                                 size = n,
                                 prob = attribution,
                                 replace = TRUE),
                          levels = pops))

    ## Sample an MLST type for each human cases:
    seqs <- unlist(lapply(seq_len(x$parameters$NPOP), function(x) {
        sample(df[, "seqID"],
               n_pop[x],
               replace = TRUE,
               prob = df[, x])
    }))
    df$Pop_human <- as.numeric(table(factor(seqs, levels = df[, "seqID"])))
    x$population <- df
    x
}
