##' Path to isource
##' @noRd
path_to_isource <- function(subpath = "") {
    file.path(.bacmeta_env$isource, subpath)
}

##' Check if isource is Compiled
##'
##' If isource is compiled, each C++ script in its src directory should have
##' a corresponding .o file, and there should be a "isource" file in isource
##' root.
##'
##' @return \code{TRUE} if bacmeta is compiled, otherwise \code{FALSE}.
is_isource_compiled <- function() {
    all(length(path_to_isource()) > 0,
        file.exists(file.path(
            path_to_isource("isource"), "isource"
        )))
}

##' Compile isource
##'
##' Copies isource source files to an environment-specific temporary directory
##' and compiles isource there. The path to the compilation directory is saved
##' in the \code{.bacmeta_env} object and is accessible via
##' \code{path_to_isource()}. If isource is already compiled, does nothing.
##'
##' @param quiet If \code{TRUE}, prints no compilation progress to console.
##'        Default is \code{FALSE}.
##'
compile_isource <- function(quiet = FALSE) {
  if (isFALSE(quiet)) cat("Compiling isource...\n")
    if (is_isource_compiled()) {
      if (isFALSE(quiet)) cat("isource already compiled.\n")
    } else {
        assign("isource", tempdir(), envir = .bacmeta_env)
        file.copy(from = system.file("isource",
                                     package = "sourceSim"),
                  recursive = TRUE,
                  to = path_to_isource())
        wd <- setwd(path_to_isource("isource"))
        on.exit(setwd(wd))
        system("make", ignore.stdout = quiet)
        if (isFALSE(quiet)) cat("Compilation successful.\n")
    }
}

##' isource functions
##' @param ... other arguments
##' @param x An object to run the isource method on
##' @export
##' @return A character vector of the proportions
isource <- function(x, ...) UseMethod("isource")

##' isource.sourceSim_result
##'
##' Runs isource asymmetric island model:
##'
##' @export
##' @param x The result of a simulation of data
##' @param iter The number of iterations to run
##' @param burnin The burin length
##' @param thinning The thinning rate
##' @param dirichlet_param The parameter on the dirichlet
##' @param group_var The variable to group the results by
isource.sourceSim_result <- function(x = NULL,
                                     iter = 20000,
                                     burnin = 1000,
                                     thinning = 50,
                                     dirichlet_param = 1,
                                     group_var = "group") {
    if (!("Pop_human" %in% names(x$population))) {
        stop("The data must contain 'Pop_human'. ",
             "You need to run the 'sample_humans' ",
             "function to create this.")
    }

    pops <- seq_len(x$parameters$NPOP) - 1
    pops <- c("Pop_human", paste0("Pop_", pops))
    npop <- length(pops) - 1

    ## Expand to long form
    df <- do.call("rbind", lapply(pops, function(y) {
        df <- data.frame(pop = y,
                         count = x$population[, y],
                         MLST = x$population[, "MLST"],
                         ST = x$population[, "seqID"])
        do.call("rbind", apply(df, 1, function(z) {
            matrix(c(rep(z[1], z[2]),
                     rep(z[3], z[2]),
                     rep(z[4], z[2])), ncol = 3)
        }))
    }))

    df <- cbind(df[, c(1, 3)], do.call("rbind", strsplit(df[, 2], "-")))
    colnames(df) <- c("group",
                      "ST",
                      "ASP",
                      "GLN",
                      "GLT",
                      "GLY",
                      "PGM",
                      "TKT",
                      "UNC")

    df[, "group"] <- as.character(
        factor(df[, "group"],
               levels = pops,
               labels = c(0, seq_len(x$parameters$NPOP))))

    df <- as.data.frame(df)
    df <- df[, c(2, 3, 4, 5, 6, 7, 8, 9, 1)]

    isource_dir <- tempdir()
    wd <- setwd(isource_dir)
    on.exit(setwd(wd))

    write.table(df,
                file = "input.txt",
                quote = FALSE,
                row.names = FALSE,
                sep = "\t")

    compile_isource()

    system2(path_to_isource("isource/isource"), args = list("input.txt",
                                                            "output.txt",
                                                            iter,
                                                            thinning,
                                                            dirichlet_param,
                                                            shQuote(group_var)))
    mcmc <- read.table("output.txt", header = TRUE, comment.char = "")
    fmcmc <- read.table("f_output.txt", header = TRUE, comment.char = "")
    g <- t(matrix(scan("g_output.txt",
                       what = double(0),
                       sep = "\t"),
                  nrow = npop))
    sim <- list(mcmc = mcmc, fmcmc = fmcmc, g = g, ng = npop)

    ## Set the burnin
    gd <- sim$mcmc$iter >= burnin
    fd <- sim$fmcmc$iter >= burnin
    df <- sim$fmcmc[fd, 2:(sim$ng + 1)]
    names(df) <- pops[-1]
    pe <- apply(df, 2, function(x) {
        c("mean" = mean(x),
          "median" = stats::median(x),
          "sd" = stats::sd(x),
          stats::quantile(x, c(.025, .975)))
    })
    pe
}

##' isource.data.frame
##'
##' Runs isource asymmetric island model on a dataframe
##'
##' @export
##' @param x The result of a simulation of data
##' @param iter The number of iterations to run
##' @param burnin The burin length
##' @param thinning The thinning rate
##' @param dirichlet_param The parameter on the dirichlet
##' @param group_var The variable to group the results by
##' @return proportions of the attribution for each population
isource.data.frame <- function(x = NULL,
                               iter = 20000,
                               burnin = 1000,
                               thinning = 50,
                               dirichlet_param = 1,
                               group_var = "group") {
## Not implemented

}
##' sample_humans
##'
##' Sample the human cases from the bacterial sources given an
##' attribution fraction.
##'
##' @param x A result of a simu() function
##' @param attribution the attribution fraction for each of the source
##'        populations
##' @param n The number of human cases
##' @export
##' @return A sourceSim_result object with added humans
sample_humans <- function(x,
                          attribution = 1,
                          n = 100) {
    stopifnot("sourceSim_result" %in% class(x))
    df <- x$population

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
