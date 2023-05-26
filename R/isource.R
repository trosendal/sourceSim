##' Path to isource
##' @noRd
path_to_isource <- function(subpath = "") {
    file.path(.sim_env$isource, subpath)
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
##' in the \code{.sim_env} object and is accessible via
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
        assign("isource", tempdir(), envir = .sim_env)
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
##' @param x The result of a simulation of data, a \code{sourceSim_result}
##'          object
##' @param iter The number of iterations to run
##' @param burnin The burin length
##' @param thinning The thinning rate
##' @param dirichlet_param The parameter on the dirichlet
##' @param group_var The variable to group the results by
##' @param simplify When TRUE the return is a matrix of
##'     proportions. When FALSE the return is a list with the same
##'     matrix and the mcmc results.
##' @param ... other arguments
##' @return The result of running the island model, an \code{isource_output}
##' object
isource.sourceSim_result <- function(x = NULL,
                                     iter = 20000,
                                     burnin = 1000,
                                     thinning = 50,
                                     dirichlet_param = 1,
                                     group_var = "group",
                                     simplify = TRUE,
                                     ...) {

    if (!("Pop_human" %in% names(x$population))) {
        stop("The data must contain 'Pop_human'. ",
             "You need to run the 'sample_humans' ",
             "function to create this.")
    }

    pops <- seq_len(x$parameters$NPOP) - 1
    pops <- c("Pop_human", paste0("Pop_", pops))

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
    colnames(df) <- c(group_var,
                      "ST",
                      "ASP",
                      "GLN",
                      "GLT",
                      "GLY",
                      "PGM",
                      "TKT",
                      "UNC")

    df[, group_var] <- as.character(
        factor(df[, group_var],
               levels = pops,
               labels = c(0, seq_len(x$parameters$NPOP))))

    df <- as.data.frame(df)
    df <- df[, c(2, 3, 4, 5, 6, 7, 8, 9, 1)]

    isource(
        df,
        pops = pops,
        iter = iter,
        burnin = burnin,
        thinning = thinning,
        dirichlet_param = dirichlet_param,
        group_var = group_var,
        simplify = simplify,
        ...)
}

##' isource.data.frame
##'
##' Runs isource asymmetric island model on a dataframe
##'
##' @export
##' @param x The result of a simulation of data, a \code{data.frame}
##' @param iter The number of iterations to run
##' @param pops the names of the populations
##' @param burnin The burin length
##' @param thinning The thinning rate
##' @param dirichlet_param The parameter on the dirichlet
##' @param group_var The variable to group the results by
##' @param simplify When TRUE the return is a matrix of
##'     proportions. When FALSE the return is a list with the same
##'     matrix and the mcmc results.
##' @param ... other arguments
##' @return proportions of the attribution for each population
isource.data.frame <- function(x = NULL,
                               pops,
                               iter = 20000,
                               burnin = 1000,
                               thinning = 50,
                               dirichlet_param = 1,
                               group_var = "group",
                               simplify = TRUE,
                               ...) {

    stopifnot(is.logical(simplify),
              identical(length(simplify), 1L))

    npop <- length(pops) - 1

    isource_dir <- tempdir()
    wd <- setwd(isource_dir)
    on.exit(setwd(wd))

    utils::write.table(
        x,
        file = "input.txt",
        quote = FALSE,
        row.names = FALSE,
        sep = "\t"
    )

    compile_isource()

    system2(
        path_to_isource("isource/isource"),
        args = list(
            "input.txt",
            "output.txt",
            iter,
            thinning,
            dirichlet_param,
            shQuote(group_var)
        )
    )

    mcmc <- utils::read.table("output.txt", header = TRUE, comment.char = "")
    fmcmc <- utils::read.table("f_output.txt", header = TRUE, comment.char = "")

    g <- t(matrix(scan(
        "g_output.txt",
        what = double(0),
        sep = "\t"),
        nrow = npop))

    sim <- list(
        mcmc = mcmc,
        fmcmc = fmcmc,
        g = g,
        ng = npop
    )

    ## Set the burnin
    fd <- sim$fmcmc$iter >= burnin
    x <- sim$fmcmc[fd, 2:(sim$ng + 1)]
    names(x) <- pops[-1]

    pe <- apply(x, 2, function(x) {
        c(
            "mean" = mean(x),
            "median" = stats::median(x),
            "sd" = stats::sd(x),
            stats::quantile(x, c(.025, .975))
        )
    })

    class(pe) <- c("isource_result_table", class(pe))

    if (isFALSE(simplify)) {
        pe <- list(pe = pe,
                   sim = sim,
                   x = x,
                   pops = pops)
        class(pe) <- c("isource_output", class(pe))
    }

    pe
}

##' plot.isource_output
##'
##' Plots described in the isource documentation
##'
##' @param x The result of the isource method
##' @param type One of c("mcmc", "hist", "bar", "evol", "pie-evol", "area")
##' @param COL A colour vector
##' @param cod The order of the groups
##' @param sc The names of the sources
##' @param ... other arguments
##' @export
plot.isource_output <- function(x,
                                type = c("mcmc",
                                         "hist",
                                         "bar",
                                         "evol",
                                         "pie-evol",
                                         "area"),
                                COL = NULL,
                                cod = NULL,
                                sc = NULL,
                                ...) {

    type <- match.arg(type)

    if (is.null(COL)) {
        COL <- rainbow(x$sim$ng)
    }
    stopifnot(length(COL) == x$sim$ng)

    if (is.null(cod)) {
        cod <- seq_len(x$sim$ng)
    }
    stopifnot(length(cod) == x$sim$ng)

    if (is.null(sc)) {
        sc <- letters[seq_len(x$sim$ng)]
    }
    stopifnot(length(sc) == x$sim$ng)

    ## SET THE BURN-IN
    gd <- x$sim$mcmc$iter >= 1000
    fd <- x$sim$fmcmc$iter >= 500

    ##
    ## PLOT 1
    ## VISUALISE DIRECTLY THE MCMC OUTPUT FOR PARAMETER F
    ## (THE PROPORTION OF ISOLATES ATTRIBUTABLE TO EACH SOURCE)
    ##
    if (type == "mcmc") {
        plot(x$sim$fmcmc$f0[fd],
             type="l",
             ylim=c(0,1),
             col=COL[1],
             ylab="Proportion")
        for(i in 2:x$sim$ng) {
            lines(x$sim$fmcmc[fd, (1 + i)], col = COL[i])
        }
    }

    ##
    ## PLOT 2 HISTOGRAMS OF THE MARGINAL DISTRIBUTIONS OF F[i] (THE
    ## PROPORTION OF ISOLATES ATTRIBUTABLE TO SOURCE i)
    ##
    if (type == "hist") {
        par(mfrow=c(3,3))
        ## PLOT THE HISTOGRAMS
        for (i in 1:x$sim$ng) {
            hist(x$sim$fmcmc[fd, (1 + i)],
                 30,
                 col = COL[i],
                 main = sc[i],
                 prob = TRUE,
                 xlim = c(0, 1),
                 xlab = "Proportion")
        }
    }

    ## TABLE 1 SUMMARIES OF THE POSTERIOR DISTRIBUTIONS OF F[i]
    ##
    ## PLOT 3 BARCHART OF THE ESTIMATED PROPORTION OF CASES
    ## ATTRIBUTABLE TO EACH SOURCE
    ##
    if (type == "bar") {

        df <- x$sim$fmcmc[fd, 2:(x$sim$ng + 1)]
        names(df) <- sc
        pe <- apply(df, 2, function(x) {
            c("mean" = mean(x),
              "median" = median(x),
              "sd" = sd(x),
              quantile(x, c(0.025, 0.975)))
            })

        mp <- barplot(pe[1, ],
                      col = COL,
                      ylim = c(0, 1),
                      ylab = "Proportion of human cases")
        segments(mp,
                 pe[4, ],
                 mp,
                 pe[5, ],
                 lwd=2)
    }

    ##
    ## PLOT 4 MCMC TRACE OF THE EVOLUTIONARY PARAMETERS
    ##
    if (type == "evol") {

        par(mfrow = c(3, 3))
        COLR <- c(COL, "black")
        for (i in 1:x$sim$ng - 1) {
            plot(x$sim$mcmc$iter[gd],
                 x$sim$mcmc[[paste("A", i, 0, "", sep = ".")]][gd],
                 type = "l",
                 col = COLR[1],
                 ylim = c(0, 1),
                 xlab = "iter",
                 ylab = "M,R",
                 main = sc[i + 1])
        }
        for (j in 2:(x$sim$ng+1)) {
            lines(x$sim$mcmc$iter[gd],
                  x$sim$mcmc[[paste("A", i, j - 1, "", sep = ".")]][gd],
                  col=COLR[j])
        }
        lines(x$sim$iter[gd],
              x$sim$mcmc[[paste("r", i, sep="")]][gd],
              col="grey")
    }

    ##
    ## PLOT 5 PIE CHARTS OF THE EVOLUTIONARY PARAMETERS
    ##
    if (type == "pie-evol") {

        COLR <- c(COL, "black")
        for(i in 0:(x$sim$ng - 1)) {
            wh0 <- which(names(x$sim$mcmc) == paste("A", i, 0, "", sep = "."))
            whng <- which(names(x$sim$mcmc) == paste("A", i, x$sim$ng, "", sep = "."))
            pie(apply(x$sim$mcmc[gd, wh0:whng],
                      2,
                      mean),
                col = COLR,
                labels = "",
                main = sc[i + 1],
                col.main = COLR[i + 1],
                radius = 1.,
                border = "white")
        }
    }

    ##
    ## PLOT 6 POSTERIOR PROBABILITY OF SOURCE FOR EACH ISOLATE
    ##
    ## Plotting order for the groups. By varying this you can improve
    ## presentation of the final image.
    ##
    if (type == "area") {
        G <- x$sim$g[, cod]
        od <- order(G[,1 ] + G[, 2], G[, 3])
        res <- 1000
        tp <- apply(G,
                    1,
                    function(x) {
                        sort(sample(1:ncol(G),
                                    res,
                                    replace = TRUE,
                                    prob = x))
                    })

        ## DO THE PLOT
        image(1:nrow(G),
              seq(0, 1, len = res),
              t(tp[, od]),
              col = COL[cod],
              ylab = "Source probability",
              xlab = "Human cases",
              bty = "n")
    }
}
