##' Generate a Parameter File for Running Bacmeta Simulation
##'
##' Loads a predefined parameter template file with pre-filled Bacmeta
##' simulation parameters (see README for details), replaces the parameter
##' values as supplied in \code{params}, and writes the modified parameter file
##' to a file named simu\[\code{suffix}\].input in the directory \code{out_path}.
##'
##' @param params A named list of parameters to be written into the parameter
##'        file. The names of all values in \code{params} must match a
##'        in the template file. If \code{params} is \code{NULL} (default),
##'        an unmodified template file will be written.
##' @param out_path The destination to which the parameter file will be
##'        written. Must be an existing directory.
##' @param suffix A suffix that will be used for the param file name. If
##'        \code{suffix} is \code{NULL} (default), the filename will simply be
##'        "simu.input", otherwise the suffix will be pasted in between "simu"
##'        and ".input". E.g. if \code{suffix} is 123, the filename will be
##'        "simu123.input". All symbols in \code{suffix} must be alphanumeric
##'        (A-Z, a-z, 0-9).
##' @return invisible \code{NULL}
##' @author Wiktor Gustafsson
##' @export
create_simu.input <- function(params = NULL,
                              out_path = path_to_bacmeta(),
                              suffix = NULL) {

    if (!dir.exists(out_path)) {
        stop(ifelse(file.exists(out_path),
                    "'out_path' must be a directory (not a file).",
                    "'out_path' is invalid (no such directory)."))
    }

    if (!is.null(suffix) && !is_alphanumeric(suffix))
            stop(paste0("All symbols in 'sufffix' must be alphanumeric ",
                        "(A-Z, a-z or 0-9)."))

    out_path <- file.path(normalizePath(out_path),
                          paste0("simu.", suffix, ".input"))

    template <- read.table(system.file("bacmeta/simu.input.default",
                                       package = "sourceSim"),
                           sep = ":")

    if (!is.null(params)) {
        stopifnot(is.list(params))

        if (!all(names(params) %in% template[, 1]))
            stop(paste0(
                "Some parameters supplied in 'params' have invalid names. ",
                "See sourceSim README for list of legal parameters."))

        if (!all(sapply(params, is.numeric)))
            stop("All supplied parameters must be numeric")

        template[match(names(params), template[, 1]), 2] <- unlist(params)
    }

    write.table(
        template,
        file = out_path,
        sep = ": \t",
        row.names = F,
        col.names = F,
        quote = F
    )

    invisible(NULL)
}

##' Generate a Migration Rate Matrix File for Running Bacmeta Simulation
##'
##' Creates a  square matrix of dimensions \code{n_populations} x
##' \code{n_populations} where each element represents a migration rate as
##' defined in the \code{rates} parameter. The rate at row i, column j
##' represents migration from population i to population j. The matrix is
##' written to a tab-separated file named migration\[\code{suffix}\].input in
##' the directory \code{out_path}.
##'
##' This migration file can then be used in a Bacmeta simulation with
##' \code{n_populations} populations, in place of a single migration rate.
##' See Bacmeta readme for details.
##'
##' @param n_populations an integer defining the number of populations (and
##'        thereby the dimensions of the migration matrix)
##' @param rates the migration rates to associate with each population pair.
##'        If \code{NULL} (default), the resulting matrix will contain only
##'        zeros. Also accepted is a matrix of size \code{n_populations}^2,
##'        or a character vector of length \code{n_populations}^2, in which
##'        case the first \code{n_populations} values are used for row 1
##'        (representing the migration rates from population 1), and so on.
##' @param out_path The destination to which the migration file will be
##'        written. Must be an existing directory.
##' @param suffix A suffix that will be used for the migration file name. If
##'        \code{suffix} is \code{NULL} (default), the filename will simply be
##'        "migration.input", otherwise the suffix will be pasted in between
##'        "migration" and ".input". E.g. if \code{suffix} is 123, the filename
##'        will be "migration123.input". All symbols in \code{suffix} must be
##'        alphanumeric (A-Z, a-z, 0-9).
##' @return invisible \code{NULL}
##' @author Wiktor Gustafsson
##' @export
create_migration.input <- function(n_populations,
                                   rates = NULL,
                                   out_path = path_to_bacmeta(),
                                   suffix = NULL) {

    stopifnot(is.numeric(n_populations) && n_population %% 1 == 0)

    if (!dir.exists(out_path)) {
        stop(ifelse(file.exists(out_path),
                    "'out_path' must be a directory (not a file).",
                    "'out_path' is invalid (no such directory)."))
    }

    if (!is.null(suffix) && !is_alphanumeric(suffix))
        stop(paste0("All symbols in 'sufffix' must be alphanumeric ",
                    "(A-Z, a-z or 0-9)."))

    out_path <- file.path(normalizePath(out_path),
                          paste0("simu.", suffix, ".input"))

    if (is.null(rates))
        return(matrix(0, nrow = n_populations, ncol = n_populations))

    stopifnot(is.numeric(rates))
    rates[is.na(rates)] <- 0

    if (is.matrix(rates)) {
        if (nrow(rates) != n_populations || ncol(rates) !=  n_populations)
            stop(paste0("'rates' matrix must be 'n_populations' x ",
                        "'n_populations' in dimensions"))
    } else {
        if (length(rates) != n_populations ** 2)
            stop("'rates' vector must be 'n_populations' ^ 2 in length.")

        rates <- t(matrix(rates, nrow = n_populations, ncol = n_populations))
    }

    write.table(rates,
                out_path,
                col.names = FALSE,
                row.names = FALSE,
                sep = "\t")

}
