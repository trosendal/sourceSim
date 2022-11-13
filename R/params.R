##' Check if a paramfile is valid
##' @noRd
valid_paramfile <- function(path) {
    path <- normalizePath(path, mustWork = TRUE)

    if (!grepl("^simu[0-9]*\\.input$", basename(path)))
        return(FALSE)

    params <- read.table(path, sep = ":")
    template <- read_paramfile()

    if (!setequal(params[, 1], template[, 1]))
        return(FALSE)

    TRUE
}

##' Read (custom or default) bacmeta simulation parameter file.
##'
##' @param path A path to a valid parameter file. The default is NULL
##'     which fetches a default parameter file from the package.
##' @param as_list Return the results as a list? If \code{FALSE}
##'     (default), simply returns the parameters as a two-column
##'     data.frame, where column 1 is the parameter names and column 2
##'     is their values. If \code{TRUE}, instead returns a name list
##'     where each element is a parameter.
##' @return a data.frame or list (depending on \code{as_list})
##'     containing the parameters.
##' @importFrom utils read.table
##' @export
##'
read_paramfile <- function(path = NULL, as_list = FALSE) {
    if (is.null(path)) {
        path <- file.path(system.file("bacmeta",
                                      package = "sourceSim"),
                          "default.params")
    } else {
        stopifnot(valid_paramfile(path))
    }

    params <- read.table(path, sep = ":")

    if (isTRUE(as_list)) {
        plist <- as.list(params[, 2])
        names(plist) <- params[, 1]
        return(plist)
    }

    names(params) <- c("param", "value")
    params
}

##' copy_paramfile
##'
##' Copies a bacmeta parameter file (default or custom) to a given destination.
##'
##' @param from Either \code{NULL} (default) which means the default param file
##'        stored in the package, or a path to a valid bacmeta simulation file.
##' @param to A directory in which to copy the parameter file.
##' @param default_params if \code{from} is \code{NULL} and this is
##'        \code{TRUE}, the paramfile will have the name "default.params" in
##'        the destination directory. Otherwise if \code{FALSE} (default),
##'        it will be name "simu.input". If \code{from} is not \code{NULL},
##'        parameter is ignored.
##' @return The path to the copied parameter file, including file name.
##'
##' @export
copy_paramfile <- function(from = NULL, to = getwd(), default_params = FALSE) {
    to <- normalizePath(to, mustWork = TRUE)
    if (!dir.exists(to))
        stop("'to' must be a directory")
    if (is.null(from)) {
        from <- file.path(system.file("bacmeta",
                                      package = "sourceSim"),
                          "default.params")

        filename <- if (isTRUE(default_params))
            "default.params"
        else
            "simu.input"

        to <- file.path(to, filename)
    } else {
        stopifnot(valid_paramfile(from))
        to <- file.path(to, basename(from))
    }

    if (!file.copy(from, to, overwrite = TRUE))
        stop("Copy of 'input' file to simulation directory failed")

    to
}

##' Check if a migration file is valid
##' @noRd
valid_migrationfile <- function(path,
                                n_populations) {
    stopifnot(file.exists(path))

    if (!grepl("^migration[0-9]*\\.input$", basename(path)))
        return(FALSE)

    migration <- as.matrix(read.table(path, sep = "\t"))

    if (is.character(migration)) {
        if (!all(migration[, ncol(migration)] == "*"))
            return(FALSE)
        migration <- migration[, -ncol(migration)]
        class(migration) <- "numeric"
    }

    if (any(is.na(migration)))
        return(FALSE)

    if (nrow(migration) != ncol(migration) ||
        nrow(migration) != n_populations ||
        !(all(sapply(migration, is.numeric)))) {
        return(FALSE)
    }

    TRUE
}

##' Generate a Parameter File for Running Bacmeta Simulation
##'
##' Loads a predefined parameter template file with pre-filled Bacmeta
##' simulation parameters (see README for details), replaces the
##' parameter values as supplied in \code{params}, and writes the
##' modified parameter file to a file named
##' simu\[\code{suffix}\].input in the directory \code{out_path}.
##'
##' @param params A named list of parameters to be written into the
##'     parameter file. The names of all values in \code{params} must
##'     match a in the template file. If \code{params} is \code{NULL}
##'     (default), an unmodified template file will be written.
##' @param out_path The destination to which the parameter file will
##'     be written. Must be an existing directory.
##' @param suffix A suffix that will be used for the param file
##'     name. If \code{suffix} is \code{NULL} (default), the filename
##'     will simply be "simu.input", otherwise the suffix will be
##'     pasted in between "simu" and ".input". E.g. if \code{suffix}
##'     is 123, the filename will be "simu123.input". All symbols in
##'     \code{suffix} must be alphanumeric (A-Z, a-z, 0-9).
##' @return The path to the generated parameter file.
##' @importFrom utils write.table
##' @author Wiktor Gustafsson
##' @export
create_simu_input <- function(params = NULL,
                              out_path = getwd(),
                              suffix = NULL) {

    if (!dir.exists(out_path)) {
        stop(ifelse(file.exists(out_path),
                    "'out_path' must be a directory (not a file).",
                    "'out_path' is invalid (no such directory)."))
    }

    filename <- "simu"

    if (!is.null(suffix)) {
        suffix <- as.character(suffix)
        if (!is_alphanumeric(suffix))
            stop(paste0("All symbols in 'sufffix' must be alphanumeric ",
                        "(A-Z, a-z or 0-9)."))
        filename <- paste0(filename, suffix)
    }

    out_path <- file.path(normalizePath(out_path),
                          paste0(filename, ".input"))

    template <- read_paramfile()

    if (!is.null(params)) {
        stopifnot(is.list(params))

        if (!all(names(params) %in% template[, 1]))
            stop(paste0(
                "Some parameters supplied in 'params' have invalid names. ",
                "See sourceSim README for list of legal parameters."))

        if (!all(sapply(params, is.numeric)))
            stop("All supplied parameters must be numeric")

        ## Due to a (likely) bug in bacmeta, a constant migration rate
        ## still needs to be given if a migration matrix is
        ## supplied. These can be anything non-zero and don't affect
        ## the simulation:
        if (!is.null(params$MIGI)) {
            if (params$MIGI == 1) {
                params$MIGR <- 0.01
            }
        }

        template[match(names(params), template[, 1]), 2] <- unlist(params)
    }

    write.table(
        template,
        file = out_path,
        sep = ": \t",
        row.names = FALSE,
        col.names = FALSE,
        quote = FALSE
    )

    out_path
}

##' Generate a Migration Rate Matrix File for Running Bacmeta Simulation
##'
##' Creates a square matrix of dimensions \code{n_populations} x
##' \code{n_populations} where each element represents a migration
##' rate as defined in the \code{rates} parameter. The rate at row i,
##' column j represents migration from population i to population
##' j. The matrix is written to a tab-separated file named
##' migration\[\code{suffix}\].input in the directory \code{out_path}.
##'
##' This migration file can then be used in a Bacmeta simulation with
##' \code{n_populations} populations, in place of a single migration
##' rate.  See Bacmeta readme for details.
##'
##' @param rates the migration rates to associate with each population
##'     pair.  If \code{NULL} (default), the resulting matrix will
##'     contain only zeros. Otherwise, a square matrix is expected.
##' @param out_path The destination to which the migration file will
##'     be written. Must be an existing directory.
##' @param suffix A suffix that will be used for the migration file
##'     name. If \code{suffix} is \code{NULL} (default), the filename
##'     will simply be "migration.input", otherwise the suffix will be
##'     pasted in between "migration" and ".input". E.g. if
##'     \code{suffix} is 123, the filename will be
##'     "migration123.input". All symbols in \code{suffix} must be
##'     alphanumeric (A-Z, a-z, 0-9).
##' @importFrom utils write.table
##' @return \code{out_path}
##' @author Wiktor Gustafsson
##' @export
create_migration_input <- function(rates = NULL,
                                   out_path = getwd(),
                                   suffix = NULL) {

    if (!dir.exists(out_path)) {
        stop(ifelse(file.exists(out_path),
                    "'out_path' must be a directory (not a file).",
                    "'out_path' is invalid (no such directory)."))
    }

    filename <- "migration"

    if (!is.null(suffix)) {
        suffix <- as.character(suffix)
        if (!is_alphanumeric(suffix))
            stop(paste0("All symbols in 'sufffix' must be alphanumeric ",
                        "(A-Z, a-z or 0-9)."))
        filename <- paste0(filename, suffix)
    }

    out_path <- file.path(normalizePath(out_path),
                          paste0(filename, ".input"))

    stopifnot(is.matrix(rates))
    rates[is.na(rates)] <- 0

    if (nrow(rates) !=  ncol(rates))
        stop(paste0("'rates' matrix must be square (equal amount of rows and ",
                    "columns)."))

    ## A bacmeta bug requires an extra column for some reason
    rates <- cbind(rates, rep("*", nrow(rates)))

    write.table(rates,
                out_path,
                col.names = FALSE,
                row.names = FALSE,
                sep = "\t",
                quote = FALSE)

    out_path
}
