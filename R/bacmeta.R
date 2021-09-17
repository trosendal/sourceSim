## Create a package environment to store path to compiled bacmeta
.bacmeta_env <- new.env()

##' Path to bacmeta
##' @noRd
path_to_bacmeta <- function(subpath = "")
    file.path(.bacmeta_env$bacmeta, subpath)

##' These files exist if bacmeta is compiled
##' @noRd
bacmeta_binaries <- function()
    c("simu",
      paste0(
          "src/",
          c("Bac.o",
            "Gene.o",
            "Para.o",
            "Pop.o",
            "Simu.o",
            "Superpop.o")
      ))

##' Check if Bacmeta is Compiled
##'
##' If bacmeta is compiled, each C++ script in its src directory should have
##' a corresponding .o file, and there should be a "Simu" file in bacmeta root.
##'
##' @return \code{TRUE} if bacmeta is compiled, otherwise \code{FALSE}.
is_bacmeta_compiled <- function() {
    all(length(path_to_bacmeta()) > 0,
        file.exists(file.path(
            path_to_bacmeta(), bacmeta_binaries()
        )))
}

##' Compile bacmeta
##'
##' Copies bacmeta source files to an environment-specific temporary directory
##' and compiles bacmeta there. The path to the compilation directory is saved
##' in the \code{.bacmeta_env} object and is accessible via
##' \code{path_to_bacmeta()}. If bacmeta is already compiled, does nothing.
##'
##' @param quiet If \code{TRUE}, prints no compilation progress to console.
##'        Default is \code{FALSE}.
##'
##' @export
compile_bacmeta <- function(quiet = FALSE) {
  if (isFALSE(quiet)) cat("Compiling bacmeta...\n")
    if (is_bacmeta_compiled()) {
      if (isFALSE(quiet)) cat("Bacmeta already compiled.\n")
    } else {
        assign("bacmeta", tempdir(), envir = .bacmeta_env)
        file.copy(from = system.file("bacmeta/src",
                                     package = "sourceSim"),
                  recursive = TRUE,
                  to = path_to_bacmeta())
        wd <- setwd(path_to_bacmeta("src"))
        on.exit(setwd(wd))
        system("make", ignore.stdout = quiet)
        if (isFALSE(quiet)) cat("Compilation successful.\n")
    }
}

##' simu
##'
##' Runs bacmeta simulation:
##' \enumerate{
##'   \item Verifies that the supplied parameter files are in a valid format,
##'   \item copies the supplied bacmeta parameter files to a simulation
##'         directory,
##'   \item compiles bacmeta if necessary,
##'   \item runs bacmeta,
##'   \item copies outputs from the simulation to a specified output directory.
##' }
##'
##' @param input The parameters for running bacmeta.
##'        Accepts one of three options:
##'        \itemize{
##'          \item \code{NULL} (default), which uses the \code{default.params}
##'          which is included in the package;
##'          \item a named list of parameters, where each name is a bacmeta-
##'          accepted parameter, and the values are numeric (note: parameters
##'          affecting file outputs from the simulation will be ignored to
##'          ensure the desired results can be captured);
##'          \item a path to an existing, correctly formatted bacmeta parameter
##'          file.
##'        }
##' @param migration Migration matrix for bacmeta simulation. Only used if
##'        the MIGI parameter in the main bacmeta parameter file is set to 1,
##'        otherwise ignored. Accepts one of three options:
##'        \itemize{
##'          \item \code{NULL} (default), which creates a zero matrix
##'          meaning there is no migration between any populations;
##'          \item a numeric matrix/vector of dimensions/length n x n, where n
##'          is the number of populations (as defined by the NPOP parameter in
##'          the main parameter file);
##'          \item a path to an existing, correctly formatted bacmeta migration
##'          file.
##'        }
##' @param suffix A suffix to append to the file name of potential parameter
##'        and/or migration files. This is ignored if \code{input} and
##'        \code{migration} are both paths to existing files.
##' @param simu_dir A directory in which the simulation should be run. Default
##'        is a new temporary directory.
##' @param plot Produce a phylogenetic plot of the results? Plots to
##'        default graphic device.
##' @return A \code{data.frame} with the resulting DNA sequences from running
##'         the simulation (see \code{read_results() function for details}).
##'
##' @export
simu <- function(input = NULL,
                 migration = NULL,
                 suffix = NULL,
                 simu_dir = tempdir(),
                 plot = TRUE) {

    simu_dir <- normalizePath(simu_dir, mustWork = TRUE)
    stopifnot(dir.exists(simu_dir))

    wd <- setwd(simu_dir)
    on.exit(setwd(wd))

    if (is.null(input) || is.list(input)) {
        if (is.null(input)) input <- list()

        input$SEQS <- ifelse(is.null(input$GENR), 10000, input$GENR)
        input$ISEQ <- 1

        paramfile <- create_simu.input(params = input,
                                       out_path = simu_dir,
                                       suffix = suffix)
    } else if (is.character(input) && length(input) == 1 && file.exists(input))
        paramfile <- copy_paramfile(from = input, to = simu_dir)
    else
        stop("Invalid value for parameter 'input'")

    default_file <- copy_paramfile(from = NULL,
                                   to = simu_dir,
                                   default.params = TRUE)

    params <- read_paramfile(path = paramfile, as_list = TRUE)

    if (params$MIGI == 1) {
        if (is.null(migration) || is.numeric(migration))
            migration <- create_migration.input(n_populations = params$NPOP,
                                                rates = migration,
                                                out_path = simu_dir,
                                                suffix = suffix)
        else if (file.exists(migration)) {
            if (!valid_migrationfile(
                migration,
                n_populations = params$NPOP))
                stop("Supplied 'migration' file has invalid name or format " ,
                     "(see bacmeta readme)")

            if (!file.copy(migration, simu_dir, overwrite = T))
                stop("Copy of 'migration' file to simulation directory failed")

        } else stop("Invalud value for parameter 'migration'")
    }

    compile_bacmeta()

    cmd <- path_to_bacmeta("simu")

    if (is.null(input))
        input_suffix <- ""
    else
        simu_suffix <- regmatches(basename(input),
                              regexec(
                                  "^simu([a-zA-Z0-9]*)\\.input$",
                                  basename(input)))[[1]][[2]]

    if (is.null(migration))
        migration_suffix <- ""
    else
        migration_suffix <- regmatches(basename(migration),
                                       regexec(
                                           "^migration([a-zA-Z0-9]*)\\.input$",
                                           basename(migration)))[[1]][[2]]

    if (nchar(simu_suffix) || nchar(migration_suffix)) {
        if (identical(simu_suffix, migration_suffix))
            cmd <- paste(cmd, "-b", simu_suffix)
        else {
            if (nchar(simu_suffix))
                cmd <- paste(cmd, "-p", simu_suffix)
            if (nchar(migration_suffix))
                cmd <- paste(cmd, "-m", migration_suffix)
        }
    }

    system(cmd)

    sequences <- file.path(simu_dir,
                                  "outputs",
                                  paste0("Sequences",
                                         params$GENR,
                                         "-",
                                         params$OPFN,
                                         ".csv"))

    initial_sequences <- file.path(simu_dir,
                                   "outputs",
                                   paste0("InitialSequences",
                                          params$OPFN,
                                          ".txt"))

    results <- read_results(sequences = sequences,
                            initial = initial_sequences)

    if (isTRUE(plot))
      plot_results(results)

    results
}
