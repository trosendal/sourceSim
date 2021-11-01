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
##' @param input The parameters for running bacmeta.  Accepts one of
##'     two options:
##'     \itemize{
##'       \item \code{NULL} (default), which uses the
##'       \code{default.params} which is included in the package;
##'
##'       \item a named list of parameters, where each name is a
##'       bacmeta-accepted parameter, and the values are numeric
##'       (note: parameters affecting file outputs from the
##'       simulation, i.e. \code{SEQI} and \code{ISEQ}, will be
##'       ignored to ensure the desired results can be captured).
##'     }
##' @param migration Migration matrix for bacmeta simulation. Only
##'     used if the \code{MIGI} parameter in the main bacmeta
##'     parameter file is set to 1, otherwise ignored. Accepts one of
##'     two options:
##'     \itemize{
##'       \item \code{NULL} (default), which creates a zero matrix
##'       meaning there is no migration between any populations. This
##'       requires the 'NPOP' parameter to be set in the main parameter file.
##'
##'       \item a numeric matrix of dimensions n x n,
##'       where n is the number of populations (as defined by the 'NPOP'
##'       parameter in the main parameter file).
##' }
##' @param plot Produce a phylogenetic plot of the results? Plots to
##'     default graphic device. Default is \code{FALSE}.
##' @return A \code{data.frame} with the resulting DNA sequences from
##'     running the simulation (see \code{read_results() function for
##'     details}).
##' @export
simu <- function(input = NULL,
                 migration = NULL,
                 plot = FALSE) {

    if (is.null(input))
        input <- list()
    stopifnot(is.list(input))

    simu_dir <- tempdir()
    wd <- setwd(simu_dir)
    on.exit(setwd(wd))

    if (is.null(input$GENR))
        input$SEQI <- 10000
    else
        input$SEQI <- input$GENR

    input$ISEQ <- 1

    if(!is.null(migration))
        input$MIGI <- 1

    paramfile <- create_simu.input(params = input,
                                   out_path = simu_dir)

    default_file <- copy_paramfile(from = NULL,
                                   to = simu_dir,
                                   default.params = TRUE)

    params <- read_paramfile(path = paramfile, as_list = TRUE)

    if (params$MIGI == 1) {
        if (is.null(migration)) {
            if (is.null(params$NPOP))
                stop(paste0("Migration matrix expected, but no matrix ",
                            "supplied and none can be generated since ",
                            "NPOP param is not set."))
            stopifnot(identical(params$NPOP %% 1, 0))
            migration <- matrix(nrow = params$NPOP, ncol = params$NPOP)
        }
        migration <- create_migration.input(rates = migration,
                                            out_path = simu_dir)
    }

    compile_bacmeta()

    system(path_to_bacmeta("simu"))

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

    results <- list(population = results,
                    parameters = params)

    class(results) <- c("sourceSim_result", class(results))

    if (isTRUE(plot))
      plot(results)

    results
}

##' Output data from a simulation
##'
##' The output of a simulation of 20000 generation of *C. jejuni* from
##' 3 metapopulations with variable overlap
##'
##' @name result
##' @docType data
##' @source From a simulation
##' @keywords datasets
##'
NULL
