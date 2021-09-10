##' Path to bacmeta
##' @noRd
path_to_bacmeta <- function(subpath = "")
    file.path(bacmeta_path_object()$bacmeta, subpath)

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
##' in the \code{.sourceSim_env} object and is accessible via
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
        assign("bacmeta", tempdir(), envir = bacmeta_path_object())
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
##'
##' \enumerate{
##'   \item Verifies that the supplied parameter files are in a valid format,
##'   \item copies the supplied bacmeta parameter files to a simulation
##'         directory,
##'   \item compiles bacmeta if necessary,
##'   \item runs bacmeta,
##'   \item copies outputs from the simulation to a specified output directory.
##' }
##'
##' @param input A path to a valid bacmeta parameter file, or \code{NULL}
##'        (default) which uses the default parameter file stored in the
##'        package.
##' @param migration A path to a valid bacmeta migration file. Only used if
##'        the "MIGI" is 1 in the parameter file supplied in \code{input}.
##'        Then, if \code{migration} is \code{NULL} (default), uses a
##'        zero-filled matrix with NPOP x NPOP dimensions (as specified in
##'        \code{input} file).
##' @param out_path An existing directory to which the simulation outputs
##'        should be copied. Default is current directory.
##' @param simu_dir A directory in which the simulation should be run. Default
##'        is a new temporary directory.
##' @return \code{out_path}
##'
##' @export
simu <- function(input = NULL,
                 migration = NULL,
                 out_path = getwd(),
                 simu_dir = tempdir()) {

    simu_dir <- normalizePath(simu_dir, mustWork = TRUE)
    out_path <- normalizePath(out_path, mustWork = TRUE)
    stopifnot(dir.exists(simu_dir), dir.exists(out_path))

    wd <- setwd(simu_dir)
    on.exit(setwd(wd))

    paramfile <- copy_paramfile(from = input, to = simu_dir)
    default_file <- copy_paramfile(from = NULL,
                                   to = simu_dir,
                                   default.params = TRUE)
    params <- read_paramfile(path = paramfile, as_list = TRUE)

    if (params$MIGI == 1) {
        if (is.null(migration)) {
          create_migration.input(n_populations = params$NPOP,
                                 )
        }


        if (!valid_migrationfile(
            migration,
            n_populations = params$NPOP))
            stop("Supplied 'migration' file has invalid name or format " ,
                 "(see bacmeta readme)")

        if (!file.copy(migration, simu_dir, overwrite = T))
            stop("Copy of 'migration' file to simulation directory failed")
    }

    compile_bacmeta()

    cmd <- path_to_bacmeta("simu")

    if (is.null(input))
        input <- "simu.input"

    simu_suffix <- regmatches(basename(input),
                              regexec("^simu([a-zA-Z0-9]*)\\.input$",
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

    simu_outputs <- file.path(simu_dir, "outputs")
    file.copy(simu_outputs, out_path, recursive = TRUE, overwrite = TRUE)

    out_path
}

## Create a package environment to store path to compiled bacmeta
.sourceSim_env <- new.env()

##' Return the bacmeta environment
##'
##' @return The environment that contains the path to compiled bacmeta
##' @export
bacmeta_path_object <- function()
  .sourceSim_env
