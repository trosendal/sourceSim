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
        if (isFALSE(quiet)) cat ("Compilation successful.\n")
    }
}

##' simulate
##'
##' @export
simulate <- function(input,
simulate <- function(input = NULL,
                     out_path,
                     migration = NA_character_,
                     simu_dir = tempdir(),
                     keep_simufiles = FALSE,
                     plot = FALSE) {

    simu_dir <- normalizePath(simu_dir, mustWork = TRUE)
    out_path <- normalizePath(out_path, mustWork = TRUE)

    if (!file.copy(input, simu_dir, overwrite = TRUE))
        stop("Copy of 'input' file to simulation directory failed")

    params <- read_paramfile(path = input, as_list = TRUE)

    if (params$MIGI == 1) {
        if (is.na(migration))
            stop("Simulation setup requires migration.input file, but none
                 supplied")

        if (!valid_migrationfile(
            migration,
            n_populations = params$NPOP))
            stop("Supplied 'migration' file has invalid name or format " ,
                 "(see bacmeta readme)")

        if (!file.copy(migration, simu_dir, overwrite = T))
            stop("Copy of 'migration' file to simulation directory failed")
    }

    wd <- setwd(simu_dir)
    on.exit(setwd(wd))

    compile_bacmeta(quiet = TRUE)

    cmd <- path_to_bacmeta("simu")

    simu_suffix <- regmatches(basename(input),
                              regexec("^simu([a-zA-Z0-9]*)\\.input$",
                                      basename(input)))[[1]][[2]]

    migration_suffix <- regmatches(basename(migration),
                              regexec("^simu([a-zA-Z0-9]*)\\.input$",
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
    if (isFALSE(keep_simufiles))
        on.exit(unlink(simu_outputs, recursive = TRUE), add = TRUE)

    return(out_path)
}

## Create a package environment to store path to compiled bacmeta
.sourceSim_env <- new.env()

##' Return the bacmeta environment
##'
##' @return The environment that contains the path to compiled bacmeta
##' @export
bacmeta_path_object <- function()
  .sourceSim_env
