##' Path to bacmeta
##' @noRd
path_to_bacmeta <- function()
    bacmeta_path_object()$bacmeta

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
    all(!is.null(path_to_bacmeta()),
        file.exists(file.path(
            path_to_bacmeta(), bacmeta_binaries()
        )))
}

##' Compile bacmeta
##' @export
compile_bacmeta <- function() {
    cat("Compiling bacmeta...\n")
    if (is_bacmeta_compiled()) {
        cat("Bacmeta already compiled.\n")
    } else {
        assign("bacmeta", tempdir(), envir = bacmeta_path_object())
        on.exit(setwd(wd))
        file.copy(from = system.file("bacmeta/src",
                                     package = "sourceSim"),
                  recursive = TRUE,
                  to = ".")
        wd <- setwd(file.path(bacmeta_path_object()$bacmeta, "src"))
        system("make")
        cat ("Compilation successful.\n")
    }
}

##' Remove bacmeta binaries ("uncompile" bacmeta)
##'
##' @export
clean_bacmeta <- function()
    unlink(file.path(path_to_bacmeta(), bacmeta_binaries()))


## Create a package environment to store path to compiled bacmeta
.sourceSim_env <- new.env()

##' Return the bactmeta environment
##'
##' @return The evironment that contains the path to compiled bactmeta
##' @export
bacmeta_path_object <- function()
  .sourceSim_env
