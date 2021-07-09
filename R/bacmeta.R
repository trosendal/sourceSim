##' Path to bacmeta
##' @noRd
path_to_bacmeta <- function()
    system.file("bacmeta", package = "sourceSim")

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
    all(file.exists(file.path(
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
        wd <- setwd(file.path(path_to_bacmeta(), "src"))
        on.exit(setwd(wd))
        system("make")
        cat ("Compilation successful.\n")
    }
}

##' Remove bacmeta binaries ("uncompile" bacmeta)
##'
##' @export
clean_bacmeta <- function()
    unlink(file.path(path_to_bacmeta(), bacmeta_binaries()))
