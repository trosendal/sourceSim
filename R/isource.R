##' Path to isource
##' @noRd
path_to_isource <- function(subpath = "")
    file.path(.bacmeta_env$isource, subpath)

##' Check if isource is Compiled
##'
##' If isource is compiled, each C++ script in its src directory should have
##' a corresponding .o file, and there should be a "isource" file in isource root.
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
