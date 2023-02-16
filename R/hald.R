##' hald functions
##' @param ... other arguments
##' @param x An object to run the hald method on
##' @export
##' @return A character vector of the proportions
hald <- function(x, ...) UseMethod("hald")

hald.sourceSim_result <- function() {

}

hald.data.frame <- function() {

}
