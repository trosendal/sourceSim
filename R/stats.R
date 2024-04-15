##' rdirichlet
##'
##' sample from a dirichlet distribution. Returns a matrix where each
##' row is a sample. The definition of a dirichlet is that each set of
##' n samples must sum to 1. Therefore the rowSums(rdirichlet()) is
##' always 1.
##'
##' @param n The number of samples from the dirichlet
##' @param alpha The shape of each gamma distribution. The length
##'     dictates the number of bins in the dirichlet.
##' @export
##' @return A matrix
rdirichlet <- function(n = 1, alpha = c(1, 1)) {
    g <- do.call("cbind", lapply(alpha, function(x) {
        stats::rgamma(n, shape = x)
    }))
    g / matrix(rowSums(g), nrow = n, ncol = length(alpha))
}

##' rmse
##'
##' Calculate the root mean square error of a series of observations to an
##' expected value. The formula is that of standard deviation, but variation
##' is measured around an expected (given) value rather than the mean of the
##' data.
##'
##' @param ob a numeric vector of observations
##' @param exp a numeric expected value. Must be of length 1.
##' @return a numeric value
##' @export
rmse <- function(ob, exp) {
    stopifnot(
        is.numeric(ob),
        is.numeric(exp),
        length(ob) > 0,
        length(exp) == 1
    )

    sqrt(sum((ob - exp)^2) / length(ob))
}
