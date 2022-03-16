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
##' @importFrom stats rgamma
##' @export
##' @return A matrix
rdirichlet <- function(n = 1, alpha = c(1, 1)) {
    g <- do.call("cbind", lapply(alpha, function(x) {
        rgamma(n, shape = x)
    }))
    g / matrix(rowSums(g), nrow = n, ncol = length(alpha))
}
