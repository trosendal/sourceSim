##' read_results
##'
##' Reads the sequence output of bacmeta simulation and produces a
##' summary dataset. Each row in the dataset represents a genotype
##' produced in the simulation. One genotype is one unique combination
##' of locus alleles. The dataset contains information about:
##'
##' \itemize{
##'   \item{Pop_n} The number of occurences of the genotype in each
##'   simulated population, in columns named "Pop_0", "Pop_1", etc;
##'
##'   \item{total} the total number of occurences (sum of all
##'   population columns);
##'
##'   \item{seqID} A unique numeric genotype/sequence ID;
##'
##'   \item{MLST} An "MLST" ID which is a concatenation of the unique
##'   allele IDs of all simulated loci, separated by hyphens;
##'
##'   \item{seq} A DNA sequence which is simply the sequences of all
##'   alleles pasted together.
##' }
##'
##' \code{seqID}, \code{MLST} and \code{seq} are all unique
##' identifiers.
##'
##' @param sequences A path to a "Sequences" file which has been
##'     output from running bacmeta.
##' @param initial A path to an "InitialSequences" file which has been
##'     output from running bacmeta. This is only used to match
##'     against the data in \code{sequences} and exclude the initial
##'     sequences from the dataset.
##' @return a \code{data.frame} with the information described above.
##' @export
read_results <- function(sequences,
                         initial) {

    ob <- readLines(sequences)
    initial <- readLines(initial)
    initial <- paste0(trimws(initial[trimws(initial) != ""]), collapse = "")

    names <- unlist(strsplit(ob[1], ";"))

    ob <- ob[-c(1, grep("^;", ob))]

    ob <- do.call("rbind", lapply(ob, function(x) {
        df <- as.data.frame(matrix(unlist(strsplit(x, ";")), nrow = 1))
        names(df) <- names
        df
    }))

    genes <- grep("^Gene", names(ob), value = T)

    for (x in genes)
        ob[, x] <- as.factor(ob[, x])

    ob$MLST <- apply(as.data.frame(lapply(ob[, genes], as.numeric)), 1,
                     function(x) paste(x, collapse = "-"))

    ob <- ob[, c("Population", "Amount", "MLST", genes)]

    ob$seq <- sapply(seq_len(nrow(ob)), function(i) {
        paste(unlist(ob[i, genes]), collapse = "")
    })

    ob$Amount <- as.numeric(ob$Amount)
    ob$Population <- as.factor(ob$Population)
    ob$ID <- seq_len(nrow(ob))
    ob$seqID <- as.numeric(as.factor(ob$seq))

    df <- do.call("rbind", lapply(unique(ob$seqID), function(x) {
        ob_inner <- ob[ob$seqID == x, c("Population", "Amount", "MLST", "seq")]
        a <- tapply(ob_inner$Amount, ob_inner$Population, "sum")
        names <- paste0("Pop_", levels(ob$Population))
        a <- as.list(a); names(a) <- names
        a[is.na(a)] <- 0
        df <- do.call(data.frame, a)
        df$total <- sum(unlist(a))
        df$seqID <- x
        df$MLST <- ob_inner$MLST[1]
        df$seq <- toupper(ob_inner$seq[1])
        df
    }))

    df[(tolower(df$seq) != tolower(initial)), ]
}

##' plot
##'
##' Plots a phylogenetic tree of sequences that are the result of a bacmeta
##' simulation
##'
##' @param x A \code{data.frame} with bacmeta simulation result sequences,
##'        as produced by \code{read_results()}.
##' @param piecol A vector of colours
##' @param legend A logical if a legend shoulöd be included
##' @param ... Other arguments
##'
##' @return A "\code{recordedplot}" object with the phylogenetic tree.
##' @import ape
##' @importFrom stats dist hclust
##' @importFrom grDevices rainbow
##' @export
plot.sourceSim_result <- function(x, piecol = NULL, legend = FALSE, ...) {
    df <- x$population

    a <- df$seq
    names(a) <- df$seqID
    a <- strsplit(a, "")
    a <- as.DNAbin(a)
    dis <- dist(a)
    hc <- hclust(dis, "average")
    hc <- as.phylo(hc)
    ## Check that we can match by the order:
    stopifnot(identical(as.numeric(hc$tip.label), df$seqID))

    labs <- as.matrix(df[, grep("^Pop_", names(df), value = T)])
    rownames(labs) <- df$seqID
    if (is.null(piecol))
        piecol <- rainbow(ncol(labs))

    plot(hc, type = "fan", show.tip.label = FALSE)
    tiplabels(pie = labs, cex = (df$total ^ 0.3) / 5, piecol = piecol)

    if (legend)
        legend("topleft", legend = colnames(labs), col = piecol, pch = 20)
}

##' print.sourceSim_result
##'
##' @param x A bacmeta simulation result
##' @param ... Other arguments
##' @importFrom utils head
##' @export
print.sourceSim_result <- function(x, ...) {

    x$population$seq <- paste0(substring(x$population$seq,
                                         first = 1,
                                         last = 5), "...")
    cat("\nFirst few rows of simulated Data:\n\n")
    print(head(x$population))
    maxl <- max(sapply(x$parameters, nchar))
    cat("\nSimulation Parameters:\n\n| param | value",
        rep(" ", maxl - 5, sep = ""),
        "|\n|-------|",
        rep("-", maxl + 1), "|\n", sep = "")
    mapply(function(x, y) {
        cat("| ", y, "  | ", x, rep(" ", maxl - nchar(x)), "|\n", sep = "")
        invisible()
    }, x = x$parameters,
    y = names(x$parameters))
    cat("|-------|", rep("-", maxl + 1), "|\n", sep = "")
}
