##' Extract Sequences from FASTA Files
##'
##' Extract all sequence data from all FASTA files in a given directory.
##'
##' @param path the path to the directory containing the FASTA files
##' @return a named character vector where each element is a concatenation of
##'         all sequences found in each file in \code{path}, i.e., every non-
##'         header line of the file. The name of each element is the file name.
##' @author Wiktor Gustafsson
##' @export
get_sequences <- function(path = system.file("extdata/MLST_alleles",
                                           package = "sourceSim")) {

    if (!dir.exists(path))
        stop(sprintf("Path '%s' does not exist or is not a directory"))

    seq_files <-
        list.files(path, pattern = "^.*\\.(fas|fasta)$", full.names = T)

    if (!length(seq_files))
        stop(sprintf("No files with file extension .fas or .fasta found at %s",
                     path))

    ret <- sapply(seq_files, function(x) {
        seq <- readLines(x)
        paste(seq[!grepl("^>", seq)], collapse = "")
    })

    names(ret) <- basename(seq_files)
    ret
}

##' Calculate Nucleotide Frequencies in a Genetic Sequence
##'
##' Given a genetic sequence \code{seq}, calculates the relative frequencies of
##' each symbol as defined in \code{bases}.
##'
##' @param seq A character vector of length 1 containing a sequences of
##'        characters representing nucleotide bases. The symbols must all be
##'        IUPAC-accepted nucleotide symbols (see README for details).
##' @param bases A length 1 character vector defining which symbols to
##'        calculate the frequencies for. Like \code{seq}, all symbols in
##'        \code{bases} must be from the IUPAC list of nucleotide codes (see
##'        README for details). Calulations will be made on each unique symbol
##'        in \code{bases}; duplicates will be ignored.
##' @return a named vector of relative frequencies, each value and its name
##'         corresponding to a symbol in \code{bases}.
##' @author Wiktor Gustafsson
##' @export
base_frequencies <- function(seq, bases = "ATCG") {
    stopifnot(is.character(seq) && length(seq) == 1 &&
                  is.character(bases) && length(bases) == 1)

    if (!grepl("^[ATCGURYSWKMBDHVN\\.-]+$", bases))
        stop(paste0(sprintf("The 'bases' string %s is invalid. Must be", bases),
                    " a non-empty, uppercase string containing a selection of ",
                    "the IUPAC nucleotide codes \n",
                    "(A, C, G, T, U, R, Y, S, W, K, M, B, D, H, V, N, ., -),\n",
                    "without separators."))

    seq <- toupper(seq)

    if (grepl(paste0("[^", bases, "]+"), seq))
        stop(paste0("Non-DNA characters found in sequence. All characters ",
                    "must be one of [", bases, "]."))

    counts <- sapply(unique(charToRaw(bases)),
        function(x) lengths(regmatches(seq, gregexpr(rawToChar(x), seq))))

    names(counts) <- unique(sapply(charToRaw(bases), rawToChar))

    counts / nchar(seq)
}

##' Generate a Parameter File for Running Bacmeta Simulation
##'
##' Loads a predefined parameter template file with pre-filled Bacmeta
##' simulation parameters (see README for details), replaces the parameter
##' values as supplied in \code{params}, and writes the modified parameter file
##' to \code{out_path}.
##'
##' @param params A named list of parameters to be written into the parameter
##'        file. The names of all values in \code{params} must match a
##'        in the template file. If \code{params} is \code{NULL} (default),
##'        an unmodified template file will be written.
##' @param out_path The destination to which the parameter file will be
##'        written. If the path is a directory, the file will be written to
##'        "params.txt" under that directory. Otherwise it will be written to
##'        a file with the specified name (and overwrite any already existing
##'        file with that name).
##' @return invisible \code{NULL}
##' @author Wiktor Gustafsson
##' @export
create_paramfile <- function(params = NULL, out_path = ".") {
    stopifnot(is.list(params) || is.null(params))

    if (dir.exists(out_path)) {
        filename <- "params.txt"
    } else {
        filename <- basename(out_path)
        out_path <- dirname(out_path)

        if (!dir.exists(out_path) || identical(out_path, filename))
            stop("Invalid param out path")
    }

    out_path <- file.path(normalizePath(out_path), filename)


    template <- read.table(system.file("params.txt",
                                       package = "sourceSim"),
                           sep = ":")

    if (!is.null(params)) {
        if (!all(names(params) %in% template[, 1]))
            stop(paste0(
                "Some parameters supplied in 'params' have invalid names. ",
                "See sourceSim README for list of legal parameters."))

        if (!all(sapply(params, is.numeric)))
            stop("All supplied parameters must be numeric")

        template[match(names(params), template[, 1]), 2] <- unlist(params)
    }

    write.table(
        template,
        file = out_path,
        sep = ": \t",
        row.names = F,
        col.names = F,
        quote = F
    )

    invisible(NULL)
}
