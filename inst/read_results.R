ob <- readLines("outputs/Sequences20000-2.csv")
init <- readLines("outputs/InitialSequences2.txt")
init <- paste0(trimws(init[trimws(init) != ""]), collapse = "")
## Just the names
names <- unlist(strsplit(ob[1], ";"))
## remove the extra lines
ob <- ob[!grepl("^;Gene0", ob)]
ob <- ob[-1]

## Now cleanup
ob <- do.call("rbind", lapply(ob, function(x) {
    df <- as.data.frame(matrix(unlist(strsplit(x, ";")), nrow = 1))
    names(df) <- names
    df
}))

ob$MLST <- paste(as.numeric(as.factor(ob$Gene0)),
                 as.numeric(as.factor(ob$Gene1)),
                 as.numeric(as.factor(ob$Gene2)),
                 as.numeric(as.factor(ob$Gene3)),
                 as.numeric(as.factor(ob$Gene4)),
                 as.numeric(as.factor(ob$Gene5)),
                 as.numeric(as.factor(ob$Gene6)), sep = "-")

ob <- ob[, c("Population", "Amount", "MLST", "Gene0", "Gene1",
             "Gene2", "Gene3", "Gene4", "Gene5", "Gene6")]
ob$Amount <- as.numeric(ob$Amount)
ob$Population <- as.factor(ob$Population)
ob$ID <- seq_len(nrow(ob))
ob$seq <- paste0(ob$Gene0, ob$Gene1, ob$Gene2, ob$Gene3, ob$Gene4, ob$Gene5, ob$Gene6)
ob$seqID <- as.numeric(as.factor(ob$seq))

## Generate a unique sequence level dataset to facilitate tip labeling
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
    df$seq <- ob_inner$seq[1]
    df
}))
View(df)
## Drop the initial sequence because it is over represented and in
## every population:
df <- df[(df$seq != init), ]
## Generate a labeled tree with unique seqeunces
library(ape)
a <- df$seq
names(a) <- df$seqID
a <- strsplit(a, "")
a <- as.DNAbin(a)
dis <- dist(a)
hc <- hclust(dis, "average")
hc <- as.phylo(hc)
## Check that we can match by the order:
stopifnot(identical(as.numeric(hc$tip.label),
          df$seqID))
labs <- as.matrix(df[, c("Pop_0", "Pop_1", "Pop_2")])
rownames(labs) <- df$seqID
plot(hc, type = "fan", show.tip.label = FALSE)
tiplabels(pie = labs, cex = (df$total^0.3)/10)
