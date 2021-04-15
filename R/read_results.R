ob <- readLines("outputs/Sequences20000-123.csv")
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
ob$ID <- seq_len(nrow(ob))

## Generate a labeled tree with all the individuals
library(ape)
a <- paste0(ob$Gene0, ob$Gene1, ob$Gene2, ob$Gene3, ob$Gene4, ob$Gene5, ob$Gene6)
names(a) <- ob$ID
a <- strsplit(a, "")
a <- as.DNAbin(a)
dis <- dist(a)
hc <- hclust(dis, "average")
hc <- as.phylo(hc)
plot(hc, type = "unroot")
