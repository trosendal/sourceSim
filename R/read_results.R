ob <- readLines("../data/Sequences20000-123.csv")
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

ob <- ob[,c("Population", "Amount", "MLST")]
