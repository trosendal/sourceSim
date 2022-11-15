library(sourceSim)
files <- list.files(c("results/experiment6"), pattern = ".Rds", full.names = TRUE)

df <- do.call("rbind", lapply(seq_len(length(files)), function(i) {
    ob <- readRDS(files[i])
    do.call("rbind", lapply(ob, function(j) {
        k <- j$result$population
        overlap <- sum(k[, 1] > 0 & k[, 2] > 0) / sum(k[, 1] > 0 | k[,2] > 0)
        data.frame(ob0 = j$attribution[1, 1],
                   ob1 = j$attribution[1, 2],
                   ob2 = j$attribution[1, 3],
                   mig01 = j$migration[2, 1],
                   mig02 = j$migration[3, 1],
                   mig12 = j$migration[3, 2],
                   mig10 = j$migration[1, 2],
                   mig20 = j$migration[1, 3],
                   mig21 = j$migration[2, 3],
                   ex0 = j$sampling[1],
                   ex1 = j$sampling[2],
                   ex2 = j$sampling[3],
                   overlap = overlap)
    }))
}))

df$err0 <- abs(df$ob0 - df$ex0)
df$err1 <- abs(df$ob1 - df$ex1)
df$err2 <- abs(df$ob2 - df$ex2)
df$migtot <- df$mig01 + df$mig10

## The relationship between migration and error
cairo_pdf("plots/large_mig_error.pdf", width = 12, height = 4)
par(mfrow = c(1, 3))
boxplot(df$err0 ~ cut(df$migtot, 40), xlab = "sum migration rate (A\u2194B)", ylab = "Error of attribution (A)", ylim = c(0,0.50))
boxplot(df$err1 ~ cut(df$migtot, 40), xlab = "sum migration rate (A\u2194B)", ylab = "Error of attribution (B)", ylim = c(0,0.50))
boxplot(df$err2 ~ cut(df$migtot, 40), xlab = "sum migration rate (A\u2194B)", ylab = "Error of attribution (C)", ylim = c(0,0.50))
dev.off()

## The relationship between migration and variance of error (Remember variance is on the squared scale when compared to error)
cairo_pdf("plots/large_mig_variance.pdf", width = 12, height = 4)
par(mfrow = c(1, 3))
plot(tapply(df$err0, cut(df$migtot, 40), var), xaxt="n", xlab = "sum migration rate (A\u2194B)", ylab = "Variance of error (A)", ylim = c(0, 0.05))
axis(side = 1, at = 1:40, labels = levels(cut(df$migtot, 40)))
plot(tapply(df$err1, cut(df$migtot, 40), var), xaxt="n", xlab = "sum migration rate (A\u2194B)", ylab = "Variance of error (B)", ylim = c(0, 0.05))
axis(side = 1, at = 1:40, labels = levels(cut(df$migtot, 40)))
plot(tapply(df$err2, cut(df$migtot, 40), var), xaxt="n", xlab = "sum migration rate (A\u2194B)", ylab = "Variance of error (C)", ylim = c(0, 0.05))
axis(side = 1, at = 1:40, labels = levels(cut(df$migtot, 40)))
dev.off()

## Relationship between overlap and error
## OBS: The first 2 bins only have 1 observation
table(cut(df$overlap, 20))
cairo_pdf("plots/large_mig_error_vs_overlap.pdf", width = 12, height = 4)
par(mfrow = c(1, 3))
boxplot(df$err0 ~ cut(df$overlap, 20), xlab = "relative overlap (A\u2194B)", ylab = "Error of attribution (A)", ylim = c(0,0.50))
boxplot(df$err1 ~ cut(df$overlap, 20), xlab = "relative overlap (A\u2194B)", ylab = "Error of attribution (B)", ylim = c(0,0.50))
boxplot(df$err2 ~ cut(df$overlap, 20), xlab = "relative overlap (A\u2194B)", ylab = "Error of attribution (C)", ylim = c(0,0.50))
dev.off()

## The relationship between migration and variance of error (Remember variance is on the squared scale when compared to error)
cairo_pdf("plots/large_mig_variance_vs_overlap.pdf", width = 12, height = 4)
par(mfrow = c(1, 3))
plot(tapply(df$err0, cut(df$overlap, 20), var), xaxt="n", xlab = "relative overlap (A\u2194B)", ylab = "Variance of error (A)", ylim = c(0, 0.05))
axis(side = 1, at = 1:20, labels = levels(cut(df$overlap, 20)))
plot(tapply(df$err1, cut(df$overlap, 20), var), xaxt="n", xlab = "relative overlap (A\u2194B)", ylab = "Variance of error (B)", ylim = c(0, 0.05))
axis(side = 1, at = 1:20, labels = levels(cut(df$overlap, 20)))
plot(tapply(df$err2, cut(df$overlap, 20), var), xaxt="n", xlab = "relative overlap (A\u2194B)", ylab = "Variance of error (C)", ylim = c(0, 0.05))
axis(side = 1, at = 1:20, labels = levels(cut(df$overlap, 20)))
dev.off()

## Migration versus overlap
cairo_pdf("plots/large_migration_vs_overlap.pdf")
boxplot(df$overlap ~ cut(df$migtot, 40), xlab = "sum migration rate (A\u2194B)", ylab = "relative overlap (A\u2194B)", ylim = c(0, 0.9))
dev.off()
