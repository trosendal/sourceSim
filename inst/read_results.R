library(sourceSim)
files <- list.files("results/", full.names = TRUE)

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
                   ex0 = j$sampling[1],
                   ex1 = j$sampling[2],
                   ex2 = j$sampling[3],
                   overlap = overlap)
    }))
}))

pops <- do.call("c", lapply(seq_len(length(files)), function(i) {
    readRDS(files[i])
}))

## Little overlap
plot(pops[[1]]$result, legend = TRUE)
## much overlap
plot(pops[[32]]$result, legend = TRUE)


## Since the migration happens between two populations we can just
## look at the first two and without squaring because it is more
## difficult to interpret.
df$err12 <- abs(df$ob0 - df$ex0) + abs(df$ob1 - df$ex1)

## The sum of the magration rates
df$migtot <- df$mig01 + df$mig02 + df$mig12

## Looks like the small migration rates result in lower squared error
boxplot(df$err12 ~ cut(df$migtot, 40), xlab = "Migration rate", ylab = "Sum squared error of attribution")

## Looks like the small migration rates result in lower squared error
boxplot(df$err12 ~ cut(df$mig01, 40),
        xlab = "Migration rate", ylab = "",
        ylim = c(0, 0.1))

boxplot(df$err12 ~ cut(df$overlap, 40),
        xlab = "overlap", ylab = "",
        ylim = c(0, 0.1))

plot(df$err12 ~ df$overlap)

## Little overlap
plot(pops[[1]]$result, legend = TRUE)
## much overlap
plot(pops[[32]]$result, legend = TRUE)

hist(df$err12)
i <- which.max(df$err12)
df[i,]
i <- which.min(df$err12)
df[i,]

plot(df$err12 ~ df$mig01)
i <- which.min(df$err12)
points(df$err12[i] ~ df$mig01[i], col = "red", pch = 20)
i <- which.max(df$err12)
points(df$err12[i] ~ df$mig01[i], col = "red", pch = 20)

## Summary of the degree of crossover in the populations

length(pops)
class(pops[i])
plot(pops[[i]]$result, legend = TRUE)

model <- lm(err12 ~ migtot, data = df)
newdata <- data.frame(migtot = seq(0, 3, by = 0.01))
newdata$err12 <- predict(model, newdata = newdata)
plot(df$err12 ~ df$migtot)
lines(newdata$err12 ~ newdata$migtot, col = "red")
summary(model)
## Look at just up to and including those with a sum of the migration
## rates of 0.00003.

## Compare low and high migration
hist(df$migtot)
df$binmig <- df$migtot > 0.002
boxplot(df$err12 ~ df$binmig, xlab = "Migration rate", ylab = "Sum error of attribution")
