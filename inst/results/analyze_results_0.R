attr_err <- function(ob, ex) {
    abs(ob - ex)
}

experiment <- "experiment0"

experiment_dir <- file.path("inst/results", experiment)

stopifnot(dir.exists(experiment_dir))

result_files <- list.files(experiment_dir, full.names = TRUE)
result_files <- result_files[tools::file_ext(result_files) == "Rds"]

stopifnot(length(result_files) > 0)

results <- lapply(result_files, readRDS)

res_df <- do.call("rbind", lapply(results, function(x) {
    attr_island <- data.frame(
        pop = colnames(x$attribution_island),
        model = "island",
        migration = unique(as.numeric(x$migration[x$migration > 0])),
        attr_exp = as.numeric(x$sampling),
        attr_ob = x$attribution_island[1, ]
    )

    attr_hald <- data.frame(
        pop = colnames(x$attribution_hald),
        model = "hald",
        migration = unique(as.numeric(x$migration[x$migration > 0])),
        attr_exp = as.numeric(x$sampling),
        attr_ob = x$attribution_hald[1, ]
    )

    attr <- rbind(attr_island, attr_hald)

    rownames(attr) <- NULL

    attr
}))

res_df$error <- mapply(
    attr_err, ob = res_df$attr_ob, ex = res_df$attr_exp
)

res_df$migration <- cut(res_df$migration, 20)

library(ggplot2)

p <- ggplot(res_df, aes(x = migration, y = error, group = migration)) +
    geom_boxplot(fill = "gray") +
    facet_grid(model ~ pop)

ggsave(file.path(experiment_dir, "boxplot.png"), p, width = 15)
