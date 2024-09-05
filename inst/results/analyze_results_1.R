library(ggplot2)

experiment <- "experiment1"

experiment_dir <- file.path("inst/results", experiment)

stopifnot(dir.exists(experiment_dir))

result_files <- list.files(experiment_dir, full.names = TRUE)
result_files <- result_files[tools::file_ext(result_files) == "Rds"]

stopifnot(length(result_files) > 0)

results <- lapply(result_files, readRDS)

res_df <- do.call("rbind", lapply(results, function(x) {
    df <- data.frame(
        mig = x$migration,
        overlap = x$overlap
    )

    rownames(df) <- NULL

    df
}))


res_df$migration <- cut(res_df$mig, 20, dig.lab = 3)
migration_bins <- levels(res_df$migration)

# res_df$migration <- cut(res_df$mig, 20, dig.lab = 3)
# migration_bins <- levels(res_df$migration)

p <- ggplot(res_df, aes(x = mig, y = overlap)) +
    geom_point() +
    ggtitle("Migration vs overlap (migp = 0)")

ggsave(file.path(
    "inst/plots/experiments", paste0(experiment, "_boxplot.pdf")
), p, width = 15)
