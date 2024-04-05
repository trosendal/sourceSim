library(ggplot2)

attr_err <- function(ob, ex) {
    abs(ob - ex)
}

experiment <- "experiment5"

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
        mig = unique(as.numeric(x$migration[x$migration > 0])),
        ov = x$overlap,
        attr_exp = as.numeric(x$sampling),
        attr_ob = x$attribution_island[1, ]
    )

    attr_hald <- data.frame(
        pop = colnames(x$attribution_hald),
        model = "hald",
        mig = unique(as.numeric(x$migration[x$migration > 0])),
        ov = x$overlap,
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

res_df$migration <- cut(res_df$mig, 20, dig.lab = 3)
migration_bins <- levels(res_df$migration)

res_df$overlap <- cut(res_df$ov, 20, dig.lab = 3)
overlap_bins <- levels(res_df$overlap)

p_mig <- ggplot(res_df, aes(x = migration, y = error, group = migration)) +
    geom_boxplot(fill = "gray") +
    facet_grid(model ~ pop) +
    ggtitle("Error vs migration (migp = 0.01, migr = (0, 0.05))") +
    scale_x_discrete(
        limits = unique(migration_bins),
        breaks = migration_bins[seq(1, length(migration_bins), by = 6)]
    ) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

ggsave(file.path("inst/plots/experiments", paste0(experiment, "_boxplot_mig.pdf")), p_mig, width = 15)

p_ov <- ggplot(res_df, aes(x = overlap, y = error, group = overlap)) +
    geom_boxplot(fill = "gray") +
    facet_grid(model ~ pop) +
    ggtitle("Error vs overlap (migp = 0.01, migr = (0, 0.05))") +
    scale_x_discrete(
        limits = unique(overlap_bins),
        breaks = overlap_bins[seq(1, length(overlap_bins), by = 6)]
    ) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

ggsave(file.path("inst/plots/experiments", paste0(experiment, "_boxplot_ov.pdf")), p_ov, width = 15)
