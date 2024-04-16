library(ggplot2)

attr_err <- function(ob, ex) {
    abs(ob - ex)
}

experiment <- "experiment7"

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
        attr_ob = x$attribution_island[1, ],
        rmse = x$rmse_island
    )

    attr_hald <- data.frame(
        pop = colnames(x$attribution_hald),
        model = "hald",
        mig = unique(as.numeric(x$migration[x$migration > 0])),
        ov = x$overlap,
        attr_exp = as.numeric(x$sampling),
        attr_ob = x$attribution_hald[1, ],
        rmse = x$rmse_hald
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

model_labs <- c("Hald model", "Asymmetric island model")
names(model_labs) <- c("hald", "island")
pop_labs <- paste("Population", c("A", "B", "C"))
names(pop_labs) <- c("Pop_0", "Pop_1", "Pop_2")

p_mig <- ggplot(res_df, aes(x = migration, y = rmse, group = migration)) +
    geom_boxplot(fill = "gray") +
    facet_grid(model ~ pop, labeller = labeller(model = model_labs, pop = pop_labs)) +
    scale_x_discrete(
        limits = unique(migration_bins),
        breaks = migration_bins[seq(1, length(migration_bins), by = 6)]
    ) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

ggsave(file.path("inst/plots/experiments", paste0(experiment, "_boxplot_mig.png")), p_mig, width = 15)

p_ov <- ggplot(res_df, aes(x = overlap, y = rmse, group = overlap)) +
    geom_boxplot(fill = "gray") +
    facet_grid(model ~ pop, labeller = labeller(model = model_labs, pop = pop_labs)) +
    scale_x_discrete(
        limits = unique(overlap_bins),
        breaks = overlap_bins[seq(1, length(overlap_bins), by = 6)]
    ) +
    labs(x = "Overlap", y = "Attribution RMSE") +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

ggsave(file.path("inst/plots/experiments", paste0(experiment, "_boxplot_ov.png")), p_ov, width = 15)
