library(ggplot2)
library(reshape2)

# attr_err <- function(ob, ex) {
#     stopifnot(identical(length(ob), length(ex)))

#     m_ob <- prod(ob) ^ (1 / length(ob))
#     m_ex <- prod(ex) ^ (1 / length(ex))

#     sqrt(
#         (1 / length(ob)) * sum(mapply(function(o, e, m_o, m_e) {
#             ((log(o) / m_ob) - (log(e) / m_ex))^2
#         }, o = ob, e = ex, MoreArgs = list(m_o = m_ob, m_e = m_ex)))
#     )
# }

attr_err <- function(ob, ex) {
    stopifnot(identical(length(ob), length(ex)))

    sum(ob * log(ob / ex))
}

experiment <- "experiment7"

experiment_dir <- file.path("inst/results", experiment)

stopifnot(dir.exists(experiment_dir))

result_files <- list.files(experiment_dir, full.names = TRUE)
result_files <- result_files[tools::file_ext(result_files) == "Rds"]

stopifnot(length(result_files) > 0)

results <- lapply(result_files, readRDS)

res_hald <- do.call("rbind", lapply(seq_len(length(results)), function(i) {
    data.frame(
        run_id = i,
        model = "hald",
        overlap = results[[i]]$overlap,
        err = attr_err(
            ob = results[[i]]$attribution_hald[1, ],
            ex = as.numeric(results[[i]]$sampling)
        )
    )
}))

res_island <- do.call("rbind", lapply(seq_len(length(results)), function(i) {
    data.frame(
        run_id = i,
        model = "island",
        overlap = results[[i]]$overlap,
        err = attr_err(
            ob = results[[i]]$attribution_island[1, ],
            ex = as.numeric(results[[i]]$sampling)
        )
    )
}))

res_df <- rbind(res_hald, res_island)

rownames(res_df) <- NULL

library(data.table)

setDT(res_df)

res_df <- res_df[order(run_id, model)]

res_df$overlap <- cut(res_df$overlap, 20, dig.lab = 3)
overlap_bins <- levels(res_df$overlap)

p <- ggplot(res_df, aes(x = overlap, y = err, group = overlap)) +
    geom_boxplot(fill = "gray") +
    facet_grid(rows = vars(model), labeller = labeller(pop = pop_labs)) +
    scale_x_discrete(
        limits = unique(overlap_bins),
        breaks = overlap_bins[seq(1, length(overlap_bins), by = 6)]
    ) +
    theme(axis.text.x = element_text(angle = 45, vjust = 0.5))

ggsave(file.path("inst/plots/experiments", paste0(experiment, "_boxplot_divergence.png")), p, width = 15)

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
