library(sourceSim)

migration <- matrix(c(0,0,0,0,0,0,0,0,0), nrow = 3)
df <- sourceSim::simu(input = list(NPOP = 3, MUTR = 4e-6),
                      migration = migration)

ob <- sample_humans(df, c(1, 0, 0))
saveRDS(ob$population, file = "sample.RDS")

sourceSim::hald(ob)
