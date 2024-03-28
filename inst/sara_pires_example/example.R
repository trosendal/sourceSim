library(sourceSim)

migration <- matrix(c(0,0,0,0,0,0,0,0,0), nrow = 3)
df <- sourceSim::simu(input = list(NPOP = 3, MUTR = 4e-6),
                      migration = migration)

## 1. Bug occurs when the source populations have only "others" and
## the human population has 0 in at least 1 source population for the
## other category. We get a dimesion problem in the result data where
## the hald model run in Bugs seems to report back the wrong number of
## dimensions. When there are 0's in the input data after tabulation
## by source and serovar you get a dimension problem in the structure
## of the output data. We need to look at model_output$sims.array and
## check which indices of lambdaij are present in the result and fill
## in those absent with 0's. Then if we can fix the
## source_output$sims.list then we can do the further calculations the
## same way are currently.
##
ob <- sample_humans(df, c(1, 0, 0))
debugonce(sourceSim:::hald.list)
sourceSim::hald(ob)

## 2. Another bug occurs
ob <- sample_humans(df, c(1, 0.5, 0.5))
debugonce(sourceSim:::hald.list)
debugonce(sourceSim:::hald.sourceSim_result)
sourceSim::hald(ob)
