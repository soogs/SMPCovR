# simulation 24 nov 2022 - compiling the results #
# soogeun park #

setwd( "C:\\Users\\park\\Desktop\\research\\project_discovR\\industry\\simulation study\\simulation 24_nov_2022\\simulation 24_nov_2022\\run_and_results\\")

load("./results1_24_nov_2022.Rdata")

results1 <- results

load("./results2_24_nov_2022.Rdata")

results2 <- results

load("./results3_24_nov_2022.Rdata")
  
results3 <- results

load("./results4_24_nov_2022.Rdata")

results4 <- results

load("./results5_24_nov_2022.Rdata")

results5 <- results

load("./results6_24_nov_2022.Rdata")

results6 <- results

load("./results7_24_nov_2022.Rdata")

results7 <- results

load("./results8_24_nov_2022.Rdata")

results8 <- results

load("./results9_24_nov_2022.Rdata")

results9 <- results


# compiling altogether #
results <- results1

results[complete.cases(results2),] <- results2[complete.cases(results2),]
results[complete.cases(results3),] <- results3[complete.cases(results3),]

results[complete.cases(results4),] <- results4[complete.cases(results4),]
results[complete.cases(results5),] <- results5[complete.cases(results5),]
results[complete.cases(results6),] <- results6[complete.cases(results6),]

results[complete.cases(results7),] <- results7[complete.cases(results7),]
results[complete.cases(results8),] <- results8[complete.cases(results8),]

results[complete.cases(results9),] <- results9[complete.cases(results9),]


anyNA(results)

sum(complete.cases(results1))
sum(complete.cases(results2))
sum(complete.cases(results3))
sum(complete.cases(results4))
sum(complete.cases(results5))
sum(complete.cases(results6))
sum(complete.cases(results7))
sum(complete.cases(results8))
sum(complete.cases(results9))

sum(complete.cases(results)) # we missed a couple of replicates for results11

boxplot(results$disc_pred, as.numeric(results$lasso0_pred), results$diacon_pred)

boxplot(results$disc_pred0, as.numeric(results$lasso0_pred0), results$diacon_pred0)

boxplot(results$disc_pred_nonzero, as.numeric(results$lasso0_pred_nonzero), results$diacon_pred_nonzero)

boxplot(results$disc_nonzero, as.numeric(results$lasso0_nonzero), results$diacon_nonzero)

boxplot(results$lasso0_zerocoefs_firstfive, results$diacon_zerocoefs_firstfive)

boxplot(results$lasso0_zerocoefs, results$diacon_zerocoefs)

# save(results, file = "../../results_combined_24_nov_2022.Rdata")
