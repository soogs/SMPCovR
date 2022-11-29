
library(ggplot2)
library(ggh4x)

load("C:\\Users\\park\\Desktop\\research\\project_discovR\\industry\\simulation study\\simulation 24_nov_2022\\results_combined_24_nov_2022.Rdata")

results <- results[complete.cases(results),]


# we are only using the design factor: number of relevant of covariates

# also, I am removing the design factor: 3 relevant covariates

results <- subset(results, py_pattern == "notall3")


pred_df <- results

relevant_factor <- pred_df$py_pattern == "notall3"
relevant_factor[relevant_factor] <- "Relev 3"
relevant_factor[relevant_factor == "FALSE"] <- "Relev 2"
relevant_factor[pred_df$py_pattern == "notall2"]

relevant_factor2 <- factor(relevant_factor, levels = c("Relev 3"))


VAF_factor <- factor(paste("VAF ", pred_df$signal_level, sep = ""), levels = c("VAF 0.9", "VAF 0.5"))

dimension_factor <- factor(paste(pred_df$dimension, " dim", sep = ""), levels = c("low dim", "high dim"))

Jy_factor <- factor(paste("Jy ",pred_df$Jy, sep = ""), levels = c("Jy 5", "Jy 20"))


# editing the label names #
dimension_factor2 <- factor(dimension_factor, levels = c("low dim", "high dim", "Jk = 100", "Jk = 15"))
dimension_factor2[dimension_factor2 == "low dim"] <- "Jk = 15"
dimension_factor2[dimension_factor2 == "high dim"] <- "Jk = 100"
dimension_factor2 <- factor(dimension_factor2, levels = c("Jk = 15", "Jk = 100"))


pred_df$dimensions <- dimension_factor2
pred_df$signal_level <- VAF_factor
pred_df$Jy <- Jy_factor

pred_df$relevant <- relevant_factor2


disc_pred_index <- which(colnames(pred_df) == "disc_pred")

lasso0_pred_index <- which(colnames(pred_df) == "lasso0_pred")

diacon_pred_index <- which(colnames(pred_df) == "diacon_pred")

diacor_pred_index <- which(colnames(pred_df) == "diacor_pred")

colnames(pred_df)[c(disc_pred_index,
                    lasso0_pred_index,
                    diacon_pred_index)] <- c("SMPCovR", "SPCovR",  "sPLS")


pred_df <- pred_df[,c(which(colnames(pred_df) == "dimensions"), which(colnames(pred_df) == "signal_level"),
                      which(colnames(pred_df) == "Jy"), which(colnames(pred_df) == "relevant"),
                      which(colnames(pred_df) == "outcome"), 
                      disc_pred_index, lasso0_pred_index, diacon_pred_index)]

pred_long <- tidyr::gather(data = pred_df, method, pred,  c("SMPCovR", "SPCovR",  "sPLS")) 

pred_long_method <- factor(pred_long$method, levels = c("SMPCovR", "SPCovR",  "sPLS"))

pred_long$method <- pred_long_method


# plot 2: conditions #
plot2_pred <- ggplot(pred_long, aes(x = method, y = pred, fill = method)) +
  theme_bw() +
  facet_nested(Jy ~ dimensions + signal_level) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, 
                                   margin = margin(0,10,0,0, "pt")),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-12,0,0,0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 7, 
                                  margin = margin(0.05,0,0.05,0, "cm"))) +
  geom_boxplot(width = 0.8, fatten = 0.9, outlier.size = 0.07, lwd = 0.3) 
# + stat_summary(fun.y = mean, geom = "point", size = 0.2, shape = 20, color = "red", position = position_dodge(0.6)) <<MEAN>>


# Prediction quality only oncerning the non-zero outcome variables
pred_nonzero_df <- results


# we are only using the design factor: number of relevant of covariates

relevant_factor <- pred_nonzero_df$py_pattern == "notall3"
relevant_factor[relevant_factor] <- "Relev 3"
relevant_factor[relevant_factor == "FALSE"] <- "Relev 2"
relevant_factor[pred_nonzero_df$py_pattern == "notall2"]

relevant_factor2 <- factor(relevant_factor, levels = c("Relev 3"))


VAF_factor <- factor(paste("VAF ", pred_nonzero_df$signal_level, sep = ""), levels = c("VAF 0.9", "VAF 0.5"))

dimension_factor <- factor(paste(pred_nonzero_df$dimension, " dim", sep = ""), levels = c("low dim", "high dim"))

Jy_factor <- factor(paste("Jy ",pred_nonzero_df$Jy, sep = ""), levels = c("Jy 5", "Jy 20"))


# editing the label names #
dimension_factor2 <- factor(dimension_factor, levels = c("low dim", "high dim", "Jk = 100", "Jk = 15"))
dimension_factor2[dimension_factor2 == "low dim"] <- "Jk = 15"
dimension_factor2[dimension_factor2 == "high dim"] <- "Jk = 100"
dimension_factor2 <- factor(dimension_factor2, levels = c("Jk = 15", "Jk = 100"))


pred_nonzero_df$dimensions <- dimension_factor2
pred_nonzero_df$signal_level <- VAF_factor
pred_nonzero_df$Jy <- Jy_factor

pred_nonzero_df$relevant <- relevant_factor2


disc_pred_nonzero_index <- which(colnames(pred_nonzero_df) == "disc_pred_nonzero")

lasso0_pred_nonzero_index <- which(colnames(pred_nonzero_df) == "lasso0_pred_nonzero")

diacon_pred_nonzero_index <- which(colnames(pred_nonzero_df) == "diacon_pred_nonzero")

diacor_pred_nonzero_index <- which(colnames(pred_nonzero_df) == "diacor_pred_nonzero")

colnames(pred_nonzero_df)[c(disc_pred_nonzero_index,
                    lasso0_pred_nonzero_index,
                    diacon_pred_nonzero_index)] <- c("SMPCovR", "SPCovR",  "sPLS")


pred_nonzero_df <- pred_nonzero_df[,c(which(colnames(pred_nonzero_df) == "dimensions"), which(colnames(pred_nonzero_df) == "signal_level"),
                      which(colnames(pred_nonzero_df) == "Jy"), which(colnames(pred_nonzero_df) == "relevant"),
                      which(colnames(pred_nonzero_df) == "outcome"), 
                      disc_pred_nonzero_index, lasso0_pred_nonzero_index, diacon_pred_nonzero_index)]

pred_nonzero_long <- tidyr::gather(data = pred_nonzero_df, method, pred_nonzero,  c("SMPCovR", "SPCovR",  "sPLS")) 

pred_nonzero_long_method <- factor(pred_nonzero_long$method, levels = c("SMPCovR", "SPCovR",  "sPLS"))

pred_nonzero_long$method <- pred_nonzero_long_method


# plot 2: conditions #
plot2_pred_nonzero <- ggplot(pred_nonzero_long, aes(x = method, y = pred_nonzero, fill = method)) +
  theme_bw() +
  facet_nested(Jy ~ dimensions + signal_level) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, 
                                   margin = margin(0,10,0,0, "pt")),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-12,0,0,0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 7, 
                                  margin = margin(0.05,0,0.05,0, "cm"))) +
  geom_boxplot(width = 0.8, fatten = 0.9, outlier.size = 0.07, lwd = 0.3) 
# + stat_summary(fun.y = mean, geom = "point", size = 0.2, shape = 20, color = "red", position = position_dodge(0.6)) <<MEAN>>



# correct classification of W #
corrects_df <- results

relevant_factor <- corrects_df$py_pattern == "notall3"
relevant_factor[relevant_factor] <- "Relev 3"
relevant_factor[relevant_factor == "FALSE"] <- "Relev 2"
relevant_factor[corrects_df$py_pattern == "notall2"]

relevant_factor2 <- factor(relevant_factor, levels = c("Relev 3"))


VAF_factor <- factor(paste("VAF ", corrects_df$signal_level, sep = ""), levels = c("VAF 0.9", "VAF 0.5"))

dimension_factor <- factor(paste(corrects_df$dimension, " dim", sep = ""), levels = c("low dim", "high dim"))

Jy_factor <- factor(paste("Jy ",corrects_df$Jy, sep = ""), levels = c("Jy 5", "Jy 20"))


# editing the label names #
dimension_factor2 <- factor(dimension_factor, levels = c("low dim", "high dim", "Jk = 100", "Jk = 15"))
dimension_factor2[dimension_factor2 == "low dim"] <- "Jk = 15"
dimension_factor2[dimension_factor2 == "high dim"] <- "Jk = 100"
dimension_factor2 <- factor(dimension_factor2, levels = c("Jk = 15", "Jk = 100"))


corrects_df$dimensions <- dimension_factor2
corrects_df$signal_level <- VAF_factor
corrects_df$Jy <- Jy_factor

corrects_df$relevant <- relevant_factor2


disc_corrects_index <- which(colnames(corrects_df) == "disc_correct")

lasso0_corrects_index <- which(colnames(corrects_df) == "lasso0_correct")

diacon_corrects_index <- which(colnames(corrects_df) == "diacon_correct")

diacor_corrects_index <- which(colnames(corrects_df) == "diacor_correct")

colnames(corrects_df)[c(disc_corrects_index,
                        lasso0_corrects_index,
                        diacon_corrects_index)] <- c("SMPCovR", "SPCovR",  "sPLS")


corrects_df <- corrects_df[,c(which(colnames(corrects_df) == "dimensions"), which(colnames(corrects_df) == "signal_level"),
                              which(colnames(corrects_df) == "Jy"), which(colnames(corrects_df) == "relevant"),
                              which(colnames(corrects_df) == "outcome"), 
                              disc_corrects_index, lasso0_corrects_index, diacon_corrects_index)]

corrects_long <- tidyr::gather(data = corrects_df, method, corrects,  c("SMPCovR", "SPCovR",  "sPLS")) 

corrects_long_method <- factor(corrects_long$method, levels = c("SMPCovR", "SPCovR",  "sPLS"))

corrects_long$method <- corrects_long_method



# plot 2: conditions #
plot2_corrects <- ggplot(corrects_long, aes(x = method, y = corrects, fill = method)) +
  theme_bw() +
  facet_nested(Jy ~ dimensions + signal_level) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, 
                                   margin = margin(0,10,0,0, "pt")),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-12,0,0,0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 7, 
                                  margin = margin(0.05,0,0.05,0, "cm"))) +
  geom_boxplot(width = 0.8, fatten = 0.9, outlier.size = 0.07, lwd = 0.3) 
1# + stat_summary(fun.y = mean, geom = "point", size = 0.2, shape = 20, color = "red", position = position_dodge(0.6)) <<MEAN>>


# corrects Py #
Py_corrects_df <- results

relevant_factor <- Py_corrects_df$py_pattern == "notall3"
relevant_factor[relevant_factor] <- "Relev 3"
relevant_factor[relevant_factor == "FALSE"] <- "Relev 2"
relevant_factor[Py_corrects_df$py_pattern == "notall2"]

relevant_factor2 <- factor(relevant_factor, levels = c("Relev 3"))


VAF_factor <- factor(paste("VAF ", Py_corrects_df$signal_level, sep = ""), levels = c("VAF 0.9", "VAF 0.5"))

dimension_factor <- factor(paste(Py_corrects_df$dimension, " dim", sep = ""), levels = c("low dim", "high dim"))

Jy_factor <- factor(paste("Jy ",Py_corrects_df$Jy, sep = ""), levels = c("Jy 5", "Jy 20"))


# editing the label names #
dimension_factor2 <- factor(dimension_factor, levels = c("low dim", "high dim", "Jk = 100", "Jk = 15"))
dimension_factor2[dimension_factor2 == "low dim"] <- "Jk = 15"
dimension_factor2[dimension_factor2 == "high dim"] <- "Jk = 100"
dimension_factor2 <- factor(dimension_factor2, levels = c("Jk = 15", "Jk = 100"))


Py_corrects_df$dimensions <- dimension_factor2
Py_corrects_df$signal_level <- VAF_factor
Py_corrects_df$Jy <- Jy_factor

Py_corrects_df$relevant <- relevant_factor2


disc_Py_corrects_index <- which(colnames(Py_corrects_df) == "disc_py_correct")

lasso0_Py_corrects_index <- which(colnames(Py_corrects_df) == "lasso0_py_correct")

diacon_Py_corrects_index <- which(colnames(Py_corrects_df) == "diacon_py_correct")

diacor_Py_corrects_index <- which(colnames(Py_corrects_df) == "diacor_py_correct")

colnames(Py_corrects_df)[c(disc_Py_corrects_index,
                           lasso0_Py_corrects_index,
                           diacon_Py_corrects_index)] <- c("SMPCovR", "SPCovR",  "sPLS")


Py_corrects_df <- Py_corrects_df[,c(which(colnames(Py_corrects_df) == "dimensions"), which(colnames(Py_corrects_df) == "signal_level"),
                                    which(colnames(Py_corrects_df) == "Jy"), which(colnames(Py_corrects_df) == "relevant"),
                                    which(colnames(Py_corrects_df) == "outcome"), 
                                    disc_Py_corrects_index, lasso0_Py_corrects_index, diacon_Py_corrects_index)]


Py_corrects_long <- tidyr::gather(data = Py_corrects_df, method, Py, c("SMPCovR", "SPCovR",  "sPLS")) 

Py_corrects_long_method <- factor(Py_corrects_long$method, levels = c("SMPCovR", "SPCovR",  "sPLS"))

Py_corrects_long$method <- Py_corrects_long_method

# plot 2: conditions #
plot2_Py_corrects <- ggplot(Py_corrects_long, aes(x = method, y = Py, fill = method)) +
  theme_bw() +
  facet_nested(Jy ~ dimensions + signal_level) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, 
                                   margin = margin(0,10,0,0, "pt")),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-12,0,0,0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 7, 
                                  margin = margin(0.05,0,0.05,0, "cm"))) +
  geom_boxplot(width = 0.8, fatten = 0.9, outlier.size = 0.07, lwd = 0.3) 
# + stat_summary(fun.y = mean, geom = "point", size = 0.2, shape = 20, color = "red", position = position_dodge(0.6)) <<MEAN>>




# zerocoefs for SPCovR and sPLS#
zerocoefs_df <- results

relevant_factor <- zerocoefs_df$py_pattern == "notall3"
relevant_factor[relevant_factor] <- "Relev 3"
relevant_factor[relevant_factor == "FALSE"] <- "Relev 2"
relevant_factor[zerocoefs_df$py_pattern == "notall2"]

relevant_factor2 <- factor(relevant_factor, levels = c("Relev 3"))


VAF_factor <- factor(paste("VAF ", zerocoefs_df$signal_level, sep = ""), levels = c("VAF 0.9", "VAF 0.5"))

dimension_factor <- factor(paste(zerocoefs_df$dimension, " dim", sep = ""), levels = c("low dim", "high dim"))

Jy_factor <- factor(paste("Jy ",zerocoefs_df$Jy, sep = ""), levels = c("Jy 5", "Jy 20"))


# editing the label names #
dimension_factor2 <- factor(dimension_factor, levels = c("low dim", "high dim", "Jk = 100", "Jk = 15"))
dimension_factor2[dimension_factor2 == "low dim"] <- "Jk = 15"
dimension_factor2[dimension_factor2 == "high dim"] <- "Jk = 100"
dimension_factor2 <- factor(dimension_factor2, levels = c("Jk = 15", "Jk = 100"))


zerocoefs_df$dimensions <- dimension_factor2
zerocoefs_df$signal_level <- VAF_factor
zerocoefs_df$Jy <- Jy_factor

zerocoefs_df$relevant <- relevant_factor2


lasso0_zerocoefs_index <- which(colnames(zerocoefs_df) == "lasso0_zerocoefs")

diacon_zerocoefs_index <- which(colnames(zerocoefs_df) == "diacon_zerocoefs")

colnames(zerocoefs_df)[c(lasso0_zerocoefs_index,
                           diacon_zerocoefs_index)] <- c("SPCovR", "sPLS")


zerocoefs_df <- zerocoefs_df[,c(which(colnames(zerocoefs_df) == "dimensions"), which(colnames(zerocoefs_df) == "signal_level"),
                                    which(colnames(zerocoefs_df) == "Jy"), which(colnames(zerocoefs_df) == "relevant"),
                                    which(colnames(zerocoefs_df) == "outcome"), 
                                    lasso0_zerocoefs_index, diacon_zerocoefs_index)]


zerocoefs_long <- tidyr::gather(data = zerocoefs_df, method, zerocoefs, c("SPCovR", "sPLS")) 

zerocoefs_long_method <- factor(zerocoefs_long$method, levels = c("SPCovR", "sPLS"))

zerocoefs_long$method <- zerocoefs_long_method

# plot 2: conditions #
plot2_zerocoefs <- ggplot(zerocoefs_long, aes(x = method, y = zerocoefs, fill = method)) +
  theme_bw() +
  facet_nested(Jy ~ dimensions + signal_level) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, 
                                   margin = margin(0,10,0,0, "pt")),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-12,0,0,0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 7, 
                                  margin = margin(0.05,0,0.05,0, "cm"))) +
  geom_boxplot(width = 0.8, fatten = 0.9, outlier.size = 0.07, lwd = 0.3) 
# + stat_summary(fun.y = mean, geom = "point", size = 0.2, shape = 20, color = "red", position = position_dodge(0.6)) <<MEAN>>




# correct inclusion of outcome variables #
Y_correct_df <- results

relevant_factor <- Y_correct_df$py_pattern == "notall3"
relevant_factor[relevant_factor] <- "Relev 3"
relevant_factor[relevant_factor == "FALSE"] <- "Relev 2"
relevant_factor[Y_correct_df$py_pattern == "notall2"]

relevant_factor2 <- factor(relevant_factor, levels = c("Relev 3"))


VAF_factor <- factor(paste("VAF ", Y_correct_df$signal_level, sep = ""), levels = c("VAF 0.9", "VAF 0.5"))

dimension_factor <- factor(paste(Y_correct_df$dimension, " dim", sep = ""), levels = c("low dim", "high dim"))

Jy_factor <- factor(paste("Jy ",Y_correct_df$Jy, sep = ""), levels = c("Jy 5", "Jy 20"))


# editing the label names #
dimension_factor2 <- factor(dimension_factor, levels = c("low dim", "high dim", "Jk = 100", "Jk = 15"))
dimension_factor2[dimension_factor2 == "low dim"] <- "Jk = 15"
dimension_factor2[dimension_factor2 == "high dim"] <- "Jk = 100"
dimension_factor2 <- factor(dimension_factor2, levels = c("Jk = 15", "Jk = 100"))


Y_correct_df$dimensions <- dimension_factor2
Y_correct_df$signal_level <- VAF_factor
Y_correct_df$Jy <- Jy_factor

Y_correct_df$relevant <- relevant_factor2


out_Y_correct_index <- which(colnames(Y_correct_df) == "disc_Y_correct")

CV_Y_correct_index <- which(colnames(Y_correct_df) == "disc_CV_Y_correct")

colnames(Y_correct_df)[c(out_Y_correct_index,
                         CV_Y_correct_index)] <- c("out-sample", "CV folds")


Y_correct_df <- Y_correct_df[,c(which(colnames(Y_correct_df) == "dimensions"), which(colnames(Y_correct_df) == "signal_level"),
                                which(colnames(Y_correct_df) == "Jy"), which(colnames(Y_correct_df) == "relevant"),
                                which(colnames(Y_correct_df) == "outcome"), 
                                out_Y_correct_index, CV_Y_correct_index)]


Y_correct_long <- tidyr::gather(data = Y_correct_df, method, Y_correct, c("out-sample", "CV folds")) 

Y_correct_long_method <- factor(Y_correct_long$method, levels = c("out-sample", "CV folds"))

Y_correct_long$method <- Y_correct_long_method

# only the out-sample results
Y_correct_long <- subset(Y_correct_long, method == "out-sample")

# plot 2: conditions #
plot2_Y_correct <- ggplot(Y_correct_long, aes(x = method, y = Y_correct, fill = method)) +
  theme_bw() +
  facet_nested(Jy ~ dimensions + signal_level) +
  theme(panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 7, 
                                   margin = margin(0,10,0,0, "pt")),
        legend.margin=margin(0,0,0,0),
        legend.box.margin=margin(-12,0,0,0),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 7, 
                                  margin = margin(0.05,0,0.05,0, "cm"))) +
  geom_boxplot(width = 0.8, fatten = 0.9, outlier.size = 0.07, lwd = 0.3) 
# + stat_summary(fun.y = mean, geom = "point", size = 0.2, shape = 20, color = "red", position = position_dodge(0.6)) <<MEAN>>

