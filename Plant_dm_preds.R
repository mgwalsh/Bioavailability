#' 50% biomass end-point predictions of greenhouse grown wheat plants from soil data
#' Soil and wheat plant wet chemistry data courtesy of FAO (doc @ https://www.dropbox.com/s/gwk07tanhu86tqj/Silanpaa%20Report.pdf?dl=0)
#' MIR soil data courtesy of ICRAF (2016)
#' M. Walsh, October 2016

# Required packages
# install.packages(c("devtools","caret","plyr","doParallel" ...)), dependencies=TRUE)
suppressPackageStartupMessages({
  require(devtools)
  require(caret)
  require(plyr)
  require(doParallel)
})

# Data setup --------------------------------------------------------------
# Run this first: https://github.com/mgwalsh/Bioavailability/blob/master/FAO_micro_setup.R
# or run ...
# SourceURL <- "https://raw.githubusercontent.com/mgwalsh/Bioavailability/master/FAO_micro_setup.R"
# source_url(SourceURL)

# Labels
HLt <- fao_cal$HL
HLv <- fao_val$HL

# Features
wett <- fao_cal[c(4, 7:9, 12:24)] ## soil wetchem
wetv <- fao_val[c(4, 7:9, 12:24)]
mirt <- fao_cal[43:1806] ## soil MIR
mirv <- fao_val[43:1806]

# bartMachine models ------------------------------------------------------
options(java.parameters = "-Xmx8000m")
library(bartMachine)

# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = TRUE,
                   summaryFunction = twoClassSummary, allowParallel = T)

# Wet chemistry features
y <- relevel(HLt, "L")
HL_wet.bar <- train(wett, y,
                    method = "bartMachine", 
                    preProc = c("center", "scale"), 
                    trControl = tc,
                    tuneLength = 3,
                    verbose = FALSE,
                    metric = "ROC",
                    seed = 1)
print(HL_wet.bar)
bar_wet <- predict(HL_wet.bar, wetv, type = "prob")
rm("HL_wet.bar")

# MIR features
set.seed(1385321)
HL_mir.bar <- train(mirt, y,
                    method = "bartMachine", 
                    preProc = c("pca"), 
                    trControl = tc,
                    tuneLength = 2,
                    verbose = FALSE,
                    metric = "ROC"
                    seed = 1)
print(HL_mir.bar)
bar_mir <- predict(HL_mir.bar, mirv, type = "prob")
rm("HL_mir.bar")

stopCluster(mc)
detach("package:bartMachine", unload=TRUE)

# RF models ---------------------------------------------------------------
require(randomForest)

# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = TRUE,
                   summaryFunction = twoClassSummary, allowParallel = T)

# Wet chemistry covariates
tg <- expand.grid(mtry=seq(2, 10, by=1))
HL_wet.rfo <- train(wett, HLt,
                    preProc = c("center","scale"),
                    method = "rf",
                    ntree = 501,
                    metric = "ROC",
                    tuneGrid = tg,
                    trControl = tc)
print(HL_wet.rfo)
rfo_wet <- predict(HL_wet.rfo, wetv, type = "prob")
rm("HL_wet.rfo")

# MIR features
set.seed(1385321)
tg <- expand.grid(mtry=seq(5, 30, by=5))
HL_mir.rfo <- train(mirt, HLt,
                    preProc = c("pca"),
                    method = "rf",
                    ntree = 501,
                    metric = "ROC",
                    tuneGrid = tg,
                    trControl = tc)
print(HL_mir.rfo)
rfo_mir <- predict(HL_mir.rfo, mirv, type = "prob")
rm("HL_mir.rfo")

stopCluster(mc)

# GBM models --------------------------------------------------------------
require(gbm)

# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "repeatedcv", repeats=5, classProbs = TRUE, summaryFunction = twoClassSummary,
                   allowParallel = T)

# Wet chemistry features
tg <- expand.grid(.n.trees=seq(10, 100, by=5), 
                  .interaction.depth = 5,
                  .shrinkage = 0.1,
                  .n.minobsinnode = 10)

HL_wet.gbm <- train(wett, HLt, 
                    method = "gbm", 
                    preProc = c("center", "scale"),
                    trControl = tc,
                    metric = "ROC",
                    tuneGrid = tg)
print(HL_wet.gbm)
gbm_wet <- predict(HL_wet.gbm, wetv, type = "prob")
rm("HL_wet.gbm")

# MIR features
tg <- expand.grid(.n.trees=seq(10, 100, by=10), 
                  .interaction.depth = 10,
                  .shrinkage = 0.1,
                  .n.minobsinnode = 10)

HL_mir.gbm <- train(mirt, HLt, 
                    method = "gbm", 
                    preProc = c("pca"),
                    trControl = tc,
                    tuneGrid = tg)
print(HL_mir.gbm)
gbm_mir <- predict(HL_mir.gbm, mirv, type = "prob")
rm("HL_mir.gbm")

stopCluster(mc)

# KNN models ------------------------------------------------------------
# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = TRUE,
                   summaryFunction = twoClassSummary, allowParallel = T)

# Wet chemistry features
HL_wet.knn <- train(wett, HLt, 
                    method = "knn",
                    tuneLength = 11,
                    preProc = c("center","scale"), 
                    trControl = tc,
                    metric ="ROC")
print(HL_wet.knn)
knn_wet <- predict(HL_wet.knn, wetv, type = "prob")
rm("HL_wet.knn")

# MIR features
HL_mir.knn <- train(mirt, HLt, 
                    method = "knn",
                    tuneLength = 11,
                    preProc = c("pca"), 
                    trControl = tc,
                    metric = "ROC")
print(HL_mir.knn)
knn_mir <- predict(HL_mir.knn, mirv, type = "prob")
rm("HL_mir.knn")

stopCluster(mc)

# Model stacking setup ----------------------------------------------------
pwetv <- as.data.frame(cbind(HLv, rfo_wet$H, gbm_wet$H, knn_wet$H, bar_wet$L))
names(pwetv) <- c("HL", "RFO", "GBM", "KNN", "BART")
pwetv$HL <- as.factor(ifelse(pwetv$HL == 1, "H", "L"))
pmirv <- as.data.frame(cbind(HLv, rfo_mir$H, gbm_mir$H, knn_mir$H, bar_mir$L))
names(pmirv) <- c("HL", "RFO", "GBM", "KNN", "BART")
pmirv$HL <- as.factor(ifelse(pmirv$HL == 1, "H", "L"))

# Write data files --------------------------------------------------------
write.csv(pwetv, "pwetv.csv", row.names=F)
write.csv(pmirv, "pmirv.csv", row.names=F)

# Remove extraneous objects from memory -----------------------------------
# rm(list=setdiff(ls(), c("pwetv", "pmirv")))

# Model stacking ----------------------------------------------------------
require(glmnet)

# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "repeatedcv", repeats = 5, classProbs = TRUE, summaryFunction = twoClassSummary,
                   allowParallel = T)

# Wet chemistry model stack
wet.ens <- train(HL ~ ., data = pwetv,
                 method = "glmnet",
                 family = "binomial",
                 metric = "ROC",
                 trControl = tc)
print(wet.ens)
wet.imp <- varImp(wet.ens)
plot(wet.imp, top=4, col="black", cex=1.2, xlab="Model importance in ensemble prediction")
ens_wet <- predict(wet.ens, pwetv, type = "prob")

# MIR model stack
set.seed(1385321)
mir.ens <- train(HL ~ ., data = pmirv,
                 method = "glmnet",
                 family = "binomial",
                 metric = "ROC",
                 trControl = tc)
print(mir.ens)
mir.imp <- varImp(mir.ens)
plot(mir.imp, top=4, col="black", cex=1.2, xlab="Model importance in ensemble prediction")
ens_mir <- predict(mir.ens, pmirv, type = "prob")

stopCluster(mc)

# Receiver operator characteristics ---------------------------------------
par(mfrow=c(2,2), mar=c(2.5,2.5,1,1))

# Wetchem predictions
boxplot(RFO~HL, notch=T, pwetv, ylim=c(0,1))
boxplot(GBM~HL, notch=T, pwetv, ylim=c(0,1))
boxplot(KNN~HL, notch=T, pwetv, ylim=c(0,1))
boxplot(BART~HL, notch=T, pwetv, ylim=c(0,1))

# MIR predictions
boxplot(RFO~HL, notch=T, pmirv, ylim=c(0,1))
boxplot(GBM~HL, notch=T, pmirv, ylim=c(0,1))
boxplot(KNN~HL, notch=T, pmirv, ylim=c(0,1))
boxplot(BART~HL, notch=T, pmirv, ylim=c(0,1))
dev.off()

# ROC curves
require(pROC)
ens_pre <- as.data.frame(cbind(HLv, ens_wet$H, ens_mir$H))
names(ens_pre) <- c("HL", "WET", "MIR")
ens_pre$HL <- as.factor(ifelse(ens_pre$HL == 1, "H", "L"))
ROCw <- roc(ens_pre$HL, ens_pre$WET)
plot(ROCw, xlim=c(1,0), ylim=c(0,1))
ROCm <- roc(ens_pre$HL, ens_pre$MIR)
plot(ROCm, xlim=c(1,0), ylim=c(0,1))

require(RcolorBrewer)
k <- adjustcolor(brewer.pal(3, "Set1")[ens_pre$HL], alpha=.8)
plot(WET~MIR, ens_pre, pch = 20, col = k, xlim= c(0,1), ylim= c(0,1))
abline(c(0,1))
