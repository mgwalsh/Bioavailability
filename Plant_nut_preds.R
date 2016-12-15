#' Plant nutrient predictions of greenhouse grown wheat plants from soil data
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
  require(randomForest)
  require(gbm)
  require(deepnet)
  require(bartMachine)
  require(glmnet)
})

# Data setup --------------------------------------------------------------
# Run this first: https://github.com/mgwalsh/Bioavailability/blob/master/FAO_micro_setup.R
# or run ...
# SourceURL <- "https://raw.githubusercontent.com/mgwalsh/Bioavailability/master/FAO_micro_setup.R"
# source_url(SourceURL)

# Plant labels
lt <- fao_cal$pZn
lv <- fao_val$pZn

# Soil features
wett <- fao_cal[c(4, 8:9, 33:42)] ## soil wetchem
mirt <- fao_cal[43:1806] # soil MIR
wetv <- fao_val[c(4, 8:9, 33:42)] ## soil wetchem
mirv <- fao_val[43:1806] # soil MIR

# RF models ---------------------------------------------------------------
# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "cv")

# Wet chemistry covariates
tg <- expand.grid(mtry=seq(2, 10, by=1))
wet.rfo <- train(wett, lt,
                 preProc = c("center", "scale"),
                 method = "rf",
                 ntree = 501,
                 tuneGrid = tg,
                 trControl = tc)
print(wet.rfo)
rfo_wet <- predict(wet.rfo, wetv)
rm("wet.rfo")

# MIR features
set.seed(1385321)
tg <- expand.grid(mtry=seq(10, 150, by=10))
mir.rfo <- train(mirt, lt,
                 preProc = c("center", "scale"),
                 method = "rf",
                 ntree = 501,
                 tuneGrid = tg,
                 trControl = tc)
print(mir.rfo)
rfo_mir <- predict(mir.rfo, mirv)
rm("mir.rfo")

stopCluster(mc)

# GBM models --------------------------------------------------------------
# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "repeatedcv", repeats=5)

# Wet chemistry features
tg <- expand.grid(.n.trees=seq(10, 200, by=5), 
                  .interaction.depth = 5,
                  .shrinkage = 0.1,
                  .n.minobsinnode = 10)

wet.gbm <- train(wett, lt, 
                 method = "gbm", 
                 preProc = c("center", "scale"),
                 trControl = tc,
                 tuneGrid = tg)
print(wet.gbm)
gbm_wet <- predict(wet.gbm, wetv)
rm("wet.gbm")

# MIR features
tg <- expand.grid(.n.trees=seq(10, 200, by=10), 
                  .interaction.depth = 10,
                  .shrinkage = 0.1,
                  .n.minobsinnode = 10)

mir.gbm <- train(mirt, lt, 
                 method = "gbm", 
                 preProc = c("center", "scale"),
                 trControl = tc,
                 tuneGrid = tg)
print(mir.gbm)
gbm_mir <- predict(mir.gbm, mirv)
rm("mir.gbm")

stopCluster(mc)

# DNN models --------------------------------------------------------------
# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "cv")

# Wet chemistry features
tg <- expand.grid(layer1 = 2:10,
                  layer2 = 0,
                  layer3 = 0,
                  hidden_dropout = 0,
                  visible_dropout = 0)

wet.dnn <- train(wett, lt, 
                 method = "dnn", 
                 preProc = c("center", "scale"), 
                 trControl = tc,
                 tuneGrid = tg)
print(wet.dnn)
dnn_wet <- predict(wet.dnn, wetv)
rm("wet.dnn")

# MIR features
set.seed(1385321)
tg <- expand.grid(layer1 = seq(5, 20, by=5),
                  layer2 = 0:2,
                  layer3 = 0,
                  hidden_dropout = 0,
                  visible_dropout = 0)

mir.dnn <- train(mirt, lt, 
                 method = "dnn", 
                 preProc = c("center", "scale"), 
                 trControl = tc,
                 tuneGrid = tg)
print(mir.dnn)
dnn_mir <- predict(mir.dnn, mirv)
rm("mir.dnn")

stopCluster(mc)

# bartMachine models ------------------------------------------------------
# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "cv")

# Wet chemistry features
wet.bar <- train(wett, lt,
                 method = "bartMachine", 
                 preProc = c("center", "scale"), 
                 trControl = tc,
                 tuneLength = 2,
                 verbose = FALSE,
                 seed = 1)
print(wet.bar)
bar_wet <- predict(wet.bar, wetv)
rm("wet.bar")

# MIR features
set.seed(1385321)
mir.bar <- train(mirt, lt,
                 method = "bartMachine", 
                 trControl = tc,
                 tuneLength = 2,
                 verbose = FALSE,
                 seed = 1)
print(mir.bar)
bar_mir <- predict(mir.bar, mirv)
rm("mir.bar")

stopCluster(mc)

# Model stacking setup ----------------------------------------------------
pwetv <- as.data.frame(cbind(lv, rfo_wet, gbm_wet, dnn_wet, bar_wet))
names(pwetv) <- c("L", "RFO", "GBM", "DNN", "BART")
pmirv <- as.data.frame(cbind(lv, rfo_mir$H, gbm_mir$H, dnn_mir$H, bar_mir$H))
names(pmirv) <- c("L", "RFO", "GBM", "DNN", "BART")

# Write data files --------------------------------------------------------
write.csv(pwetv, "pwetv.csv", row.names=F)
write.csv(pmirv, "pmirv.csv", row.names=F)

# Remove extraneous objects from memory -----------------------------------
# rm(list=setdiff(ls(), c("pwetv", "pmirv")))

# Model stacking ----------------------------------------------------------
# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "cv", classProbs = TRUE, summaryFunction = twoClassSummary,
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
boxplot(DNN~HL, notch=T, pwetv, ylim=c(0,1))
boxplot(BART~HL, notch=T, pwetv, ylim=c(0,1))

# MIR predictions
boxplot(RFO~HL, notch=T, pmirv, ylim=c(0,1))
boxplot(GBM~HL, notch=T, pmirv, ylim=c(0,1))
boxplot(DNN~HL, notch=T, pmirv, ylim=c(0,1))
boxplot(BART~HL, notch=T, pmirv, ylim=c(0,1))
dev.off()

require(RcolorBrewer)
k <- adjustcolor(brewer.pal(3, "Set1")[ens_pre$HL], alpha=.8)
plot(WET~MIR, ens_pre, pch = 20, col = k, xlim= c(0,1), ylim= c(0,1))
abline(c(0,1))
