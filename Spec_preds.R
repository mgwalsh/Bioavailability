# FAO data MIR spectral predictions
# Soil and wheat plant wet chemistry data courtesy of FAO (doc @ https://www.dropbox.com/s/gwk07tanhu86tqj/Silanpaa%20Report.pdf?dl=0)
# MIR soil data courtesy of ICRAF (2016)
# M. Walsh, March 2019

# Data setup --------------------------------------------------------------
# Run this first: https://github.com/mgwalsh/Bioavailability/blob/master/FAO_micro_setup.R
# or run ...
# SourceURL <- "https://raw.githubusercontent.com/mgwalsh/Bioavailability/master/FAO_micro_setup.R"
# source_url(SourceURL)

# Labels ... insert the relevant label
# str(fao_cal) ## check potential labels
lt <- log(fao_cal$Fe+1) ## variables prefaced by "p" are potential plant labels
lv <- log(fao_val$Fe+1) ## ensure that validation and training labels are the same

# Soil spectral features
mirt <- fao_cal[33:1796] # soil MIR features
mirv <- fao_val[33:1796] # ensure that validation features are the same

# RF models ---------------------------------------------------------------
library(doParallel)
library(randomForest)

# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "cv", allowParallel = T)

# Tuning parameters
tg <- expand.grid(mtry=seq(10, 150, by=10))

# Fit model
mir.rfo <- train(mirt, lt,
                 preProc = c("center", "scale"),
                 method = "rf",
                 ntree = 501,
                 tuneGrid = tg,
                 trControl = tc)
print(mir.rfo)
rf_mirt <- predict(mir.rfo, mirt) ## predict training set
rf_mirv <- predict(mir.rfo, mirv) ## predict validation set

stopCluster(mc)
detach("package:randomForest", unload=TRUE)

# GBM models --------------------------------------------------------------
library(plyr)
library(gbm)

# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "repeatedcv", repeats=5, allowParallel = T)

# Tuning parameters
tg <- expand.grid(.n.trees=seq(10, 200, by=10), 
                  .interaction.depth = 10,
                  .shrinkage = 0.1,
                  .n.minobsinnode = 10)

# Fit model
mir.gbm <- train(mirt, lt, 
                 method = "gbm",
                 trControl = tc,
                 tuneGrid = tg)
print(mir.gbm)
gb_mirt <- predict(mir.gbm, mirt) ## predict training set
gb_mirv <- predict(mir.gbm, mirv) ## predict validation set

stopCluster(mc)
detach("package:gbm", unload=TRUE)

# PLS models --------------------------------------------------------------
library(pls)

# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "repeatedcv", repeats = 5, allowParallel = TRUE)

# Fit models
mir.pls <- train(mirt, lt,
                 preProc = c("center", "scale"),
                 method = "pls",
                 tuneGrid = expand.grid(ncomp=seq(2, 20, by=1)),
                 trControl = tc)
print(mir.pls)
pls_mirt <- predict(mir.pls, mirt) ## predict validation set
pls_mirv <- predict(mir.pls, mirv) ## predict validation set

stopCluster(mc)
detach("package:pls", unload=TRUE)

# bartMachine models ------------------------------------------------------
options(java.parameters = "-Xmx8000m")
library(bartMachine)

# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
tc <- trainControl(method = "cv", returnResamp = "all", allowParallel = T)

# Fit model
mir.bar <- train(mirt, lt,
                 method = "bartMachine", 
                 preProc = c("center", "scale"),
                 trControl = tc,
                 tuneLength = 2,
                 seed = 123)
print(mir.bar)
bar_mirt <- predict(mir.bar, mirt) ## predict training set
bar_mirv <- predict(mir.bar, mirv) ## predict validation set

stopCluster(mc)
detach("package:bartMachine", unload=TRUE)

# Model stacking setup ----------------------------------------------------
pmirt <- as.data.frame(cbind(lt, rf_mirt, gb_mirt, pls_mirt, bar_mirt))
names(pmirt) <- c("L", "RFOm", "GBMm", "PLSm", "BARTm")
pmirv <- as.data.frame(cbind(lv, rf_mirv, gb_mirv, pls_mirv, bar_mirv))
names(pmirv) <- c("L", "RFOm", "GBMm", "PLSm", "BARTm")

# Remove extraneous objects from memory -----------------------------------
# rm(list=setdiff(ls(), pmirv"))

# Model stacking ----------------------------------------------------------
library(MASS)

# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "cv", allowParallel = T)

# MIR model stack
set.seed(1385321)
mir.ens <- train(L ~ ., data = pmirv,
                 method = "glm",
                 family = "gaussian",
                 trControl = tc)

print(mir.ens)
summary(mir.ens)

# Predictions -------------------------------------------------------------
ens_mirt <- as.data.frame(predict(mir.ens, pmirt))
names(ens_mirt) <- c("ENSm")
pmirt <- cbind(pmirt, ens_mirt)
ens_mirv <- as.data.frame(predict(mir.ens, pmirv))
names(ens_mirv) <- c("ENSm")
pmirv <- cbind(pmirv, ens_mirv)
pmira <- rbind(pmirt, pmirv) ## combined predictions based on ensemble model

stopCluster(mc)

# Write data files --------------------------------------------------------
write.csv(pmirt, "Fe_pmirt.csv", row.names=F)
write.csv(pmirv, "Fe_pmirv.csv", row.names=F)
write.csv(pmira, "Fe_pmira.csv", row.names=F)

# Plot ensemble predictions -----------------------------------------------

