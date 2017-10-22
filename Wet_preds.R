# FAO data soil wet chemistry plant nutrient uptake predictions
# Soil and wheat plant wet chemistry data courtesy of FAO (doc @ https://www.dropbox.com/s/gwk07tanhu86tqj/Silanpaa%20Report.pdf?dl=0)
# MIR soil data courtesy of ICRAF (2016)
# M. Walsh, February 2017

# Data setup --------------------------------------------------------------
# Run this first: https://github.com/mgwalsh/Bioavailability/blob/master/FAO_micro_setup.R
# or run ...
# SourceURL <- "https://raw.githubusercontent.com/mgwalsh/Bioavailability/master/FAO_micro_setup.R"
# source_url(SourceURL)
rm(mirdat)

# Labels ... insert the relevant label
# str(fao_cal) ## check potential labels
lt <- fao_cal$uZn ## variables prefaced by "u" are uptake labels
lv <- fao_val$uZn ## ensure that validation and training labels are the same

# Soil spectral features
wett <- fao_cal[c(4:24)] ## soil wet chem features
wetv <- fao_val[c(4:24)] ## ensure that validation features are the same

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
tg <- expand.grid(mtry=seq(1, 15, by=1))

# Fit model
wet.rfo <- train(wett, lt,
                 preProc = c("center", "scale"),
                 method = "rf",
                 ntree = 501,
                 tuneGrid = tg,
                 trControl = tc)
print(wet.rfo)
rfo_wet <- predict(wet.rfo, wetv) ## predict validation set

stopCluster(mc)

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
wet.gbm <- train(wett, lt, 
                 method = "gbm",
                 trControl = tc,
                 tuneGrid = tg)
print(wet.gbm)
plot(varImp(wet.gbm))
gbm_wet <- predict(wet.gbm, wetv) ## predict validation set

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
wet.pls <- train(wett, lt,
                 preProc = c("center", "scale"),
                 method = "pls",
                 tuneGrid = expand.grid(ncomp=seq(2, 10, by=1)),
                 trControl = tc)
print(wet.pls)
plot(varImp(wet.pls))
pls_wet <- predict(wet.pls, wetv) ## predict validation set

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
wet.bar <- train(wett, lt,
                 method = "bartMachine", 
                 preProc = c("center", "scale"),
                 trControl = tc,
                 tuneLength = 2,
                 seed = 123)
print(wet.bar)
plot(varImp(wet.bar))
bar_wet <- predict(wet.bar, wetv)

stopCluster(mc)
detach("package:bartMachine", unload=TRUE)

# Model stacking setup ----------------------------------------------------
pwetv <- as.data.frame(cbind(lv, rfo_wet, gbm_wet, pls_wet, bar_wet))
names(pwetv) <- c("L", "RFO", "GBM", "PLS", "BART")

# Remove extraneous objects from memory -----------------------------------
# rm(list=setdiff(ls(), pwetv"))

# Model stacking ----------------------------------------------------------
library(glmnet)

# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "cv", allowParallel = T)

# wet chem model stack
set.seed(1385321)
wet.ens <- train(L ~ ., data = pwetv,
                 method = "glmnet",
                 family = "gaussian",
                 trControl = tc)
print(wet.ens)
plot(varImp(wet.ens))
ens_wet <- as.data.frame(predict(wet.ens, pwetv))
names(ens_wet) <- c("ENS")
pwetv <- cbind(pwetv, ens_wet)

stopCluster(mc)

# Write data files --------------------------------------------------------
write.csv(pwetv, "Zn_pwetv.csv", row.names=F) ## adjust output name

# Prediction plots --------------------------------------------------------
# Plot individual model predictions
par(mfrow=c(2,2), mar=c(5,4.5,1,1))
lmin <- 0
lmax <- max(pwetv$L)
plot(L ~ RFO, pwetv, cex=1.2, xlim=c(lmin, lmax), ylim=c(lmin, lmax), xlab = "RFO prediction", ylab = "Observed uptake", cex.lab=1.3)
abline(c(0,1), col="red")
plot(L ~ GBM, pwetv, cex=1.2, xlim=c(lmin, lmax), ylim=c(lmin, lmax), xlab = "GBM prediction", ylab = "Observed uptake", cex.lab=1.3)
abline(c(0,1), col="red")
plot(L ~ PLS, pwetv, cex=1.2, xlim=c(lmin, lmax), ylim=c(lmin, lmax), xlab = "PLS prediction", ylab = "Observed uptake", cex.lab=1.3)
abline(c(0,1), col="red")
plot(L ~ BART, pwetv, cex=1.2, xlim=c(lmin, lmax), ylim=c(lmin, lmax), xlab = "BART prediction", ylab = "Observed uptake", cex.lab=1.3)
abline(c(0,1), col="red")

# Ensemble predictions 
par(mfrow=c(1,1), mar=c(5,4.5,1,1))
plot(L ~ ENS, pwetv, cex=1.2, xlim=c(lmin, lmax), ylim=c(lmin, lmax), xlab = "Model ensemble prediction", ylab = "Observed uptake", cex.lab=1.3)
abline(c(0,1), col="red")
