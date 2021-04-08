# Malawi Maize grain Zn predictions
# Data courtesy of the GeoNutrion project (http://www.geonutrition.com/)
# M. Walsh, April 2020

# Install & load packages -------------------------------------------------
# package names
packages <- c("downloader", "devtools", "caret", "plyr", "MASS", "randomForest", "gbm", "Cubist", "quantreg", "doParallel")

# install packages
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# load packages
invisible(lapply(packages, library, character.only = TRUE))

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("MWI_Zn_data", showWarnings=F)
setwd("./MWI_Zn_data")
dir.create("./Results")

samp <- read.table("MWI_Zn_data.csv", header = T, sep = ",")
samp <- samp[complete.cases(samp[,c(7:18)]), ] ## removes incomplete cases
samp <- within(samp, rm(Oxalates)) # I drop the sum of oxalates and fit these individually

# set calibration/validation set randomization seed
seed <- 12358
set.seed(seed)

# split data into calibration and validation sets
gsIndex <- createDataPartition(samp$Zn, p = 4/5, list = F, times = 1)
gs_cal <- samp[ gsIndex,]
gs_val <- samp[-gsIndex,]

# Soil calibration labels
labs <- c("Zn")
lcal <- as.vector(t(gs_cal[labs]))

# soil calibration features
fcal <- gs_cal[,7:17]

# Training models with cross-validation -----------------------------------
# Random forest
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", allowParallel = TRUE)
tg <- expand.grid(mtry = seq(1,10, by=1)) ## model tuning steps

# model training
rf <- train(fcal, lcal,
            preProc = c("center","scale"),
            method = "rf",
            ntree = 301,
            metric = "RMSE",
            importance = TRUE,
            tuneGrid = tg,
            trControl = tc)

rfImp <- varImp(rf)
plot(rfImp)
gs_val$rf <- predict(rf, gs_val)
stopCluster(mc)
fname <- paste("./Results/", labs, "_rf.rds", sep = "")
saveRDS(rf, fname)

# GBM
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(1385321)
tc <- trainControl(method = "cv", importance = TRUE, allowParallel = TRUE)

## for initial <gbm> tuning guidelines see @ https://stats.stackexchange.com/questions/25748/what-are-some-useful-guidelines-for-gbm-parameters
tg <- expand.grid(interaction.depth = seq(2,5, by=1), shrinkage = 0.01, n.trees = seq(101,501, by=50),
                  n.minobsinnode = 50) ## model tuning steps

# model training
gb <- train(fcal, lcal, 
            method = "gbm", 
            preProc = c("center", "scale"),
            trControl = tc,
            tuneGrid = tg,
            metric = "RMSE")

gbImp <- varImp(gb)
plot(gbImp)
gs_val$gb <- predict(gb, gs_val)
stopCluster(mc)
fname <- paste("./Results/", labs, "_gb.rds", sep = "")
saveRDS(gb, fname)

# Cubist
# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# control setup
set.seed(seed)
tc <- trainControl(method="repeatedcv", number=10, repeats=3, allowParallel = T)
# tg <- needs tuning

cu <- train(fcal, lcal, 
            method = "cubist", 
            trControl = tc,
            importance = TRUE,
            metric = "RMSE")

cuImp <- varImp(cu)
plot(cuImp)
gs_val$cu <- predict(cu, gs_val)
stopCluster(mc)
fname <- paste("./Results/", labs, "_cu.rds", sep = "")
saveRDS(cu, fname)

# Stacking model on validation set ----------------------------------------
lval <- as.vector(t(gs_val[labs]))
fval <- gs_val[,19:21] ## subset validation features

# start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# model setup
set.seed(seed)
tc <- trainControl(method="repeatedcv", number=10, repeats=3, allowParallel=T)

st <- train(fval, lval,
            method = "glmStepAIC",
            trControl = tc,
            metric = "RMSE")

gs_val$st <- predict(st, gs_val) ## stacked predictions
stopCluster(mc)
fname <- paste("./Results/", labs, "_st.rds", sep = "")
saveRDS(st, fname)

# Prediction check --------------------------------------------------------
samp$rf <- predict(rf, samp)
samp$gb <- predict(gb, samp)
samp$cu <- predict(cu, samp)
samp$st <- predict(st, samp)

par(pty="s")
plot(Zn~st, xlab="Ensemble Zn prediction (mg/kg)", ylab="Measured grain Zn (mg/kg)", cex.lab=1.2, 
     xlim=c(10,50), ylim=c(10,50), samp)
stQ <- rq(Zn~st, tau=c(0.025,0.5,0.975), data=samp) ## quantile regression fit
print(stQ)
curve(stQ$coefficients[4]*x+stQ$coefficients[3], add=T, from=10, to=50, col="red", lwd=2)
abline(c(0,1), col="grey", lwd=1)

