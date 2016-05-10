#' PLS regression predictions of wheat plant micronutrient contents
#' data courtesy FAO (1982) & ICRAF (2016)
#' M. Walsh, May 2016

# Required packages
# install.packages(c("devtools","caret","doParallel","pls")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(devtools)
  require(caret)
  require(doParallel)
  require(pls)
})

# Data setup --------------------------------------------------------------
# Run this first: https://github.com/mgwalsh/Bioavailability/blob/master/FAO_micro_setup.R
# or run ...
# SourceURL <- "https://raw.githubusercontent.com/mgwalsh/Bioavailability/master/FAO_micro_setup.R"
# source_url(SourceURL)

# Target variables: wheat dry matter yield (pDM, mg) & micro-nutrient concentrations (ppm)
pB  <- fao_cal$pB
pCu <- fao_cal$pCu
pMn <- fao_cal$pMn
pMo <- fao_cal$pMo
pZn <- fao_cal$pZn
pFe <- fao_cal$pFe

# Covariates
wetc <- fao_cal[c(4:24)] ## Wet chemistry calibration data
wetv <- fao_val[c(4:24)] ## Wet chemistry validation data from 8 randomly selected countries
mirc <- fao_cal[c(32:1795)] ## MIR calibration data
mirv <- fao_val[c(32:1795)] ## MIR validation data for 8 randomly selected countries

# PLS models --------------------------------------------------------------
# Start doParallel to parallelize model fitting
mc <- makeCluster(detectCores())
registerDoParallel(mc)

# Control setup
set.seed(1385321)
tc <- trainControl(method = "repeatedcv", number = 10, repeats = 3, allowParallel = TRUE)

# Plant Boron content (ppm) -----------------------------------------------
# Wet chemistry covariates
pB_wet.pls <- train(wetc, log(pB),
                    preProc = c("center", "scale"),
                    method = "pls",
                    tuneGrid = expand.grid(ncomp=seq(2,20,by=1)),
                    trControl = tc)
print(pB_wet.pls)
pB_wet.imp <- varImp(pB_wet.pls)
plot(pB_wet.imp, top=20)

# MIR covariates
pB_mir.pls <- train(mirc, log(pB),
                    preProc = c("center", "scale"),
                    method = "pls",
                    tuneGrid = expand.grid(ncomp=seq(2,40,by=1)),
                    trControl = tc)
print(pB_mir.pls)
pB_mir.imp <- varImp(pB_mir.pls)
plot(pB_mir.imp, top=20)

# Plant Copper content (ppm) ----------------------------------------------
# Wet chemistry covariates
pCu_wet.pls <- train(wetc, log(pCu),
                     preProc = c("center", "scale"),
                     method = "pls",
                     tuneGrid = expand.grid(ncomp=seq(2,20,by=1)),
                     trControl = tc)
print(pCu_wet.pls)
pCu_wet.imp <- varImp(pCu_wet.pls)
plot(pCu_wet.imp, top=20)

# MIR covariates
pCu_mir.pls <- train(mirc, log(pCu),
                     preProc = c("center", "scale"),
                     method = "pls",
                     tuneGrid = expand.grid(ncomp=seq(2,40,by=1)),
                     trControl = tc)
print(pCu_mir.pls)
pCu_mir.imp <- varImp(pCu_mir.pls)
plot(pCu_mir.imp, top=20)

# Plant Manganese content (ppm) -------------------------------------------
# Wet chemistry covariates
pMn_wet.pls <- train(wetc, log(pMn),
                     preProc = c("center", "scale"),
                     method = "pls",
                     tuneGrid = expand.grid(ncomp=seq(2,20,by=1)),
                     trControl = tc)
print(pMn_wet.pls)
pMn_wet.imp <- varImp(pMn_wet.pls)
plot(pMn_wet.imp, top=20)

# MIR covariates
pMn_mir.pls <- train(mirc, log(pMn),
                     preProc = c("center", "scale"),
                     method = "pls",
                     tuneGrid = expand.grid(ncomp=seq(2,40,by=1)),
                     trControl = tc)
print(pMn_mir.pls)
pMn_mir.imp <- varImp(pMn_mir.pls)
plot(pMn_mir.imp, top=20)

# Plant Molybdenum content (ppm) ------------------------------------------
# Wet chemistry covariates
pMo_wet.pls <- train(wetc, log(pMo),
                     preProc = c("center", "scale"),
                     method = "pls",
                     tuneGrid = expand.grid(ncomp=seq(2,20,by=1)),
                     trControl = tc)
print(pMo_wet.pls)
pMo_wet.imp <- varImp(pMo_wet.pls)
plot(pMo_wet.imp, top=20)

# MIR covariates
pMo_mir.pls <- train(mirc, log(pMo),
                     preProc = c("center", "scale"),
                     method = "pls",
                     tuneGrid = expand.grid(ncomp=seq(2,40,by=1)),
                     trControl = tc)
print(pMo_mir.pls)
pMo_mir.imp <- varImp(pMo_mir.pls)
plot(pMo_mir.imp, top=20)

# Plant Zinc content (ppm) ------------------------------------------------
# Wet chemistry covariates
pZn_wet.pls <- train(wetc, log(pZn),
                     preProc = c("center", "scale"),
                     method = "pls",
                     tuneGrid = expand.grid(ncomp=seq(2,20,by=1)),
                     trControl = tc)
print(pZn_wet.pls)
pZn_wet.imp <- varImp(pZn_wet.pls)
plot(pZn_wet.imp, top=20)

# MIR covariates
pZn_mir.pls <- train(mirc, log(pZn),
                     preProc = c("center", "scale"),
                     method = "pls",
                     tuneGrid = expand.grid(ncomp=seq(2,40,by=1)),
                     trControl = tc)
print(pZn_mir.pls)
pZn_mir.imp <- varImp(pZn_mir.pls)
plot(pZn_mir.imp, top=20)

# Plant Iron content (ppm) ------------------------------------------------
# Wet chemistry covariates
pFe_wet.pls <- train(wetc, log(pFe),
                     preProc = c("center", "scale"),
                     method = "pls",
                     tuneGrid = expand.grid(ncomp=seq(2,20,by=1)),
                     trControl = tc)
print(pFe_wet.pls)
pFe_wet.imp <- varImp(pFe_wet.pls)
plot(pFe_wet.imp, top=20)

# MIR covariates
pFe_mir.pls <- train(mirc, log(pFe),
                     preProc = c("center", "scale"),
                     method = "pls",
                     tuneGrid = expand.grid(ncomp=seq(2,40,by=1)),
                     trControl = tc)
print(pFe_mir.pls)
pFe_mir.imp <- varImp(pFe_mir.pls)
plot(pFe_mir.imp, top=20)

