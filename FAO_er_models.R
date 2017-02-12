#' Trace mineral soil-plant enrichment ratio models
#' Soil and wheat plant wet chemistry data courtesy of FAO (doc @ https://www.dropbox.com/s/gwk07tanhu86tqj/Silanpaa%20Report.pdf?dl=0)
#' MIR soil data courtesy of ICRAF (2016)
#' M. Walsh, February 2017

# Required packages
# install.packages(c("devtools","quantreg")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(devtools)
  require(quantreg)
})

# Data setup --------------------------------------------------------------
# Run this first: https://github.com/mgwalsh/Bioavailability/blob/master/FAO_micro_setup.R
# or run ...
# SourceURL <- "https://raw.githubusercontent.com/mgwalsh/Bioavailability/master/FAO_micro_setup.R"
# source_url(SourceURL)

# Basic allometry ---------------------------------------------------------
hidm <- mirdat[ which(mirdat$HL=='H'), ] ## select high dry matter production reference group

# Boron
bref <- quantile(hidm$pB/hidm$B, p=1/2)
mirdat$dB <- (mirdat$pB/mirdat$B)/bref
plot(log(dB)~log(B), mirdat)
dB.lm <- lm(log(dB)~log(B), mirdat)

# Copper
curef <- quantile(hidm$pCu/hidm$Cu, p=1/2)
mirdat$dCu <- (mirdat$pCu/mirdat$Cu)/curef
plot(log(dCu)~log(Cu), mirdat)
dCu.lm <- lm(log(dCu)~log(Cu), mirdat)

# Manganese
mnref <- quantile(hidm$pMn/hidm$Mn, p=1/2)
mirdat$dMn <- (mirdat$pMn/mirdat$Mn)/mnref
plot(log(dMn)~log(Mn), mirdat)
dMn.lm <- lm(log(dMn)~log(Mn), mirdat)

# Molybdenum
moref <- quantile(hidm$pMo/hidm$Mo, p=1/2)
mirdat$dMo <- (mirdat$pMo/mirdat$Mo)/moref
plot(log(dMo)~log(Mo), mirdat)
dMo.lm <- lm(log(dMo)~log(Mo), mirdat)

# Zinc
znref <- quantile(hidm$pZn/hidm$Zn, p=1/2)
mirdat$dZn <- (mirdat$pZn/mirdat$Zn)/znref
plot(log(dZn)~log(Zn), mirdat)
dZn.lm <- lm(log(dZn)~log(Zn), mirdat)

# Iron
feref <- quantile(hidm$pFe/hidm$Fe, p=1/2)
mirdat$dFe <- (mirdat$pFe/mirdat$Fe)/feref
plot(log(dFe)~log(Fe), mirdat)
dFe.lm <- lm(log(dFe)~log(Fe), mirdat)
