#' Trace mineral soil-plant enrichment ratio models
#' Soil and wheat plant wet chemistry data courtesy of FAO (doc @ https://www.dropbox.com/s/gwk07tanhu86tqj/Silanpaa%20Report.pdf?dl=0)
#' MIR soil data courtesy of ICRAF (2016)
#' M. Walsh, February 2017

# Required packages
# install.packages(c("devtools","arm")), dependencies=TRUE)
suppressPackageStartupMessages({
  require(devtools)
  require(quantreg)
})

# Data setup --------------------------------------------------------------
# Run this first: https://github.com/mgwalsh/Bioavailability/blob/master/FAO_micro_setup.R
# or run ...
# SourceURL <- "https://raw.githubusercontent.com/mgwalsh/Bioavailability/master/FAO_micro_setup.R"
# source_url(SourceURL)

# Allometric plots --------------------------------------------------------
hidm <- mirdat[ which(mirdat$HL=='H'), ] ## select high dry matter production reference group
par(mfrow=c(2,3), mar=c(5,4.5,1,1))

# Boron
bref <- quantile(hidm$pB/hidm$B, p=1/2)
mirdat$dB <- (mirdat$pB/mirdat$B)/bref
plot(log(dB)~log(B), mirdat, xlab="log[Soil B]", cex=1.2, cex.lab=1.5)
abline(h=0)
abline(rq(log(dB)~log(B), tau=0.5, mirdat), col="red", lwd=2)

# Copper
curef <- quantile(hidm$pCu/hidm$Cu, p=1/2)
mirdat$dCu <- (mirdat$pCu/mirdat$Cu)/curef
plot(log(dCu)~log(Cu), mirdat, xlab="log[Soil Cu]", cex=1.2, cex.lab=1.5)
abline(h=0)
abline(rq(log(dCu)~log(Cu), tau=0.5, mirdat), col="red", lwd=2)

# Manganese
mnref <- quantile(hidm$pMn/hidm$Mn, p=1/2)
mirdat$dMn <- (mirdat$pMn/mirdat$Mn)/mnref
plot(log(dMn)~log(Mn), mirdat, xlab="log[Soil Mn]", cex=1.2, cex.lab=1.5)
abline(h=0)
abline(rq(log(dMn)~log(Mn), tau=0.5, mirdat), col="red", lwd=2)

# Molybdenum
moref <- quantile(hidm$pMo/hidm$Mo, p=1/2)
mirdat$dMo <- (mirdat$pMo/mirdat$Mo)/moref
plot(log(dMo)~log(Mo), mirdat, xlab="log[Soil Mo]", cex=1.2, cex.lab=1.5)
abline(h=0)
abline(rq(log(dMo)~log(Mo), tau=0.5, mirdat), col="red", lwd=2)

# Zinc
znref <- quantile(hidm$pZn/hidm$Zn, p=1/2)
mirdat$dZn <- (mirdat$pZn/mirdat$Zn)/znref
plot(log(dZn)~log(Zn), mirdat, xlab="log[Soil Zn]", cex=1.2, cex.lab=1.5)
abline(h=0)
abline(rq(log(dZn)~log(Zn), tau=0.5, mirdat), col="red", lwd=2)

# Iron
feref <- quantile(hidm$pFe/hidm$Fe, p=1/2)
mirdat$dFe <- (mirdat$pFe/mirdat$Fe)/feref
plot(log(dFe)~log(Fe), mirdat, xlab="log[Soil Fe]", cex=1.2, cex.lab=1.5)
abline(h=0)
abline(rq(log(dFe)~log(Fe), tau=0.5, mirdat), col="red", lwd=2)

# dev.off()
