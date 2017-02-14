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
par(mfrow=c(2,3), mar=c(5,5,1,1))

# Boron
bref <- quantile(hidm$pB/hidm$B, p=c(0.05,0.5,0.95))
mirdat$dB <- (mirdat$pB/mirdat$B)
plot(log(dB)~log(B), mirdat, xlab=expression(paste("log"[e], " (", "B"[s], ")")), ylab=expression(paste("log"[e], " (", "B"[p], " / ", "B"[s], ")")), cex=1.2, cex.lab=1.5)
abline(h=log(bref[1]), lty=2)
abline(h=log(bref[2]), lwd=2)
abline(h=log(bref[3]), lty=2)
abline(rq(log(dB)~log(B), tau=0.5, mirdat), col="red", lwd=2)

# Copper
curef <- quantile(hidm$pCu/hidm$Cu, p=c(0.05,0.5,0.95))
mirdat$dCu <- (mirdat$pCu/mirdat$Cu)
plot(log(dCu)~log(Cu), mirdat, xlab=expression(paste("log"[e], " (", "Cu"[s], ")")), ylab=expression(paste("log"[e], " (", "Cu"[p], " / ", "Cu"[s], ")")), cex=1.2, cex.lab=1.5)
abline(h=log(curef[1]), lty=2)
abline(h=log(curef[2]), lwd=2)
abline(h=log(curef[3]), lty=2)
abline(rq(log(dCu)~log(Cu), tau=0.5, mirdat), col="red", lwd=2)

# Manganese
mnref <- quantile(hidm$pMn/hidm$Mn, p=c(0.05,0.5,0.95))
mirdat$dMn <- (mirdat$pMn/mirdat$Mn)
plot(log(dMn)~log(Mn), mirdat, xlab=expression(paste("log"[e], " (", "Mn"[s], ")")), ylab=expression(paste("log"[e], " (", "Mn"[p], " / ", "Mn"[s], ")")), cex=1.2, cex.lab=1.5)
abline(h=log(mnref[1]), lty=2)
abline(h=log(mnref[2]), lwd=2)
abline(h=log(mnref[3]), lty=2)
abline(rq(log(dMn)~log(Mn), tau=0.5, mirdat), col="red", lwd=2)

# Molybdenum
moref <- quantile(hidm$pMo/hidm$Mo, p=c(0.05,0.5,0.95))
mirdat$dMo <- (mirdat$pMo/mirdat$Mo)
plot(log(dMo)~log(Mo), mirdat, xlab=expression(paste("log"[e], " (", "Mo"[s], ")")), ylab=expression(paste("log"[e], " (", "Mo"[p], " / ", "Mo"[s], ")")), cex=1.2, cex.lab=1.5)
abline(h=log(moref[1]), lty=2)
abline(h=log(moref[2]), lwd=2)
abline(h=log(moref[3]), lty=2)
abline(rq(log(dMo)~log(Mo), tau=0.5, mirdat), col="red", lwd=2)

# Zinc
znref <- quantile(hidm$pZn/hidm$Zn, p=c(0.05,0.5,0.95))
mirdat$dZn <- (mirdat$pZn/mirdat$Zn)
plot(log(dZn)~log(Zn), mirdat, xlab=expression(paste("log"[e], " (", "Zn"[s], ")")), ylab=expression(paste("log"[e], " (", "Zn"[p], " / ", "Zn"[s], ")")), cex=1.2, cex.lab=1.5)
abline(h=log(znref[1]), lty=2)
abline(h=log(znref[2]), lwd=2)
abline(h=log(znref[3]), lty=2)
abline(rq(log(dZn)~log(Zn), tau=0.5, mirdat), col="red", lwd=2)

# Iron
feref <- quantile(hidm$pFe/hidm$Fe, p=c(0.05,0.5,0.95))
mirdat$dFe <- (mirdat$pFe/mirdat$Fe)
plot(log(dFe)~log(Fe), mirdat, xlab=expression(paste("log"[e], " (", "Fe"[s], ")")), ylab=expression(paste("log"[e], " (", "Fe"[p], " / ", "Fe"[s], ")")), cex=1.2, cex.lab=1.5)
abline(h=log(feref[1]), lty=2)
abline(h=log(feref[2]), lwd=2)
abline(h=log(feref[3]), lty=2)
abline(rq(log(dFe)~log(Fe), tau=0.5, mirdat), col="red", lwd=2)

dev.off()
