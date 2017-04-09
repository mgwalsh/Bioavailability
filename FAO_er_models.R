#' Trace mineral soil-plant uptake/bioavailability models
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

# Allometric plots --------------------------------------------------------
hidm <- wetdat[ which(wetdat$HL=='H'), ] ## select high dry matter production reference group
par(mfrow=c(2,3), mar=c(5,5,1,1))

# Boron
pref <- quantile(hidm$pB/hidm$B, p=c(0.05,0.5,0.95))
sref <- quantile(hidm$B, p=c(0.05,0.5,0.95))
wetdat$dB <- (wetdat$pB/wetdat$B)
plot(log(dB)~log(B), wetdat, xlab=expression(paste("log"[e], " (", "B"[s], ")")),
     ylab=expression(paste("log"[e], " (", "B"[p], " / ", "B"[s], ")")), cex=1.2, cex.lab=1.5)
abline(h=log(pref[1]), col="blue") ## 5th percentile reference value
abline(h=log(pref[2]), lwd=2, col="red") ## median
abline(h=log(pref[3]), col="blue") ## 95th precentile
abline(v=log(sref[1]), col="blue") ## 5th percentile reference value
abline(v=log(sref[2]), lwd=2, col="red") ## median
abline(v=log(sref[3]), col="blue") ## 95th precentile
# abline(rq(log(dB)~log(B), tau=0.5, wetdat), col="red", lwd=2)

# Copper
pref <- quantile(hidm$pCu/hidm$Cu, p=c(0.05,0.5,0.95))
sref <- quantile(hidm$Cu, p=c(0.05,0.5,0.95))
wetdat$dCu <- (wetdat$pCu/wetdat$Cu)
plot(log(dCu)~log(Cu), wetdat, xlab=expression(paste("log"[e], " (", "Cu"[s], ")")),
     ylab=expression(paste("log"[e], " (", "Cu"[p], " / ", "Cu"[s], ")")), cex=1.2, cex.lab=1.5)
abline(h=log(pref[1]), col="blue") ## 5th percentile reference value
abline(h=log(pref[2]), lwd=2, col="red") ## median
abline(h=log(pref[3]), col="blue") ## 95th precentile
abline(v=log(sref[1]), col="blue") ## 5th percentile reference value
abline(v=log(sref[2]), lwd=2, col="red") ## median
abline(v=log(sref[3]), col="blue") ## 95th precentile
# abline(rq(log(dCu)~log(Cu), tau=0.5, wetdat), col="red", lwd=2)

# Manganese
pref <- quantile(hidm$pMn/hidm$Mn, p=c(0.05,0.5,0.95))
sref <- quantile(hidm$Mn, p=c(0.05,0.5,0.95))
wetdat$dMn <- (wetdat$pMn/wetdat$Mn)
plot(log(dMn)~log(Mn), wetdat, xlab=expression(paste("log"[e], " (", "Mn"[s], ")")),
     ylab=expression(paste("log"[e], " (", "Mn"[p], " / ", "Mn"[s], ")")), cex=1.2, cex.lab=1.5)
abline(h=log(pref[1]), col="blue") ## 5th percentile reference value
abline(h=log(pref[2]), lwd=2, col="red") ## median
abline(h=log(pref[3]), col="blue") ## 95th precentile
abline(v=log(sref[1]), col="blue") ## 5th percentile reference value
abline(v=log(sref[2]), lwd=2, col="red") ## median
abline(v=log(sref[3]), col="blue") ## 95th precentile
# abline(rq(log(dMn)~log(Mn), tau=0.5, wetdat), col="red", lwd=2)

# Molybdenum
pref <- quantile(hidm$pMo/hidm$Mo, p=c(0.05,0.5,0.95))
sref <- quantile(hidm$Mo, p=c(0.05,0.5,0.95))
wetdat$dMo <- (wetdat$pMo/wetdat$Mo)
plot(log(dMo)~log(Mo), wetdat, xlab=expression(paste("log"[e], " (", "Mo"[s], ")")),
     ylab=expression(paste("log"[e], " (", "Mo"[p], " / ", "Mo"[s], ")")), cex=1.2, cex.lab=1.5)
abline(h=log(pref[1]), col="blue") ## 5th percentile reference value
abline(h=log(pref[2]), lwd=2, col="red") ## median
abline(h=log(pref[3]), col="blue") ## 95th precentile
abline(v=log(sref[1]), col="blue") ## 5th percentile reference value
abline(v=log(sref[2]), lwd=2, col="red") ## median
abline(v=log(sref[3]), col="blue") ## 95th precentile
# abline(rq(log(dMo)~log(Mo), tau=0.5, wetdat), col="red", lwd=2)

# Zinc
# x11()
# par(mfrow=c(1,1), mar=c(5,5,1,1))
pref <- quantile(hidm$pZn/hidm$Zn, p=c(0.05,0.5,0.95))
sref <- quantile(hidm$Zn, p=c(0.05,0.5,0.95))
wetdat$dZn <- (wetdat$pZn/wetdat$Zn)
plot(log(dZn)~log(Zn), wetdat, xlab=expression(paste("log"[e], " (", "Zn"[s], ")")),
     ylab=expression(paste("log"[e], " (", "Zn"[p], " / ", "Zn"[s], ")")), cex=1.2, cex.lab=1.5)
abline(h=log(pref[1]), col="blue") ## 5th percentile reference value
abline(h=log(pref[2]), lwd=2, col="red") ## median
abline(h=log(pref[3]), col="blue") ## 95th precentile
abline(v=log(sref[1]), col="blue") ## 5th percentile reference value
abline(v=log(sref[2]), lwd=2, col="red") ## median
abline(v=log(sref[3]), col="blue") ## 95th precentile
# abline(rq(log(dZn)~log(Zn), tau=0.05, hidm), col="blue", lwd=1)
# abline(rq(log(dZn)~log(Zn), tau=0.5, hidm), col="red", lwd=2)
# abline(rq(log(dZn)~log(Zn), tau=0.95, hidm), col="blue", lwd=1)
# ZnQ <- rq(log(dZn)~log(Zn), tau=0.5, wetdat)
# summary(ZnQ)
# dev.copy(pdf, 'Zn_ER.pdf')
# dev.off()

# Iron
pref <- quantile(hidm$pFe/hidm$Fe, p=c(0.05,0.5,0.95))
sref <- quantile(hidm$Fe, p=c(0.05,0.5,0.95))
wetdat$dFe <- (wetdat$pFe/wetdat$Fe)
plot(log(dFe)~log(Fe), wetdat, xlab=expression(paste("log"[e], " (", "Fe"[s], ")")),
     ylab=expression(paste("log"[e], " (", "Fe"[p], " / ", "Fe"[s], ")")), cex=1.2, cex.lab=1.5)
# abline(h=log(pref[1]), col="blue") ## 5th percentile reference value
abline(h=log(pref[2]), lwd=2, col="red") ## median
abline(h=log(pref[3]), col="blue") ## 95th precentile
abline(v=log(sref[1]), col="blue") ## 5th percentile reference value
abline(v=log(sref[2]), lwd=2, col="red") ## median
abline(v=log(sref[3]), col="blue") ## 95th precentile
# abline(rq(log(dFe)~log(Fe), tau=0.5, wetdat), col="red", lwd=2)
# dev.off() ## remember to turn graphics device off!

# Quantreg models ---------------------------------------------------------
# Zn uptake as f(DM yield)
wetdat$uZn <- wetdat$pDM/10000*wetdat$pZn ## Zn uptake/pot
ZnY <- rq(log(uZn)~log(pDM), tau=c(0.025,0.5,0.975), data=wetdat)
print(ZnY)
summary(ZnY)
par(mfrow=c(1,3), mar=c(5,5,1,1))
plot(uZn~pDM, wetdat, xlab="Dry matter (mg / pot)", xlim=c(0,2500), ylab="Zn uptake (mg / pot)", ylim=c(0,25), cex=1.1, cex.lab=1.3)
curve(exp(ZnY$coefficients[1])*x^ZnY$coefficients[2], add=T, from=0, to=2500, col="blue", lwd=2)
curve(exp(ZnY$coefficients[5])*x^ZnY$coefficients[6], add=T, from=0, to=2500, col="blue", lwd=2)
curve(exp(ZnY$coefficients[3])*x^ZnY$coefficients[4], add=T, from=0, to=2500, col="red", lwd=2)

# QUEFTS-type fit
# ZnQ <- rq(log(pDM)~log(uZn), tau=c(0.025,0.5,0.975), data=wetdat)
# print(ZnQ)
# summary(ZnQ)
# plot(pDM~uZn, wetdat, xlab="Zn uptake (mg / pot)", xlim=c(0,25), ylab="Dry matter (mg / pot)", ylim=c(0,2500), cex=1.1, cex.lab=1.3)
# curve(exp(ZnQ$coefficients[1])*x^ZnQ$coefficients[2], add=T, from=0, to=25, col="blue", lwd=2)
# curve(exp(ZnQ$coefficients[5])*x^ZnQ$coefficients[6], add=T, from=0, to=25, col="blue", lwd=2)
# curve(exp(ZnQ$coefficients[3])*x^ZnQ$coefficients[4], add=T, from=0, to=25, col="red", lwd=2)

# Zn uptake as f(Zn concentration in soil)
ZnS <- rq(log(uZn)~log(Zn), tau=c(0.025,0.5,0.975), data=wetdat) ## Zn uptake/pot
print(ZnS)
summary(ZnS)
plot(uZn~Zn, wetdat, xlab="Soil Zn (ppm)", xlim=c(0,200), ylim=c(0,25), ylab="Zn uptake (mg / pot)", cex=1.1, cex.lab=1.3)
curve(exp(ZnS$coefficients[1])*x^ZnS$coefficients[2], add=T, from=0, to=200, col="blue", lwd=2)
curve(exp(ZnS$coefficients[5])*x^ZnS$coefficients[6], add=T, from=0, to=200, col="blue", lwd=2)
curve(exp(ZnS$coefficients[3])*x^ZnS$coefficients[4], add=T, from=0, to=200, col="red", lwd=2)

# Zn uptake index as f(Zn concentration in soil & DM yield)
ZnI <- rq(log(uZn)~log(Zn)*log(pDM), tau=c(0.025,0.5,0.975), data=wetdat) ## Zn uptake/pot
print(ZnI)
summary(ZnI)
