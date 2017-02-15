#' Wheat plant micro-nutrient bioavailability quantiles
#' Soil and wheat plant wet chemistry data courtesy of FAO (doc @ https://www.dropbox.com/s/gwk07tanhu86tqj/Silanpaa%20Report.pdf?dl=0)
#' MIR soil data courtesy of ICRAF (2016)
#' M. Walsh, May 2016

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

# Quantile regressions ----------------------------------------------------
par(mfrow=c(2,3), mar=c(5,5,1,1))

# Boron
plot(pB~B, xlab=expression(paste("B"[s], " (ppm)")), ylab=expression(paste("B"[p], " (ppm)")), cex.lab=1.5, 
     xlim=c(-0.1,10.1), cex=1.2, mirdat)
BQ <- rq(log(pB)~log(B), tau=c(0.05,0.25,0.5,0.75,0.95), data=mirdat)
print(BQ)
curve(exp(BQ$coefficients[1])*x^BQ$coefficients[2], add=T, from=0, to=9.5, col="blue", lwd=1)
curve(exp(BQ$coefficients[3])*x^BQ$coefficients[4], add=T, from=0, to=9.5, col="blue", lty=2)
curve(exp(BQ$coefficients[5])*x^BQ$coefficients[6], add=T, from=0, to=9.5, col="red", lwd=2)
curve(exp(BQ$coefficients[7])*x^BQ$coefficients[8], add=T, from=0, to=9.5, col="blue", lty=2)
curve(exp(BQ$coefficients[9])*x^BQ$coefficients[10], add=T, from=0, to=9.5, col="blue", lwd=1)

# Copper
plot(pCu~Cu, xlab=expression(paste("Cu"[s], " (ppm)")), ylab=expression(paste("Cu"[p], " (ppm)")), cex.lab=1.5,
     xlim=c(-0.1,60.1), cex=1.2, mirdat)
CuQ <- rq(log(pCu)~log(Cu), tau=c(0.05,0.25,0.5,0.75,0.95), data=mirdat)
print(CuQ)
curve(exp(CuQ$coefficients[1])*x^CuQ$coefficients[2], add=T, from=0, to=58, col="blue", lwd=1)
curve(exp(CuQ$coefficients[3])*x^CuQ$coefficients[4], add=T, from=0, to=58, col="blue", lty=2)
curve(exp(CuQ$coefficients[5])*x^CuQ$coefficients[6], add=T, from=0, to=58, col="red", lwd=2)
curve(exp(CuQ$coefficients[7])*x^CuQ$coefficients[8], add=T, from=0, to=58, col="blue", lty=2)
curve(exp(CuQ$coefficients[9])*x^CuQ$coefficients[10], add=T, from=0, to=58, col="blue", lwd=1)

# Manganese
plot(pMn~Mn, xlab=expression(paste("Mn"[s], " (ppm)")), ylab=expression(paste("Mn"[p], " (ppm)")), cex.lab=1.5, 
     xlim=c(-0.1,300.1), cex=1.2, mirdat)
MnQ <- rq(log(pMn)~log(Mn), tau=c(0.05,0.25,0.5,0.75,0.95), data=mirdat)
print(MnQ)
curve(exp(MnQ$coefficients[1])*x^MnQ$coefficients[2], add=T, from=0, to=290, col="blue", lwd=1)
curve(exp(MnQ$coefficients[3])*x^MnQ$coefficients[4], add=T, from=0, to=290, col="blue", lty=2)
curve(exp(MnQ$coefficients[5])*x^MnQ$coefficients[6], add=T, from=0, to=290, col="red", lwd=2)
curve(exp(MnQ$coefficients[7])*x^MnQ$coefficients[8], add=T, from=0, to=290, col="blue", lty=2)
curve(exp(MnQ$coefficients[9])*x^MnQ$coefficients[10], add=T, from=0, to=290, col="blue", lwd=1)

# Molybdenum
plot(pMo~Mo, xlab=expression(paste("Mo"[s], " (ppm)")), ylab=expression(paste("Mo"[p], " (ppm)")), cex.lab=1.5,
     xlim=c(-0.1,3.1), ylim=c(-0.1,4.1), cex=1.2, mirdat)
MoQ <- rq(log(pMo)~log(Mo), tau=c(0.05,0.25,0.5,0.75,0.95), data=mirdat)
print(MoQ)
curve(exp(MoQ$coefficients[1])*x^MoQ$coefficients[2], add=T, from=0, to=2.9, col="blue", lwd=1)
curve(exp(MoQ$coefficients[3])*x^MoQ$coefficients[4], add=T, from=0, to=2.9, col="blue", lty=2)
curve(exp(MoQ$coefficients[5])*x^MoQ$coefficients[6], add=T, from=0, to=2.9, col="red", lwd=2)
curve(exp(MoQ$coefficients[7])*x^MoQ$coefficients[8], add=T, from=0, to=2.9, col="blue", lty=2)
curve(exp(MoQ$coefficients[9])*x^MoQ$coefficients[10], add=T, from=0, to=2.9, col="blue", lwd=1)

# Zinc
# x11()
# par(mfrow=c(1,1), mar=c(5,5,1,1))
plot(pZn~Zn, xlab=expression(paste("Zn"[s], " (ppm)")), ylab=expression(paste("Zn"[p], " (ppm)")), cex.lab=1.5,
     xlim=c(-0.1,200.1), cex=1.2, mirdat)
ZnQ <- rq(log(pZn)~log(Zn), tau=c(0.05,0.25,0.5,0.75,0.95), data=mirdat)
print(ZnQ)
curve(exp(ZnQ$coefficients[1])*x^ZnQ$coefficients[2], add=T, from=0, to=190, col="blue", lwd=1)
curve(exp(ZnQ$coefficients[3])*x^ZnQ$coefficients[4], add=T, from=0, to=190, col="blue", lty=2)
curve(exp(ZnQ$coefficients[5])*x^ZnQ$coefficients[6], add=T, from=0, to=190, col="red", lwd=2)
curve(exp(ZnQ$coefficients[7])*x^ZnQ$coefficients[8], add=T, from=0, to=190, col="blue", lty=2)
curve(exp(ZnQ$coefficients[9])*x^ZnQ$coefficients[10], add=T, from=0, to=190, col="blue", lwd=1)
# dev.copy(pdf, 'Zn_quant.pdf')
# dev.off()

# Iron
plot(pFe~Fe, xlab=expression(paste("Fe"[s], " (ppm)")), ylab=expression(paste("Fe"[p], " (ppm)")), cex.lab=1.5,
     xlim=c(-0.1,500.1), cex=1.2, mirdat)
FeQ <- rq(log(pFe)~log(Fe), tau=c(0.05,0.25,0.5,0.75,0.95), data=mirdat)
print(FeQ)
curve(exp(FeQ$coefficients[1])*x^FeQ$coefficients[2], add=T, from=0, to=490, col="blue", lwd=1)
curve(exp(FeQ$coefficients[3])*x^FeQ$coefficients[4], add=T, from=0, to=490, col="blue", lty=2)
curve(exp(FeQ$coefficients[5])*x^FeQ$coefficients[6], add=T, from=0, to=490, col="red", lwd=2)
curve(exp(FeQ$coefficients[7])*x^FeQ$coefficients[8], add=T, from=0, to=490, col="blue", lty=2)
curve(exp(FeQ$coefficients[9])*x^FeQ$coefficients[10], add=T, from=0, to=490, col="blue", lwd=1)
