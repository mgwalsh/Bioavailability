#' Wheat plant dry matter micro-nutrient bioavailability prediction setup
#' Soil and wheat plant wet chemistry data courtesy of FAO, 1982 (doc @ https://www.dropbox.com/s/gwk07tanhu86tqj/Silanpaa%20Report.pdf?dl=0)
#' MIR soil data courtesy of ICRAF
#' M. Walsh, May 2016

# install.packages(c("downloader","compositions","MASS","RColorBrewer","caret"), dependencies=T)
suppressPackageStartupMessages({
  require(downloader)
  require(compositions)
  require(MASS)
  require(RColorBrewer)
  require(caret)
 })

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("FAO_data", showWarnings=F)
setwd("./FAO_data")

# Download 
download("https://www.dropbox.com/s/hhdoxswpfb9vlz1/FAO_micro_bioavailability.zip?dl=0", "FAO_micro_bioavailability.zip", mode="wb")
unzip("FAO_micro_bioavailability.zip", overwrite=T)
cid <- read.table("countries.csv", header=T, sep=",") ## country ID's
soils <- read.table("soils.csv", header=T, sep=",") ## FAO soil chem data (in ppm, except for pH, EC, CaCO3, CEC, texture, Volwt)
soils <- merge(cid, soils, by="CC")
mir <- read.table("mir.csv", header=T, sep=",") ## ICRAF MIR data
mir <- mir[!duplicated(mir[,1]), ]
plant <- read.table("plants.csv", header=T, sep=",") ## FAO plant dry matter yield (pDM, mg) and micronutrient data (ppm)
wetdat <- merge(soils, plant, by="SSID")
wetdat <- wetdat[!duplicated(wetdat), ]

# High/Low biomass label (mg DM / pot)
qlevel <- quantile(wetdat$pDM, p=1/2)
wetdat$HL <- as.factor(ifelse(wetdat$pDM > qlevel, "H", "L"))

# Soil nutrient profile setup ---------------------------------------------
fpart <- names(wetdat[c(13:17, 19:21, 23:24)])
cdata <- wetdat[fpart]
cdata$Fv <- 1000000-rowSums(cdata[fpart]) ## calculates "fill value" (Fv), in mg/kg soil
scoda <- acomp(cdata)

# Parallel coordinates plot of soil nutrient profiles with High/Low biomass labels
cdata$HL <- as.factor(wetdat$HL)
k <- adjustcolor(brewer.pal(3, "Set1")[cdata$HL], alpha=.8)
parcoord(cdata[,1:11], col = k, var.label=T)

# Parallel coordinates plot of plant nutrients with High/Low biomass labels
pdata <- wetdat[,26:31]
names(pdata) <- c("B","Cu","Mn","Mo","Zn","Fe")
parcoord(pdata, col = k, var.label=T)

# Plot of MIR spectra with High/Low biomass labels
w <- read.csv("wavelengths.csv", header=F)
spec <- as.matrix(mir[2:1765])
matplot(w, t(spec), type="l", col=k, lty=1, xlim=c(4000,500), ylim=c(-0.02,0.02), 
        xlab="Wavenumber", ylab="1rst derivative log(1/R)", cex.lab=1.3)

# Sequential binary partion & isometric log ratio (ilr) transform
bpart <- t(matrix(c( 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,-1,
                    -1,-1,-1,-1,-1, 1, 1, 1, 1, 1, 0,
                    -1,-1, 1, 1, 1, 0, 0, 0, 0, 0, 0,
                     0, 0, 1,-1,-1, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 1,-1, 0, 0, 0, 0, 0, 0,
                     1,-1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
                     0, 0, 0, 0, 0, 1, 1,-1, 1,-1, 0,
                     0, 0, 0, 0, 0,-1, 1, 0,-1, 0, 0,
                     0, 0, 0, 0, 0, 1, 0, 0,-1, 0, 0,
                     0, 0, 0, 0, 0, 0, 0, 1, 0,-1, 0), ncol=11, nrow=10, byrow=T))
CoDaDendrogram(X=acomp(scoda), signary=bpart, type="lines") ## compositional balance mobile graph				
idata <- as.data.frame(ilr(scoda, V=bpart))
parcoord(idata, col = k, var.label=T)

# Assemble dataframes
wetdat <- cbind(wetdat, idata)
mirdat <- merge(wetdat, mir, by="SSID")

# Plot soil-plant nutrient scatter plots
par(mfrow=c(3,2), mar=c(4.5,4.5,1,1))
plot(pB~B, cex=1.2, xlab="Soil B (ppm)", ylab="Plant B (ppm)", cex.lab=1.5, wetdat)
plot(pCu~Cu, cex=1.2, xlab="Soil Cu (ppm)", ylab="Plant Cu (ppm)", cex.lab=1.5, wetdat)
plot(pMn~Mn, cex=1.2, xlab="Soil Mn (ppm)", ylab="Plant Mn (ppm)", cex.lab=1.5, wetdat)
plot(pMo~Mo, cex=1.2, xlab="Soil Mo (ppm)", ylab="Plant Mo (ppm)", cex.lab=1.5, wetdat)
plot(pZn~Zn, cex=1.2, xlab="Soil Zn (ppm)", ylab="Plant Zn (ppm)", cex.lab=1.5, wetdat)
plot(pFe~Fe, cex=1.2, xlab="Soil Fe (ppm)", ylab="Plant Fe (ppm)", cex.lab=1.5, wetdat)

# Train/Test set partition ------------------------------------------------
set.seed(1385321)
faoIndex <- createDataPartition(mirdat$SSID, p = 3/4, list = FALSE, times = 1)
fao_cal <- mirdat[ faoIndex,] ## random 75% for calibration
fao_val <- mirdat[-faoIndex,] ## random 25% for validation

# Write data files --------------------------------------------------------
write.csv(fao_cal, "fao_cal.csv", row.names=F)
write.csv(fao_val, "fao_val.csv", row.names=F)

# Remove extraneous objects from memory -----------------------------------
rm(list=setdiff(ls(), c("mirdat", "fao_cal", "fao_val")))
dev.off()

