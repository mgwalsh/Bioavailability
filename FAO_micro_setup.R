# Wheat plant dry matter yield and micro-nutrient uptake prediction setup
# Soil and wheat plant wet chemistry data courtesy of FAO, 1982 (doc @ https://www.dropbox.com/s/gwk07tanhu86tqj/Silanpaa%20Report.pdf?dl=0)
# MIR soil data courtesy of ICRAF
# M. Walsh, March 2019

# install.packages(c("downloader","MASS","RColorBrewer","caret"), dependencies=T)
suppressPackageStartupMessages({
  require(downloader)
  require(MASS)
  require(caret)
 })

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("FAO_data", showWarnings=F)
setwd("./FAO_data")

# Download 
download("https://www.dropbox.com/s/hhdoxswpfb9vlz1/FAO_micro_bioavailability.zip?raw=1", "FAO_micro_bioavailability.zip", mode="wb")
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

# Assemble dataframes
mirdat <- merge(wetdat, mir, by="SSID")

# Train/Test set partition ------------------------------------------------
set.seed(1385321)
faoIndex <- createDataPartition(mirdat$SSID, p = 4/5, list = FALSE, times = 1)
fao_cal <- mirdat[ faoIndex,] ## random 80% for calibration
fao_val <- mirdat[-faoIndex,] ## random 20% for validation

# Write data files --------------------------------------------------------
write.csv(fao_cal, "fao_cal.csv", row.names=F)
write.csv(fao_val, "fao_val.csv", row.names=F)

# Remove extraneous objects from memory -----------------------------------
rm(list=setdiff(ls(), c("wetdat", "mirdat", "fao_cal", "fao_val")))
