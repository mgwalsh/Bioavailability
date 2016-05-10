#' Wheat plant micro-nutrient bioavailability prediction setup
#' Soil and wheat plant wet chemistry data courtesy of FAO, 1982 (doc @ https://www.dropbox.com/s/gwk07tanhu86tqj/Silanpaa%20Report.pdf?dl=0)
#' MIR soil data courtesy of ICRAF
#' M. Walsh, May 2016

# install.packages(c("downloader","caret"), dependencies=T)
suppressPackageStartupMessages({
  require(downloader)
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
soils <- read.table("soils.csv", header=T, sep=",") ## FAO soil chem data (in ppm, except for pH(water), EC, CaCO3, CEC, texture, Volwt)
soils <- merge(cid, soils, by="CC")
mir <- read.table("mir.csv", header=T, sep=",") ## ICRAF MIR data
plant <- read.table("plants.csv", header=T, sep=",") ## FAO plant dry matter yield (pDM, mg) and micronutrient data (ppm)

# Assemble dataframes
wetdat <- merge(soils, plant, by="SSID")
mirdat <- merge(wetdat, mir, by="SSID")

# Train/Test set partition ------------------------------------------------
set.seed(1385321)
faoIndex <- createDataPartition(mirdat$SSID, p = 4/5, list = FALSE, times = 1)
fao_cal <- mirdat[ faoIndex,]
fao_val <- mirdat[-faoIndex,]

# Write data files --------------------------------------------------------
write.csv(fao_cal, "fao_cal.csv", row.names=F)
write.csv(fao_val, "fao_val.csv", row.names=F)
