#' Wheat plant micro-nutrient bioavailability prediction setup
#' Soil and wheat plant data courtesy FAO, 1982 (doc @ https://www.dropbox.com/s/gwk07tanhu86tqj/Silanpaa%20Report.pdf?dl=0)
#' MIR soil data courtesy ICRAF
#' M. Walsh, May 2016

# install.packages(c("downloader"), dependencies=T)
suppressPackageStartupMessages({
  require(downloader)
 })

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("FAO_data", showWarnings=F)
setwd("./FAO_data")

# Download 
download("https://www.dropbox.com/s/hhdoxswpfb9vlz1/FAO_micro_bioavailability.zip?dl=0", "FAO_micro_bioavailability.zip", mode="wb")
unzip("FAO_micro_bioavailability.zip", overwrite=T)
cid <- read.table("countries.csv", header=T, sep=",") ## country ID's
soils <- read.table("soils.csv", header=T, sep=",") ## FAO soil chemistry data
soils <- merge(cid, soils, by="CC")
mir <- read.table("mir.csv", header=T, sep=",") ## ICRAF MIR data
plant <- read.table("plants.csv", header=T, sep=",") ## FAO plant DM yield and micronutrient data

# Assemble dataframes
wetdat <- merge(soils, plant, by="SSID")
mirdat <- merge(wetdat, mir, by="SSID")

# Train/Test set partition ------------------------------------------------
country <- names(table(mirdat$Country))
set.seed(1385321)
train <- sample(country, 0.8*length(country)) ## sample 80% of of countries for calibration
fao_cal <- mirdat[ mirdat$Country%in%train, ] ## calibration data
fao_val <- mirdat[!mirdat$Country%in%train, ] ## validation data

# Write data files --------------------------------------------------------
write.csv(fao_cal, "fao_cal.csv", row.names=F)
write.csv(fao_val, "fao_val.csv", row.names=F)
