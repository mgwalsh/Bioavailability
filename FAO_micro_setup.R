#' Plant micro-nutrient bioavailability prediction setup
#' Soil and plant data courtesy FAO, 1982 (doc @ https://www.dropbox.com/s/gwk07tanhu86tqj/Silanpaa%20Report.pdf?dl=0)
#' MIR soil data courtesy ICRAF
#' M. Walsh, May 2016

# install.packages(c("downloader","compositions","colorRamps","RColorBrewer","MASS"), dependencies=T)
suppressPackageStartupMessages({
  require(downloader)
  require(compositions)
})

# Data setup --------------------------------------------------------------
# Create a data folder in your current working directory
dir.create("FAO_data", showWarnings=F)
setwd("./FAO_data")

# Download
download("https://www.dropbox.com/s/hhdoxswpfb9vlz1/FAO_micro_bioavailability.zip?dl=0", "FAO_micro_bioavailability.zip", mode="wb")
unzip("FAO_micro_bioavailability.zip", overwrite=T)
cid <- read.table("countries.csv", header=T, sep=",") ## country ID's
soil <- read.table("soil.csv", header=T, sep=",") ## FAO soil chemistry data
mir <- read.table("mir.csv", header=T, sep=",") ## ICRAF MIR data
plant <- read.table("plant.csv", header=T, sep=",") ## FAO plant yield and micronutrient data
dat <- merge(cid, soil, mir, plant, by="SSID")
