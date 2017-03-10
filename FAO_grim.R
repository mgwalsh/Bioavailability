#' Graphical interaction models of greenhouse grown wheat plants from soil data
#' Soil and wheat plant wet chemistry data courtesy of FAO (doc @ https://www.dropbox.com/s/gwk07tanhu86tqj/Silanpaa%20Report.pdf?dl=0)
#' M. Walsh, March 2017

# install.packages(c("downloader","gRim","statnet"), dependencies=T)
suppressPackageStartupMessages({
  require(downloader)
  require(gRim)
  require(statnet)
})

# Data setup --------------------------------------------------------------
# Run this first: https://github.com/mgwalsh/Bioavailability/blob/master/FAO_micro_setup.R
# or run ...
# SourceURL <- "https://raw.githubusercontent.com/mgwalsh/Bioavailability/master/FAO_micro_setup.R"
# source_url(SourceURL)

# Select variables
psm <- wetdat[,26:42]
