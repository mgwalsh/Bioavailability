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
spx <- wetdat[c(32,4,9,7,33:42)]

# Graphical models --------------------------------------------------------
# Fit interaction model
gm0 <- mmod(~.^., data = spx)
gm1 <- stepwise(gm0)

# Plot interaction graph
dag <- ugList(terms(gm1), result="matrix")
net <- as.network(x = dag, directed = FALSE, loops = FALSE, matrix.type = "adjacency")
plot.network(net, vertex.col = "white", vertex.cex = 5, displaylabels = T, label.pos = 5, mode = "circle")

# Fit naive Bayes model
gm2 <- mmod(~pH:HL+CEC:HL+Sand:HL+V1:HL+V2:HL+V3:HL+V4:HL+V5:HL+V6:HL+V7:HL+V8:HL+V9:HL+V10:HL, data=spx)
gm3 <- stepwise(gm2)
dag <- ugList(terms(gm3), result="matrix")
net <- as.network(x = dag, directed = FALSE, loops = FALSE, matrix.type = "adjacency")
plot.network(net, vertex.col = "white", vertex.cex = 6, displaylabels = T, label.pos = 5, mode = "circle")

