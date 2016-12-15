require(compositions)

# Compositional analysis setup --------------------------------------------
fpart <- names(wetdat[c(13:17, 19:21, 23:24)]) 
cdata <- wetdat[fpart]
cdata$Fv <- 1000000-rowSums(cdata[fpart]) ## calculates "fill value" (Fv), in mg/kg soil
scoda <- acomp(cdata)

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
wetdat <- cbind(wetdat, idata)

