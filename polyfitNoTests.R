#POLYFITNOTESTS Fit polynomial to data.
#
# This is based on the Matlab function polyfit but tests have been stripped
# out to speed it up. Therefore it is only suitable for use where there are
# unlikely to be problems with the conditioning.
#
# edited Joanne Williams 2018 joll@noc.ac.uk
#
# ---------------
# Code history: 
#
# polyfitNoTests.m Matlab script Converted to R in November 2023 by 
# Nikos Kalligeris

polyfitNoTests <- function(x,y,n) {
  
  # library(matlib)
  
  # compute Vandermonde matrix.
  V <- outer(x, seq(0, n), `^`)
  # and flip 
  N<-V[,(n+1):1] 
  
  # Solve least squares problem through orthogonal-triangular decomposition
  res <- QR(N)
  tmp <- t(res$Q) %*% matrix(t(y))
  # and by solving the system of linear equations
  p <- solve(res$R, tmp) # vector p contains the polynomial coefficients
  
  return(p) 
}