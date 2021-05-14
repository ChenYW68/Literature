# Rcpp::sourceCpp("./src/nonPara.cpp")
source("E:/Literature/semiBase/R/util.R")
source("./R/PSTVB_Packages.R")
n <- 100
delta <- 0.1
x.coords <- 0 + (1:(n - 1))*delta
y.coords <- 0 + (1:(n - 1))*delta


Coords <- matrix(NA, ncol = 1, nrow = 1)
while (nrow(Coords)!=n) {
  Coords <- cbind(sample(x.coords, n, replace = T), 
                  sample(y.coords, n, replace = T)
  )
  Coords <- unique(Coords)
}

D <- fields::rdist(Coords, Coords)
range(D)
sigmaMat <- Matern(D, range = 0.5, nu = 0.5)
w0 <- as.vector(mvnfast::rmvn(1, rep(0, nrow(D)), 
                              sigma = sigmaMat, 
                              ncores = 10))
C <- semiCovs( w0, 
                     Coords,
                     Kernel = 0,
                     h = 3e-1, 
                     nuUnifb = 1,
                     nu = 0,
                     nThreads = 10)

x <- sigmaMat[upper.tri(sigmaMat)]
y <- C$C
# y[which(y < 1e-3)] <- 0

plot(x, y)

d <- fields::rdist(Coords, Coords)
d <- as.vector(d[upper.tri(d)])
plot(d, x, ylim = c(min(x, y), max(x, y)))
points(d, y, col = "red")
