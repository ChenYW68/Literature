source("./R/PSTVB_Packages.R")
source("./R/SetSimuModel.R")
source("E:/Literature/semiBase/R/util.R")
n = 10
Nt = 20
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
simDa <- siMuIncF(n = n, Nt = Nt, 
                  x.0 = c(0),
                  y.0 = c(0), 
                  delta = 0.1,
                  para = list(
                    Phis = as.numeric(0.2), 
                    nu = c(0.5, 0.5), 
                    sigma.sq.s = c(5e-1, 1e-1),
                    sigma.sq.t = c(1, 1), 
                    # Phit = as.numeric(Phit),
                    # rho = 0.1,
                    # tau.sq = 0, 
                    beta = c(1, 5)))
Vs <- Matern(d = simDa$D, range = 0.5,
            smoothness = 0.5,
            phi = 1)
Vt <- diag(Nt)
phi.t <- 0.3
Vt <- phi.t^abs(row(Vt) - col(Vt))/(1 - phi.t^2)

## Cross covariance
Vc <- kronecker(Vs, Vt)
Time <- 1:Nt/Nt

W_ts <- t(matrix(Matrix::crossprod(Matrix::chol(Vc), 
                                   rep(rnorm(n * Nt))),
                 # + rnorm(n * Nt, 0, sd = sqrt(0.1))
                 nrow = Nt, ncol = n))

source("E:/Literature/semiBase/R/util.R")
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
start_time <- Sys.time()
Kernel = 0 
prob = c(0.8, 0.1)
h <- 1e-1
C <- semiCovst(y = W_ts, Time = Time, 
                Coord = simDa$loc,
                Kernel = Kernel, h = h, 
                prob = prob, nThreads = 10)
Ct <- semiCovt(y = W_ts, Time, Kernel = Kernel, h = h,
               prob = prob[2], nThreads = 10)
Cs <- semiCovs(y= W_ts, Coord= simDa$loc, Kernel = Kernel,
               h = h,prob = prob[1],
               nThreads = 10)
end_time <- Sys.time()
end_time - start_time
Vc.est1 <- kronecker(Cs$ModCov$Cov, Ct$ModCov$Cov)#
Vc.est2 <- kronecker(Cs$Cov, Ct$Cov)#ModCov$Cov
M1 <- min(C$C, Vc, Vc.est1, Vc.est2)
M2 <- max(C$C, Vc, Vc.est1, Vc.est2)

bk = 100
library(plot.matrix)
pdf(file = paste0("./figure", "/CovTest", ".pdf"),
    width = 9, height = 7)
par(mfrow = c(3, 2))
plot(Vc, border = NA, main = "Vst", 
     breaks= seq(M1, M2,, bk), 
     col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
     key=list(side=4, cex.axis=size), fmt.key="%.2f")
plot(Vc.est1, border = NA, main = "Vst.sep.Mod", 
     breaks= seq(M1, M2,, bk),  
     col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
     key=list(side=4, cex.axis=size), fmt.key="%.2f")
plot(Vc.est2, border = NA, main = "Vst.sep.taper", 
     breaks= seq(M1, M2,, bk), 
     col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
     key=list(side=4, cex.axis=size), fmt.key="%.2f")

plot(Vc.est2, border = NA, main = "Vst.est", 
     breaks= seq(M1, M2,, bk), 
     col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
     key=list(side=4, cex.axis=size), fmt.key="%.2f")
plot(C$ModCov$Cov, border = NA, main = "Vst.Mod", 
     breaks= seq(M1, M2,, bk), 
     col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
     key=list(side=4, cex.axis=size), fmt.key="%.2f")
plot(C$Cov, border = NA, main = "Vst.taper", 
     breaks= seq(M1, M2,, bk), 
     col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
     key=list(side=4, cex.axis=size), fmt.key="%.2f")

# 
# image.plot(Vc, zlim = c(M1, M2), main = "Vst")
# image.plot(Vc.est1, zlim = c(M1, M2), main = "Vst.sep.Mod")
# image.plot(Vc.est2, zlim = c(M1, M2), main = "Vst.sep.taper")
# image.plot(C$C, zlim = c(M1, M2), main = "Vst.est")
# image.plot(C$ModCov$Cov, zlim = c(M1, M2), main = "Vst.Mod")
# image.plot(C$Cov, zlim = c(M1, M2), main = "Vst.taper")
dev.off()