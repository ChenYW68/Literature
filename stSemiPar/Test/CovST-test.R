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
                    sigma.sq.s = c(1e0, 1e-1),
                    sigma.sq.t = c(1, 1), 
                    # Phit = as.numeric(Phit),
                    # rho = 0.1,
                    # tau.sq = 0, 
                    beta = c(0, 0)))
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
range(simDa$D)
source("E:/Literature/semiBase/R/util.R")
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
start_time <- Sys.time()
Kernel = 0 
prob = c(2, 2)
h <- 1e-1
C <- semiCovst(y = simDa$W_ts, Time = simDa$time, 
                Coord = simDa$loc,
                Kernel = Kernel, h = h, 
                prob = prob, nThreads = 10)
Ct <- semiCovt(y = simDa$W_ts, simDa$time, Kernel = Kernel, h = h,
               prob = prob[2], nThreads = 10)
Cs <- semiCovs(y = simDa$W_ts, Coord = simDa$loc, Kernel = Kernel,
               h = h,prob = prob[1],
               nThreads = 10)
end_time <- Sys.time()
end_time - start_time
Vc.est1 <- kronecker(Cs$ModCov$Cov, Ct$ModCov$Cov)#
Vc.est2 <- kronecker(Cs$Cov, Ct$ModCov$Cov)#ModCov$Cov
M1 <- min(C$C, simDa$Vc, Vc.est1, Vc.est2)
M2 <- max(C$C, simDa$Vc, Vc.est1, Vc.est2)
# range(simDa$Vc)
# range(kronecker(Vs, simDa$Vt))
# range(simDa$Vs)
# range(simDa$Vt)
bk = 100
library(plot.matrix)
pdf(file = paste0("./figure", "/CovTest", ".pdf"),
    width = 10, height = 10)
par(mfrow = c(4, 3))
{
plot(simDa$Vc, border = NA, main = "Vst", 
     breaks= NULL, # breaks= seq(M1, M2,, bk), 
     col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
     key=list(side=4, cex.axis=size), 
     # axis.col = list(side = 3), 
     fmt.key="%.2f")
plot(Vc.est1, border = NA, main = "Vst.sep.Mod", 
     breaks= NULL, # breaks= seq(M1, M2,, bk),  
     col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
     key=list(side=4, cex.axis=size), fmt.key="%.2f")
plot(Vc.est2, border = NA, main = "Vst.sep.taper", 
     breaks = NULL, #breaks= seq(M1, M2,, bk), 
     col = topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
     key = list(side=4, cex.axis=size), fmt.key="%.2f")

plot(C$C, border = NA, main = "Vst.est", 
     breaks= NULL, #breaks= seq(M1, M2,, bk), 
     col = topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
     key = list(side=4, cex.axis=size), fmt.key="%.2f")
plot(C$ModCov$Cov, border = NA, main = "Vst.Mod", 
     breaks = NULL, #breaks= seq(M1, M2,, bk), 
     col = topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
     key = list(side=4, cex.axis=size), fmt.key="%.2f")
plot(C$Cov, border = NA, main = "Vst.taper", 
     breaks= NULL, #breaks= seq(M1, M2,, bk), 
     col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
     key=list(side=4, cex.axis=size), fmt.key="%.2f")
}
# 
# image.plot(Vc, zlim = c(M1, M2), main = "Vst")
# image.plot(Vc.est1, zlim = c(M1, M2), main = "Vst.sep.Mod")
# image.plot(Vc.est2, zlim = c(M1, M2), main = "Vst.sep.taper")
# image.plot(C$C, zlim = c(M1, M2), main = "Vst.est")
# image.plot(C$ModCov$Cov, zlim = c(M1, M2), main = "Vst.Mod")
# image.plot(C$Cov, zlim = c(M1, M2), main = "Vst.taper")
# dev.off()

{
  # pdf(file = paste0("./figure", "/CovTest", ".pdf"),
  #     width = 9, height = 7)
  # par(mfrow = c(2, 4))
  size = 0.8
  library(plot.matrix)
  brk <- 10
    M1 <- min(simDa$Vs, Cs$ModCov$Cov, Cs$Cov)
    M2 <- max(simDa$Vs, Cs$ModCov$Cov, Cs$Cov)
    plot(simDa$Vs, border = NA, main = "Vs", breaks=NULL, 
         # breaks= seq(M1, M2,, bk), 
         col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
         key=list(side=4, cex.axis=size), fmt.key="%.2f")
    
    plot(Cs$ModCov$Cov, border = NA, main = "Vs.est",
         # breaks= seq(M1, M2,, bk),  #
         breaks=NULL, 
         col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
         key=list(side=4, cex.axis=size), fmt.key="%.2f")  
    plot(Cs$Cov, border = NA, main = "Vs.taper", 
         # breaks= seq(M1, M2,, bk),  #
         breaks = NULL, 
         col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
         key=list(side=4, cex.axis=size), fmt.key="%.2f") 
    
    # image.plot(simDa$Vs, main = "Vs")
    # image.plot(fit$Cs$ModCov$Cov, main = "Vs.est")
    # image.plot(fit$Cs$Cov, main = "Vs.taper")

  brk <- 100

    M1 <- min(simDa$Vt, Ct$ModCov$Cov, Ct$Cov)
    M2 <- max(simDa$Vt, Ct$ModCov$Cov, Ct$Cov)
    plot(simDa$Vt, border = NA, main = "Vt", 
         # breaks= seq(M1, M2,, bk),  #
         breaks=NULL, 
         col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
         key=list(side=4, cex.axis=size), 
         fmt.key="%.2f")
    plot(Ct$ModCov$Cov, border = NA, main = "Vt.est",
         # breaks= seq(M1, M2,, bk),  #
         breaks = NULL,  
         col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
         key = list(side = 4, cex.axis=size), fmt.key="%.2f")
    plot(Ct$Cov, border = NA, main = "Vt.taper",
         # breaks= seq(M1, M2,, bk),  #
         breaks=NULL, 
         col= topo.colors(brk), #brewer.pal(5,"Spectral"), #tim.colors(10), 
         key=list(side=4, cex.axis=size), fmt.key="%.2f")
    # image.plot(simDa$Vt, main = "Vt")
    # image.plot(fit$Ct$ModCov$Cov, main = "Vt.est")
    # image.plot(fit$Ct$Cov, main = "Vt.taper")

}
dev.off()