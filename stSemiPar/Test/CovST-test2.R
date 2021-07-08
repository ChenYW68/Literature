rm(list=ls())
source("./R/PSTVB_Packages.R")
source("./R/SetSimuModel.R")
source("./R/spT.validation.R")
source("./R/stSemiPar.R")
source("./R/stGCVfun.R")
source("./R/stSemiPar_WLS.R")
# Y_ts <- siMuIncF(n = 100, Nt = 10,
#                  x.0 = c(0),
#                  y.0 = c(0),
#                  delta = 0.1,
#                  para = list(tau.sq = 0.001, Phis = 0.3,
#                              nu = 1, sigma.sq.s = 1,
#                              sigma.sq.t = 1,
#                              Phit = 0.8,
#                              rho = 0.1,
#                              beta = c(1, 5)),
#                  nRatio = 0.8)
# range(Y_ts$D)
source("E:/Literature/semiBase/R/util.R")
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
# install.packages("E:/Literature/semiBase_1.0.zip", repos = NULL, type = "win.binary")
# library(semiBase)
DSN_01 <- odbcConnect(
  "DSN_01",
  uid = "myname",
  pwd = "mypwd",
  believeNRows = FALSE,
  case = "toupper"
)

##############################################################
seed <- 1:50
n <- 50
Nt <- 20
M <- c("WI", "WEC_t", "WEC_tw", "WEC_st", "WEC_stw", "WLS")
##############################################################
Phis <- "0.5"
sigma.sq.s <- 0.5
##############################################################
prob = c(1.0e0, 1e0)
Kernel <- c(0, 0)
nIter = 20
simDa <- siMuIncF(n = n, Nt = Nt, 
                  x.0 = c(0),
                  y.0 = c(0), 
                  delta = 0.1,
                  para = list(
                    Phis = as.numeric(Phis), 
                    nu = c(0.5, 0.5), 
                    sigma.sq.s = c(sigma.sq.s),
                    sigma.sq.t = c(1), 
                    # Phit = as.numeric(Phit),
                    # rho = 0.1,
                    # tau.sq = 0, 
                    beta = c(1, 5)))
Train <- list(Y_ts = simDa$Y_ts,
              X_ts = simDa$X_ts,#rep(1, nrow(TrainData)),simDa$train.X
              Z_ts = simDa$Z_ts,
              loc = simDa$loc,
              theta = simDa$theta,
              time = simDa$time
)
data <- Train





source("E:/Literature/semiBase/R/util.R")
# fit <- semiCovt(y = data$Y_ts, 
#          Time = data$time, 
#           Kernel = 0, 
#           h = 1e-1, prob = 5, nuUnifb = 1,
#           nu = 0.5, nThreads = 10)
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
start_time <- Sys.time()
C <- semiCovst(y = data$Y_ts, 
          Time = data$time, 
          Coord = data$loc, 
          Kernel = 0,
          h = 1e-2, 
          prob = c(0.1, 1.5), nuUnifb = 1,
          nu = 0.5, nThreads = 10)
end_time <- Sys.time()
end_time - start_time

# all.equal(as.vector(C$Cov2), as.vector(C$Cov))
range(eigen(C$ModCov$Cov)$value)
range(eigen(C$Cov.taper)$value)
range(eigen(C$Cov)$value)
range(eigen(simDa$Vc)$value)
M1 <- min(simDa$Vc,  as.matrix(C$Cov.taper),as.matrix(C$ModCov$Cov), as.matrix(C$Cov))
M2 <- max(simDa$Vc,  as.matrix(C$Cov.taper),as.matrix(C$ModCov$Cov), as.matrix(C$Cov))
par(mfrow = c(2, 2))
image.plot(simDa$Vc, zlim = c(M1, M2))
image.plot(as.matrix(C$Cov), zlim = c(M1, M2))
image.plot(as.matrix(C$Cov.taper), zlim = c(M1, M2))
image.plot(C$ModCov$Cov, zlim = c(M1, M2))

# i <- sample(1:n, 1)
# plot(simDa$Vc[i,], C$Cov[i, ])
# plot(simDa$Vc[i,], C$ModCov$Cov[i, ])