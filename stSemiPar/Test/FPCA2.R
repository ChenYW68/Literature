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
n <- 100
Nt <- 40
M <- c("WI", "WEC_t", "WEC_tw", "WEC_st", "WEC_stw", "WLS")
##############################################################
Phis <- c("0.05", "0.01", "0.1")
sigma.sq.s <- c(5, 2, 1)
##############################################################
prob = c(1.0e0, 1e0)
Kernel <- c(0, 0)
nIter = 30

data <- siMuIncF(n = n, Nt = Nt, 
                 x.0 = c(0),
                 y.0 = c(0), 
                 delta = 0.1,
                 para = list(
                   Phis = as.numeric(Phis), 
                   nu = c(0.5, 0.5, 0.5), 
                   sigma.sq.s = c(sigma.sq.s),
                   nugget = c(0.1), 
                   # Phit = as.numeric(Phit),
                   # rho = 0.1,
                   # tau.sq = 0, 
                   beta = c(0, 0)))

dim(data$Psi_mat)
dim(data$Cs)



library(fda)
Nt <- ncol(data$Y_ts)
n <- nrow(data$Y_ts)

source("E:/Literature/semiBase/R/util.R")
fd <- sp.FPCA(data$Y_ts, data$time, nharm = 3, lambda = 1e-7)
plot.sp.fd(data$time, data$Vt, PCA = fd$PCA)
source("E:/Literature/semiBase/R/util.R")
Psi <- build_Psi_fun(Psi = fd$PCA, n = n, Nt = nrow(fd$PCA))
Psi

n = 50
D <- data$D[1:n, 1:n]
D1 <- build_spCov_fun(n = n, J = 3, D, range = c(rep(0.1, 3)),
                            sigma.sq.s = c(rep(1, 3)),
                            nu = c(rep(0.5, 3)))

D1 <- as.matrix(D1)
source("E:/Literature/semiBase/R/util.R")
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
D2 <- build_spCov(Coord = data$loc[1:n, ], J = 3, 
            range = c(rep(0.1, 3)),
            sigma.sq.s = c(rep(1, 3)),
            nu = c(rep(0.5, 3)), 
            CovModel = 0)
# all.equal(D1, D2$C)
all.equal(as.vector(D1), 
          as.vector(D2))
J = 3
source("E:/Literature/semiBase/R/util.R")
# Psi.crossprod <- crossprod(data$Psi_mat)
system.time(
r2 <- optim(rep(-2, J + 1), optim_semPara, 
               method = "L-BFGS-B",#L-BFGS-B
              y = data$Y_ts,
              Coord = data$loc, 
              Psi = data$Psi_mat, 
              # Phi.tem = Psi.crossprod,
              J = 3, sigma.sq.s = sigma.sq.s,
              CovModel = 0,
              lower = c(rep(log(0.01), 3), log(0.05)), 
              upper = c(rep(log(0.5), 3), log(0.5)),
              control = list(REPORT = 1,
                           trace = 1))
)
exp(r2$par)
# lower = lower, #upper = 6,
# control = list(REPORT = 1,
#                trace = 0)
# , hessian = F
# )





diff(data$time)[1]
daybasis365 = create.fourier.basis(c(0, 1),Nt, 
                                   period = 1)
plot(daybasis365)
harmLfd = vec2Lfd(c(0,(2*pi/Nt)^2,0), c(0, 1))

tempfdPar = fdPar(daybasis365, harmLfd, 1e-3)
tempfd = smooth.basis(data$time, t(data$Y_ts),tempfdPar)

# plot(tempfd$fd)


ptemppca = pca.fd(tempfd$fd, nharm = 3, 
                  harmfdPar = tempfdPar)
#temppca$values are the eigenvalues
ptemppca$varprop
ptemppca$values
var(ptemppca$scores[, 1])
# Get FPCs
pharmfd = ptemppca$harmonics
pharmvals = eval.fd(data$time, pharmfd)
pharmvals[, 1]%*% pharmvals[, 1]
# pharmvals[, 1] <- pharmvals[, 1]/sqrt(as.vector(pharmvals[, 1]%*% pharmvals[, 1]))

# quartz()
par(mfrow=c(2, 2))
for(i in 1:3){
  plot(data$time, data$Vt[, i])
}
matplot(1:Nt,pharmvals,xlab='day',ylab='PCs',
        lwd=2,lty=1,cex.lab=1.5,cex.axis=1.5,type='l')
# legend(0.5,-0.05,c('PC1','PC2'),col=1:3,lty=0.1,lwd=1)
title('Temperature Principle Component Functions')


# Show the mean curves - temppca$meanfdr(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(tempfd$fd,xlab='time',ylab='simulation',
     cex.lab=1.5,cex.axis=1.5,col=4)
lines(ptemppca$meanfd,lwd=2.5,col=2)

# functional principal components
harmfd = ptemppca$harmonics
harmvals = eval.fd(data$time,harmfd)
dim(harmvals) # The top 4 FPCs

par(mfrow=c(1,1), mar = c(8, 8, 4, 2))
plot(data$time, harmvals[, 1], xlab='time', ylab='PCs',
     lwd=4, lty=1, cex.lab=2, cex.axis=2, type='l')


# Now let's do a bit of reconstruction for Montreal

# Do FPCA on the temperature curves without the Montreal curve
s <- sample(1:50, 1)
# s <- 20
Stemppca = pca.fd(tempfd$fd[-s],nharm=1,harmfdPar=tempfdPar)
# Get the FPC
harms = Stemppca$harmonics
# Get the mean curve
meanfd = Stemppca$meanfd

# Montreal data
Mdat = data$Y_ts[s, ]
time <- data$time
tt <- 10
length(Mdat)
# evaluate FPC in the days [1:132]
Stempvals = eval.fd(time[1:tt],harms)
# evaluate the mean curve in the days [1:132]
mtempvals = eval.fd(time[1:tt],meanfd)

# Remove Mean curve from the Montreal data
Mdat2 = Mdat[1:tt]-mtempvals

# Obtain the FPC scores of the Montreal curve
coef = lm(Mdat2~Stempvals-1)$coef
coef

# Prediction for the Montreal curve
Rfd = coef[1]*harms[1]+#coef[2]*harms[2]+
  #   coef[3]*harms[3]+coef[4]*harms[4]+
  Stemppca$meanfd


Rvals = eval.fd(time,Rfd)

par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(time,Rvals,type='l',lwd=1,col=4,xlab='time',ylab='y',
     cex.lab=2.5,cex.axis=2.5)
points(time,Mdat,col=2,lwd=5)
points(time[1:tt],Mdat[1:tt],col=3)
# lines(Stemppca$meanfd,lty=2,col=5,lwd=2)
