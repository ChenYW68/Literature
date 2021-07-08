rm(list=ls())
source("./R/PSTVB_Packages.R")
source("./R/SetSimuModel.R")
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
Nt <- 20
M <- c("WI", "WEC_t", "WEC_tw", "WEC_st", "WEC_stw", "WLS")
##############################################################
Phis <- c("0.05", "0.01", "0.1")
sigma.sq.s <- c(5, 2, 0.5)
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
# Train <- list(Y_ts = simDa$Y_ts,
#               X_ts = simDa$X_ts,#rep(1, nrow(TrainData)),simDa$train.X
#               Z_ts = simDa$Z_ts,
#               loc = simDa$loc,
#               theta = simDa$theta,
#               time = simDa$time#,
#               # Q_temp = Matrix::solve(simDa$Vc)
# )
# data <- Train

library(fda)
Nt <- ncol(data$Y_ts)
n <- nrow(data$Y_ts)
daybasis365 = create.fourier.basis(c(0,Nt), Nt)
# plot(daybasis365)
harmLfd = vec2Lfd(c(0,(2*pi/(10*Nt))^2, 0), c(0, Nt))

tempfdPar = fdPar(daybasis365,harmLfd, 2e3)
tempfd = smooth.basis(1:Nt, t(data$Y_ts),tempfdPar)

# plot(tempfd$fd)


ptemppca = pca.fd(tempfd$fd,nharm=3,
                  harmfdPar=tempfdPar, 
                  centerfns = T)
#temppca$values are the eigenvalues
ptemppca$varprop
ptemppca$values[1:5]
# Get FPCs
pharmfd = ptemppca$harmonics
pharmvals = eval.fd(1:Nt,pharmfd)
pharmvals[, 1]%*% pharmvals[, 1]



# quartz()
par(mfrow=c(2, 2))
for(i in 1:3){
plot(data$time, data$Vt[, i])
}
matplot(data$time,pharmvals,xlab='time',ylab='PCs',
        lwd=2,lty=1,cex.lab=1.5,cex.axis=1.5,type='l')
# legend(0.5,-0.05,c('PC1','PC2'),col=1:3,lty=0.1,lwd=1)
title('Temperature Principle Component Functions')


# Show the mean curves - temppca$meanfdr(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(tempfd$fd,xlab='time',ylab='simulation',cex.lab=1.5,cex.axis=1.5,col=4)
lines(ptemppca$meanfd,lwd=2.5,col=2)

# functional principal components
harmfd = ptemppca$harmonics
harmvals = eval.fd(1:Nt,harmfd)
dim(harmvals) # The top 4 FPCs

par(mfrow=c(1,1),mar = c(8, 8, 4, 2))
plot(1:Nt,harmvals[,1],xlab='time',ylab='PCs',
     lwd=4,lty=1,cex.lab=2,cex.axis=2,type='l')


# Now let's do a bit of reconstruction for Montreal

# Do FPCA on the temperature curves without the Montreal curve
# s <- sample(1:50, 1)
s <- 20
Stemppca = pca.fd(tempfd$fd[-s],nharm=1,harmfdPar=tempfdPar)
# Get the FPC
harms = Stemppca$harmonics
# Get the mean curve
meanfd = Stemppca$meanfd

# Montreal data
Mdat = data$Y_ts[s, ]
time <- 1:Nt
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

par(mfrow=c(1,1),mar = c(8, 8, 4, 2)) #
plot(time,Rvals,type='l',lwd=1,col=4,xlab='time',ylab='y',
     cex.lab=2.5,cex.axis=2.5)
points(time,Mdat,col=2,lwd=5)
points(time[1:tt],Mdat[1:tt],col=3)
# lines(Stemppca$meanfd,lty=2,col=5,lwd=2)
