source("./R/PSTVB_Packages.R")
source("./R/SetSimuModel.R")
source("E:/Literature/semiBase/R/util.R")
n = 10
Nt = 15
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

# semiCovst <- function(y, Time, Coord,
#                       Kernel, h, d = NULL,
#                       prob = 1, nuUnifb,
#                       nu, nThreads){
  
  y = simDa$Y_ts
  Time <- simDa$time
  Coord <- simDa$loc
  Kernel <- 0
  h <- 0.1
  d = NULL
  prob = 1
  nuUnifb = 1
  nu = 0.5
  nThreads = 1
  
  y <- as.matrix(y)
  n <- nrow(y)
  Nt <- ncol(y)
  
  thresh.s <- 1
  if(is.null(d)){
    d <- fields::rdist(Coord, Coord)
    thresh.s <- max(d)*prob
    rho.space <- fields::Wendland(d, theta = thresh.s, dimension = 1, k =1)
    # rho.space[which(rho.space > 0)] <- 1
    lower.s <- which(lower.tri(d), arr.ind = T)
    d <- d[lower.tri(d, diag = F)]
  }
  m1 <- length(d)
  # thresh <- #quantile(d, probs = prob)
  index.upp <- which(d > thresh.s)
  index.low <- which(d <= thresh.s)
  
  d <- c(d[index.low])
  m <- length(d)
  
  Vt <- diag(Nt)
  Vt <- abs(row(Vt) - col(Vt))
  thresh.t <- max(Vt)*prob
  rho.time <- fields::Wendland(Vt, theta = thresh.t, dimension = 1, k =1)
  
  if(length(h)!=4){h <- rep(h, 4)}
  
  # cv_index <- (which(lower.tri(Cmat), arr.ind = T))
  # data.table::setorderv(cv_index, c("row", "col"))
  # method 1
  lower <- which(lower.tri(rho.time), arr.ind = T)
  index <- data.frame(as.data.frame(lower), value = rho.time[lower])
  cv_index.0 <- CvIndex <- as.matrix(index[index$value!=0, 1:2])
  index.upp <- as.matrix(index[index$value==0, 1:2])
  
  # method 2
  # cv_index.0 <- cv_index <- which(lower.tri(Cmat), arr.ind = T)
  
  
  cv_len <- nrow(CvIndex)
  
  storage.mode(y) <- "double"
  storage.mode(Time) <- "double"
  storage.mode(Coord) <- "double"
  
  storage.mode(CvIndex) <- "integer"
  storage.mode(n) <- "integer"
  storage.mode(Nt) <- "integer"
  storage.mode(cv_len) <- "integer"
  storage.mode(Kernel) <- "integer"
  
  storage.mode(h) <- "double"
  storage.mode(d) <- "double"
  storage.mode(m) <- "integer"
  storage.mode(nuUnifb) <- "integer"
  storage.mode(nu) <- "double"
  storage.mode(nThreads) <- "integer"
  Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
  Cst <- semiCov_Cst( y, Time, Coord, CvIndex, 
                      n, Nt, cv_len, Kernel, h, 
                      d, m, nuUnifb, nu, nThreads)
  
  Ct <- semiCovt(y, Time, Kernel, h,
                 prob = 1, nuUnifb, 
                 nu, nThreads)
    
  Cmat <- matrix(0, nrow = n*Nt, ncol= n*Nt)
  for(i in 1:nrow(lower.s)){
    temp <- matrix(0, nrow = Nt, ncol= Nt)
    temp[cv_index.0] <- Cst$Cov[, i]
    temp <- temp + t(temp)
    diag(temp) <- Cst$Var[, i]
   
    Cmat[((lower.s[i, 1] - 1)*Nt + 1):(lower.s[i, 1]*Nt), 
         ((lower.s[i, 2] - 1)*Nt + 1):(lower.s[i, 2]*Nt)] <- temp
   
  }
  Cmat <- t(Cmat) + Cmat
  range(Cmat)
  
  for(i in 1:n){
    Cmat[((i - 1)*Nt + 1):(i*Nt), ((i - 1)*Nt + 1):(i*Nt)] <- Ct$Cov
  }
  

  Var <- diag(Cmat)
  diag(Cmat) <- 0
  Var <- as.double(Var)
  Cmat1 <- as.double(Cmat)
  # a <- ModCov(Cmat, Var, n)
  ModCov.0 <- ModCov(Cmat1, Var, n*Nt)
  Cov <- ModCov.0$Cov * kronecker(rho.space, rho.time)
  # return(list(ModCov = ModCov.0, Cov = Cov, Var.st = Var))
  # }
  
#   return(Cst)
# }