source("./R/PSTVB_Packages.R")
source("./R/SetSimuModel.R")
source("E:/Literature/semiBase/R/util.R")
simDa <- siMuIncF(n = 100, Nt = 20, 
                  x.0 = c(0),
                  y.0 = c(0), 
                  delta = 0.1,
                  para = list(tau.sq = 0.5, 
                              Phis = 0.8, 
                              nu = 1, 
                              sigma.sq.s = 10,
                              sigma.sq.t = 1, 
                              Phit = 0.5,
                              rho = 0.5,
                              beta = c(1, 5)))
simDa$Vt
simDa$Vc
simDa$Vs

dd <- simDa$D[upper.tri(simDa$D)]
Time <- simDa$time
t1 = Time[1]
t2 = Time[2]
y = simDa$W_ts
coord = simDa$loc

CovEst <- function(t1, t2, y = simDa$W_ts,
                   coord = simDa$loc, d, 
                   Time, h = c(5e-1, 5e-1)){
  # y = Y_ts$W_ts
  Nt <- ncol(y)
  n <- nrow(y)
  N <- n*Nt*(Nt - 1)
  # h <- 1
  # t1 <- 1
  # t2 <- 2
  # Time <- 1:(Nt)
  d1 <- Time - t1
  d2 <- Time - t2
  K1 <- exp(-abs(d1)/h[2])
  
  K2 <- exp(-abs(d2)/h[2])
  
  A1 <- A2 <- 0
  R1 <- matrix(c(1, NA, NA), ncol = 1)
 
  # for(k in 1:m){
    for(s1 in 1:n){
      for(s2 in 1:n){
        cs <- exp(-abs(as.vector(Rdist(coord[s1,], coord[s2,])$Dist) - d)/h[1])
        for(i in 1:Nt){
          for(j in 1:Nt){
            if(j!=i){
              
              R1[2, 1] <- (Time[i] - t1)/h[2]
              R1[3, 1] <- (Time[j] - t2)/h[2]
              # cat(R1[3, 1], "\n")
              A1 <- A1 + (R1%*% t(R1))*K1[i]*K2[j]*cs
              A2 <- A2 + R1*K1[i]*K2[j]*cs*y[s1, i]*y[s2, j]
              #print(round(c(i - 1, j -1,  R1[3, 1], as.vector(A2)), 3))
            }
          }
        }
      }
    # }
  alpha <- solve(A1) %*% A2
}
return(alpha[1])
}
h <- c(5e-1, 5e-1)
Nt <- ncol(simDa$W_ts)
Cov <- matrix(NA, nrow = Nt, ncol = Nt)
Time <- simDa$time
Vc <- simDa$Vc[(Nt + 1):(2*Nt), 1:Nt]
for(k in 1:1){
for(t1 in 1:(Nt - 1)){
  for(t2 in (t1 + 1):Nt){
    Cov[t1, t2] <- CovEst(Time[t1], Time[t2], y = simDa$W_ts, 
                          coord = simDa$loc, d = dd[k],
                          Time = Time, h = h)
    cat("t = ", t1, "Cov = ",  Cov[t1, t2],
        "; True = ", Vc[t1, t2],"\n")
    # if(t == (s + 1)){
    plot(Vc[1:t1, ], Cov[1:t1, ])
    # }
    # points(simDa$train.sigmaMat[s, (s + 1):t],
    #      Cov[s, (s + 1):t])
  }
 }
}


source("./R/PSTVB_Packages.R")
source("./R/SetSimuModel.R")
source("E:/Literature/semiBase/R/util.R")
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
simDa <- siMuIncF(n = 5, Nt = 6, 
                  x.0 = c(0),
                  y.0 = c(0), 
                  delta = 0.1,
                  para = list(tau.sq = 1, 
                              Phis = 0.1, 
                              nu = 1, 
                              sigma.sq.s = 1,
                              sigma.sq.t = 1, 
                              Phit = 0.5,
                              rho = 0.5,
                              beta = c(0, 0)))
simDa$Vt
simDa$Vc
simDa$Vs

dd <- simDa$D[upper.tri(simDa$D)]
Time <- simDa$time
t1 = Time[1]
t2 = Time[2]
y = simDa$W_ts
coord = simDa$loc

CovEst <- function(y = simDa$W_ts,
                   coord = simDa$loc,
                   prob = 1,
                   h = c(1e-1)){
  # y = Y_ts$W_ts
  Nt <- ncol(y)
  n <- nrow(y)
  # h <- 1
  # t1 <- 1
  # t2 <- 2
  # Time <- 1:(Nt)
  d <- Rdist(coord, coord)$Dist
  d <- d[upper.tri(d)]
  m1 <- length(d)
  thresh <- quantile(d, probs = prob)
  index.upp <- which(d > thresh)
  index.low <- which(d <= thresh)
  
  d <- c(0, d[index.low])
  
 
  R1 <- matrix(c(1, 0), ncol = 1)
  l = 1
  alpha <- vector()
  for(k in 1:length(d)){
    A1 <- A2 <- 0
  for(s1 in 1:n){
    for(s2 in 1:n){
      ds <- as.vector(Rdist(coord[s1,], coord[s2,])$Dist) - d[k]
      cs <- exp(-abs(ds)/h)/h
      for(t in 1:Nt){
          # if(j!=i){
            R1[2, 1] <- ds/h
            # cat(R1[3, 1], "\n")
            A1 <- A1 + (R1%*% t(R1))*cs
            A2 <- A2 + R1*cs*y[s1, t]*y[s2, t]
            #print(round(c(i - 1, j -1,  R1[3, 1], as.vector(A2)), 3))
          # }
        }
    }
  }
    alpha[l] <- (solve(A1) %*% A2)[1]
    l <- l + 1
  }
  Var <- alpha[1]
  Cov <- rep(0, m1)
  Cov[index.upp] <- 0
  Cov[index.low] <- alpha[-1]
  
  Cmat <- matrix(0, nrow = n, ncol= n)
  cv_index <- as.data.frame(which(upper.tri(Cmat), arr.ind = T))
  setorderv(cv_index, c("row", "col"))
  
  cv_index.0 <- as.matrix(cv_index)
  Cmat[cv_index.0] <- Cov
  
  Cmat <- (t(Cmat) + Cmat)
  diag(Cmat) <- Var
  
  return(Cmat)
}

h = 2e-1
alpha <- CovEst(y = simDa$W_ts, coord = simDa$loc, 
                prob = 1, h = h)
alpha
# all.equal(Cs, alpha)
h = 2e-1
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
start_time <- Sys.time()
source("E:/Literature/semiBase/R/util.R")
Cs <- semiCovs(y = simDa$W_ts, 
                Coord = simDa$loc,     
                Kernel = 0,
                h = h, 
                prob = 1,
                nuUnifb = 0,
                nu = 0,
                nThreads = 15)
Cs
end_time <- Sys.time()
end_time - start_time
all.equal(Cs, Cs)
library(reticulate)
scipy <- import("scipy")


simDa$Vs
Cs
plot(simDa$Vs, Cs$Vec)

Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
start_time <- Sys.time()
Ct1<- semiCovt(simDa$W_ts,
             Time = simDa$time, 
             Kernel = 0,
             h = 5e-1,
             nuUnifb = 0,
             nu = 0,
             nThreads = 10)
eigen(Ct$Cmat)$value


end_time <- Sys.time()
end_time - start_time
all.equal(Ct$Cmat, Ct$Cov)

Ct$Cmat
simDa$Vt
plot(simDa$Vt, Ct$Cov)

Vc.est <- kronecker(Cs$Cov, Ct$Cov)
plot(simDa$Vc[upper.tri(simDa$Vc)], Vc.est[upper.tri(Vc.est)])

theta <- theta_Cov_fun(y_ts = simDa$Y_ts, 
                        z_ts = simDa$Z_ts,
                        Time = simDa$time,
                        Q = kronecker(solve(Cs), 
                                      solve(Ct$Cmat)),
                        h = 5e-1)
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
h = 2e-1
Cs <- semiCovs(y = simDa$W_ts, 
               Coord = simDa$loc,     
               Kernel = 0,
               h = h, 
               prob = 1,
               nuUnifb = 0,
               nu = 0,
               nThreads = 15)
Cs0 <- as.matrix(Cs$Cmat)
n <- nrow(Cs0)
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")

Cs1 <- Cs0
Cs1 <- t(Cs1) + Cs1
diag(Cs1) <- 0 #Cs$Var
eigen(Cs1)

Cs2 <- Cs0
# diag(Cs2) <- Cs$Var
Var <- rep(Cs$Var, n)
Var <- as.double(Var)
Cs2 <- as.double(Cs2)
n <- as.integer(n)
a <- ModCov(Cs2, Var, n)
all.equal(a$Cov, Cs$Mcov)

M <- matrix(a$Vec, nrow = 5, ncol = 5)
M

modify.Cov(Cs$Cmat, Cs$Var)

M%*% diag(a$Eig)%*% t(M)



theta_est <- function(y_ts = simDa$Y_ts, 
                      z_ts = simDa$Z_ts,
                      Time = simDa$time,
                      Q = solve(simDa$Vc),
                      Kernel = 0,
                      h = 5e-1,
                      nuUnifb = 0,
                      nu = 0,
                      nThreads = 1){
  Nt <- ncol(y_ts)
  n <- nrow(y_ts)
  
  y <- as.vector(t(y_ts))
  Pz <- dim(z_ts)[1]
  Z <- matrix(0, nrow = n*Nt, ncol = 2*Pz)
  dt = 1
  for(j in 1:Pz){
    Z[, (Pz + j)] <-  Z[, j] <- as.vector(t(z_ts[j, , ])) 
  }
  
  storage.mode(y) <- "double"
  storage.mode(Z) <- "double"
  storage.mode(Time) <- "double"
  storage.mode(Q) <- "double"
  
  storage.mode(n) <- "integer"
  storage.mode(Nt) <- "integer"
  storage.mode(Pz) <- "integer"
  storage.mode(Kernel) <- "integer"
  storage.mode(h) <- "double"
  
  storage.mode(nuUnifb) <- "integer"
  storage.mode(nu) <- "double"
  storage.mode(nThreads) <- "integer"
  
  fit <- theta_Est_Ct(y, Z, Time, Q, n, Nt, Pz, 
               Kernel, h, nuUnifb, nu, nThreads)
  
  return(fit)
}
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
FIT <- theta_est(y_ts = simDa$Y_ts, 
                      z_ts = simDa$Z_ts,
                      Time = simDa$time,
                      Q = solve(simDa$Vc),
                      Kernel = 0,
                      h = 5e-1,
                      nuUnifb = 0,
                      nu = 0,
                      nThreads = 1)


source("./R/PSTVB_Packages.R")
source("./R/SetSimuModel.R")
source("E:/Literature/semiBase/R/util.R")
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
simDa <- siMuIncF(n = 5, Nt = 6, 
                  x.0 = c(0),
                  y.0 = c(0), 
                  delta = 0.1,
                  para = list(tau.sq = 1, 
                              Phis = 0.1, 
                              nu = 1, 
                              sigma.sq.s = 1,
                              sigma.sq.t = 1, 
                              Phit = 0.5,
                              rho = 0.5,
                              beta = c(0, 0)))

theta_est <- function(y_ts = simDa$Y_ts, 
                      z_ts = simDa$Z_ts,
                      Time = simDa$time,
                      Q = solve(simDa$Vc),
                      Kernel = 0,
                      h = 5e-1,
                      nuUnifb = 0,
                      nu = 0,
                      nThreads = 1){
  Nt <- ncol(y_ts)
  n <- nrow(y_ts)
  
  y <- as.vector(t(y_ts))
  Pz <- dim(z_ts)[1]
  Z <- matrix(0, nrow = n*Nt, ncol = 2*Pz)
  dt = 1
  for(j in 1:Pz){
    Z[, (Pz + j)] <-  Z[, j] <- as.vector(t(z_ts[j, , ])) 
  }
  
  storage.mode(y) <- "double"
  storage.mode(Z) <- "double"
  storage.mode(Time) <- "double"
  storage.mode(Q) <- "double"
  
  storage.mode(n) <- "integer"
  storage.mode(Nt) <- "integer"
  storage.mode(Pz) <- "integer"
  storage.mode(Kernel) <- "integer"
  storage.mode(h) <- "double"
  
  storage.mode(nuUnifb) <- "integer"
  storage.mode(nu) <- "double"
  storage.mode(nThreads) <- "integer"
  
  fit <- theta_Est_Ct(y, Z, Time, Q, n, Nt, Pz, 
                      Kernel, h, nuUnifb, nu, nThreads)
  
  return(fit)
}
Rcpp::sourceCpp("E:/Literature/semiBase/src/nonPara.cpp")
FIT <- theta_est(y_ts = simDa$Y_ts, 
                 z_ts = simDa$Z_ts,
                 Time = simDa$time,
                 Q = solve(simDa$Vc),
                 Kernel = 0,
                 h = 5e-1,
                 nuUnifb = 0,
                 nu = 0,
                 nThreads = 1)
FIT$alpha
# FIT$S[1, 1:30]
source("E:/Literature/semiBase/R/util.R")
fit <- theta_Cov_fun(y_ts = simDa$Y_ts, 
                     z_ts = simDa$Z_ts,
                     Time = simDa$time,
                     Q = solve(simDa$Vc),
                     h = 5e-1)
fit$alpha
FIT$S[1:5, 1:5]
fit$S[1:5, 1:5]
all.equal(FIT$S, fit$S)
