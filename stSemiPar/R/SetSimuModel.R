f1 <- function(t){
  z = sin(t*pi) #exp(-((x- 0.5)^2 + (y- 0.5)^2))#+f2d(y)#sin(((x^2 + y^2))*pi)
  return(z)
  # return(sqrt(x^2+y^2))sin(((x^2 + y^2))*pi)+f2d(y)#
}
# t = seq(0, 1,, 10)
# plot(t, f1(t))

f2 <- function(t){
  # z = cos(2*x) + sin(6*y) +
  #   (x - pi*y + 2*pi^2/9)
  z = (1 - cos(t*pi))#cos((sqrt((x - 0.5)^2 + (y - 0.5)^2))*pi)#+f2d(x)#exp(-((x - 0.5)^2+ (y - 0.5)^2)/100)
  return(z)
}
# t = seq(0, 1,, 10)
# plot(t, f2(t))
siMuIncF <- function(n = 1e2, Nt = 10, 
                     x.0 = c(0),
                     y.0 = c(0), 
                     delta = 0.1,
                     para = list(tau.sq = 0.5, Phis = 0.3, 
                                 nu = 1, sigma.sq.s = 1,
                                 sigma.sq.t = 1, 
                                 Phit = 0.8,
                                 rho = 0.1,
                                 beta = c(1, 1)),
                     nRatio = 0.8){
  x.coords <- seq(0, 1, , n)#x.0 + (1:(n - 1))*delta
  y.coords <- seq(0, 1, , n)#y.0 + (1:(n - 1))*delta
  
  Coords <- matrix(NA, ncol = 1, nrow = 1)
  
  while (nrow(Coords)!=n) {
    Coords <- cbind(sample(x.coords, n, replace = T), 
                    sample(y.coords, n, replace = T)
    )
    Coords <- unique(Coords)
  }
  D <- fields::rdist(Coords, Coords)
  time <- seq(0, 1,, Nt)
  
  mu <- matrix(NA, nrow = n, ncol = Nt)
  X_ts <-  array(0, dim = c(2,  nrow(Coords), Nt),
                  dimnames = list(c("Intercept", "BetaX"),
                                  c(1:n),
                                  as.character(1:Nt)              
                ))
  Z_ts <- thetaF <- array(0, dim = c(2, nrow(Coords), Nt),
                        dimnames = list(c("Intercept", "ThetaZ"), 
                                        c(1:n), 
                                        as.character(1:Nt)
                                        
                        ))
  x.phi <- 0.2*max(D)
  
  L <- Matrix::chol(Matern(d = D, range = x.phi, 
                           smoothness = 0.5, phi = 1))

  for(t in 1:Nt){
    X_ts[1, , t] <- rep(n, 1)#rgamma(n, 10, 1)
   
    X_ts[2, , t] <- Matrix::crossprod(L, rep(rnorm(n)))
      # mvnfast::rmvn(1, rep(0, n), 
      #             sigma = exp(- D/(x.phi)),  #quantile(D, probs = 0.2)
      #             ncores = 10)
    mu[, t] <- as.vector(t(X_ts[1:2, , t]) %*% para$beta)
    
    Z_ts[1, , t] <- rnorm(n, 0, 1)
    Z_ts[2, , t] <- rnorm(n, 0, 1)
    
    thetaF[1, , t] <- Z_ts[1, , t]*f1(time[t])
    thetaF[2, , t] <- Z_ts[2, , t]*f2(time[t])
 
  }
  
  library(RandomFields)
  # if(para$Phis == 0){
  #   Vs <-  diag(n) 
  # }else{
  #   Vs <- Matern(d = D, range = para$Phis,
  #                smoothness = para$nu,
  #                phi = para$sigma.sq.s) 
  # }
  # 
  # ## Time covariance -- phi = 0.8
  # Vt <- diag(Nt) 
  # Vt <- para$sigma.sq.t * para$Phit^abs(row(Vt) - col(Vt))/(1 - para$Phit^2)
  # 
  # ## Cross covariance
  # Vc <- kronecker(Vs, Vt)
  # 
  # W_ts <- t(matrix(Matrix::crossprod(Matrix::chol(Vc), rep(rnorm(n * Nt))) + 
  #                    rnorm(n * Nt, 0, sd = sqrt(para$tau.sq)),
  #                  nrow = Nt, ncol = n))

  #K-L expansion
  W_st <- 0
  p <- 2
  # Vs <- list(p)
  for(i in 1:p){
    Vs <- Matern(d = D, range = para$Phis,
                 smoothness = para$nu[i],
                 phi = para$sigma.sq.s[i]) 
    # Vs.2 <- Matern(d = D, range = para$Phis[2],
    #                smoothness = para$nu[2],
    #                phi = para$sigma.sq.s[2]) 
    Ws <- Matrix::crossprod(Matrix::chol(Vs), rep(rnorm(n)))
    # Ws.2 <- Matrix::crossprod(Matrix::chol(Vs.2), rep(rnorm(n)))
    if(i == 1){
      Wt <- rep(1, Nt)
    }else{
      Wt <- sqrt(2)*cos(2*pi*(1:Nt))
    }
    
   W_st <- W_st + tensor::tensor(as.vector(Ws), Wt)
  }
  
  
  
  # model <- RMnsst(phi = RMexp(scale = para$Phis, 
  #                             var = para$sigma.sq.s),
  #                 psi = RMfbm(alpha = 1, 
  #                             var = para$sigma.sq.t), delta = 2)# +
  #           RMnugget(var = para$tau.sq)# nugget
  # model <- RMiaco(nu=1, lambda=1.5, delta=0.5) +
  #   RMnugget(var = para$tau.sq)
  
  
  # plot(model, dim = 2)
  # W_ts <- RFsimulate(model, x = Coords[, 1],
  #                y = Coords[, 2], T = 1:Nt)
  
  y <-  thetaF[1,,] + thetaF[2,,] + mu + W_st + 
        matrix(rnorm(n*Nt, sd = sqrt(para$tau.sq)), 
               nrow = n, ncol = Nt)
  
  theta <- matrix(0, nrow = Nt, ncol = 2)
  for (i in 1:2) {
    if(i == 1){
      theta[, i] <- f1(time)
    }
    if(i == 2){
      theta[, i] <- f2(time)
    }
  }
  
  
  # y <- mu + t(matrix(W_ts@data[["variable1"]], nrow = n, ncol = Nt)) + 
  #      Z_ts[1,,] + Z_ts[2,,] + 
  #             matrix(rnorm(Nt*n, 0, sd = sqrt(para$tau.sq)),
  #                          nrow = Nt, ncol = n)

  colnames(y) <- 1:Nt
  rownames(y) <- 1:n
  Yts_Xts <- list(Y_ts = y, X_ts = X_ts, 
                  Z_ts = Z_ts, 
                  time = time,
                  theta = theta, 
                  D = D,
                  loc = Coords,
                  W_ts = W_st,
                  Vt = NULL,
                  Vs = NULL,
                  Vc = NULL)
  return(Yts_Xts)
}
