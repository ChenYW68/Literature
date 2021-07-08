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
f3 <- function(x)
{
  g <- function(z)
    return(exp( - (z - 1)^2) + exp( -0.8 * (z + 1)^2)
           - 0.05 * sin(8 * (z + 0.1)))
  g(x)
}
# t = seq(0, 1,, 10)
# plot(t, f2(t))
siMuIncF <- function(n = 1e2, Nt = 10, 
                     x.0 = c(0),
                     y.0 = c(0), 
                     delta = 0.1,
                     para = list(Phis = 0.3, 
                                 nu = 1, 
                                 sigma.sq.s = 1,
                                 sigma.sq.t = 1, 
                                 # Phit = 0.8,
                                 # rho = 0.1,
                                 tau.sq = 0.5, 
                                 beta = c(1, 1)),
                     nRatio = 0.8){
  x.coords <- seq(0, sqrt(n)/2, , n)#x.0 + (1:(n - 1))*delta
  y.coords <- seq(0, sqrt(n)/2, , n)#y.0 + (1:(n - 1))*delta
  
  Coords <- matrix(NA, ncol = 1, nrow = 1)
  
  while (nrow(Coords)!=n) {
    Coords <- cbind(sample(x.coords, n, replace = T), 
                    sample(y.coords, n, replace = T)
    )
    Coords <- unique(Coords)
  }
  D <- fields::rdist(Coords, Coords)
  # time <- seq(0, 1,, Nt)
  time <- (2*(1:Nt) - 1)/(2*Nt)
  
  mu <- matrix(NA, nrow = n, ncol = Nt)
  Px <- 1
  X_ts <-  array(0, dim = c(Px, nrow(Coords), Nt),
                 dimnames = list(paste0("X", 1:Px),
                                 c(1:n), as.character(1:Nt)              
                 ))
  Pz <- 1
  Z_ts <- thetaF <- array(0, dim = c(Pz, nrow(Coords), Nt),
                          dimnames = list(paste0("Z", 1:Pz), 
                                          c(1:n), as.character(1:Nt)
                                          
                          ))
  

  
  for(t in 1:Nt){
    X_ts[1, , t] <- rnorm(n, 0, 0.5)#10*time[t]^2 + rnorm(n, 0, 0.5)#rnorm(n, 0, 1)#
    if(Px > 1){
      x.phi <- 0.1*max(D)
      L <- Matrix::chol(Matern(d = D, range = x.phi, 
                               smoothness = 0.5, phi = 1))
      
      X_ts[2, , t] <- Matrix::crossprod(L, rep(rnorm(n)))
      # mvnfast::rmvn(1, rep(0, n), 
      #             sigma = exp(- D/(x.phi)),  #quantile(D, probs = 0.2)
      #             ncores = 10)
      mu[, t] <- as.vector(t(X_ts[1:2, , t]) %*% para$beta)
    }else{
      mu[, t] <-  X_ts[1, , t]*para$beta[1]
    }
   
    
    Z_ts[1, , t] <- rep(1, n)#rnorm(n, 0, 1)#rnorm(n, 0, 1)
    # Z_ts[2, , t] <- Matrix::crossprod(L, rep(rnorm(n)))
    
    thetaF[1, , t] <- Z_ts[1, , t]*f1(time[t])
    # thetaF[2, , t] <- Z_ts[2, , t]*f2(time[t])
    
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
  J <- para$J
  
  Psi <- matrix(0, nrow = Nt, ncol = J)
  
  Psi[, 1] <- cos(2*pi*time)
  c1 <- Psi[, 1] %*% Psi[, 1]
  Psi[, 1] <-  Psi[, 1]/sqrt(as.vector(c1))
  if(J > 1){
    for(j in 2:J)
    Psi[, j] <- sin(2*pi*time)
    c2 <- Psi[, j] %*% Psi[, j]
    Psi[, j] <-  Psi[, j]/sqrt(as.vector(c2)) 
  }
 
  # 
  
  # 
  # 
 
  # Psi[, 3] <- sqrt(2)*cos(4*pi*time)
  Xi <- matrix(0, nrow = n, ncol = J)
  # Vs <- list(p)
  for(i in 1:J){
    if(para$Phis[i] !=0){
      Cs <- Matern(d = D, range = para$Phis[i],
                   smoothness = para$nu[i],
                   phi = para$sigma.sq.s[i]) 
      # Vs.2 <- Matern(d = D, range = para$Phis[2],
      #                smoothness = para$nu[2],
      #                phi = para$sigma.sq.s[2]) 
      Xi[, i] <- Matrix::crossprod(Matrix::chol(Cs), rep(rnorm(n)))
    }else{
      Cs <- NULL
      Xi <- rnorm(n, 0 , sd = sqrt(para$sigma.sq.s[i]))
    }
    
    # Ws.2 <- Matrix::crossprod(Matrix::chol(Vs.2), rep(rnorm(n)))
    # if(i == 1){
    #   Wt <- rep(1, Nt)
    # }else{
    # Wt <- sqrt(2)*cos(2*pi*time)
    # }
    W_st <- W_st + tensor::tensor(as.vector(Xi[, i]), Psi[, i])# + para$nugget*rnorm(n*Nt)
  }
  # W_st <-  W_st
  Psi_mat <- build_Psi_fun(Psi, n = n, Nt = Nt)
 
  # Cs1 <- build_spCov_fun(n, J, D, range = para$Phis,
  #                        sigma.sq.s = para$sigma.sq.s,
  #                        nu = para$nu)
  
  Cs <- build_spCov(Coord = Coords, J = J, 
                    range = para$Phis,
                    sigma.sq.s = para$sigma.sq.s,
                    nu = para$nu, 
                    CovModel = 0)
  
  # Cst <- Psi_mat %*% Cs %*% t(Psi_mat) + para$nugget *diag(n*Nt)
  
  
  Qst <- Solve_Cst(Cs, Psi_mat, para$nugget)$inv.Cst
  # image(Qst)
  # Vt =tensor::tensor(as.vector(Psi[, i]), Psi[, i])
  # Vc <- kronecker(Cs, Vt)
  # diag(Vc) <-  diag(Vc) + para$nugget
  # Vc[1:10]
  # Cst[1:10]
  # all.equal(as.vector(Vc), as.vector(Cst))
  # 
  # solve(Cst)[1:10]
  # solve(Vc)[1:10]
  # Qst[1:10]
  # Cst <- 2*para$sigma.sq.s[1]*tensor::tensor(cos(2*pi*time), cos(2*pi*time)) #para$sigma.sq.s[1]*
  # if(para$Phis !=0){
  #   Vc <- kronecker(Vs, Vt)/para$sigma.sq.s[1] #tensor::tensor(Wt, Wt))
  # }else{
  #   # Vc <- kronecker(diag(n), Vt)
  #   Vc <- Vt/para$sigma.sq.s[1]
  # }
  # diag(Vc) <-  diag(Vc) + para$sigma.sq.t[1]
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
  
  y <- thetaF[1,,] + mu + W_st  + matrix(rnorm(n*Nt, sd = sqrt(para$nugget)),
                                         nrow = n, ncol = Nt)
  p.theta = 1
  theta <- matrix(0, nrow = Nt, ncol = p.theta)
  for (i in 1:p.theta) {
    if(i == 1){
      theta[, i] <- f1(time)
    }
    if(i == 2){
      theta[, i] <- f3(time)
    }
    if(i == 3){
      theta[, i] <- f1(time)
    }
  }
  # Cst = Psi_mat %*% Cs %*% Matrix::t(Psi_mat)
  # y <- mu + t(matrix(W_ts@data[["variable1"]], nrow = n, ncol = Nt)) + 
  #      Z_ts[1,,] + Z_ts[2,,] + 
  #             matrix(rnorm(Nt*n, 0, sd = sqrt(para$tau.sq)),
  #                          nrow = Nt, ncol = n)
  
  colnames(y) <- 1:Nt
  rownames(y) <- 1:n
  Yts_Xts <- list(Y_ts = y, 
                  X_ts = X_ts, 
                  Z_ts = Z_ts, 
                  time = time,
                  theta = theta,
                  thetaF = thetaF,
                  D = D,
                  loc = Coords,
                  W_ts = W_st,
                  Psi_mat = Psi_mat,
                  Vt = Psi,
                  Vs = Xi,
                  Cs = Cs,
                  # Cst = Cst,
                  Qst = Qst)
  return(Yts_Xts)
}
