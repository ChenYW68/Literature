stSemiPar <- function(y_ts = Y_ts$Y_ts, 
                      x_ts = Y_ts$X_ts, 
                      z_ts = Y_ts$Z_ts,
                      loc = NULL,
                      Vc = NULL,
                      prob = 1,
                      time = Y_ts$time,
                      Inde = TRUE,
                      method = 1,
                      h = c(1e-1, 1e-1, 1e-1),
                      nIter = 10){
#######################################################
#######################################################
X_ts_Transf <- function(Nt, X_ts, beta)
{
  t = seq_len(Nt)
  x_ts = sapply(t, function(t) t(X_ts[, , t]) %*% beta
                , simplify = "matrix")
  return(x_ts)
}
#######################################################
#######################################################
# source("E:/Literature/semiBase/R/util.R")
Nt <- ncol(y_ts)
n <- nrow(y_ts)
y <- as.vector(y_ts)
Px <- dim(x_ts)[1]
Pz <- dim(z_ts)[1]
lm.da <- data.frame(y = y)
if(!is.null(Px)){
  for (i in 1:Px) {
    lm.da <- cbind(lm.da, as.vector(x_ts[i,,])) 
  }
  colnames(lm.da) <- c("y", paste0("X", 1:Px))
  fmla <- as.formula(paste("y ~ -1 + ", paste(paste0("X", 1:Px),
                                          collapse = "+")))
  fit <- lm(fmla, data = lm.da)
  beta <- fit$coefficients
  curr.beta <- beta
  print(beta)
}else{
  beta <- 0
}

iter <- 0
# beta 
if(!is.null(Px)){
  X <- X_ts_Transf(Nt, x_ts, beta)
  fix.residuals <- y_ts - X
}else{
  fix.residuals <- y_ts
}
#theta
theta <- theta_fun(y_ts = fix.residuals, 
                   z_ts = z_ts,
                   Time = time,
                   Q = diag(Nt), #solve(Y_ts$Vt)；diag(ncol(Y_ts$Y_ts))
                   h = h[3])
theta0 <- theta
fix.semi.residuals <- fix.residuals - theta$y.fit
varH <- 1
error <- 1
curr.alpha <- theta$alpha[, 1:2]
if(!is.null(Vc)){
  Q <- solve(Vc)
}
# DQ = gpuR::eigen(Q)
# L = DQ$vectors %*% sqrt(diag(DQ$values)) %*% t(DQ$vectors)
# L <- solve(t(chol(Vc)))
while((iter < nIter) & (error > 1e-3)){
  S <- theta$S
  #1. Update covariance matrix
  # if(!Inde){
  varZ <- time
  residual.sq <- as.vector(colMeans(fix.semi.residuals)^2)
  # varH <- KernSmooth::dpill(varZ, residual.sq)
  if(is.null(Vc)){
  # if(is.na(varH)){
   varH <- lokern::glkerns(varZ, residual.sq)$bandwidth
  # }

  # Cov = semiCov(y = fix.semi.residuals, 
  #                Time = time, 
  #                Kernel = 0,
  #                h = c(h[1], varH),
  #                nuUnifb = 0,
  #                nu = 0,
  #                nThreads = 10)
  # Q <- solve(Cov$Cmat) #Cmatsolve(Y_ts$Vt)#diag(Nt)#solve(Cov$Cmat)
  
  
  
  # }else{
  #   Q <- diag(Nt)
  # }
   
  Cs <- semiCovs(fix.semi.residuals, 
                 loc,     
                 Kernel = 0,
                 h = h[2], 
                 prob = prob,
                 nuUnifb = 1,
                 nu = 0.8,
                 nThreads = 10)
  # h[1] <- 1E-3
  Ct <- semiCovt(fix.semi.residuals,
                 Time = time, 
                 Kernel = 0,
                 h = c(h[3], varH),
                 nuUnifb = 1,
                 nu = 0.8,
                 nThreads = 10)
    # Cov.t <- Ct$ModCov$Cov  
  # Matrix::diag(Ct$ModCov$Cov) <- Matrix::diag(Ct$ModCov$Cov)/ Cs$Var.s
  
   Q <- kronecker(Matrix::solve(Cs$ModCov$Cov),
                  Matrix::solve(Ct$ModCov$Cov))
   #Q = kronecker(diag(n), solve(Ct$Cmat))
  }
    
  if(Inde){
    Q <- as.vector(diag(Q))*diag(Nt) 
  }
  
    
  #2. Update beta
  if(!is.null(Px)){
    # if(is.null(Vc)){
    # XX <- XY <- 0
    # for(s in 1:n){
    #  XX <- XX + x_ts[, s, ] %*% Q %*% t(x_ts[, s, ])
    #  XY <- XY + x_ts[, s, ] %*% Q %*% (y_ts[s, ] - theta$y.fit[s, ])
    # }
    # beta <- solve(XX) %*% XY
    # }else{
      XX <- NULL
      for(i in 1:Px){
        XX <- cbind(XX, as.vector(t(x_ts[i, , ])))
      }
      beta <- solve(t(XX) %*% Q %*% XX) %*% t(XX) %*% Q %*% as.vector(t((y_ts - theta$y.fit)))
      # beta <- as.vector(beta)
    # }
    # beta[1:2, 1] <- c(1, 5) 
    fix.effect.fit <- X_ts_Transf(Nt, x_ts, beta)
    fix.residuals <- y_ts - fix.effect.fit
  }else{
    fix.effect.fit <- 0
  }
    
    # DQ = eigen(Q)
    # L = DQ$vectors %*% (diag(DQ$values))^(0.5) %*% t(DQ$vectors)
    # theta <- theta_fun(y_ts = fix.residuals, 
    #                    z_ts = z_ts,
    #                    Time = time,
    #                    Q = diag(Nt), #solve(Y_ts$Vt)；diag(ncol(Y_ts$Y_ts))
    #                    h = h[2]
    #                    )
    # plot(fix.residuals, theta$y.fit)
    # Z.r <- t(matrix(L %*% (as.vector(t(fix.residuals - theta$y.fit))),
    #         nrow = Nt, ncol = n)) + theta$y.fit
    
    
  #3. Update theta
  if(method == 1){
    # theta <- theta_fun(y_ts = fix.residuals,
    #                    z_ts = z_ts,
    #                    Time = time,
    #                    Q = Q, #diag(Nt)
    #                    h = h[2]
    #                   )
    theta <- theta_est(y_ts = fix.residuals, 
               z_ts = z_ts,
               Time = time,
               Q =  Q,
               Kernel = 0,
               h = h[1],
               nuUnifb = 0,
               nu = 0,
               nThreads = 10)
    
    # theta <- theta_Cov_fun(y_ts = fix.residuals,
    #               z_ts = z_ts,
    #               Time = time,
    #               Q = Q, #diag(Nt)
    #               h = h[1]
    # )
    
    # theta <- theta_fun(y_ts = Z.r, 
    #                    z_ts = z_ts,
    #                    Time = time,
    #                    Q = diag(Nt), #solve(Y_ts$Vt)；diag(ncol(Y_ts$Y_ts))
    #                    h = h[2]*(n)^(1e-1))#*(n*Nt)^(1e-2)
    
  }else{
    theta <- theta_Wang_fun(y_ts = fix.residuals,
                            z_ts = z_ts,
                            Time = time,
                            Q = Q, #diag(Nt)
                            S0 = S,
                            pre_theta = theta$alpha,
                            h = h[1]
    )
  }
  
  fix.semi.residuals <- fix.residuals - theta$y.fit
  if(!is.null(Px)){
    error <- mean(abs(curr.beta - beta)/abs(beta))
    curr.beta <- beta
  }else{
    error <- mean(abs(theta$alpha[, 1:2] - curr.alpha)/max(abs(curr.alpha), 1e-5))
    curr.alpha <- theta$alpha[, 1:2]
  }
 
  iter <- iter + 1
  cat("***************************************\n")
  cat(
    " iter: ",
    iter,
    "\n",
    "beta: ",
    Round(beta, 4),
    "\n",
    "varH: ",
    Round(varH, 4),
    "\n",
    "error: ",
    Round(error, 4),
    "\n"
  )
  cat("--------------------------------------\n")
  # plot(y_ts, fix.effect.fit + theta$y.fit)
}
# XY <- matrix(NA, nrow = Px, ncol = n*Nt)
# XX <- 0
# for(s in 1:n){
#   XX <- XX + (x_ts[, s, ] %*% Q %*% t(x_ts[, s, ]))
#   XY[, ((s - 1)*Nt + 1):(s*Nt)] <- x_ts[, s, ] %*% Q
# }
# A = solve(XX) %*% XY
# S.beta <- matrix(0, nrow = n*Nt, ncol = n*Nt)
# k <- 1
# for(t in 1:Nt){
#   for (s in 1:n){
#     for(l in 1:Px){
#       S.beta[k, ] <- S.beta[k, ] + x_ts[l, s, t]*A[l, ]
#     }
#     k <- k + 1; 
#   }
# }
# fit <- matrix(S.beta %*% as.vector(t(y_ts)), nrow = n, ncol = Nt)
# plot(fit, y_ts)
# plot(fix.effect.fit, y_ts)
# all.equal(fix.effect.fit, fit)
# 

return(list(beta = beta, theta = theta, #C = kronecker(Cs$ModCov$Cov,
                                                     # Ct$ModCov$Cov), 
            # S = theta$S + S.beta - theta$S %*% S.beta,
            fit.value = fix.effect.fit + theta$y.fit,
            y_ts = y_ts, fix.effect.fit = fix.effect.fit,
            time = time, nIter = iter))
}