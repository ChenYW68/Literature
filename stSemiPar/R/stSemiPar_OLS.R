stSemi_OLS <- function(y_ts = Y_ts$Y_ts, 
                      x_ts = Y_ts$X_ts, 
                      z_ts = Y_ts$Z_ts,
                      loc = NULL,
                      Vc = NULL,
                      prob = c(1, 1),
                      time = Y_ts$time,
                      Inde = TRUE,
                      method = 1,
                      h = c(1e-1, 1e-1, 1e-1),
                      nIter = 10){

  #######################################################
  #######################################################
  # source("E:/Literature/semiBase/R/util.R")
  Nt <- ncol(y_ts)
  n <- nrow(y_ts)
  y <- as.vector(y_ts)
  Px <- dim(x_ts)[1]
  Pz <- dim(z_ts)[1]
  
  if(Pz == 1){
    Z <- as.vector(z_ts[,, 1])
    for (t in 2:Nt) {
      Z <- Matrix::bdiag(Z, as.vector(z_ts[,, t])) 
    }
  }else{
    Z <- cbind(t(z_ts[,, 1]))
    for (t in 2:Nt) {
      Z <- Matrix::bdiag(Z, t(z_ts[,, t])) 
    }
  }
 
  X <- NULL
  if(!is.null(Px)){
    for(i in 1:Px){
      X <- cbind(X, as.vector((x_ts[i, , ])))
    }
  }
  
  X.Z <- cbind(X, Z) 
  
  theta <- theta.WLS(y = y_ts, X = X.Z, Q = NULL, 
                     n = n, Px = Px, Pz = Pz)
  
  # theta$alpha
  # plot(time, theta1[, 1])
  # plot(theta1[,2], theta$alpha[, 2])
  
  
  
  # theta <- theta_fun(y_ts = fix.residuals, 
  #                    z_ts = z_ts,
  #                    Time = time,
  #                    Q = diag(Nt), #solve(Y_ts$Vt)；diag(ncol(Y_ts$Y_ts))
  #                    h = h[3])
  theta0 <- theta
  fix.semi.residuals <- y_ts - theta$y.fit #fix.residuals
  varH <- 1
  error <- 1
  iter <- 0
  curr.alpha <- theta$alpha[, 1:Pz]
  curr.beta <- theta$beta
  # if(!is.null(Vc)){
  #   Q <- solve(Vc)
  # }
  # DQ = gpuR::eigen(Q)
  # L = DQ$vectors %*% sqrt(diag(DQ$values)) %*% t(DQ$vectors)
  # L <- solve(t(chol(Vc)))
  
  
  # fit <- locpoly(x = time, y = theta$alpha[, 1], 
  #                bandwidth = h, binned = F)
  # fit <- ksmooth(x = time, y = theta$alpha[, 1],
  #                bandwidth = h, x.points = time)
  # plot(fit$x, fit$y, type = "l", col = "red")
  # lines(time, theta$alpha[, 1])
  # plot(time, fit$alpha[1, ], type = "l", col = "red")
  # lines(time, theta$alpha[, 1])
  
  
  while((iter < nIter) & (error > 1e-3)){
      # Covraiance
      Cs <- semiCovs(fix.semi.residuals, 
                     Coord = loc,     
                     Kernel = 0,
                     h = h[2], 
                     prob = prob[1],
                     nuUnifb = 1,
                     nu = 0.8,
                     nThreads = 10)
      # h[1] <- 1E-3
      Ct <- semiCovt(fix.semi.residuals,
                     Time = time, 
                     Kernel = 0,
                     h = c(h[3], h[3]),
                     prob = prob[2],
                     nuUnifb = 1,
                     nu = 0.8,
                     nThreads = 10)
      # Cov.t <- Ct$ModCov$Cov  
      # Matrix::diag(Ct$ModCov$Cov) <- Matrix::diag(Ct$ModCov$Cov)/ Cs$Var.s
      
      Q <- kronecker(Matrix::solve(Cs$Cov),
                     Matrix::solve(Ct$ModCov$Cov))#Ct$ModCov$Cov
   
    #2. Update theta
     theta <- theta.WLS(y = y_ts, X = X.Z, Q = Q, 
                         n = n, Px = Px, Pz = Pz)

    fix.semi.residuals <- y_ts - theta$y.fit #fix.residuals
    if(!is.null(Px)){
      error <- mean(abs(curr.beta - theta$beta)/abs(theta$beta))
      curr.beta <- theta$beta
    }else{
      error <- mean(abs(theta$alpha[, 1:Pz] - curr.alpha)/max(abs(curr.alpha), 1e-5))
      curr.alpha <- theta$alpha[, 1:Pz]
    }
    
    iter <- iter + 1
    
    cat("***************************************\n")
    cat(
      " iter: ",
      iter,
      "\n",
      "beta: ",
      Round(theta$beta, 4),
      "\n",
      # "varH: ",
      # Round(varH, 4),
      # "\n",
      "error: ",
      Round(error, 4),
      "\n"
    )
    cat("--------------------------------------\n")
    if(iter == 1){error <- 1}
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
  theta.loc <- NULL
  for(k in 1:Pz){
    h <- KernSmooth::dpill(x = time, y = theta$alpha[, k])
    theta.loc <- cbind(theta.loc, 
                       local_linear_kernel(y = theta$alpha[, k], 
                                           covZ = time, 
                                           covModel = 3,
                                           h = h,
                                           nu = 0, 
                                           nuUnifb = 0,
                                           nThreads = 10)$alpha[1, ])
  }
  Theta <- list()
  Theta$S <- theta$S
  Theta$alpha <- theta.loc
  return(list(beta = beta, theta = Theta,  Theta = theta,
              Cs = Cs,#NULL, #
              Ct = Ct,
              # C = kronecker(Cs$ModCov$Cov, Ct$ModCov$Cov), 
              # S = theta$S + S.beta - theta$S %*% S.beta,
              y_ts = y_ts,
              fit.value = X_ts_Transf(x_ts, theta$beta) + 
                          X_ts_Transf(z_ts, theta.loc),
              fit.value.2 = X_ts_Transf(x_ts, theta$beta) + 
                X_ts_Transf(z_ts, theta$alpha),
              fix.effect.fit = X_ts_Transf(x_ts, theta$beta),
              time = time, nIter = iter))
}