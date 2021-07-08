stSemiPar <- function(y_ts = NULL, 
                      x_ts = NULL, 
                      z_ts = NULL,
                      loc = NULL,
                      prob = c(1.5, 1.5),
                      time = NULL,
                      method = c("WI"),
                      Kernel = c(0, 0),
                      h = c(1e-1, 1e-1, 1e-1),
                      nuUnifb = 1,
                      nu = 0.8,
                      nThreads = 10,
                      nIter = 10){
  #######################################################
  # y_ts = data$Y_ts
  # x_ts = data$X_ts
  # z_ts = data$Z_ts
  # loc = data$loc
  # time = data$time
  # 
  # h = c(1e-1, 1e-1, 1e-1)
  # # method <- M[1]
  # nuUnifb = 1
  # nu = 0.8
  # nThreads = 10
  #######################################################
  #######################################################
  # source("E:/Literature/semiBase/R/util.R")
  Nt <- ncol(y_ts)
  n <- nrow(y_ts)
  y <- as.vector(y_ts)
  if("matrix" %in% class(x_ts)){
    temp <- x_ts
    x_ts <-  array(0, dim = c(1, nrow(x_ts), ncol(x_ts)),
                   dimnames = list(c("x_ts"),
                                   c(1:nrow(x_ts)),
                                   as.character(1:ncol(x_ts))              
                   ))
    x_ts[1,,] <- temp
  }
  
  if("matrix" %in% class(z_ts)){
    Z <- z_ts
    z_ts <-  array(0, dim = c(1, nrow(z_ts), ncol(z_ts)),
                   dimnames = list(c("z_ts"),
                                   c(1:nrow(z_ts)),
                                   as.character(1:ncol(z_ts))              
                   ))
    z_ts[1,,] <- Z
  }
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
  
  # beta 
  if(!is.null(Px)){
    X <- X_ts_Transf(x_ts, beta)
    fix.residuals <- y_ts - X
  }else{
    fix.residuals <- y_ts
  }
  #theta
  
  theta <- theta_WI(y_ts = fix.residuals,
                    z_ts = z_ts,
                    Time = time,
                    Q = diag(Nt), #solve(Y_ts$Vt)；diag(ncol(Y_ts$Y_ts))
                    Kernel = Kernel[1],
                    h = h[1])
  fix.semi.residuals <- fix.residuals - theta$y.fit
  fix.semi.residuals <- fix.semi.residuals - mean(fix.semi.residuals)
  curr.alpha <- theta$alpha[, 1:Pz]
  # if(!is.null(Vc)){
  #   Q <- solve(Vc)
  # }
  # DQ = gpuR::eigen(Q)
  # L = DQ$vectors %*% sqrt(diag(DQ$values)) %*% t(DQ$vectors)
  # L <- solve(t(chol(Vc)))
  
  varH <- 1
  Cs <- Ct <- NULL
  S0 <- theta$St
  iter <- 0; error <- 1
  Q_temp0 <- Q_temp
  while((iter < nIter)&(error >=1e-3) ){ #& (error >=1e-3)
    
    # 1. Update covariance matrix
    # Ct <- semiCovt(fix.semi.residuals,
    #                Time = time, 
    #                Var.s = NULL,
    #                Kernel = Kernel[1],
    #                h = c(h[3], h[3]),
    #                prob = prob[2],
    #                nuUnifb = nuUnifb,
    #                nu = nu,
    #                nThreads = nThreads)
    # Cst <- semiCovst(y = fix.semi.residuals, 
    #                Time = time, 
    #                Coord = loc, 
    #                Kernel = Kernel[2],
    #                h = h, 
    #                prob = c(0.1, 1.5), 
    #                nuUnifb = nuUnifb,
    #                nu = nu, nThreads = nThreads)
    # if(method %in% c("WI", "WEC_s")){  # only time dimension
    # Q <- Matrix::solve(Ct$ModCov$Cov)
    # if(method %in% c("WI")){
    #   Q <- diag(Nt)*diag(Q)
    # }
    # Q_temp <- kronecker(Matrix::diag(n), Q)
    # type = 0
    
    # Q_temp0 <- diag(Q_temp)*diag(n*Nt)
    # }
    if(method %in% c("WI", "WEC_t")){ 
      Q_temp0 <- diag(Nt)
    }
    
    # if(method %in% c("WEC_tw", "WEC_st", "WEC_stw")){ # space-time
    #   # h[1] <- 1E-3
    #   if(method %in% c("WEC_st", "WEC_stw")){
    #     Cs <- semiCovs(fix.semi.residuals, 
    #                    loc,     
    #                    Kernel = Kernel[2],
    #                    h = h[2], 
    #                    prob = prob[1],
    #                    nuUnifb = nuUnifb,
    #                    nu = nu,
    #                    nThreads = nThreads)
    #     # Cov.t <- Ct$ModCov$Cov  
    #     # Matrix::diag(Ct$ModCov$Cov) <- Matrix::diag(Ct$ModCov$Cov)/ Cs$Var.s
    #     cat("\n---------------\n")
    #     cat("var: ", Round(Cs$Var.s, 4))
    #     cat("\n---------------\n")
    #     Q_temp <- Q <- kronecker(Matrix::solve(Cs$Cov/Cs$Var.s), #estCov/Cs$Var.s
    #                    Matrix::solve(Ct$ModCov$Cov))#ModCov$Cov
    #   }else{
    #     Q_temp <- Q <- kronecker(Matrix::diag(n), #estCov
    #                    Matrix::solve(Ct$ModCov$Cov)) #ModCov$Cov
    #   }
    #   type = 1
    # }
    type = 1
    # Q_temp <- as.matrix(Matrix::solve(Cst$ModCov$Cov)) #ModCov$
    
    # Q_temp <- 
    #3. Update theta 
    theta <- semPar.space.time( y_ts = fix.residuals,
                                z_ts = z_ts,
                                Time = time,
                                Q = Q_temp0,
                                S = theta$St,
                                Kernel = Kernel[1],
                                h = h[1],
                                nuUnifb = nuUnifb,
                                nu = nu,
                                nThreads = nThreads,
                                method = method)
    
    
    #3. Update beta
    if(!is.null(Px)){
      if(type){
        XX <- NULL
        for(i in 1:Px){
          XX <- cbind(XX, as.vector(t(x_ts[i, , ])))
        }
        beta <- solve(t(XX) %*% Q_temp %*% XX) %*% t(XX) %*% Q_temp %*%
          as.vector(t((y_ts - theta$y.fit)))
        
        # SS <- theta$St
        # Matrix::diag(SS) <- 1 - Matrix::diag(SS)
        # # X^T %*% (I - SS)^T %*% Q %*% (I - SS)- theta$y.fit
        # X.T <- Matrix::t(SS %*% XX ) %*% Q_temp #%*% SS
        # beta <- solve(X.T %*% XX) %*% X.T %*% as.vector(t(y_ts))
        # # # beta <- as.vector(beta)
      }else{
        XX <- XY <- 0
        if(Px != 1){
          for(s in 1:n){
            XX <- XX + x_ts[, s, ] %*% Q %*% t(x_ts[, s, ])
            XY <- XY + x_ts[, s, ] %*% Q %*% (y_ts[s, ] - theta$y.fit[s, ])
          }
        }else{
          for(s in 1:n){
            XX <- XX + t(cbind(x_ts[, s, ])) %*% Q %*% cbind(x_ts[, s, ])
            XY <- XY + t(cbind(x_ts[, s, ])) %*% Q %*% (y_ts[s, ] - theta$y.fit[s, ])
          }
        }
        beta <- solve(XX) %*% XY
      }
      # beta[1:2, 1] <- c(1, 5) 
      fix.effect.fit <- X_ts_Transf(x_ts, beta)
      fix.residuals <- y_ts - fix.effect.fit
    }else{
      fix.effect.fit <- 0
      fix.residuals <- y_ts 
    }
    
    fix.semi.residuals <- fix.residuals - theta$y.fit
    fix.semi.residuals <- fix.semi.residuals - mean(fix.semi.residuals)
    # print(range(theta$alpha))
    if(!is.null(Px)){
      error <- mean(abs(curr.beta - beta)/abs(beta))
      curr.beta <- beta
    }else{
      error <- mean(abs(theta$alpha[, 1:Pz] - curr.alpha)/max(abs(curr.alpha), 1e-5))
      curr.alpha <- theta$alpha[, 1:Pz]
    }
    
    iter <- iter + 1
    cat("--------------------------------------\n")
    cat("**************************************\n")
    cat(
      " iter: ",
      iter,
      "\n",
      "beta: ",
      Round(beta, 4),
      "\n",
      # "var: ",
      # Round(Cs$Var.s, 4),
      # "\n",
      "error: ",
      Round(error, 4),
      "\n"
    )
    cat("**************************************\n")
    cat("--------------------------------------\n")
    if(iter == 1){error <- 1}
    # plot(y_ts, fix.effect.fit + theta$y.fit)
  }
  
  # if(method %in% c("WI")){
  #   C = kronecker(diag(n), diag(Ct$ModCov$Cov)*diag(Nt))
  #   Var.s <- NA
  # }
  # if(method %in% c("WEC_t", "WEC_tw")){
  #   C = kronecker(diag(n), Ct$ModCov$Cov)
  #   Var.s <- NA
  # }
  # if(method %in% c("WEC_st", "WEC_stw")){
  #   C = kronecker(Cs$Cov, Ct$ModCov$Cov)/(Cs$Var.s)^0
  #   Var.s <- Cs$Var.s
  # }
  return(list(beta = beta, theta = theta, 
              # Cs = Cs,#NULL, #
              # Ct = Ct,
              # C = C, 
              # Cst = Cst$ModCov$Cov,
              # Var.s = mean(Cst$Var.st),
              St = theta$St,
              fit.value = fix.effect.fit + theta$y.fit,
              y_ts = y_ts, fix.effect.fit = fix.effect.fit,
              time = time, nIter = iter))
} 


# if(method %in% "WDstR"){
#   W2 <- theta_stWang( y_ts = fix.residuals,
#                         z_ts = z_ts,
#                         Time = time,
#                         Q = Q, S = theta$St,
#                         Kernel = Kernel,
#                         h = h[1], nuUnifb = nuUnifb,
#                         nu = nu, nThreads = 1)
# }else{


# W1 = stWang(y_ts = fix.residuals,
#                z_ts = z_ts,
#                Time = time,
#                Q = Q,
#                theta_St = theta$St,
#                h = h[1])
# 
# all.equal(W1$alpha[, 1], W2$alpha[, 1])
# all.equal(W1$S, W2$S)
# all.equal(W1$St, W2$St)
# all.equal(W1$y.fit, W2$y.fit)

# }





# if(method == 1){
#   
#   if(Vc == 3){
#   W <- theta_stWang( y_ts = simDa$Y_ts,
#                      z_ts = simDa$Z_ts,
#                      Time = simDa$time,
#                      Q = Q,
#                      S = theta$St,
#                      Kernel = 0,
#                      h = h[1],
#                      nuUnifb = 0,
#                      nu = 0,
#                      nThreads = 1)
#   }
#   if(Vc == 2){
#     theta <- theta_est(y_ts = fix.residuals,
#                        z_ts = z_ts,
#                        Time = time,
#                        Q =  Q,
#                        Kernel = 0,
#                        h = h[1],
#                        nuUnifb = 0,
#                        nu = 0,
#                        nThreads = 10)
#   }
#   if(Vc == 1){
#     theta <- theta_fun(y_ts = fix.residuals,
#                        z_ts = z_ts,
#                        Time = time,
#                        Q = Q, #diag(Nt)
#                        h = h[1]
#       )
#   }
#   if(Vc == 0){
#     theta <- theta_fun(y_ts = fix.residuals,
#                        z_ts = z_ts,
#                        Time = time,
#                        Q = diag(Nt), #diag(Nt)
#                        h = h[1]
#     )
#   }
#   # theta <- theta_Cov_fun(y_ts = fix.residuals,
#   #               z_ts = z_ts,
#   #               Time = time,
#   #               Q = Q, #diag(Nt)
#   #               h = h[1]
#   # )
# }
# if((method == 2) & (Vc == 1)){
#   theta <- theta_Wang_fun(y_ts = fix.residuals,
#                           z_ts = z_ts,
#                           Time = time,
#                           Q = Q, #diag(Nt)
#                           S0 = theta$S,
#                           pre_theta = theta$alpha,
#                           h = h[1]
#   )
# }

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
