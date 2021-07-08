stSemiPar <- function(y_ts = NULL, 
                      x_ts = NULL, 
                      z_ts = NULL,
                      loc = NULL,
                      time = NULL,
                      method = c("WI"),
                      Kernel = c(0, 0),
                      h = 1e-1,
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
  # #
  # h = 1e-1
  # # # method <- M[1]
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
  # curr.beta <- beta
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

# fix.semi.residuals <- fix.semi.residuals - mean(fix.semi.residuals)
curr.alpha <- theta$alpha[, 1:Pz]
# if(!is.null(Vc)){
#   Q <- solve(Vc)
# }
# DQ = gpuR::eigen(Q)
# L = DQ$vectors %*% sqrt(diag(DQ$values)) %*% t(DQ$vectors)
# L <- solve(t(chol(Vc)))

#   Stage 2
S0 <- theta$St
beta <- beta_est(y_ts, x_ts, theta, Q = NULL)
curr.beta <- beta$beta
fix.residuals <- beta$fix.residuals
fix.semi.residuals <- fix.residuals - theta$y.fit

iter <- 0; error <- 1
# Q_temp0 <- Q_temp
J <- 2
D <- max(fields::rdist(loc, loc))
# Q0 <- Q

#   Stage 3
while((iter < nIter)&(error >=1e-3) ){ #& (error >=1e-3)
 
  # 1. Update covariance matrix
  if(method %in% c("WI")){ 
    Q <- diag(Nt)
    # type <- 0
  }
  if(method %in% c("WEC_t", "WEC_tw", "WEC_st", "WEC_stw")){ # space-time
    # source("E:/Literature/semiBase/R/util.R")
  fd <- sp.FPCA(fix.semi.residuals, time = 0:(Nt - 1),
                nharm = J, lambda = 1e3)
  # fd$pca_fd$values[1:3]
  # fd$PCA[, 1] %*% fd$PCA[, 1]
  # for(j in 1:J){
  #   c1 <- fd$PCA[, j] %*% fd$PCA[, j]
  #   fd$PCA[, j] <-  fd$PCA[, j]/sqrt(as.vector(c1))
  # }

  Psi <- build_Psi_fun(Psi = fd$PCA, n = n, Nt = nrow(fd$PCA))
  range_nugget <- exp(optim(rep(-2, J + 1), optim_semPara, 
              method = "L-BFGS-B",#L-BFGS-B
              y = fix.semi.residuals,
              Coord = loc, 
              Psi = Psi, 
              # Phi.tem = Psi.crossprod,
              J = J, sigma.sq.s = fd$pca_fd$values[1:J],
              CovModel = Kernel[2],
              lower = c(rep(log(0.01), J), log(0.1)), 
              upper = c(rep(log(D*0.2), J), log(1.0)),
              control = list(REPORT = 1,
                             trace = 1))$par)
  cat("**************************************\n")
  cat("sigma.sq: ", round(fd$pca_fd$values[1:J], 3), "\n")
  cat("range_nugget: ", round(range_nugget, 3), "\n")
  Cs <- build_spCov(Coord = loc, J = J, 
                    range = range_nugget[1:J],
                    sigma.sq.s = fd$pca_fd$values[1:J],
                    nu = rep(0.5, J), 
                    CovModel = Kernel[2])
  Q <- as.matrix(Solve_Cst(Cs, Psi, fd$pca_fd$values[J + 1])$inv.Cst)
 
  }
  type <- 1
  if(method %in% c("WEC_t")){
    method0 = "WI"
  }else{
    method0 = method
  }
  #3. Update theta 
  theta <- semPar.space.time( y_ts = fix.residuals,
                              z_ts = z_ts,
                              Time = time,
                              Q = Q,
                              S = theta$St,
                              Kernel = Kernel[1],
                              h = h[2],
                              nuUnifb = nuUnifb,
                              nu = nu,
                              nThreads = nThreads,
                              method = method0)
  #3. Update beta
  if(!is.null(Px)){
    # if(method %in% c("WI")){
    #   Q <- diag(n*Nt)
    # }
    # else{
    #   Q0 <- Q
    # }
    if(type){
      XX <- NULL
      for(i in 1:Px){
        XX <- cbind(XX, as.vector(t(x_ts[i, , ])))
      }
      # beta <- solve(t(XX) %*% Q %*% XX) %*% t(XX) %*% Q %*%
      #   as.vector(t((y_ts - theta$y.fit)))
      SS <- diag(n*Nt) - theta$St
      # SS <- theta$St
      # Matrix::diag(SS) <- 1 - Matrix::diag(SS)
      # X^T %*% (I - SS)^T %*% Q %*% (I - SS)- theta$y.fit
     
      # X.T <- Matrix::t(SS %*% XX ) %*% Q0 #%*% SS
      # beta <- solve(X.T %*% XX) %*% X.T %*% as.vector(t(y_ts - theta$y.fit))
      # 
      
      X.tild <- Matrix::t(SS %*% XX)
      beta <- solve(X.tild %*% Q %*% t(X.tild)) %*%
        X.tild %*% Q %*% SS %*% as.vector(t(y_ts))

    }
    # else{
    #   XX <- XY <- 0
    #   if(Px != 1){
    #     for(s in 1:n){
    #       XX <- XX + x_ts[, s, ] %*% Q %*% t(x_ts[, s, ])
    #       XY <- XY + x_ts[, s, ] %*% Q %*% (y_ts[s, ] - theta$y.fit[s, ])
    #     }
    #   }else{
    #     for(s in 1:n){
    #       XX <- XX + t(cbind(x_ts[, s, ])) %*% Q %*% cbind(x_ts[, s, ])
    #       XY <- XY + t(cbind(x_ts[, s, ])) %*% Q %*% (y_ts[s, ] - theta$y.fit[s, ])
    #     }
    #   }
    #   beta <- solve(XX) %*% XY
    # }
     # beta[1:2, 1] <- c(1, 5) 
    fix.effect.fit <- X_ts_Transf(x_ts, beta)
    fix.residuals <- y_ts - fix.effect.fit
  }else{
    fix.effect.fit <- 0
    fix.residuals <- y_ts 
  }

    fix.semi.residuals <- fix.residuals - theta$y.fit
    # fix.semi.residuals <- fix.semi.residuals - mean(fix.semi.residuals)
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
    # if(iter == 1){error <- 1}
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
if(method %in% c("WEC_t", "WEC_tw", "WEC_st", "WEC_stw")){
  return(list(beta = beta, theta = theta, 
              St = theta$St,
              Q = as(Q, "sparseMatrix"),
              Cs = Cs,
              Psi = Psi,
              Vt = fd$PCA,
              range_nugget = range_nugget,
              sigma.sq.s = fd$pca_fd$values,
              fit.value = fix.effect.fit + theta$y.fit,
              y_ts = y_ts, fix.effect.fit = fix.effect.fit,
              time = time, nIter = iter))
}else{
  return(list(beta = beta, theta = theta, 
              St = theta$St,
              Q = as(Q, "sparseMatrix"),
              Cs = NA,
              Psi = NA,
              Vt = NA,
              range_nugget = rep(NA, J + 1),
              sigma.sq.s = rep(NA, J),
              fit.value = fix.effect.fit + theta$y.fit,
              y_ts = y_ts, fix.effect.fit = fix.effect.fit,
              time = time, nIter = iter))
}
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
