CovEst <- function(t1, t2, y, Time, h = 1e-1){
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
  K1 <- exp(-abs(d1)/h)
  K2 <- exp(-abs(d2)/h)

  A1 <- A2 <- 0
  R1 <- matrix(c(1, NA, NA), ncol = 1)
  for(s in 1:n){
  for(i in 1:Nt){
    for(j in 1:Nt){
      if(j!=i){

        R1[2, 1] <- (Time[i] - t1)/h
        R1[3, 1] <- (Time[j] - t2)/h
        # cat(R1[3, 1], "\n")
        A1 <- A1 + (R1%*% t(R1))*K1[i]*K2[j]

        A2 <- A2 + R1*K1[i]*K2[j]*y[s, i]*y[s, j]
        #print(round(c(i - 1, j -1,  R1[3, 1], as.vector(A2)), 3))
      }
    }
  }
  }
  alpha <- solve(A1) %*% A2
  return(alpha[1])
}

Nt <- ncol(simDa$W_ts)
Cov <- matrix(NA, nrow = Nt, ncol = Nt)
Time <- simDa$time
for(t1 in 1:(Nt - 1)){
  for(t2 in (t1 + 1):Nt){
    Cov[t1, t2] <- CovEst(Time[t1], Time[t2], 
                          simDa$W_ts, Time = Time,
                          h = 5e-1)
    cat("t = ", t1, "Cov = ",  Cov[t1, t2],
        "; True = ", simDa$Vt[t1, t2],"\n")
    # if(t == (s + 1)){
      plot(simDa$Vt[1:t1, ], Cov[1:t1, ])
    # }
    # points(simDa$train.sigmaMat[s, (s + 1):t],
    #      Cov[s, (s + 1):t])
  }

}


all.equal(Cov1$Cmat[upper.tri(Cov1$Cmat)], Cov[upper.tri(Cov)])

















theta_Wang_fun <- function(y_ts = Y_ts$Y_ts, 
                           z_ts = Y_ts$Z_ts,
                           Time = Y_ts$time,
                           Q = diag(ncol(Y_ts$Y_ts)),
                           S0 = S,
                           pre_theta = NULL,
                           h = 5e-1){
  Nt <- ncol(y_ts)
  n <- nrow(y_ts)
  Pz <- dim(z_ts)[1]
  alpha <- matrix(NA, nrow = Nt, ncol = Pz)
  Y <- as.vector(t(y_ts))
  k <- 1
  S <- matrix(NA, nrow = n*Nt, ncol = n*Nt)
  y.fit <- matrix(0, nrow = n, ncol = Nt)
  if(is.null(pre_theta)){pre_theta = matrix(1, nrow = Nt, ncol = Pz)}
  
  for (t in 1:Nt) {
    S0.s <-  S2.s <- 0
    dt <- (Time - Time[t])
    K <- (exp(-abs(dt)/h))
    A1.s <- A2.s <-NULL
    
    A1.s <- matrix(NA, nrow = 2*Pz, ncol = Nt*n)
    for(s in 1:n){
      # Z.s <- NULL
      
      ys <- y_ts[s, ]
      index <- NULL
      for(i in 1:(Nt)){
        index <- c(index, s + n*(i - 1))
      }
      Sc <- S0[index,]   #((s - 1)*Nt + 1):(s*Nt)
      # Sc %*% Y
      # plot(Sc %*% Y, ys)
      # y.fit <- matrix(S0 %*% Y, nrow = n, ncol = Nt)
      # plot(t(y.fit), Y)
      
      S1.s <- 0 
      for(j in 1:Nt){
        Z.s <- matrix(0, nrow = 2*Pz, ncol = Nt)
        Z.s[1:Pz, j] <- z_ts[1:Pz, s, j]
        Z.s[(Pz + 1):(2*Pz), j] <-  dt[j]*z_ts[1:Pz, s, j]/h
        
        # 1
        # lam <- 0
        # for(i in 1:Pz){
        #   lam <- lam + pre_theta[, i]*z_ts[i, s, ] 
        # }
        # Yj <- y - lam
        # Yj[j] <- y[j]
        
        Sj <- Sc
        Sj[j, ] <- 0
        c0 <- K[j]*Z.s %*% Q
        S0.s <- S0.s + c0 %*% Sj
        S1.s <- S1.s + c0 #Sj %*% Y
        # S1.s <- S1.s + Z.s %*% K^(0.5)%*% Q %*% K^(0.5) %*% Yj
        # 2
        # S2.s <- S2.s + Z.s %*% K^(0.5)%*% Q %*% K^(0.5) %*% t(Z.s)
        S2.s <- S2.s + c0 %*% t(Z.s)
        
        # C0[, ((l - 1)*Nt + 1):(l*Nt)] <- c0 
        # l <- l + 1
      }
      # sec.1-y
      A1.s[, ((s - 1)*Nt + 1):(s*Nt)] <- S1.s 
      
      # A1.s <- cbind(A1.s, S1.s) 
      
    }
    
    A <- solve(S2.s) %*% (A1.s - S0.s)
    # y_s(t) = h_s(t)
    
    
    for (j in 1:n) {
      S[k, ] <- z_ts[1, j, t]*A[1, ] + z_ts[2, j, t]*A[2, ]
      k <- k + 1;
    }
    alpha[t, ] <- (A %*% Y)[1:Pz]
    # y.fit[, t] <-  z_ts[1,, t]*alpha[t, 1] + z_ts[2,, t]*alpha[t, 2]
    for(j in 1:Pz){
      y.fit[, t] <- y.fit[, t] + z_ts[j,, t]*alpha[t, j]
    }
    
    # y.fit1[, t] <- z_ts[1,,t] + z_ts[2,,t]
  }
  y.fit1 <- matrix(S %*% Y, nrow = n, ncol = Nt)
  all.equal(as.vector((y.fit1)), as.vector(y.fit))
  
  return(list(S = S, alpha = alpha, 
              y.fit = y.fit))
}


theta_fun <- function(y_ts = simDa$Y_ts, 
                      z_ts = simDa$Z_ts,
                      Time = simDa$time,
                      Q = diag(ncol(simDa$Y_ts)),
                      h = 5e-1){
    Nt <- ncol(y_ts)
    n <- nrow(y_ts)
    Pz <- dim(z_ts)[1]
    alpha <- matrix(NA, nrow = Nt, ncol = Pz)
    Y <- as.vector(t(y_ts))
    k <- 1
    S <- matrix(NA, nrow = n*Nt, ncol = n*Nt)
    y.fit <- matrix(0, nrow = n, ncol = Nt)
   for (t in 1:Nt) {
     S2.s <- 0
     S11.s <- NULL
     dt <- (Time - Time[t])
     for(s in 1:n){
       # Z.s <- NULL
       # for(i in 1:Pz){
       #   Z.s <- cbind(Z.s, z_ts[i, s, ])
       # }
       # for(i in 1:Pz){
       #   Z.s <- cbind(Z.s, dt*(z_ts[i, s, ])/h)
       # }
       
       Z.s <- matrix(0, nrow = Nt, ncol = 2*Pz)
       Z.s[, 1:Pz] <- t(z_ts[1:Pz, s, ])
       Z.s[,(Pz + 1):(2*Pz)] <-  dt*t(z_ts[1:Pz, s, ])/h
     
     
       # all.equal(as.vector(Z.s1[, 1:4]), as.vector(Z.s))
     
       K <- diag(exp(-abs(dt)/h))
       S1.s <- t(Z.s) %*% K^(0.5)%*% Q %*% K^(0.5)
       S11.s <- cbind(S11.s, S1.s)
       S2.s <- S2.s + (S1.s %*% Z.s)
     }
     
     A <- solve(S2.s) %*% S11.s
     # y_s(t) = h_s(t)Y
     
     
     for (j in 1:n) { 
       S[k, ] <- z_ts[1, j, t]*A[1, ] + z_ts[2, j, t]*A[2, ]
       k <- k + 1;
     }
     alpha[t, ] <- (A %*% Y)[1:Pz]
     # y.fit[, t] <-  z_ts[1,, t]*alpha[t, 1] + z_ts[2,, t]*alpha[t, 2]
     for(j in 1:Pz){
       y.fit[, t] <- y.fit[, t] + z_ts[j,, t]*alpha[t, j]
     }
     
     # y.fit1[, t] <- z_ts[1,,t] + z_ts[2,,t]
   }
   y.fit1 <- matrix(S %*% Y, nrow = n, ncol = Nt)
   all.equal(((y.fit1)), (y.fit))
   # plot(as.vector(t(y.fit)), (Y))
  return(list(S = S, alpha = alpha, 
              Y = Y, y.fit = y.fit))
}

theta_Cov_fun <- function(y_ts = simDa$Y_ts, 
                          z_ts = simDa$Z_ts,
                          Time = simDa$time,
                          Q = solve(simDa$Vc),
                          h = 5e-1){
  Nt <- ncol(y_ts)
  n <- nrow(y_ts)
  Pz <- dim(z_ts)[1]
  alpha <- matrix(NA, nrow = Nt, ncol = Pz)
  Y <- as.vector(t(y_ts))
  
  k <- 1
  S <- matrix(0, nrow = n*Nt, ncol = n*Nt)
  y.fit <- matrix(0, nrow = n, ncol = Nt)
  Q <- as(Q, "sparseMatrix")
  for (t in 1:Nt) {
    dt <- rep((Time - Time[t]), times = n)
    Z.s <- matrix(0, nrow = n*Nt, ncol = 2*Pz)
    for(j in 1:Pz){
      Z.s[, j] <- as.vector(t(z_ts[j, , ])) 
      Z.s[, (Pz + j)] <- dt*as.vector(t(z_ts[j, , ]))/h
    }
    K <- as(diag(exp(-abs(dt)/h)), "sparseMatrix")
    K_0 <- K %*% Q %*% K
    St <- (solve(t(Z.s)%*% K_0 %*% Z.s) %*% t(Z.s)%*%  
             K_0)[1:2, ]#%*% Y 
    
    # for(k in 1:Pz){
    #   j <- j + 1
    #   S[j, ] <- S[j, ] + as.vector(t(z_ts[k, ,]))*St[k, ]
    # }
    for (s in 1:n) { 
      for(l in 1:Pz){
        S[k, ] <- S[k, ] + z_ts[l, s, t]*St[l, ]# + z_ts[2, j, t]*A[2, ]
      }
      k <- k + 1;
    }
    
    alpha[t, ] <- (St %*% Y)[1:Pz]
    # y.fit[, t] <-  z_ts[1,, t]*alpha[t, 1] + z_ts[2,, t]*alpha[t, 2]
    for(l in 1:Pz){
      y.fit[, t] <- y.fit[, t] + z_ts[l,, t]*alpha[t, l]
    } 
  }
  
  y.fit1 <- matrix(S %*% Y, nrow = n, ncol = Nt)
   # all.equal(((y.fit1)), (y.fit))
  plot(as.vector(t(y.fit1)), (Y))


  return(list(S = S, alpha = alpha, 
              Y = Y, y.fit = y.fit))
}



theta <- theta_fun(y_ts = Y_ts$Y_ts, 
                  z_ts = Y_ts$Z_ts,
                  Time = Y_ts$time,
                  Q = diag(ncol(Y_ts$Y_ts)) #solve(Y_ts$Vt)；diag(ncol(Y_ts$Y_ts))
                  )
par(mfrow = c(2, 2))
plot(Y_ts$time,Y_ts$theta1, 
     ylim = c(min(Y_ts$theta1, theta$alpha[, 1]), 
              max(Y_ts$theta1, theta$alpha[, 1])))
lines(Y_ts$time, theta$alpha[, 1], col = "red")
plot(Y_ts$time, Y_ts$theta2, 
     ylim = c(min(Y_ts$theta2, theta$alpha[, 2]), 
              max(Y_ts$theta2, theta$alpha[, 2])))
lines(Y_ts$time, theta$alpha[, 2], col = "red")
plot(theta$Y, t(theta$y.fit))

spT.validation(Y_ts$theta1, theta$alpha[, 1])
spT.validation(Y_ts$theta2, theta$alpha[, 2])
