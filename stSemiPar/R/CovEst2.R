CovEst <- function(t1, t2, y, coord, d, Time, h = 1e-2){
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
  for(k in 1:m)
  for(s1 in 1:n){
    for(s2 in 1:n){
    for(i in 1:Nt){
      for(j in 1:Nt){
        if(j!=i){
          
          R1[2, 1] <- (Time[i] - t1)/h
          R1[3, 1] <- (Time[j] - t2)/h
          # cat(R1[3, 1], "\n")
          A1 <- A1 + (R1%*% t(R1))*K1[i]*K2[j]*exp(-abs(dist(coord[s1,], coord[s2,]) - d)/h)
          
          A2 <- A2 + R1*K1[i]*K2[j]*y[s1, i]*y[s2, j]
          #print(round(c(i - 1, j -1,  R1[3, 1], as.vector(A2)), 3))
        }
      }
    }
    }
  }
  
}
  alpha <- solve(A1) %*% A2
  return(alpha[1])
}

Nt <- ncol(Y_ts$W_ts)
Cov <- matrix(NA, nrow = Nt, ncol = Nt)
Time <- Y_ts$time
for(t1 in 1:(Nt - 1)){
  for(t2 in (t1 + 1):Nt){
    Cov[t1, t2] <- CovEst(Time[t1], Time[t2], Y_ts$W_ts, Time = Time)
    cat("t = ", t1, "Cov = ",  Cov[t1, t2],
        "; True = ", Y_ts$Vt[t1, t2],"\n")
    # if(t == (s + 1)){
    plot(Y_ts$Vt[1:t1, ], Cov[1:t1, ])
    # }
    # points(simDa$train.sigmaMat[s, (s + 1):t],
    #      Cov[s, (s + 1):t])
  }
  
}