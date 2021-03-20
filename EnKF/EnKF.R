library(mvtnorm)


#The R functions for the Ensemble Kalman filter (EnKF)
#can be downloaded from: 
source("http://www.datall-analyse.nl/R/EnKF.R")
#Have a look at the function EnKF, and notice that at the beginning of
#the script you will find a description of the function's arguments.
EnKF



##EXAMPLE 1

#In the following example, the EnKF will be used to simulate heat
#conduction in a one-dimensional bar. The heat conduction in the
#bar is governed by the partial differential equation (PDE):
# d T(x,t) / d t = k * d^2 T(x,t) / d x^2 + u(x,t) + w(x,t),
#where T(x,t) is the temperature at position x and time t, k the heat conduction
#coefficient, u(x,t) an external heat source at position x and time t,
#and w(x,t) represents Gaussian noise at position x and time t.
#We will employ a numerical method (more specifically, the "explicit method")
#for solving the above PDE. In applying this numerical method, we will devide
#the bar into 100 equally spaced nodes, which yields 100 temperature states.
#We will assume that the boundary conditions are T(0,t)=100 and T(L,t)=50 (where
#L is the length of the bar). The intial condition is set to 10 for the interior
#nodes, and 100 and 50 for the nodes at x=0 and x=L, respectively.
#Finally, note that in the example below we suppose that u(x,t) is
#zero (i.e., no heat sources act on the bar). The process noise (covariance),
#w(x,t), is set to 5e-4*I (where I is a 100x100 identity matrix).
#
#Measurements of the temperature are taken at nodes 10, 20, ..., 90. The
#noise (covariance) of these measurements is 5e-1*I (where I is
#a 9x9 identity matrix).



##TRUE STATES
#Generate the true states T(x,t).
dt <- 1/100 #time step size used in the explicit method
tt <- 1000 #upper bound of time window (that is, max number of time steps)
st <- seq(0, tt, by=dt) #lower time bounds of the integration intervals
ns <- length(st) #number of integrations
nc <- 100 #number of nodes
dx <- 1 #spacing between nodes
x <- matrix(0, ncol=nc, nrow=ns) #specify matrix for states
x[1,] <- c(100, rep(10, nc-2), 50) + rnorm(nc, 0, 1) #initial conditions (at t=0)

#Specify noises (variances)
wk <- 5e-4 #process noise
vk <- 5e-1 #measurement noise

#Parameters in the PDE
k <- 5 #heat conduction coefficient

#Construct the state transition matrix
A <- matrix(0, ncol=nc, nrow=nc)

#Entries of the interior nodes in the transition matrix
for(i in 2:(nc-1)) {
  A[i, (i-1):(i+1)] <- c(1, -2, 1)
}

#Taking into account the spacing between the nodes (=dx) and
#the heat conduction coefficient (=k)
A <- (dx)^2*k*A

#The explicit method yields stable and convergent results if (k*dt)/(dx)^2 < .5
(k*dt)/(dx)^2

#Simulate true states
for (i in 2:ns) {
  #states
  x[i, ] <- x[i-1, ] + (A%*%x[i-1, ])*dt + rnorm(nc, 0, sqrt(wk))}



#Plot of simulated states at specific times.
#The times are computed as k*dt, where k is the time step index
#and dt the size of the time steps.
#However, note that x[1, ] represents the temperature states at t=0. Thus, for
#t=5 we need to select x[501, ].
par(mfrow=c(2, 2))
plot(1:nc, x[501,], type="l", xlab="Node, i", ylab="Temperature (i)",
     main="t=5", ylim=c(0,120))
plot(1:nc, x[1001,], type="l", xlab="Node, i", ylab="Temperature (i)",
     main="t=10", ylim=c(0,120))
plot(1:nc, x[10001,], type="l", xlab="Node, i", ylab="Temperature (i)",
     main="t=100", ylim=c(0,120))
plot(1:nc, x[100001,], type="l", xlab="Node, i", ylab="Temperature (i)",
     main="t=1000", ylim=c(0,120))
par(mfrow=c(1, 1))



##MEASUREMENTS
#In this example we will create measurements for the first 500 time steps.
#More specifically, we assume that we take measurements from t=.01 to t=5.
#Again, let k be the time step index and dt the size of the time steps.
#In that way, we may compute the time (t) as k*dt. Thus, for k=1 and dt=1/100
#the corresponding time is 1*(1/100)=.01

#Extract the temperature states for t=.01 to t=5.
#Note that x[1, ] represents the temperature states at t=0. Thus, for t=.01
#to t=5 we need to select x[2:501, ]
x <- x[2:501,]
#Select temperatures at nodes 10, 20, ..., 90
np <- length(seq(10,90,10)) #number of nodes
dataEx1 <- x[, seq(10, 90, 10)]

#Add the measurement noise to the selected states
for (i in 1:500) {
  dataEx1[i, ] <- dataEx1[i, ] + rnorm(np, 0, sqrt(vk))}

#Plot the generated measurements at t=5
plot(1:nc, x[500,], type="l", xlab="Node, i", ylab="Temperature (i)",
     ylim=c(0,100), col=gray(level=.8))
points(seq(10,90,10), dataEx1[500,], col="darkgreen", cex=.5)
legend("topright", pch=c(NA, 1), lty=c(1, NA),
       col=c(gray(level=.8), "darkgreen"),
       legend=c("true", "measured"),
       bty="n", y.intersp=1.2, cex=.7)



##ENSEMBLE KALMAN FILTER (EnKF)

#Dynamic model: specifying 100 temperature states.
#Note that you may change the specified initial state estimates (at t=0)
#below for the 100 states and see how such changes influence the behavior
#of the ensemble Kalman filter.
ex1 <- list(m0=c(100, rep(10, nc-2), 50), #initial state estimates
            #error covariances of the initial state estimates:
            #this matrix reflects the uncertainty in our initial state estimates
            #you may change the values in this matrix and see how they influence
            #the behavior of the Kalman filter
            C0=diag(rep(1e+1, nc)),
            #measurement noise
            V=diag(rep(vk, 9)),
            #process noise
            W=diag(rep(wk, nc)))

#Specify the state transition function:
#WARNING: always use arguments x and k when specifying the GGfunction
GGfunction <- function (x, k){
  x + (A%*%x)*dt}

#Specify the observation/measurement function:
#WARNING: always use arguments x and k when specifying the FFfunction
FFfunction <- function (x, k){
  x[seq(10, 90, 10)]}



##Compute the filtered (a posteriori) state estimates with the EnKF,
#and employ 10 ensemble members in the EnKF 
enkf1 <- EnKF(y=dataEx1, mod=ex1, size=10,
              GGfunction=GGfunction, FFfunction=FFfunction)

#As a comparison, increase size of the ensemble to 20
enkf2 <- EnKF(y=dataEx1, mod=ex1, size=20,
              GGfunction=GGfunction, FFfunction=FFfunction)

#And, finally, an EnKF with 100 ensemble members
enkf3 <- EnKF(y=dataEx1, mod=ex1, size=100,
              GGfunction=GGfunction, FFfunction=FFfunction)


#Plot the filtered state estimates at t=5
plot(1:nc, x[500,], type="l", col=c(gray(level=.5)),
     ylim=range(c(enkf1$m[501,], enkf2$m[501,], enkf3$m[501,])),
     xlab="Node, i", ylab="Temperature (i)", main="t=5")
lines(1:nc, enkf1$m[501,], lty=2, col="blue", lwd=1)
lines(1:nc, enkf2$m[501,], lty=2, col="red", lwd=1)
lines(1:nc, enkf3$m[501,], lty=2, col="darkgreen", lwd=1)
legend("topright", lty=c(1, 2, 2, 2),
       col=c(gray(level=.5), "blue", "red", "darkgreen"),
       legend=c("true state", "EnKf with 10 members",
                "EnKf with 20 members", "EnKf with 100 members"),
       bty="n", y.intersp=1.2, cex=.7)



#Error plot
e1 <- sapply(1:nrow(x), function (i) mean((x[i,]-enkf1$m[i+1,])^2))
e2 <- sapply(1:nrow(x), function (i) mean((x[i,]-enkf2$m[i+1,])^2))
e3 <- sapply(1:nrow(x), function (i) mean((x[i,]-enkf3$m[i+1,])^2))
rangeError <- max(cbind(e1, e2, e3))
plot(1:nrow(x), e1, type="l", lty=1, col="blue",
     ylim=c(0, rangeError), ylab="Mean Squared Error", xlab="Time index, k")
lines(1:nrow(x), e2, lty=1, col="red")
lines(1:nrow(x), e3, lty=1, col="darkgreen")
legend("topright", lty=c(2, 2, 2),
       col=c("blue", "red", "darkgreen"),
       legend=c("EnKf with 10 members",
                "EnKf with 20 members", "EnKf with 100 members"),
       bty="n", y.intersp=1.2, cex=.7)



#Compute the confidence intervals for the filtered state estimates of the EnKF

#95% confidence intervals at t=5

#EnKF with 10 members
seFX10 <- sqrt(enkf1$C[[501]])
ciFX10 <- enkf1$m[501,] + qnorm(.05/2)*seFX10%o%c(1, -1)

#plot
plot(1:nc, x[500,], type="l", col=c(gray(level=.7)), lwd=1.5,
     xlab="Node, i", ylab="Temperature (i)",
     main="EnKF with 10 members", ylim=c(0,100))
lines(1:nc, ciFX10[,1], lty=2, col="blue")
lines(1:nc, ciFX10[,2], lty=2, col="blue")
legend("topright", lty=c(1, 2), 
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "filtered state"),
       bty="n", y.intersp=1.2, cex=.7)



#EnKF with 20 members
seFX20 <- sqrt(enkf2$C[[501]])
ciFX20 <- enkf2$m[501,] + qnorm(.05/2)*seFX20%o%c(1, -1)

#plot
plot(1:nc, x[500,], type="l", col=c(gray(level=.7)), lwd=1.5,
     xlab="Node, i", ylab="Temperature (i)",
     main="EnKF with 20 members", ylim=c(0,100))
lines(1:nc, ciFX20[,1], lty=2, col="blue")
lines(1:nc, ciFX20[,2], lty=2, col="blue")
legend("topright", lty=c(1, 2), 
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "filtered state"),
       bty="n", y.intersp=1.2, cex=.7)



#EnKF with 100 members
seFX100 <- sqrt(enkf3$C[[501]])
ciFX100 <- enkf3$m[501,] + qnorm(.05/2)*seFX100%o%c(1, -1)

#plot
plot(1:nc, x[500,], type="l", col=c(gray(level=.7)), lwd=1.5,
     xlab="Node, i", ylab="Temperature (i)",
     main="EnKF with 100 members", ylim=c(0,100))
lines(1:nc, ciFX100[,1], lty=2, col="blue")
lines(1:nc, ciFX100[,2], lty=2, col="blue")
legend("topright", lty=c(1, 2), 
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "filtered state"),
       bty="n", y.intersp=1.2, cex=.7)





##EXAMPLE 2
#Similar to Example 1, but in this example we will apply two external heat
#sources. The two heat sources act on the one-dimensional at nodes 33 and 67.
#These two heat sources are active from time steps 1 to 400, while at time step
#401 till the end of our heat conduction experiment they are turned off.
#In this example, the boundary conditions are T(0,t)=200 and T(L,t)=150.
#The intial condition is set to 100 for the interior nodes, and 200 and
#150 for the nodes at x=0 and x=L, respectively.
#The process noise (covariance) is set to 1e-2*I (where I is a 100x100
#identity matrix).
#
#Measurements of the temperature are taken at nodes 10, 20, ..., 90. The
#noise (covariance) of these measurements is 1*I (where I is
#a 9x9 identity matrix).



##TRUE STATES

#Generate the true states T(x,t).
dt <- 1/500 #time step size used in the explicit method
tt <- 100 #upper bound of time window (that is, max number of time steps)
st <- seq(0, tt, by=dt) #lower time bounds of the integration intervals
ns <- length(st) #number of integrations
nc <- 100 #number of nodes
x <- matrix(0, ncol=nc, nrow=ns) #specify matrix for states
dx <- 1 #spacing between nodes
x[1,] <- c(200, rep(100, nc-2), 150) + rnorm(nc, 0, 1) #initial conditions (at t=0)

#Specify noises (variances)
wk <- 1e-2 #process noise
vk <- 1 #measurement noise

#Parameters
k <- 5 #heat conduction coefficient

#State transition matrix
A <- matrix(0, ncol=nc, nrow=nc)

#Interior nodes
for(i in 2:(nc-1)) {
  A[i, (i-1):(i+1)] <- c(1, -2, 1)
}

#Taking into account the spacing between the nodes (=dx) and
#the heat conduction coefficient (=k)
A <- (dx)^2*k*A

#Stability and convergence test for the explicit method
(k*dt)/(dx)^2

#Control input matrix
B <- matrix(0, ncol=2, nrow=100)
B[33, 1] <- B[67, 2] <- 1 #heat sources act at nodes 33 and 67

#Heat source
timeHS <- 400 #time steps during which the heat sources act on the bar
hs <- c(0, c(rep(1, timeHS), rep(0, length(st)-timeHS)))

#Simulate true states
for (i in 2:ns) {
  #states
  x[i, ] <- x[i-1, ] + (A%*%x[i-1, ])*dt +
    B%*%c(2*hs[i], .5*hs[i]) +
    rnorm(nc, 0, sqrt(wk))}


#Plot of simulated states
par(mfrow=c(2, 2))
plot(1:nc, x[501,], type="l", xlab="Node, i", ylab="Temperature (i)",
     main="t=5")
plot(1:nc, x[1001,], type="l", xlab="Node, i", ylab="Temperature (i)",
     main="t=10")
plot(1:nc, x[5001,], type="l", xlab="Node, i", ylab="Temperature (i)",
     main="t=50")
plot(1:nc, x[10001,], type="l", xlab="Node, i", ylab="Temperature (i)",
     main="t=100")
par(mfrow=c(1, 1))


##MEASUREMENTS
#We will create measurements for the first 500 time steps.
#More specifically, we assume that we take measurements from t=.01 to t=5.
#Again, let k be the time step index and dt the size of the time steps.
#In that way, we may compute the time (t) as k*dt. Thus, for k=1 and dt=1/100
#the corresponding time is 1*(1/100)=.01

#Extract the temperature states for t=.01 to t=5.
#Note that x[1, ] represents the temperature states at t=0. Thus, for t=.01
#to t=5 we need to select x[2:501, ]
x <- x[2:501,]
#Select temperatures at nodes 10, 20, ..., 90
np <- length(seq(10,90,10)) #number of nodes
dataEx2 <- x[, seq(10, 90, 10)]

#Add measurement noise to the selected states
for (i in 1:500) {
  dataEx2[i, ] <- dataEx2[i, ] + rnorm(np, 0, sqrt(vk))}

#Plot the generated measurements at t=5
plot(1:nc, x[500,], type="l", xlab="Node, i", ylab="Temperature (i)",
     col=gray(level=.8))
points(seq(10,90,10), dataEx2[500,], col="darkgreen", cex=.5)
legend("topright", pch=c(NA, 1), lty=c(1, NA),
       col=c(gray(level=.8), "darkgreen"),
       legend=c("true", "measured"),
       bty="n", y.intersp=1.2, cex=.7)



##ENSEMBLE KALMAN FILTER (EnKF)

#Dynamic model: specifying 100 temperature states.
#Note that you may change the specified initial state estimates (at t=0)
#below for the 100 states and see how such changes influence the behavior
#of the ensemble Kalman filter.
ex2 <- list(m0=c(200, rep(100, nc-2), 150), #initial state estimates
            #error covariances of the initial state estimates:
            #this matrix reflects the uncertainty in our initial state estimates
            #you may change the values in this matrix and see how they influence
            #the behavior of the Kalman filter
            C0=diag(rep(1e+1, nc)),
            #measurement noise
            V=diag(rep(vk, 9)),
            #process noise
            W=diag(rep(wk, nc)))

#Specify the state transition function:
#WARNING: always use arguments x and k when specifying the GGfunction
GGfunction <- function (x, k){
  x + (A%*%x)*dt + B%*%c(2*hs[k+1], .5*hs[k+1])}

#Specify the observation/measurement function:
#WARNING: always use arguments x and k when specifying the FFfunction
FFfunction <- function (x, k){
  x[seq(10, 90, 10)]}



##Compute the filtered (a posteriori) state estimates with the EnKF
#and employ 10 ensemble members in the EnKF 
enkf1 <- EnKF(y=dataEx2, mod=ex2, size=10,
              GGfunction=GGfunction, FFfunction=FFfunction)

#As a comparison, increase size of the ensemble to 20
enkf2 <- EnKF(y=dataEx2, mod=ex2, size=20,
              GGfunction=GGfunction, FFfunction=FFfunction)

#And, finally, an EnKF with 100 ensemble members
enkf3 <- EnKF(y=dataEx2, mod=ex2, size=100,
              GGfunction=GGfunction, FFfunction=FFfunction)


#Plot the filtered state estimates at t=5
plot(1:nc, x[500,], type="l", col=c(gray(level=.5)),
     ylim=range(c(x[500,], enkf1$m[501,], enkf2$m[501,], enkf3$m[501,])),
     xlab="Node, i", ylab="Temperature (i)", main="t=5")
lines(1:nc, enkf1$m[501,], lty=2, col="blue", lwd=1)
lines(1:nc, enkf2$m[501,], lty=2, col="red", lwd=1)
lines(1:nc, enkf3$m[501,], lty=2, col="darkgreen", lwd=1)
legend("topright", lty=c(1, 2, 2, 2),
       col=c(gray(level=.5), "blue", "red", "darkgreen"),
       legend=c("true state", "EnKf with 10 members",
                "EnKf with 20 members", "EnKf with 100 members"),
       bty="n", y.intersp=1.2, cex=.7)



#Error plot
e1 <- sapply(1:nrow(x), function (i) mean((x[i,]-enkf1$m[i+1,])^2))
e2 <- sapply(1:nrow(x), function (i) mean((x[i,]-enkf2$m[i+1,])^2))
e3 <- sapply(1:nrow(x), function (i) mean((x[i,]-enkf3$m[i+1,])^2))
rangeError <- range(cbind(e1, e2, e3))
plot(1:nrow(x), e1, type="l", lty=1, col="blue", log="y",
     ylim=rangeError, ylab="Mean Squared Error", xlab="Time index, k")
lines(1:nrow(x), e2, lty=1, col="red")
lines(1:nrow(x), e3, lty=1, col="darkgreen")
legend("topright", lty=c(2, 2, 2),
       col=c("blue", "red", "darkgreen"),
       legend=c("EnKf with 10 members",
                "EnKf with 20 members", "EnKf with 100 members"),
       bty="n", y.intersp=1.2, cex=.7)



#Compute the confidence intervals for the filtered state estimates of the EnKF

#95% confidence intervals at t=5

#EnKF with 10 members
seFX10 <- sqrt(enkf1$C[[501]])
ciFX10 <- enkf1$m[501,] + qnorm(.05/2)*seFX10%o%c(1, -1)

#plot
plot(1:nc, x[500,], type="l", col=c(gray(level=.7)),
     lwd=1.5, ylim=range(ciFX10),
     xlab="Node, i", ylab="Temperature (i)", main="EnKF with 10 members")
lines(1:nc, ciFX10[,1], lty=2, col="blue")
lines(1:nc, ciFX10[,2], lty=2, col="blue")
legend("topright", lty=c(1, 2), 
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "filtered state"),
       bty="n", y.intersp=1.2, cex=.7)



#EnKF with 20 members
seFX20 <- sqrt(enkf2$C[[501]])
ciFX20 <- enkf2$m[501,] + qnorm(.05/2)*seFX20%o%c(1, -1)

#plot
plot(1:nc, x[500,], type="l", col=c(gray(level=.7)),
     lwd=1.5, ylim=range(ciFX20),
     xlab="Node, i", ylab="Temperature (i)", main="EnKF with 20 members")
lines(1:nc, ciFX20[,1], lty=2, col="blue")
lines(1:nc, ciFX20[,2], lty=2, col="blue")
legend("topright", lty=c(1, 2), 
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "filtered state"),
       bty="n", y.intersp=1.2, cex=.7)



#EnKF with 100 members
seFX100 <- sqrt(enkf3$C[[501]])
ciFX100 <- enkf3$m[501,] + qnorm(.05/2)*seFX100%o%c(1, -1)

#plot
plot(1:nc, x[500,], type="l", col=c(gray(level=.7)),
     lwd=1.5, ylim=range(ciFX100),
     xlab="Node, i", ylab="Temperature (i)", main="EnKF with 100 members")
lines(1:nc, ciFX100[,1], lty=2, col="blue")
lines(1:nc, ciFX100[,2], lty=2, col="blue")
legend("topright", lty=c(1, 2), 
       col=c(gray(level=.5), "blue"),
       legend=c("true state", "filtered state"),
       bty="n", y.intersp=1.2, cex=.7)