library(RODBC)
DSN_01 <- odbcConnect(
                  "DSN_01",
                  uid = "myname",
                  pwd = "mypwd",
                  believeNRows = FALSE,
                  case = "toupper"
                )
source("./R/PSTVB_Packages.R")
source("./R/SetSimuModel.R")
source("./R/spT.validation.R")
source("./R/stSemiPar.R")
source("./R/stGCVfun.R")
source("E:/Literature/semiBase/R/util.R")

n <- 100
Nt <- 10
Phis <- 0.2
delta <- 0.1

time <- seq(0, 1,, Nt)
theta.1 <- matrix(rep(f1(time), times = 1), 
                  nrow =  Nt, ncol = 1)
theta.2 <- matrix(rep(f1(time), times = 1), 
                  nrow =  Nt, ncol = 1)

theta.1.1 <- matrix(rep(f1(time), times = 50), 
                  nrow =  Nt, ncol = 50)
theta.2.1 <- matrix(rep(f1(time), times = 50), 
                  nrow =  Nt, ncol = 50)


Inde <- T
method <- 1
tab <- paste0(Inde, "_", method, "_", n, "_", Nt, "_", 
              substr(Phis, 3, 3)) 
Result1 = sqlQuery(DSN_01, paste0("SELECT * FROM ", tab))
load(paste0("./data/", tab, "_alpha_est.RData"))
alpha.est1 <- alpha.est

# alpha.est1[1:10, ] - theta.1
# index <- which(Result1$Y_RMSE > 5)
# if(length(index)>0){
#   Result1 <- Result1[-index, ]
#   alpha.est1 <- alpha.est1[, -index]
#   # theta.1.1 <- theta.1.1[, -index]
# }




Inde <- FALSE
tab <- paste0(Inde, "_", method, "_", n, "_", Nt, "_", 
              substr(Phis, 3, 3), "T") 
Result2 = sqlQuery(DSN_01, paste0("SELECT * FROM ", tab))
load(paste0("./data/", tab, "_alpha_est.RData"))
alpha.est2 <- alpha.est

# index <- which(Result2$Y_RMSE > 5)
# if(length(index)>0){
# Result2 <- Result2[-index, ]
# alpha.est2 <- alpha.est2[, -index]
# theta.1.2 <- theta.1.1[, -index]
# theta.2.2 <- theta.2.1[, -index]
# }
# head(Result0[, c(3:5, 8, 11:15)])
# colMeans(Result0)

Inde <- F
method <- 2
tab <- paste0(Inde, "_", method, "_", n, "_", Nt, "_", 
              substr(Phis, 3, 3)) 
Result3 = sqlQuery(DSN_01, paste0("SELECT * FROM ", tab))
load(paste0("./data/", tab, "_alpha_est.RData"))
alpha.est3 <- alpha.est

# index <- which(Result3$Y_RMSE > 5)
# if(length(index)>0){
# Result3 <- Result3[-index, ]
# alpha.est3 <- alpha.est3[, -index]
# theta.1.3 <- theta.1.1[, -index]
# theta.2.3 <- theta.2.1[, -index]
# }

# head(Result1[, c(3:5, 8, 11:15)])
# colMeans(Result1)

# (E(beta.hat) - beta) and sd(beta.hat)
round(c(mean(Result1[, 3]) - 1,
        mean(Result2[, 3]) - 1, 
        mean(Result3[, 3]) - 1, 
        sd((Result1[, 3] - 0)^1),
        sd((Result2[, 3] - 0)^1),
        sd((Result3[, 3] - 0)^1),
        mean((Result1[, 3] - 1)^2), 
        mean((Result2[, 3] - 1)^2), 
        mean((Result3[, 3] - 1)^2)), 
        4)

round(c(mean(Result1[, 4]) - 5,
        mean(Result2[, 4]) - 5, 
        mean(Result3[, 4]) - 5, 
        sd((Result1[, 4] - 0)^1),
        sd((Result2[, 4] - 0)^1),
        sd((Result3[, 4] - 0)^1),
        mean((Result1[, 4] - 5)^2), 
        mean((Result2[, 4] - 5)^2), 
        mean((Result3[, 4] - 5)^2)),
      4)


# head(Result1[, c(3:5, 8, 11:15)])
# head(Result2[, c(3:5, 8, 11:15)])
# head(Result3[, c(3:5, 8, 11:15)])

# colMeans(Result1[,c(3:5, 8, 11:15)])
# colMeans(Result2[,c(3:5, 8, 11:15)])
# colMeans(Result3[,c(3:5, 8, 11:15)])

# theta
# round(c(mean(rowMeans(alpha.est1[1:10, ])- theta.1), 
#         mean(rowMeans(alpha.est2[1:10, ])- theta.1),
#         mean(rowMeans(alpha.est3[1:10, ])- theta.1),
#         mean(sqrt(rowVar(alpha.est1[1:10, ]))),
#         mean(sqrt(rowVar(alpha.est2[1:10, ]))),
#         mean(sqrt(rowVar(alpha.est3[1:10, ]))),
#         mean((alpha.est1[1:10, ]- theta.1.1)^2), 
#         mean((alpha.est2[1:10, ]- theta.1.1)^2), 
#         mean((alpha.est3[1:10, ]- theta.1.1)^2)),
#         4)
# 
# round(c(mean(rowMeans(alpha.est1[11:(2*Nt), ])- theta.2), 
#         mean(rowMeans(alpha.est2[11:(2*Nt), ])- theta.2),
#         mean(rowMeans(alpha.est3[11:(2*Nt), ])- theta.2),
#         mean(sqrt(rowVar(alpha.est1[11:(2*Nt), ]))),
#         mean(sqrt(rowVar(alpha.est2[11:(2*Nt), ]))),
#         mean(sqrt(rowVar(alpha.est3[11:(2*Nt), ]))),
# 		    mean((alpha.est1[11:(2*Nt), ]- theta.2.1)^2), 
#         mean((alpha.est2[11:(2*Nt), ]- theta.2.1)^2), 
#         mean((alpha.est3[11:(2*Nt), ]- theta.2.1)^2)), 4)

