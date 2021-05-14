#######################################################
X_ts_Transf <- function(Nt, X_ts, beta)
{
  t = seq_len(Nt)
  x_ts = sapply(t, function(t) t(X_ts[, , t]) %*% beta
                 , simplify = "matrix")
  return(x_ts)
}