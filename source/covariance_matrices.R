# Define a function for evaluating a square-exponential covariance matrix with
# non-isotropic correlation length i.e. there is a different correlation length 
# for each input
ARD_SE_cov <- function(x, lambda, beta, delta) {
  N = nrow(x) # number of data points
  d = ncol(x) # number of input dimensions
  K = diag(rep(1/lambda, N))
  for (i in 1:N) {
    for (j in (i + seq_len(N-i))) {
      dx = x[i,] - x[j,]
      K[i, j] = exp(-sum(beta*dx^2))/lambda
      K[j, i] = K[i, j]
    }
  }
  # Add nugget term to diagonals
  K = K + diag(rep(delta,N))
  
  return(K)
}

# Define a function for evaluating a square-exponential covariance matrix with 
# non-isotropic correlation length as above, but for a non-symmetric covariance 
# matrix. This is used when evaluating the covariance between two different sets
# of points.
ARD_SE_cov_non_sym <- function(x1, x2, lambda, beta) {
  N1 = nrow(x1)   # Number of data points in x1
  N2 = nrow(x2)   # Number of data points in x2
  d = ncol(x1)    # Number of input dimensions
  K = matrix(0,N1,N2)
  for (i in 1:N1) {
    for (j in 1:N2) {
      dx = x1[i,] - x2[j,]
      K[i, j] = exp(-sum(beta*dx^2))/lambda
    }
  }
  
  return(K)
}