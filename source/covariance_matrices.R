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

full_field_cov = function(x, lambda, beta) {
  # Auto covariance matrix for a full-field emulator composed of p basis 
  # functions, for a set of N samples, resulting in a symmetric block diagonal 
  # matrix with the ith block being a covariance matrix for the ith base.
  # lambda is a vector of p emulator precisions, with ith entry for the ith base
  # beta is a vector of p*q correlation lengths, where q is the dimension of x,
  # with each q consecutive entries corresponding to a base
  
  p = length(lambda) # Determine number of bases
  N = nrow(x) # Determine number of samples
  d = ncol(x) # Determine number of dimensions of x
  K = matrix(0, N*p, N*p)
  
  for (i in 1:p) {
    # Calculate the covariance matrix for the ith base
    K[((i-1)*N+1):(i*N), ((i-1)*N+1):(i*N)] = ARD_SE_cov(x, lambda[i], beta[((i-1)*d+1):(i*d)], 0)
  }
  
  return(K)
}

full_field_cov_non_sym = function(x1, x2, lambda, beta) {
  # Cross covariance function for a full-field emulator composed of p basis 
  # functions, for two sets of sample inputs x1 and x2, resulting in a block 
  # diagonal matrix with the ith diagonal being an N1 x N2 cross covariance
  # matrix for the ith base.
  
  p = length(lambda)
  N1 = nrow(x1)
  N2 = nrow(x2)
  d = ncol(x1)
  K = matrix(0, N1*p, N2*p)
  
  for (i in 1:p) {
    K[((i-1)*N1+1):(i*N1),((i-1)*N2+1):(i*N2)] = ARD_SE_cov_non_sym(x1, x2, lambda[i], beta[((i-1)*d+1):(i*d)])
  }
  
  return(K)
}

ff_calibration_cov = function(Sigma_u, Sigma_w, Sigma_uw) {
  # Assemble joint covariance matrix for experimental and model training data 
  # used in the full-field model calibration.
  # Sigma_u = p_eta x p_eta emulator autocovariance matrix for the experimental 
  # data
  # Sigma_w = (m*p_eta)x(m*p_eta) emulator autocovariance matrix for the 
  # training data points
  # Sigma_uw = p_eta x (m*p_eta) cross-covariance matrix between the 
  # experimental data and model data
  # where m is the number of training data points, and p_eta is the number of
  # basis functions used to represent the full-field model output
  # (Note: This same code could be used in the univariate case, though not sure
  # if it would hold for multiple experimental data points)
  nu = nrow(Sigma_u)
  nw = nrow(Sigma_w)
  Sigma_z = matrix(0,nu+nw,nu+nw)
  Sigma_z[1:nu,1:nu] = sigma_u
  Sigma_z[-(1:nu),-(1:nu)] = sigma_w
  Sigma_z[1:nu,-(1:nu)] = sigma_uw
  Sigma_z[-(1:nu),1:nu] = t(sigma_uw)
  return(Sigma_z)
}
