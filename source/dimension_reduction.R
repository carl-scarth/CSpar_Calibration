# Header file used to compute reduced-dimensional quantities used in the
# calibration method published in:
# D. Higdon et al, "Computer Model Calibration Using High-Dimensional Output",
# Journal of the American Statistical Association, 2008.

reduce_dimension <- function(eta, K, a_eta, b_eta, orthog_K = TRUE) {
  
  # This portion of the code deals with the matrix algebra from Section 2.2.4 of 
  # Higdon et al., calculating the necessary quantities for passing to Stan
  
  # eta = n_eta x m matrix of simulation outputs, where n_eta is the number of
  #     output values per simulation, and m is the number of training data runs
  # K = n_eta x p_eta matrix of basis functions used to decompose eta, where
  #     p_eta is the number of basis functions
  # a_eta = Shape parameter of the gamma prior placed on the error associated 
  #         with truncating the expansion for eta
  # b_eta = Rate parameter of the gamma prior on the truncation error term
  # orthog_K = Boolean stating if the basis functions K are orthogonal
  
  # Calculate inverse of K^T*K. K is the matrix of emulator basis functions for 
  # the full model output arranged as in Section 2.2.2 of Higdon et al. The 
  # inverse is calculated via closed-form expressions for K^T*K.
  
  p_eta = ncol(K) # Number of basis functions
  m = ncol(eta) # Number of training samples

  # Populate dense matrix of products of products with KTK_dense[i,j] = k_i'*k_j
  KTK_dense = t(K)%*%K
  
  # If k_i are orthogonal then K^T*K is diagonal and the inverse is specified 
  # directly via the closed-form expression. Otherwise populate the full matrix 
  # K^T*K and calculate the inverse of this matrix in R.
  if (orthog_K){
    KTKinv = matrix(0, m*p_eta, m*p_eta)
    for (i in 1:p_eta){
      KTKinv[((i-1)*m+1):(i*m),((i-1)*m+1):(i*m)] = diag(m)*(1/KTK_dense[i,i])
    }
  } else {
    KTK = matrix(0, m*p_eta, m*p_eta)
    for (i in 1:p_eta){
      # Populate diagonal terms
      KTK[((i-1)*m+1):(i*m),((i-1)*m+1):(i*m)] = diag(m)*KTK_dense[i,i]
      # Compute off-diagonal terms. Need to use seq_len as colon operator in R 
      # doesn't assume the list must always increase, and j in 1:(i+1) throws an 
      # error when i = p_eta
      for (j in (i + seq_len(p_eta-i))){
        KTK[((i-1)*m+1):(i*m),((j-1)*m+1):(j*m)] = diag(m)*KTK_dense[i,j]
        KTK[((j-1)*m+1):(j*m),((i-1)*m+1):(i*m)] = diag(m)*KTK_dense[i,j]
      }
    }
    # Determine inverse of K'*K.
    KTKinv = solve(KTK)
  }
  
  # Calculate K^T*eta. Note that eta in the paper eta is reshaped into a 
  # m x n_eta vector, and so the resulting expression must also be reshaped.
  # Doing the matrix product first is more efficient and has the same result.
  KTeta = t(eta)%*%(K)
  KTeta = as.vector(KTeta)
  
  # Do pseudo-inverse of basis matrix to get reduced-dimensional representation 
  # of model output
  z_hat = as.vector(KTKinv%*%KTeta)
  
  # Adjust prior parameters of the expansion truncation error to account for the
  # dimension reduction(Eq. 11 Higdon et al.)
  a_eta_dash = a_eta+(0.5*(m*(n_eta-p_eta))) 
  # Stack model output eta into a single vector
  eta_vec = as.vector(eta)
  # Re-arranged b_eta_dash from Higdon et al. for computational efficiency using
  # eta'*(I - K*(K'*K)^-1*K')*eta = eta'*eta - (K'*eta)'*w_hat
  b_eta_dash = as.numeric(b_eta + (0.5*(t(eta_vec)%*%eta_vec - t(KTeta) %*% z_hat)))
  
  return(list(a_eta_dash, b_eta_dash, z_hat, KTKinv))
  
}