# Header file used to compute reduced-dimensional quantities used in the
# calibration method published in:
# D. Higdon et al, "Computer Model Calibration Using High-Dimensional Output",
# Journal of the American Statistical Association, 2008.
library(abind)

svd_basis <- function(eta, p_eta = NULL, exp_tol = NULL, print_output = FALSE, export_basis = TRUE, csv_label = NULL){
  
  # Perform an SVD upon eta to generate a basis of principal components for 
  # reducing the dimension of eta, and return the first p_eta basis functions.
  # If p_eta is not specified this is instead computed by determining how many
  # basis functions are required to capture exp_tol fraction of the variance
  # in the training data.
  
  # eta = n_eta x m matrix of simulation outputs, where n_eta is the number of
  #     output values per simulation, and m is the number of training data runs
  # p_eta = integer value of number of basis functions to be returned
  # exp_tol = desired fraction of the variance which is to be returned
  # print_output = Boolean indicating whether or not to print output to the 
  #     console, and produce convergence plots
  # export_basis = Boolean indicating whether or not to export basis functions
  #     to a csv for plotting externally
  # csv_label = String for labelling exported csv file if this is produced
  
  m = ncol(eta) # number of training samples
  n_eta = nrow(eta) # Number of output values per simulation
  
  # Perform SVD on the centred data
  dt_svd = svd(eta)
  
  # Calculate the fraction of total variance captured by each singular value if 
  # this is required to assess convergence
  if ((!is.null(exp_tol)) | print_output) {
    dr_sqr = sum(dt_svd$d^2)
    d_r = dt_svd$d
    # Fraction of the total variance captured by each singular value
    d_r_norm = d_r^2/dr_sqr
    # Determine the number of basis functions required to capture a specified
    # fraction of the variance if basis_tol has been specified
    if (! is.null(exp_tol)){
      # How many functions are needed achieve a tolerance fraction of the variance?
      basis_tol = 1:m
      basis_tol = basis_tol[cumsum(d_r_norm) > (1-exp_tol)]
      basis_tol = basis_tol[1]
      if (is.null(p_eta)){
        p_eta = basis_tol
      }
      if (print_output){
        print("Number of basis functions required to represent output within tolerance = ")
        print(basis_tol)
      }
    }
    if (print_output){
      # Sanity check that weights of SVD have zero mean and unit variance
      print("mean of reduced dimension output w = ")
      print(colMeans(dt_svd$v))
      print("standard deviations of reduced dimension output w = ")
      print(colSds(dt_svd$v*sqrt(m-1)))
    }
  }
  
  # If neither p_eta nor exp_tol has been provided, retain all basis functions
  if (is.null(p_eta)){
    print("WARNING: Neither p_eta nor exp_tol has been specified, and so all bases are retained")
    p_eta = m
  }
  
  # If needed plot magnitude of singular value with increasing number of bases
  # to indicate how much each base contributes to the output variance.
  if (print_output){
    par(mfrow = c(1,2))
    plot(1:p_eta,d_r_norm[1:p_eta],"type"="p","col"="red","pch"=4,"lwd"=3,cex=1.5,
       'xlab' = "Feature",'ylab'="normalised d_i",cex.axis=1.75,cex.lab=1.75)
    # Omit first point for greater clarity on convergence
    plot(2:p_eta,d_r_norm[2:p_eta],"type"="p","col"="red","pch"=4,"lwd"=3,cex=1.5,
       'xlab' = "Feature",'ylab'="normalised d_i",cex.axis=1.75,cex.lab=1.75)
    title("Convergence with Number of Basis Functions", line = 0, outer = TRUE, cex.main=1.75)
  }
    
  # Extract the first p_eta basis functions from the svd, and standardise so the
  # coefficients (columns of sqrt(m-1)*v) have unit variance
  K = dt_svd$u[,1:p_eta]*matrix(dt_svd$d[1:p_eta],nrow = n_eta, ncol = p_eta, byrow = TRUE)/sqrt(m-1)

  # write basis functions and mean vector to file for external plotting if needed
  if (export_basis){
    out_frame = as.data.frame(K)
    for (i in 1:p_eta){colnames(out_frame)[i] <- sprintf("K_basis_%d",i)}
    if (is.null(csv_label)){
      write.csv(out_frame, "outputs/basis.csv", row.names = FALSE)    
    } else {
      write.csv(out_frame, paste("outputs/basis_",csv_label,".csv", sep=""), row.names = FALSE)    
    }
  }
  
  # Return relevant quantities
  return(list(K, p_eta))
}

reduce_dimension_emulator <- function(eta, K, a_eta = NULL, b_eta = NULL, orthog_K = TRUE) {
  
  # Perform the matrix algebra from Section 2.2.4 of Higdon et al., calculating
  # the necessary quantities for passing to Stan. Focus on quantities relating
  # to the emulator
  
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
  
  # If required, adjust prior parameters of the expansion truncation error to 
  # account for the dimension reduction(Eq. 11 Higdon et al.)
  if (!is.null(a_eta) & !is.null(b_eta)){
    a_eta_dash = a_eta+(0.5*(m*(n_eta-p_eta))) 
    # Stack model output eta into a single vector
    eta_vec = as.vector(eta)
    # Re-arranged b_eta_dash from Higdon et al. for computational efficiency using
    # eta'*(I - K*(K'*K)^-1*K')*eta = eta'*eta - (K'*eta)'*w_hat
    b_eta_dash = as.numeric(b_eta + (0.5*(t(eta_vec)%*%eta_vec - t(KTeta) %*% z_hat)))  
    return(list(a_eta_dash, b_eta_dash, z_hat, KTKinv))
  } else {
    return(list(z_hat, KTKinv))
  }
  
}

reduce_dimension_calibration <- function(y, B, W_y, q_y = 1, a_y = NULL, b_y = NULL){
  
  # Perform the matrix algebra from Section 2.2.4 of Higdon et al., calculating
  # the necessary quantities for passing to Stan. Focus on quantities relating
  # to the experimental data
  
  # y = n_y-vector of experimental measurements, where n_y is the number of data
  #     points
  # B = n_y x p_B matrix of basis functions used to decompose eta, where p_B
  #     depends upon whether, and how the discrepancy is included. For the most 
  #     general implementation of Higdon et al. this will equal p_eta + p_delta
  # W_y = Prior precision of observation error, which can be passed as an n_y x
  #     n_y matrix, or if the precision is diagonal, as a n_y vector of 
  #     the diagonal
  # q_y = Number of output components (if output is a vector field)
  # a_y = Shape parameter of the gamma prior on the observation error
  # b_y = Rate parameter of the gamma prior on the observation error
  
  # Calculate B'*W_y*B
  p_B = ncol(B)
  n_y = length(y)/q_y
  BTWyB = matrix(0, nrow = p_B, ncol = p_B) # Total matrix product
  BTWyB_sep = array(NA, dim = c(0, p_B, p_B)) # matrix product separated into different vector components
  for (i in 1:q_y) {
    inds_i = ((i-1)*n_y+1):(i*n_y)
    if (is.vector(W_y)){
      # BTWyB = t(K_y)%*%(W_y*K_y) # For diagonal W_y 
      BTWyB_i = t(K_y[inds_i,])%*%(W_y[inds_i]*K_y[inds_i,]) # For diagonal W_y 
    } else {
      # BTWyB = t(K_y[inds_i,])%*%W_y[inds_i,inds_i]%*%K_y[inds_i,]
      BTWyB_i = t(K_y[inds_i,])%*%W_y[inds_i,inds_i]%*%K_y[inds_i,]
    }
    BTWyB = BTWyB + BTWyB_i
    BTWyB_sep = abind(BTWyB_sep, array(BTWyB_i, dim = c(1, p_B, p_B)), along=1)
  }
  BTWyBinv = solve(BTWyB)
  
  # Calculate reduced dimensional model coefficients for the experimental data.
  if (is.vector(W_y)){
    BTWyy = t(B)%*%(W_y*y) # For diagonal W_y
  } else {
    BTWyy = t(B)%*%W_y%*%y
  }
  
  u_hat = as.vector(BTWyBinv%*%BTWyy)
  
  # If required, adjust prior parameters of the observation error to account 
  # for the dimension reduction
  # I'M NOT SURE IF THESE EXPRESSIONS ARE CORRECT WITH q_y > 0
  if (!is.null(a_y) & !is.null(b_y)){
    a_y_dash = a_y+(0.5*(n_y-p_eta)) # Adjusted shape parameter for the lambda_y prior.
    # Re-arraged version of b_y_dash from that in Eq. (11) for efficiency.
    if (is.vector(W_y)){
      b_y_dash = as.numeric(b_y + (0.5*(t(y)%*%(W_y*y) - t(BTWyy) %*% u_hat))) # For diagonal W_y
    } else {
      b_y_dash = as.numeric(b_y + (0.5*(t(y)%*%W_y%*%y - t(BTWyy) %*% u_hat)))
    }
    # return(list(a_y_dash, b_y_dash, u_hat, BTWyBinv))
    return(list(a_y_dash, b_y_dash, u_hat, BTWyB_sep))
  } else {
    # return(list(u_hat, BTWyBinv))
    return(list(u_hat, BTWyB_sep))
  }
}

adjust_error_covariance <- function(Sigma_z, KTKinv, BTWyBinv, lambda_eta, lambda_y) {
  # Adust covariance matrix Sigma_z of the calibration statistical model for the
  # effect of dimension-reduction upon the model decomposition truncation error,
  # and experimental observation error, which are defined for full-field output
  # Sigma_z = joint emulator covariance matrix for training data and experiment
  # KTKinv = inverse of the inner product of the basis matrix used to decompose
  # full-field model output
  # BTWy_Binv = inverse of matrix product of basis matrix at experimental data 
  # points, and prior observation error precision
  # lambda_eta = scalar sample of truncation error precision
  # lambda_y = scalar sample of observation error precision
  n_exp = nrow(BTWyBinv)
  Sigma_z_hat = matrix(0, nrow(Sigma_z), ncol(Sigma_z))
  Sigma_z_hat[1:n_exp,1:n_exp] = BTWyBinv/lambda_y
  Sigma_z_hat[-(1:n_exp),-(1:n_exp)] = KTKinv/lambda_eta
  Sigma_z_hat = Sigma_z_hat + Sigma_z
  
  return(Sigma_z_hat)
}