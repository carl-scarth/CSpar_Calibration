# Code for making predictions using fitted Gaussian process emulators, and 
# calibrated models

library(MASS) # Needed for mvrnorm
source("source/covariance_matrices.R")
source("source/dimension_reduction.R")

zero_mean_gp_pred <- function(y, K_x, k_x_xstar, k_xstar, sam_gp = F, N_sam = 1, inv_K = F, nugget = F){
  # Make predictions of a zero-mean Gaussian process, where y is the training
  # data output, K_x is the autocovariance matrix of the training data, 
  # k_x_xstar is the crosscovariance of training data points with test data, and
  # k_xstar is the autocovariance of the test data point
  # sam_gp = Boolean which dictates whether to draw N_sam samples from the GP
  # inv_K = Boolean, which if specified, allows the inverse of K_x to be passed
  # for computational efficiency, as this quantity is independed of the test 
  # point and can often be calculated outside of loops
  
  # Calculate predictive mean and covariance. Use a different method if the 
  # inverse of the training data covariance matrix has been provided
  if (!inv_K){
    f_mu = t(k_x_xstar) %*% solve(K_x, y)
    f_sigma = k_xstar - (t(k_x_xstar) %*% solve(K_x, k_x_xstar))
  } else {
    # Add nugget if requested
    if (nugget){
      K_x = K_x + diag(rep(1e-8, nrow(K_x)))
    }
    K_x_inv = K_x # Inverse of K has been passed to the function instead
    f_mu = t(k_x_xstar) %*% K_x_inv %*% y
    f_sigma = k_x_xstar - (t(k_x_xstar) %*% K_x_inv %*% k_x_xstar)
  }

  # Sample from the Gaussian process if required
  if (sam_gp) {
    f_star = mvrnorm(n = N_sam, f_mu, f_sigma)
  } else {
    f_star = NULL
  }
  return(list(f_mu, f_sigma, f_star))
}

full_field_gp_pred <- function(N_subsam, x_train, z_train, x_pred, beta_w, lambda_w, lambda_eta, K ,KTKinv, sam_gp = F, output_coeff_sam = F, output_ff_sam = F, output_coeff_mean = T, output_ff_mean = T){
  # Make N_subsam predictions from a Gaussian process emulator fitted to full
  # field model output, using the method described by Higdon et al.
  # Sub-sampling from the posterior helps reduce memory requirement, preventing
  # errors when dealing with very high-dimensional output.
  # beta_w = n_sam_post x p*q matrix of posterior correlation length samples
  # lambda_w = n_sam_post x p*q matrix of posterior emulator precision samples
  # lambda_eta = n_sam_post vector of posterior expansion truncation error samples
  # K = emulator basis matrix
  # KTKinv = inverse of inner product of emulator basis matrix with itself
  # sam_gp = Boolean which dictates whether to sample from the GP
  # output_coeff_sam = Boolean on whether to output info for the basis coefficients
  # output_ff_sam = Boolean on whether to output info for the full-field response
  # output_coeff_mean = Boolean on whether to average across posterior predictions for the basis coefficients
  # output_ff_mean = Boolean on whether to average across posterior predictions for the full-field
  
  N_post = nrow(beta_w) # Determine overall number of posterior samples
  N_pred = nrow(x_pred) # Determine number of points at which predictions are needed
  N_eta = nrow(K) # Determine number of output values per simulation
  p = ncol(lambda_w) # Determine number of bases
  
  
  # Sub-sample from the posterior distribution
  subsam_ind = sample.int(N_post, N_subsam) # Determine index of posterior samples
  beta_w = beta_w[subsam_ind,]
  lambda_w = lambda_w[subsam_ind,]
  lambda_eta = lambda_eta[subsam_ind]
  
  # Initialise output arrays
  if (!is.null(K)){has_K = TRUE} else {has_K = FALSE}
  out_list = initialise_out_list(p=p, N_eta=N_eta, N_post=N_subsam, N_pred=N_pred, K=has_K, sam_gp=sam_gp, output_coeff_sam=output_coeff_sam, output_ff_sam=output_ff_sam, output_coeff_mean=output_coeff_mean, output_ff_mean=output_ff_mean)

  # Loop over all posterior samples
  for (i in 1:N_subsam){
    print(i) # Print count of loop
    # Construct covariance matrices for making predictions
    Sigma_z = full_field_cov(x_train, lambda_w[i,], beta_w[i,])
    # Adjust sigma_z to account for the dimension reduction
    Sigma_z_hat = Sigma_z + KTKinv/lambda_eta[i]
    # Determine covariance of emulator prediction with itself
    sigma_w_star = diag(1/lambda_w[i,])
  
    # For each point at which a prediction is needed calculate the mean and 
    # covariance matrices. Generate samples if required
    for (j in 1:N_pred){
      # Determine cross covariance of training data with the jth prediction
      Sigma_z_w_star = full_field_cov_non_sym(x_train, t(as.matrix(x_pred[j,])), lambda_w[i,], beta_w[i,])
      # Make predictions in reduced dimensional space
      pred = zero_mean_gp_pred(z_train, Sigma_z, Sigma_z_w_star, sigma_w_star, sam_gp = sam_gp)
      w_star_mu = pred[[1]]
      w_star_sigma = pred[[2]]
      w_star = pred[[3]]

      # Store emulator posterior samples of basis coefficients if requested
      if (output_coeff_sam){
        out_list$w_star_mu[,i,j] = w_star_mu
        out_list$w_star_sigma[,,i,j] = w_star_sigma
        if (sam_gp) {
          out_list$w_star[,i,j] = w_star
        }
      }      
    
      # Store full-field posterior prediction samples if requested
      eta_mu = K %*% w_star_mu
      eta_sigma = sqrt(as.matrix(rowSums((K %*% w_star_sigma) * K)))
      if (sam_gp) {
        eta_sam = K %*% w_star
      }
      if (output_ff_sam){
        out_list$eta_mu[,i,j] = eta_mu
        # Only output diagonals of covariance as full matrix is too high-
        # dimensional. Take square-root for standard
        out_list$eta_sigma[,i,j] = eta_sigma
        if (sam_gp) {
          out_list$eta_sam[,i,j] = eta_sam
        }
      }
    
      # If averages are required across posterior samples, keep a running total
      # to reduce memory requirements, compared with calculating mean at the end
      if (output_coeff_mean) {
        out_list$w_star_mu_mu[,j] = out_list$w_star_mu_mu[,j] + w_star_mu/N_subsam
        out_list$w_star_sigma_mu[,,j] = out_list$w_star_sigma_mu[,,j] + w_star_sigma/N_subsam
        if (sam_gp) {
          out_list$w_star_sam_mu[,j] = out_list$w_star_sam_mu[,j] + w_star/N_subsam
        }
      }
      # Add to average across full-field posterior samples if required
      if (output_ff_mean){
        out_list$eta_mu_mu[,j] = out_list$eta_mu_mu[,j] + eta_mu/N_subsam
        out_list$eta_sigma_mu[,j] = out_list$eta_sigma_mu[,j] + eta_sigma/N_subsam
        if (sam_gp) {
          out_list$eta_sam_mu[,j] = out_list$eta_sam_mu[,j] + eta_sam/N_subsam
        }
      }
    }
  }

  
  return(out_list)
}

full_field_calibration_pred_fixed_em <- function(N_subsam, tc, tf, z_hat, beta_w, lambda_w, lambda_eta, lambda_y, KTKinv, BTWyBinv, K = NULL, K_y = NULL, sam_gp = F, nugget = F, output_coeff_sam = F, output_ff_sam = F, output_coeff_mean = T, output_ff_mean = T, output_ff_std = F) { 
  # Make N_post_pred predictions replicating the full-field test data for a 
  # single condition using a calibrated Gaussian process, with emulator 
  # parameters fixed to constant values. Sub-sampling from the posterior helps
  # reduce memory requirement, preventing errors when dealing with very 
  # high-dimensional output
  # tc = m x q matrix of calibration input training data values
  # tf = N_sam_post x q matrix of posterior samples of calibration inputs
  # z_hat = (m+1)*p_eta vector of emulator coefficients for experimental and 
  # model training data
  # beta_w = p_eta*q vector of correlation lengths
  # lambda_w = p*q vector of emulator precisions
  # lambda_eta = scalar-valued expansion truncation error
  # lambda_y = N_sam_post vector of observation error posterior samples
  # KTKinv = inverse of inner product of emulator basis matrix with itself
  # BTWyBinv = inverse of matrix product of basis matrix of experimental data 
  # and prior observation error precision
  # sam_gp = Boolean which dictates whether to sample from the GP
  # nugget = Boolean which dictates whether to add a nugget to the covariance
  # output_coeff_sam = Boolean on whether to output info for the basis coefficients
  # output_ff_sam = Boolean on whether to output info for the full-field response
  # output_coeff_mean = Boolean on whether to average across posterior predictions for the basis coefficients
  # output_ff_mean = Boolean on whether to average across posterior predictions for the full-field
  # output_ff_std = Boolean on whether to output standard deviation across posterior predictions

  N_post = nrow(tf) # Determine overall number of posterior samples
  if (!is.null(K)){N_eta = nrow(K)} else {N_eta = NULL} # Determine number of output values per simulation
  if (!is.null(K_y)){N_y = nrow(K_y)} else {N_y = NULL} # Determine number of experimental data points
  p_eta = length(lambda_w) # Determine number of bases (this might not work for fixed point)

  # Sub-sample from the posterior distribution
  subsam_ind = sample.int(N_post, N_subsam) # Determine index of posterior samples
  tf = tf[subsam_ind, ]
  lambda_y = lambda_y[subsam_ind]
  
  # Initialise output arrays
  if (!is.null(K)){has_K = TRUE} else {has_K = FALSE}
  if (!is.null(K_y)){has_K_y = TRUE} else {has_K_y = FALSE}
  out_list = initialise_out_list(p=p_eta, N_eta = N_eta, N_y = N_y, N_post=N_subsam, K=has_K, K_y=has_K_y, sam_gp = sam_gp, output_coeff_sam = output_coeff_sam, output_ff_sam = output_ff_sam, output_coeff_mean = output_coeff_mean, output_ff_mean = output_ff_mean, output_ff_std = output_ff_std)

  # Covariance matrix for the emulator at the experimental data points. This 
  # reduces to a diagonal matrix for a single experiment
  Sigma_u = diag(1/lambda_w, nrow=p_eta, ncol=p_eta)
  # Define the emulator covariance matrix evaluated at the training data points
  Sigma_w = full_field_cov(tc, lambda_w, beta_w)
  # Loop over selected posterior samples
  for (i in 1:N_subsam){
    print(i) # Print count of loop
    # Construct covariance matrices for making predictions
    # Calculate emulator cross-covariance between training data and experimental data points
    Sigma_uw = full_field_cov_non_sym(tf[i,,drop=FALSE], tc, lambda_w, beta_w)
    # Calculate emulator auto-covariance for training data points
    Sigma_z = ff_calibration_cov(Sigma_u, Sigma_w, Sigma_uw)
    # Adjust to account for the expansion truncation and observation error
    Sigma_z_hat = adjust_error_covariance(Sigma_z, KTKinv, BTWyBinv, lambda_eta, lambda_y[i])
    # Calculate cross-covariance of prediction with (model) training data. This 
    # is given by sigma_u_w' as the prediction has the controlled input value x 
    # as the training data, though this does not hold in general
    Sigma_w_w_star = t(Sigma_uw)
    # Define cross-covariance for joint data z with prediction w_star
    Sigma_z_w_star = rbind(Sigma_u, Sigma_w_w_star)
    # Determine auto-covariance of the prediction
    Sigma_w_star = Sigma_u
    # Determine GP mean and covariance in reduced-dimensional space, and if 
    # required sample from the GP
    pred = zero_mean_gp_pred(z_hat, Sigma_z_hat, Sigma_z_w_star, Sigma_w_star, sam_gp = sam_gp, nugget = nugget)
    w_star_mu = pred[[1]]
    w_star_sigma = pred[[2]]
    w_star = pred[[3]]
    
    # Store emulator posterior samples of basis coefficients if requested
    if (output_coeff_sam){
      out_list$w_star_mu[,i] = w_star_mu
      out_list$w_star_sigma[,,i] = w_star_sigma
      if (sam_gp) {
        out_list$w_star[,i] = w_star
      }
    }
    
    # Convert from reduced-dimensional space to full-field, output at 
    # experimental data points, or both, depending upon which basis matrices
    # have been provided
    if (!is.null(K)){
      eta_mu = c(K %*% w_star_mu)
      eta_sigma = sqrt((rowSums((K %*% w_star_sigma) * K)))
      if (sam_gp) {
        eta_sam = c(K %*% w_star)
      }
    }
    if (!is.null(K_y)){
      eta_mu_y = c(K_y %*% w_star_mu)
      eta_sigma_y = sqrt((rowSums((K_y %*% w_star_sigma) * K_y)))
      if (sam_gp) {
        eta_sam_y = c(K_y %*% w_star)
      }
    }
    if (output_ff_sam & !is.null(K)){
      out_list$eta_mu[,i] = eta_mu
      # Only output diagonals of covariance as full matrix is too high-
      # dimensional. Take square-root for standard
      out_list$eta_sigma[,i] = eta_sigma
      if (sam_gp) {
        out_list$eta_sam[,i] = eta_sam
      }
    }
    if (output_ff_sam & !is.null(K_y)) {
      out_list$eta_mu_y[,i] = eta_mu_y
      # Only output diagonals of covariance as full matrix is too high-
      # dimensional. Take square-root for standard
      out_list$eta_sigma_y[,i] = eta_sigma_y
      if (sam_gp) {
        out_list$eta_sam_y[,i] = eta_sam_y
      }
    }
    
    # If averages are required across posterior samples, keep a running total
    # to reduce memory requirements, compared with calculating mean at the end
    if (output_coeff_mean) {
      out_list$w_star_mu_mu = out_list$w_star_mu_mu + w_star_mu/N_subsam
      out_list$w_star_sigma_mu = out_list$w_star_sigma_mu + w_star_sigma/N_subsam
      if (sam_gp) {
        out_list$w_star_sam_mu = out_list$w_star_sam_mu + w_star/N_subsam
      }
    }
    # Add to average across full-field posterior samples if required
    if ((output_ff_mean | output_ff_std) & !is.null(K)){
      out_list$eta_mu_mu = out_list$eta_mu_mu + eta_mu/N_subsam
      out_list$eta_sigma_mu = out_list$eta_sigma_mu + eta_sigma/N_subsam
      if (sam_gp) {
        out_list$eta_sam_mu = out_list$eta_sam_mu + eta_sam/N_subsam
      }
    }
    if ((output_ff_mean | output_ff_std) & !is.null(K_y)){
      out_list$eta_mu_y_mu = out_list$eta_mu_y_mu + eta_mu_y/N_subsam
      out_list$eta_sigma_y_mu = out_list$eta_sigma_y_mu + eta_sigma_y/N_subsam
      if (sam_gp) {
        out_list$eta_sam_y_mu = out_list$eta_sam_y_mu + eta_sam_y/N_subsam
      }
    }
    if (output_ff_std & !is.null(K)){
      out_list$eta_mu_sigma = out_list$eta_mu_sigma + (eta_mu^2)/N_subsam
      if (sam_gp) {
        out_list$eta_sam_sigma = out_list$eta_sam_sigma + (eta_sam^2)/N_subsam
      }
    }
    if (output_ff_std & !is.null(K_y)){
      out_list$eta_mu_y_sigma = out_list$eta_mu_y_sigma + (eta_mu_y^2)/N_subsam
      if (sam_gp) {
        out_list$eta_sam_y_sigma = out_list$eta_sam_y_sigma + (eta_sam_y^2)/N_subsam
      }
    }
  }

  if (output_ff_std & !is.null(K)){
    out_list$eta_mu_sigma = sqrt(out_list$eta_mu_sigma - (out_list$eta_mu_mu^2))
    if (sam_gp){
      out_list$eta_sam_sigma = sqrt(out_list$eta_sam_sigma - (out_list$eta_sam_mu^2))
    }
  }
  if (output_ff_std & !is.null(K_y)){
    out_list$eta_mu_y_sigma = sqrt(out_list$eta_mu_y_sigma - (out_list$eta_mu_y_mu^2))
    if (sam_gp){
      out_list$eta_sam_y_sigma = sqrt(out_list$eta_sam_y_sigma - (out_list$eta_sam_y_mu^2))
    }
  }
  
  return(out_list)
  
}

# The below functions deal with initialising, and updating output lists used to
# pass the requested predictions to the main R code. The contents of this list
# will depend upon values of the Boolean variables passed to the prediction code
initialise_out_list <- function(p=NULL, N_eta = NULL, N_y = NULL, N_post=1, N_pred=1, K=FALSE, K_y=FALSE, sam_gp = FALSE, output_coeff_sam = FALSE, output_ff_sam = FALSE, output_coeff_mean = FALSE, output_ff_mean = FALSE, output_ff_std = FALSE){
  # Initialise output list depended upon values of Boolean variables
  # N_post = number of posterior samples
  # N_pred = number of predictions made by the Gaussian process
  # Consider getting rid of K_y and just running the prediciton code separately
  # Alternative is to sample from w then do the inverse dimension reduction 
  # seperately
  # Possibility to automate further with a for loop?
  # Drop dimension of an array in R??
  
  out_list = list() # Create list
  # Do I want samples of the coefficients?
  if (output_coeff_sam){
    out_list[["w_star_mu"]] = array(0, c(p, N_post, N_pred))
    out_list[["w_star_sigma"]] = array(0, c(p, p, N_post, N_pred))
    if (sam_gp) {
      out_list[["w_star"]] = array(0, c(p, N_post, N_pred))
    }
  }
  
  # Do I want to output posterior samples for full-field response
  if (output_ff_sam & K) {
    out_list[["eta_mu"]] = array(0, c(N_eta, N_post, N_pred))
    out_list[["eta_sigma"]] = array(0, c(N_eta, N_post, N_pred))
    if (sam_gp) {
      out_list[["eta_sam"]] = array(0, c(N_eta, N_post, N_pred))
    }
  }
  if (output_ff_sam & K_y) {
    out_list[["eta_mu_y"]] = array(0, c(N_y, N_post, N_pred))
    out_list[["eta_sigma_y"]] = array(0, c(N_y, N_post, N_pred))
    if (sam_gp) {
      out_list[["eta_sam_y"]] = array(0, c(N_y, N_post, N_pred))
    }
  }
  # Output mean of reduced coefficients across posterior samples?
  if (output_coeff_mean) {
    out_list[["w_star_mu_mu"]] = array(0, c(p, N_pred))
    out_list[["w_star_sigma_mu"]] = array(0, c(p, p, N_pred))
    if (sam_gp) {
      out_list[["w_star_sam_mu"]] = array(0, c(p, N_pred))
    }
  }
  # Do I want to output averages across samples of full-field output?
  # Note: I need the full-field mean to calculate the full field standard 
  # deviation
  if ((output_ff_mean | output_ff_std) & K) {
    out_list[["eta_mu_mu"]] = array(0, c(N_eta, N_pred))
    out_list[["eta_sigma_mu"]] = array(0, c(N_eta, N_pred))
    if (sam_gp) {
      out_list[["eta_sam_mu"]] = array(0, c(N_eta, N_pred))
    }
  }
  if ((output_ff_mean | output_ff_std) & K_y) {
    out_list[["eta_mu_y_mu"]] = array(0, c(N_y, N_pred))
    out_list[["eta_sigma_y_mu"]] = array(0, c(N_y, N_pred))
    if (sam_gp) {
      out_list[["eta_sam_y_mu"]] = array(0, c(N_y, N_pred))
    }
  }
  # Do I want to output standard deviation across posterior distribution?
  if (output_ff_std & K) {
    out_list[["eta_mu_sigma"]] = array(0, c(N_eta, N_pred))
    if (sam_gp) {
      out_list[["eta_sam_sigma"]] = array(0, c(N_eta, N_pred))
    }
  }
  if (output_ff_std & K_y) {
    out_list[["eta_mu_y_sigma"]] = array(0, c(N_y, N_pred))
    if (sam_gp) {
      out_list[["eta_sam_y_sigma"]] = array(0, c(N_y, N_pred))
    }
  }
  
  # Collapse last dimension if not needed
  for (i in 1:length(out_list)){
    if (dim(out_list[[i]])[length(dim(out_list[[i]]))] == 1){
      dim(out_list[[i]]) <- c(dim(out_list[[i]])[-length(dim(out_list[[i]]))])
    }
  }
  
  return(out_list)
}