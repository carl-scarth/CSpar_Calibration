# Code for making predictions using fitted Gaussian process emulators, and 
# calibrated models

library(MASS) # Needed for mvrnorm
source("source/covariance_matrices.R")
source("source/dimension_reduction.R")

zero_mean_gp_pred <- function(y, K_x, k_x_xstar, k_xstar, sam_gp = FALSE, N_sam = 1, inv_K = FALSE, nugget = FALSE){
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

full_field_gp_pred <- function(N_post_pred, x_train, z_train, x_pred, beta_w, lambda_w, lambda_eta, K ,KTKinv, sam_gp = FALSE, output_coeff_sam = TRUE, output_ff_sam = TRUE, output_coeff_mean = TRUE, output_ff_mean = TRUE){
  # Make N_post_pred predictions from a Gaussian process emulator fitted to full
  # field model output, using the method described by Higdon et al.
  # Sub-sampling from the posterior can help reduce required memory, which can 
  # quickly lead to errors when dealing with very high-dimensional output
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
  
  N_sam_post = nrow(beta_w) # Determine overall number of posterior samples
  N_pred = nrow(x_pred) # Determine number of points at which predictions are needed
  N_eta = nrow(K) # Determine number of output values per simulation
  p = ncol(lambda_w) # Determine number of bases
  
  
  # Sub-sample from the posterior distribution
  subsam_ind = sample.int(N_sam_post, N_post_pred) # Determine index of posterior samples
  beta_w = beta_w[subsam_ind,]
  lambda_w = lambda_w[subsam_ind,]
  lambda_eta = lambda_eta[subsam_ind]
  
  # Initialise output arrays
  # Do I want to output posterior samples for basis coefficients?
  if (output_coeff_sam){
    w_star_mu_out = array(0, c(p, N_post_pred, N_pred))
    w_star_sigma_out = array(0, c(p, p, N_post_pred, N_pred))
    if (sam_gp) {
      w_star_out = array(0, c(p, N_post_pred, N_pred))
    } else {
      w_star_out = NULL
    }
  }
  # Do I want to output posterior samples for full-field response
  if (output_ff_sam){
    eta_mu_out = array(0, c(N_eta, N_post_pred, N_pred))
    eta_sigma_out = array(0, c(N_eta, N_post_pred, N_pred))
    if (sam_gp) {
      eta_sam_out = array(0, c(N_eta, N_post_pred, N_pred))
    } else {
      eta_sam_out = NULL
    }
  }
  # Do I want to output averages across samples of basis coefficients?
  if (output_coeff_mean) {
    w_star_mu_mu_out = matrix(0, p, N_pred)
    w_star_sigma_mu_out = array(0, c(p, p, N_pred))
    if (sam_gp) {
      w_star_sam_mu_out = matrix(0, p, N_pred)
    }
  }
  # Do I want to output averages across samples of full-field output?
  if (output_ff_mean) {
    eta_mu_mu_out = matrix(0, N_eta, N_pred)
    eta_sigma_mu_out = matrix(0, N_eta, N_pred)
    if (sam_gp) {
      eta_sam_mu_out = matrix(0, N_eta, N_pred)
    }
  }

  # Loop over all GP samples
  for (i in 1:N_post_pred){
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
        w_star_mu_out[,i,j] = w_star_mu
        w_star_sigma_out[,,i,j] = w_star_sigma
        if (sam_gp) {
          w_star_out[,i,j] = w_star
        }
      }      
    
      # Store full-field posterior prediction samples if requested
      eta_mu = K %*% w_star_mu
      eta_sigma = sqrt(as.matrix(rowSums((K %*% w_star_sigma) * K)))
      if (sam_gp) {
        eta_sam = K %*% w_star
      }
      if (output_ff_sam){
        eta_mu_out[,i,j] = eta_mu
        # Only output diagonals of covariance as full matrix is too high-
        # dimensional. Take square-root for standard
        eta_sigma_out[,i,j] = eta_sigma
        if (sam_gp) {
          eta_sam_out[,i,j] = eta_sam
        }
      }
    
      # If averages are required across posterior samples, keep a running total
      # to reduce memory requirements, compared with calculating mean at the end
      if (output_coeff_mean) {
        w_star_mu_mu_out[,j] = w_star_mu_mu_out[,j] + w_star_mu/N_post_pred
        w_star_sigma_mu_out[,,j] = w_star_sigma_mu_out[,,j] + w_star_sigma/N_post_pred # Might be a problem addding arrays, may need to use apply?
        if (sam_gp) {
          w_star_sam_mu_out[,j] = w_star_sam_mu_out[,j] + w_star/N_post_pred
        }
      }
      # Add to average across full-field posterior samples if required
      if (output_ff_mean){
        eta_mu_mu_out[,j] = eta_mu_mu_out[,j] + eta_mu/N_post_pred
        eta_sigma_mu_out[,j] = eta_sigma_mu_out[,j] + eta_sigma/N_post_pred
        if (sam_gp) {
          eta_sam_mu_out[,j] = eta_sam_mu_out[,j] + eta_sam/N_post_pred
        }
      }
    }
  }

  # create list of the requested outputs
  out_list = list()
  # Add items to list
  if (output_coeff_sam){
    out_list[["w_star_mu"]] = w_star_mu_out
    out_list[["w_star_sigma"]] = w_star_sigma_out
    if (sam_gp){
      out_list[["w_star"]] = w_star_out
    }
  }
  if (output_ff_sam){
    out_list[["eta_mu"]] = eta_mu_out
    out_list[["eta_sigma"]] = eta_sigma_out
    if (sam_gp){
      out_list[["eta_sam"]] = eta_sam_out
    }
  }
  if (output_coeff_mean){
    out_list[["w_star_mu_mu"]] = w_star_mu_mu_out
    out_list[["w_star_sigma_mu"]] = w_star_sigma_mu_out
    if (sam_gp){
      out_list[["w_star_sam_mu_out"]] = w_star_sam_mu_out
    }
  }
  if (output_ff_mean){
    out_list[["eta_mu_mu"]] = eta_mu_mu_out
    out_list[["eta_sigma_mu"]] = eta_sigma_mu_out
    if (sam_gp){
      out_list[["eta_sam_mu"]] = eta_sam_mu_out
    }
  }
  
  return(out_list)
}

# Make a prediction replicating the full-field test data for a single condition
# using a calibrated Gaussian process, with emulator parameters fixed to 
# constant values
full_field_calibration_pred_fixed_em <- function(N_post_pred, tc, tf, z_hat, beta_w, lambda_w, lambda_eta, lambda_y, KTKinv, BTWyBinv, K = NULL, K_y = NULL, sam_GP = FALSE, nugget = FALSE, output_coeff_sam = TRUE, output_ff_sam = TRUE){
  # N_post_pred = 
  # tf = 
  # WRITE LIST OF INPUTS
  # COPY AND PASTE CODE HERE. WHERE THERE IS REPETIION WITH THE ABOVE CODE,
  # PACKAGE INTO SEPARATE FUNCTIONS
  N_sam_post = nrow(tf) # Determine overall number of posterior samples
  N_eta = nrow(KTKinv) # Determine number of output values per simulation
  N_y = nrow(BTWyBinv) # Determine number of experimental data points
  p_eta = ncol(lambda_w) # Determine number of bases (this might not work for fixed point)
  print(p_eta)
  
  # Sub-sample from the posterior distribution
  subsam_ind = sample.int(N_sam_post, N_post_pred) # Determine index of posterior samples
  # do this 
  tf = tf[sam_ind, ]
  lambda_y = lambda_y[subsam_ind]
  
  # Initialise output arrays
  # Do I want to output posterior samples for basis coefficients?
  if (output_coeff_sam){
    w_star_mu_out = matrix(0, p_eta, N_post_pred)
    w_star_sigma_out = array(0, c(p_eta, p_eta, N_post_pred))
    if (sam_gp) {
      w_star_out = matrix(0, p_eta, N_post_pred)
    } else {
      w_star_out = NULL
    }
  }

  # Do I want to output posterior samples for full-field response
  if (output_ff_sam){
    eta_mu_out = matrix(0, N_eta, N_post_pred)
    eta_sigma_out = matrix(0, N_eta, N_post_pred)
    if (sam_gp) {
      eta_sam_out = matrix(0, N_eta, N_post_pred)
      eta_y_out = matrix(0, n_y, N_post_pred)
    } else {
      eta_sam_out = NULL
    }
  }
  # Commented for now until I figure out what I've already implemented
  # Do I want to output averages across samples of basis coefficients?
  # if (output_coeff_mean) {
  #   w_star_mu_mu_out = matrix(0, p, N_pred)
  #   w_star_sigma_mu_out = array(0, c(p, p, N_pred))
  #  if (sam_gp) {
  #    w_star_sam_mu_out = matrix(0, p, N_pred)
  #  }
  #}
  # Do I want to output averages across samples of full-field output?
  #if (output_ff_mean) {
  #  eta_mu_mu_out = matrix(0, N_eta, N_pred)
  #  eta_sigma_mu_out = matrix(0, N_eta, N_pred)
  #  if (sam_gp) {
  #    eta_sam_mu_out = matrix(0, N_eta, N_pred)
  #  }
  #}
  
  # Covariance matrix for the emulator at the experimental data points. This 
  # reduces to a diagonal matrix for a single experiment
  Sigma_u = diag(1/lambda_w, nrow=p_eta, ncol=p_eta)
  # Define the covariance matrix for the emulator weights, evaluated at the
  # training data points.
  Sigma_w = full_field_cov(tc, lambda_w, beta_w)
  View(Sigma_w)
  # Loop over selected posterior samples
  for (i in 1:N_post_pred){
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
    Sigma_z_w_star = rbind(sigma_u, sigma_w_w_star)
    # Determine auto-covariance of the prediction
    Sigma_w_star = Sigma_u
    # Determine GP mean and covariance in reduced-dimensional space, and if 
    # required sample from the GP
    pred = zero_mean_gp_pred(z_hat, Sigma_z_hat, Sigma_z_w_star, sigma_w_star, sam_gp = sam_gp, nugget = nugget)
    w_star_mu = pred[[1]]
    w_star_sigma = pred[[2]]
    w_star = pred[[3]]
    View(w_star_mu)
    View(w_star_sigma)
    View(w_star)
    
    # Store emulator posterior samples of basis coefficients if requested
    if (output_coeff_sam){
      w_star_mu_out[,i] = w_star_mu
      w_star_sigma_out[,,i] = w_star_sigma
      if (sam_gp) {
        w_star_out[,i,j] = w_star
      }
    }
    
    # Convert from reduced-dimensional space to full-field, output at 
    # experimental data points, or both, depending upon which basis matrices
    # have been provided
    if (!is.null(K)){
      eta_mu = K %*% w_star_mu
      eta_sigma = sqrt(as.matrix(rowSums((K %*% w_star_sigma) * K)))
      if (sam_gp) {
        eta_sam = K %*% w_star
      }
    }
    if (!is.null(K_y)){
      eta_mu_y = K_y %*% w_star_mu
      eta_sigma_y = sqrt(as.matrix(rowSums((K_y %*% w_star_sigma) * K_y)))
      if (sam_GP) {
        eta_sam_y = K_y %*% w_star
      }
    }
    # SADSADSDSDSADSADSADSDSADSADSADSADSADSADSADSADSA
    # Continue below and check with above, will need to uncomment some stuff above
    # add in K_y
    # Then check if it works under certain conditions
    # will need to correct below indices as well
    # Consider reinstalling stan to get round import issues
    # Store full-field samples if requested
    if (output_ff_sam){
      eta_mu_out[,i,j] = eta_mu
      # Only output diagonals of covariance as full matrix is too high-
      # dimensional. Take square-root for standard
      eta_sigma_out[,i,j] = eta_sigma
      if (sam_gp) {
        eta_sam_out[,i,j] = eta_sam
      }
    }
    sddfsadsdsad
    
    # If averages are required across posterior samples, keep a running total
    # to reduce memory requirements, compared with calculating mean at the end
    if (output_coeff_mean) {
      w_star_mu_mu_out[,j] = w_star_mu_mu_out[,j] + w_star_mu/N_post_pred
      w_star_sigma_mu_out[,,j] = w_star_sigma_mu_out[,,j] + w_star_sigma/N_post_pred # Might be a problem addding arrays, may need to use apply?
      if (sam_gp) {
        w_star_sam_mu_out[,j] = w_star_sam_mu_out[,j] + w_star/N_post_pred
      }
    }
    # Add to average across full-field posterior samples if required
    if (output_ff_mean){
      eta_mu_mu_out[,j] = eta_mu_mu_out[,j] + eta_mu/N_post_pred
      eta_sigma_mu_out[,j] = eta_sigma_mu_out[,j] + eta_sigma/N_post_pred
      if (sam_gp) {
        eta_sam_mu_out[,j] = eta_sam_mu_out[,j] + eta_sam/N_post_pred
      }
    }
    
    # Copied from calibration code...
    eta_sam[,i] = (K_eta %*% w_star[,i])*sd_dt + mu_dt
    eta_y[,i] = (K_y %*% w_star[,i])*sd_dt + mu_y
  }
  
  # Think about whether it's necessary to produce outputs for y in this code, 
  # or if it's better to just have a separate prediction code using K_y in 
  # place of K_eta. This will depend on if I need both to make predicitons, 
  # which I can't remember   
}