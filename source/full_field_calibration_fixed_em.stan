// Bayesian model calibration code using full-field model output and expeimental
// data using the method from:
// D. Higdon et al, "Computer Model Calibration Using High-Dimensional Output",
// Journal of the American Statistical Association,2008.
// This is a simplified code wherein there is no discrepancy, only one 
// experiment, and the emulator parameters are fixed.

#include covariance_matrices.stan

data {
  int<lower=0> m;               // number of computer simulations
  int<lower=0> q;               // number of uncontrolled inputs, t
  int<lower=0> n_eta;           // number of output values per simulation
  int<lower=0> n_y;             // number of experimental observations (per test)
  int<lower=0> p_eta;           // number of basis functions used for emulator
  real<lower=0> a_y_dash;       // Adjusted value of the shape parameter for the lambda_y prior in accordance with the normal-gamma model of Higdon et al.
  real<lower=0> b_y_dash;       // Adjusted value of the rate parameter for the lambda_y prior in accordance with the normal-gamma model of Higdon et al.
  real<lower=0> lambda_eta;     // Precision parameter of truncation error
  vector[(m+1)*p_eta] z_hat; // Vector of regression weights, representing both model output for the experimental data and emulator training data, as defined in Higdon et al.
  vector[q] tf_param_1;              // prior means of calibration parameters
  vector[q] tf_param_2;           // prior standard deviations of calibration parameters
  vector<lower=0,upper=1>[q*p_eta] rho_w; // Correlation length of emulator
  vector<lower=0>[p_eta] lambda_w; // Precision of emulator
  // real<lower=0> lambda_w; // Precision of emulator
  matrix[m, q] tc;              // Matrix of uncontrolled inputs at the training data points
  matrix[m*p_eta,m*p_eta] KTKinv; // Inverse of the inner product of the matrix of basis functions used to approximate model output, K, with themselves. See Higdon et al. for details on K.
  matrix[p_eta,p_eta] BTWyBinv; // Inverse of product B'*W_y*B, where B is a matrix containing emualtor basis terms, as outlined in Higdon et al., and W_y is the observation error precision matrix
}

transformed data {
  vector[(m+1)*p_eta] mu_z_hat;                // prior mean vector
  row_vector[q*p_eta] beta_w; // Transformed correlation lengths

  // Assemble joint vector of experimental and model outputs
  mu_z_hat = rep_vector(0, (m+1)*p_eta);  // set prior mean vector to zero
  beta_w = -4.0 * log(rho_w');
}

parameters {
  // tf: calibration parameters
  // lambda_y: precision parameter associated with experimental error
  // bounding tf prevents the sampler from extrapolating a long way outside of 
  // the training data, which causes convergence problems. For a uniform
  // -distributed prior this may make sense, however, reasonable bounds should 
  // be chosen for a normal prior (say [-0.1,1.1]), in conjunction with a well-
  // characterised prior
  row_vector<lower=-0.1,upper=1.1>[q] tf_gauss; // Gaussian-distributed priors
  // row_vector<lower=0,upper=1>[q-1] tf_gauss; // Gaussian-distributed priors
  real<lower=0> lambda_y;
  // real<lower=tf_param_1[q],upper=tf_param_2[q]> tf_unif; // uniform-distributed prior
}

transformed parameters {
  // cacluate correlation length from transformed parameter rho_eta
  row_vector[q] tf; // Vector containing both uniform and Gaussian distributed priors
  
  // tf = append_col(tf_gauss, tf_unif);
  tf = tf_gauss;
}

model {
  // declare variables
  matrix[p_eta, p_eta] sigma_u;     // Covariance matrix for emulator weights, evaluated at experimental data points
  matrix[m*p_eta, m*p_eta] sigma_w; // Covariance matrix for emulator weights, evaluated at training data points
  matrix[p_eta, m*p_eta] sigma_uw;  // Cross-covariance matrix for emulator weights between training data and experimental data points
  matrix[(m+1)*p_eta, (m+1)*p_eta] sigma_z; // Covariance matrix for joint experimental and modelling data
  matrix[(m+1)*p_eta, (m+1)*p_eta] sigma_z_hat; // Covariance matrix for joint experimental and model data, adjusted with model and emulator error terms transformed into low-dimensional space
  matrix[(m+1)*p_eta, (m+1)*p_eta] L_z_hat;             // Cholesky decomposition of covariance matrix
    
  // Define the covariance matrix for the emulator weights, evaluated at the experimental data points.
  // As we only have one experiment, and R(theta,theta) = 1, this also reduces to a diagonal matrix.
  sigma_u = diag_matrix(rep_vector(1,p_eta)./lambda_w);
  // print(lambda_w);
  // print("");
  // print(sigma_u);
  // print("");
  // print(beta_w);
  // print("");
  
  // Define the covariance matrix for the emulator weights, evaluated at the training data points. Note
  // that this works differently to my full-field emulator code, which defined the relationship between
  // inputs and outputs within a loop. This won't be possible here as the outputs are provided in their
  // high-dimensional format. It may be possible to slightly speed this code up by at least breaking up the
  // relationship between zs and u, v, and w into smaller chunks according to the non-zero regions of the
  // covariance matrix
  sigma_w = rep_matrix(0,m*p_eta,m*p_eta);
  for (i in 1:p_eta) {
    // Calculate the covariance matrix of training data points for the ith emulator weight
    sigma_w[(i-1)*m+1:i*m,(i-1)*m+1:i*m] = ARD_SE_cov(tc, lambda_w[i], beta_w[(i-1)*q+1:i*q], 0);
    //sigma_w[(i-1)*m+1:i*m,(i-1)*m+1:i*m] = ARD_SE_cov(tc, lambda_w, beta_w[(i-1)*q+1:i*q], 0);
  }
  // print(sigma_w);
  // print("");
  
  // Evaluate cross-covariance of emulator weights between training data and experimental data points
  sigma_uw = rep_matrix(0,p_eta,m*p_eta);
  // print(tf);
  // print("");
  for (i in 1:p_eta) {
    sigma_uw[i,(i-1)*m+1:i*m] = to_row_vector(ARD_SE_cov_non_sym(to_matrix(tf), tc, lambda_w[i], beta_w[(i-1)*q+1:i*q]));
    //sigma_uw[i,(i-1)*m+1:i*m] = to_row_vector(ARD_SE_cov_non_sym(tf, tc, lambda_w, beta_w[(i-1)*q+1:i*q]));
  }
  // print(sigma_uw);
  // print("");
  
  // Assemble covariance matrix for joint experimental and model data
  sigma_z = rep_matrix(0,(m+1)*p_eta, (m+1)*p_eta);
  sigma_z[1:p_eta,1:p_eta] = sigma_u;
  sigma_z[p_eta+1:,p_eta+1:] = sigma_w;
  sigma_z[1:p_eta,p_eta+1:] = sigma_uw;
  sigma_z[p_eta+1:,1:p_eta] = sigma_uw';

  // Adjust the covariance matrix sigma_z to include transformed emulator and experimental error terms
  sigma_z_hat = rep_matrix(0,(m+1)*p_eta, (m+1)*p_eta);
  sigma_z_hat[1:p_eta,1:p_eta] = BTWyBinv/lambda_y;
  sigma_z_hat[p_eta+1:,p_eta+1:] = KTKinv/lambda_eta;
  
  // print(BTWyBinv);
  // print("");
  // print(lambda_y);
  // print("");
  // print(diagonal(sigma_z));
  // print("");
  // print(diagonal(sigma_z_hat));
  // print("");
  
  sigma_z_hat = sigma_z_hat + sigma_z;
  // print(diagonal(sigma_z_hat));
  // print("");
  // Add small nugget
  sigma_z_hat = sigma_z_hat + diag_matrix(rep_vector(1e-6,(m+1)*p_eta));
  // print(diagonal(sigma_z_hat));
  // print("");
  
  // print(sigma_z_hat);
  // print("");
  
  // Perform Cholesky decomposition of sigma_z_hat for efficiency
  // either force symmetry somehow using cholesky_decompose function,
  L_z_hat = cholesky_decompose(sigma_z_hat);
  
  // Specify prior distribution of z_hat
  z_hat ~ multi_normal_cholesky(mu_z_hat, L_z_hat);
  
  // Specify priors here
  lambda_y ~ gamma(a_y_dash, b_y_dash);               // Precision parameter for observation error, gamma (shape, rate)
  // define priors on calibration parameters
  for (i in 1:q){
    tf_gauss[i] ~ normal(tf_param_1[i], tf_param_2[i]);
  }
  // Remaining calibration parameters are uniformly distributed by default, between previously defined bounds
}
