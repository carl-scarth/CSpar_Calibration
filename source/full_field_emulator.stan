// Code for fitting a Gaussian process emulator to full-field model output using
// Bayesian inferrence using the method from:
// D. Higdon et al, "Computer Model Calibration Using High-Dimensional Output",
// Journal of the American Statistical Association,2008.
// This code contains only the emulator component of the framework.

#include covariance_matrices.stan

data {
  int<lower=0> m;                 // number of computer simulations
  int<lower=0> q;                 // number of inputs, t
  int<lower=0> n_eta;             // number of output values per simulation
  int<lower=0> p_eta;             // number of basis functions
  real<lower=0> a_eta_dash;       // Adjusted value of the shape parameter for lambda_eta prior
  real<lower=0> b_eta_dash;       // Adjusted value of the rate parameter for the lambda_eta prior
  row_vector[m*p_eta] z_hat;      // Vector of regression weights representing model output
  matrix[m, q] tc;                // Matrix of inputs at the training data points
  matrix[m*p_eta,m*p_eta] KTKinv; // Inverse of the inner product of the emulator basis matrix
}

transformed data {
  vector[m*p_eta] mu_z_hat; // prior mean vector

  mu_z_hat = rep_vector(0, m*p_eta);  // set prior mean to zero
}

parameters {
  row_vector<lower=0,upper=1>[q*p_eta] rho_w; // reparameterisation of correlation length, beta_w
  vector<lower=0>[p_eta] lambda_w; // precision for emulator
  real<lower=0> lambda_eta; // precision for series truncation error
}

transformed parameters {
  row_vector[q*p_eta] beta_w;
  
  beta_w = -4.0 * log(rho_w); // cacluate correlation length from rho_w
}

model {
  matrix[m*p_eta, m*p_eta] sigma_z;     // Covariance matrix for emulator training data
  matrix[m*p_eta, m*p_eta] sigma_z_hat; // Covariance matrix adjusted for dimension reduction
  matrix[m*p_eta, m*p_eta] L_z_hat;     // Cholesky decomposition of covariance matrix
    
  // Populate the covariance matrix
  sigma_z = rep_matrix(0,m*p_eta,m*p_eta);
  for (i in 1:p_eta) {
    // Calculate the covariance matrix for the ith basis coefficient
    sigma_z[(i-1)*m+1:i*m,(i-1)*m+1:i*m] = ARD_SE_cov(tc, lambda_w[i], beta_w[(i-1)*q+1:i*q], 0);
  }
  
  // Adjust the covariance matrix to account for the dimension reduction
  sigma_z_hat = rep_matrix(0,m*p_eta, m*p_eta);
  sigma_z_hat = KTKinv/lambda_eta;
  sigma_z_hat = sigma_z_hat + sigma_z;
  
  // Specify prior distribution of z_hat
  L_z_hat = cholesky_decompose(sigma_z_hat);
  z_hat ~ multi_normal_cholesky(mu_z_hat, L_z_hat);
  
  // Specify priors on emulator hyperparameters
  rho_w[1:p_eta*q] ~ beta(1.0, 0.1);
  lambda_w[1:p_eta] ~ gamma(5.0, 5.0);        // gamma (shape, rate) (from Higdon)
  //lambda_w[1:p_eta] ~ gamma(10.0, 10.0);    // gamma (shape, rate) (from Chong)
  lambda_eta ~ gamma(a_eta_dash, b_eta_dash); // gamma (shape, rate}.
}
