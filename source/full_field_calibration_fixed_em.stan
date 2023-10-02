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
  int<lower=0> p_eta;           // number of emulator basis functions
  real<lower=0> a_y_dash;       // Adjusted value of the shape parameter for the lambda_y prior
  real<lower=0> b_y_dash;       // Adjusted value of the rate parameter for the lambda_y prior
  real<lower=0> lambda_eta;     // Precision of truncation error
  vector[(m+1)*p_eta] z_hat;    // Vector of regression weights representing model output for the experimental data and emulator training data
  vector[q] tf_param_1;         // First parameter of calibration parameter priors
  vector[q] tf_param_2;         // Second parameter of calibration parameter priors
  vector<lower=0,upper=1>[q*p_eta] rho_w; // Emulator correlation parameters
  vector<lower=0>[p_eta] lambda_w;        // Emulator precision parameters
  matrix[m, q] tc;                        // Matrix of uncontrolled inputs at the training data points
  matrix[m*p_eta,m*p_eta] KTKinv;         // Inverse of the inner product of the emulator basis matrix
  matrix[p_eta,p_eta] BTWyBinv;           // Inverse of B'*W_y*B, where B is the emulator basis matrix interpolated to the experimental data points, and W_y is the observation error precision
}

transformed data {
  vector[(m+1)*p_eta] mu_z_hat; // prior mean vector
  row_vector[q*p_eta] beta_w;   // Transformed correlation parameters

  mu_z_hat = rep_vector(0, (m+1)*p_eta); // set prior mean vector to zero
  beta_w = -4.0 * log(rho_w');
}

parameters {
  // lambda_y: precision parameter associated with experimental error
  real<lower=0> lambda_y;
  // For now, define calibration parameters via an external stan file, which
  // can be written in R depending on the required priors. An alternative method
  // would be via optional variables, though this seems too rigid
  #include calibration_parameters.stan
}

transformed parameters {
  // cacluate correlation length from transformed parameter rho_eta
  // stan file containing instructions on how to concatenate the different 
  // calibration parameters defined in parameters, into single vector tf. Will
  // depend upon the chosen priors
  row_vector[q] tf; // Vector containing both uniform and Gaussian distributed priors
  #include concat_parameters.stan
}

model {
  matrix[p_eta, p_eta] sigma_u;                 // Covariance matrix for emulator at experimental data points
  matrix[m*p_eta, m*p_eta] sigma_w;             // Covariance matrix for emulator at training data points
  matrix[p_eta, m*p_eta] sigma_uw;              // Emulator cross-covariance  training data and experimental data points
  matrix[(m+1)*p_eta, (m+1)*p_eta] sigma_z;     // Joint covariance matrix for joint experimental and model data
  matrix[(m+1)*p_eta, (m+1)*p_eta] sigma_z_hat; // Covariance matrix for joint experimental and model data, adjusted with model and emulator error terms transformed into low-dimensional space
  matrix[(m+1)*p_eta, (m+1)*p_eta] L_z_hat;     // Cholesky decomposition of covariance
  
  // I'm happy with my comments up to this point
  // Define the covariance matrix for the emulator weights, evaluated at the experimental data points.
  // As we only have one experiment, and R(theta,theta) = 1, this also reduces to a diagonal matrix.
  sigma_u = diag_matrix(rep_vector(1,p_eta)./lambda_w);
  
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
  
  // Evaluate cross-covariance of emulator weights between training data and experimental data points
  sigma_uw = rep_matrix(0,p_eta,m*p_eta);
  // print(tf);
  // print("");
  for (i in 1:p_eta) {
    sigma_uw[i,(i-1)*m+1:i*m] = to_row_vector(ARD_SE_cov_vect_mat(tf, tc, lambda_w[i], beta_w[(i-1)*q+1:i*q]));
  }
  
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
