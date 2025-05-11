// Bayesian model calibration code using full-field model output and expeimental
// data using the method from:
// D. Higdon et al, "Computer Model Calibration Using High-Dimensional Output",
// Journal of the American Statistical Association,2008.
// This is a simplified code wherein there is no discrepancy, only one 
// experiment, and the emulator parameters are fixed.

#include /covariance_matrices.stan
data {
  int<lower=0> m;                // number of computer simulations
  int<lower=0> q;                // number of uncontrolled inputs, t
  int<lower=0> n_eta;            // number of output values per simulation
  int<lower=0> n_y;              // number of experimental observations (per test)
  int<lower=0> p_eta;            // number of emulator basis functions
  int<lower=0> q_y;              // number of components of output vector
  real<lower=0> a_y_dash[q_y];   // Adjusted value of the shape parameter for the lambda_y prior
  real<lower=0> b_y_dash[q_y];   // Adjusted value of the rate parameter for the lambda_y prior
  real<lower=0> lambda_eta;      // Precision of truncation error
  vector[m*p_eta] w_hat;     // Vector of regression weights representing model output for the experimental data and emulator training data
  row_vector[q] tf_param_1;          // First parameter of calibration parameter priors
  row_vector[q] tf_param_2;          // Second parameter of calibration parameter priors
  vector<lower=0,upper=1>[q*p_eta] rho_w; // Emulator correlation parameters
  vector<lower=0>[p_eta] lambda_w;        // Emulator precision parameters
  //real<lower=0> lambda_w;        // Emulator precision parameters
  matrix[m, q] tc;                        // Matrix of uncontrolled inputs at the training data points
  matrix[m*p_eta,m*p_eta] KTKinv;         // Inverse of the inner product of the emulator basis matrix
  matrix[p_eta,p_eta] BTB[q_y];           // Product B'*B, where B is the emulator basis matrix interpolated to the experimental data points
  matrix[p_eta,q_y] BTy;                  // Product of B'*y, where y is the experimental data
}

transformed data {
  vector[(m+1)*p_eta] mu_z_hat; // prior mean vector
  row_vector[q*p_eta] beta_w;   // Transformed correlation parameters

  mu_z_hat = rep_vector(0, (m+1)*p_eta); // set prior mean vector to zero
  beta_w = -4.0 * log(rho_w');
}

parameters {
  // lambda_y: precision parameter associated with experimental error
  vector<lower=0>[q_y] lambda_y;
  // For now, define calibration parameters via an external stan file, which
  // can be written in R depending on the required priors. An alternative method
  // would be via optional variables, though this seems too rigid
#include /calibration_parameters.stan
}

transformed parameters {
  // cacluate correlation length from transformed parameter rho_eta
  // stan file containing instructions on how to concatenate the different 
  // calibration parameters defined in parameters, into single vector tf. Will
  // depend upon the chosen priors
  row_vector[q] tf; // Vector containing both uniform and Gaussian distributed priors
#include /concat_parameters.stan
}

model {
  vector[p_eta] BTeps_yy;                       // Sum of inner prodcts of B_iTy_i and lambda_y
  vector[p_eta] u_hat;                          // Reduced dimensional representation of experimental data
  vector[(m+1)*p_eta] z_hat;                    // Combined experimental and model data
  matrix[p_eta, p_eta] sigma_u;                 // Covariance matrix for emulator at experimental data points
  matrix[m*p_eta, m*p_eta] sigma_w;             // Covariance matrix for emulator at training data points
  matrix[p_eta, m*p_eta] sigma_uw;              // Emulator cross-covariance  training data and experimental data points
  matrix[(m+1)*p_eta, (m+1)*p_eta] sigma_z;     // Joint covariance matrix for joint experimental and model data
  matrix[(m+1)*p_eta, (m+1)*p_eta] sigma_z_hat; // Covariance matrix for joint experimental and model data, adjusted with model and emulator error terms transformed into low-dimensional space
  matrix[p_eta, p_eta] prec_y;                  // Observation error precision in the reduced dimensional space
  matrix[(m+1)*p_eta, (m+1)*p_eta] L_z_hat;     // Cholesky decomposition of covariance
  
  // Calculate the emulator auto-covariance for the experimental data.
  // This is diagonal for a single experimental data point
  sigma_u = diag_matrix(rep_vector(1,p_eta)./lambda_w);
  
  // Emulator auto-covariance for the training data points
  sigma_w = rep_matrix(0,m*p_eta,m*p_eta);
  for (i in 1:p_eta) {
    // Calculate the covariance matrix for the ith basis coefficient
    sigma_w[(i-1)*m+1:i*m,(i-1)*m+1:i*m] = ARD_SE_cov(tc, lambda_w[i], beta_w[(i-1)*q+1:i*q], 0);
    //sigma_w[(i-1)*m+1:i*m,(i-1)*m+1:i*m] = ARD_SE_cov(tc, lambda_w, beta_w[(i-1)*q+1:i*q], 0);
  }
  
  // Emulator cross covariance between training data and experimental data.
  sigma_uw = rep_matrix(0,p_eta,m*p_eta);
  for (i in 1:p_eta) {
    sigma_uw[i,(i-1)*m+1:i*m] = to_row_vector(ARD_SE_cov_vect_mat(tf, tc, lambda_w[i], beta_w[(i-1)*q+1:i*q]));
    //sigma_uw[i,(i-1)*m+1:i*m] = to_row_vector(ARD_SE_cov_vect_mat(tf, tc, lambda_w, beta_w[(i-1)*q+1:i*q]));
  }
  
  // Assemble covariance matrix for joint experimental and model data
  sigma_z = rep_matrix(0,(m+1)*p_eta, (m+1)*p_eta);
  sigma_z[1:p_eta,1:p_eta] = sigma_u;
  sigma_z[p_eta+1:,p_eta+1:] = sigma_w;
  sigma_z[1:p_eta,p_eta+1:] = sigma_uw;
  sigma_z[p_eta+1:,1:p_eta] = sigma_uw';

  // Adjust the covariance matrix to include transformed emulator and 
  // experimental error terms
  BTeps_yy = rep_vector(0,p_eta);
  prec_y = rep_matrix(0,p_eta,p_eta);
  sigma_z_hat = rep_matrix(0,(m+1)*p_eta, (m+1)*p_eta);
  // Calculate inner product sum of B_iT*B_i*lambda_i and take inverse
  // Also sum of B_iT*y_i*lambda_i
  for (i in 1:q_y){
    BTeps_yy = BTeps_yy + lambda_y[i]*BTy[:,i];
    prec_y = prec_y + lambda_y[i]*BTB[i];
  }
  prec_y = inverse_spd(prec_y);
  // Add contribution to observation error for the ith vector component
  sigma_z_hat[1:p_eta,1:p_eta] = prec_y;
  sigma_z_hat[p_eta+1:,p_eta+1:] = KTKinv/lambda_eta;
  sigma_z_hat = sigma_z_hat + sigma_z;
  
  // Calculate reduced dimensional term for experimental data
  u_hat = prec_y*BTeps_yy;
  z_hat = append_row(u_hat, w_hat);
  
  // Add small nugget
  sigma_z_hat = sigma_z_hat + diag_matrix(rep_vector(1e-8,(m+1)*p_eta));
  
  // Specify prior distribution of z_hat
  L_z_hat = cholesky_decompose(sigma_z_hat);
  z_hat ~ multi_normal_cholesky(mu_z_hat, L_z_hat);
  
  // Specify priors on model hyperparameters
  // Precision of observation error, gamma (shape, rate)
  for (i in 1:q_y){
    lambda_y[i] ~ gamma(a_y_dash[i], b_y_dash[i]);
  }
  // define priors on calibration parameters
#include calibration_priors.stan
}
