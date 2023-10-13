// Code for fitting a Gaussian process emulator using Bayesian inferrence, and making subsequent predictions
// Modified from: https://github.com/adChong/bc-stan, [1]
// prediction code modified from: 
// https://betanalpha.github.io/assets/case_studies/gaussian_processes.html#21_Simulating_From_A_Gaussian_Process [2]
// The main modification is to strip out the experimental data from [1], leaving only the emulator.
// efficient posterior predictions are made using the closed-form Gaussian process expressions, with code taken from [2].
// For expressions, see Chapter 2, "Gaussian Processes for Machine Learning", Rasmussen and Williams, MIT Press, 2006, ISBN 0-262-18253-X.
// http://www.gaussianprocess.org/gpml/

functions {
  // Define a function for evaluating a square-exponential covariance matrix with non-isotropic correlation length i.e. there 
  // is a different correlation length for each input (this is necessary as the in-built Stan SE covariance is isotropic)
  matrix ARD_SE_cov(matrix x, real lambda, row_vector beta, real delta) {
    int N = rows(x);    // number of training data points
    int d = cols(x);    // number of inputs
    matrix[N, N] K;     // Covariance matrix
    {
    // Main body of the function
    // Declare variables
    row_vector[d] dx;   // variable used to determine difference between inputs
    K = diag_matrix(rep_vector((1/lambda), N));
    for (i in 1:N) {
      for (j in (i+1):N) {
        dx = x[i] - x[j];
        K[i, j] = (beta .* dx) * dx';
        K[i, j] = exp(-K[i, j])/lambda;
        K[j, i] = K[i, j];
      }
    }
    // add nugget term to diagonals
    K = K + diag_matrix(rep_vector(delta,N));
    }
    return K;
  }
  
  // Define a function for evaluating a square-exponential covariance matrix with non-isotropic correlation length as above,
  // but for a non-symmetric covariance matrix. This is used when evaulating the covariance between two different sets of points.
  matrix ARD_SE_cov_non_sym(matrix x1, matrix x2, real lambda, row_vector beta) {
    int N1 = rows(x1);  // number of training data points
    int N2 = rows(x2);  // number of points at which predictions are required
    int d = cols(x1);   // number of inputs
    matrix[N1, N2] K;   // Declare covariance matrix
    {
    // Main body of the function
    // Declare variables
    row_vector[d] dx;   // variable used to determine difference between inputs
    K = rep_matrix(0,N1,N2);
    for (i in 1:N1) {
      for (j in 1:N2) {
        dx = x1[i] - x2[j];
        K[i, j] = (beta .* dx) * dx';
        K[i, j] = exp(-K[i, j])/lambda;
      }
    }
    }
    return K;
  }
  
  // Define a function for analytical expressions for GP posterior predictions for a zero-mean prior, and squared-exponential covariance
  vector gp_pred_rng(matrix x_t, vector y, matrix x, real lambda, row_vector beta, real delta) {
    int N1 = rows(y);   // number of training data points
    int N2 = rows(x_t); // number of test points at which predictions are required
    int d = cols(x);    // number of inputs
    vector[N2] f_pred;  // declare function output, a vector of predictions
    {
    // Main body of the function
    // Declare variables
    matrix[N1, N1] K;       // Covariance matrix of training data points with themselves
    matrix[N1, N1] L_K;     // Cholesky decomposition of covariance matrix
    vector[N1] L_K_div_y;
    vector[N1] K_div_y;
    matrix[N1,N2] k_x_xt;   // Covariance matrix of training data points with test data points
    vector[N2] f_mu;        // Posterior predictive GP mean
    matrix[N1, N2] v_pred;
    matrix[N2, N2] C_xt_xt; // Covariance matrix of test data points with themselves
    matrix[N2, N2] f_cov;   // Posterior predictive GP covariance
    // Define covariance matrix combining training and test data
    K = ARD_SE_cov(x, lambda, beta, 0);
    // Perform Cholesky decomposition of covariance matrix for subsequent efficient matrix algebra
    L_K  = cholesky_decompose(K);
    // Define covariance matrix of training data with test data. Note that this matrix is not symmetric
    k_x_xt = ARD_SE_cov_non_sym(x, x_t, lambda, beta);
    // Determine covariance of test points in x_t with themselves
    // note: a small nugget term, delta, is added to ensure covariance is positive-semi-definite when many predictions are required
    C_xt_xt = ARD_SE_cov(x_t, lambda, beta, delta);
    // Efficient code for calculating emulator predictive mean, k_x_xt^T * K^-1 * y
    L_K_div_y = mdivide_left_tri_low(L_K, y);
    K_div_y = mdivide_right_tri_low(L_K_div_y', L_K)'; // K^-1 * y
    // Evaluate posterior predictive mean
    f_mu = (k_x_xt' * K_div_y); // posterior mean
    // Now cacluate posterior predictive covariance, C_xt_xt - k_x_xt^T * K^-1 * k_x_xt
    v_pred = mdivide_left_tri_low(L_K, k_x_xt);
    // calcluate posterior predictive covariance
    f_cov = C_xt_xt - v_pred' * v_pred;
    // Sample from Gaussian process
    f_pred = multi_normal_rng(f_mu, f_cov);
    }
    return f_pred;
  }
}

data {
  int<lower=0> m;               // number of computer simulations
  int<lower=0> q;               // number of uncontrolled inputs, t
  int<lower=0> n_eta;           // number of output values per simulation
  int<lower=0> n_y;             // number of experimental observations (per test)
  int<lower=0> p_eta;           // number of basis functions used for emulator
  int<lower=0> p_delta;         // number of basis functions used for discrepancy
  int<lower=0> rank_B;          // Rank of matrix combining basis functions for emulator and discrepancy
  //int<lower=1> n_pred;        // number of points at which predictions are required
  real<lower=0> a_eta_dash;     // Adjusted value of the shape parameter for the lambda_eta prior in accordance with the normal-gamma model of Higdon et al.
  real<lower=0> a_y_dash;       // Adjusted value of the shape parameter for the lambda_y prior in accordance with the normal-gamma model of Higdon et al.
  real<lower=0> b_eta_dash;     // Adjusted value of the rate parameter for the lambda_eta prior in accordance with the normal-gamma model of Higdon et al.
  real<lower=0> b_y_dash;       // Adjusted value of the rate parameter for the lambda_y prior in accordance with the normal-gamma model of Higdon et al.
  vector[rank_B + m*p_eta] z_hat; // Vector of regression weights, representing both model output and discrepancy as defined in Higdon et al.
  vector[q] tf_param_1;         // Parameter 1 of prior distribution
  vector[q] tf_param_2;         // Parameter 2 of prior distribution
  matrix[m, q] tc;              // Matrix of uncontrolled inputs at the training data points
  matrix[m*p_eta,m*p_eta] KTKinv; // Inverse of the inner product of the matrix of basis functions used to approximate model output, K, with themselves. See Higdon et al. for details on K.
  matrix[rank_B,rank_B] BTWyBinv; // Inverse of product B'*W_y*B, where B is a matrix containing emualtor and discrepancy basis terms, as outlined in Higdon et al., and W_y is the observation error precision matrix
  matrix[rank_B,(p_eta+p_delta)] L; // Matrix used to map between combined matrix of emulator and discrepancy basis terms B, and full rank equivalent B_tilde
  
  //matrix[n_pred, q] t_pred; // "uncontrolled" input values at which predictions are required
}

transformed data {
  vector[rank_B+m*p_eta] mu_z_hat;                // prior mean vector

  // Assemble joint vector of experimental and model outputs
  mu_z_hat = rep_vector(0, rank_B+m*p_eta);  // set prior mean vector to zero
}

parameters {
  // tf: calibration parameters
  // rho_eta: reparameterisation of correlation length, beta_eta
  // lambda_eta: precision parameter for eta
  // lambda_delta: precision parameter for delta
  // lambda_y: precision parameter associated with experimental error
  // bounding tf prevents the sampler from extrapolating a long way outside of 
  // the training data, which causes convergence problems. For a uniform
  // -distributed prior this may make sense, however, reasonable bounds should 
  // be chosen for a normal prior (say [-0.1,1.1]), in conjunction with a well-
  // characterised prior
  row_vector<lower=0,upper=1>[q-1] tf_gauss; // Gaussian-distributed priors
  row_vector<lower=0,upper=1>[q*p_eta] rho_w;
  vector<lower=0>[p_eta] lambda_w; 
  vector<lower=0>[p_delta] lambda_v;
  real<lower=0> lambda_y;
  real<lower=0> lambda_eta;
  real<lower=tf_param_1[q],upper=tf_param_2[q]> tf_unif; // uniform-distributed prior
}

transformed parameters {
  // cacluate correlation length from transformed parameter rho_eta
  row_vector[q*p_eta] beta_w;
  row_vector[q] tf; // Vector containing both uniform and Gaussian distributed priors
  
  beta_w = -4.0 * log(rho_w);
  tf = append_col(tf_gauss, tf_unif);
}

model {
  // declare variables
  matrix[p_delta, p_delta] sigma_v; // Covariance matrix for discrepancy weights, evaluated at experimental data points
  matrix[p_eta, p_eta] sigma_u;     // Covariance matrix for emulator weights, evaluated at experimental data points
  matrix[m*p_eta, m*p_eta] sigma_w; // Covariance matrix for emulator weights, evaluated at training data points
  matrix[p_eta, m*p_eta] sigma_uv;  // Cross-covariance matrix for emulator weights between training data and experimental data points
  matrix[p_eta+p_delta,p_eta+p_delta] sigma_vu_join; // joined matrix combining sigma_u and sigma_v for subsequent matrix algebra
  matrix[rank_B, m*p_eta] Lsigma_uv; // Contributions from cross-covariance matrix sigma_uv, transformed to reduced-rank equivalent
  matrix[rank_B+m*p_eta, rank_B+m*p_eta] sigma_z;       // Covariance matrix for joint experimental and modelling data
  matrix[rank_B+m*p_eta, rank_B+m*p_eta] sigma_z_hat;   // Covariance matrix for joint experimental and model data, adjusted with model and emulator error terms transformed into low-dimensional space
  matrix[rank_B+m*p_eta, rank_B+m*p_eta] L_z_hat;             // Cholesky decomposition of covariance matrix
    
  // Define the covariance for the discrepancy weights. Here I assume that each weight is an independent 
  // Gaussian process (i.e. F = p_delta, |G_i| = 1 for all i = 1,...p_delta) following the notation of 
  // Higdon. Note that this is just one option from the more general approach outlined in Higdon. 
  sigma_v = diag_matrix(rep_vector(1,p_delta)./lambda_v);
  
  // Define the covariance matrix for the emulator weights, evaluated at the experimental data points.
  // As we only have one experiment, and R(theta,theta) = 1, this also reduces to a diagonal matrix.
  sigma_u = diag_matrix(rep_vector(1,p_eta)./lambda_w);
  
  // Combine contributions from sigma_u and sigma_v into a single block-diagonal matrix for subsequent transformation
  // onto their full-rank equivalent combining terms from both
  sigma_vu_join = rep_matrix(0,p_eta+p_delta,p_eta+p_delta);
  sigma_vu_join[1:p_delta,1:p_delta] = sigma_v;
  sigma_vu_join[p_delta+1:,p_delta+1:] = sigma_u;
  
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
  }
  
  // Evaluate cross-covariance of emulator weights between training data and experimental data points
  sigma_uv = rep_matrix(0,p_eta,m*p_eta);
  for (i in 1:p_eta) {
    sigma_uv[i,(i-1)*m+1:i*m] = to_row_vector(ARD_SE_cov_non_sym(to_matrix(tf), tc, lambda_w[i], beta_w[(i-1)*q+1:i*q]));
  }
  
  // Transform contributions from sigma uv to into their full rank equivalent
  Lsigma_uv = L*append_row(rep_matrix(0,p_delta,m*p_eta),sigma_uv);
  
  // Assemble covariance matrix for joint experimental and model data
  sigma_z = rep_matrix(0,rank_B+m*p_eta, rank_B+m*p_eta);
  sigma_z[1:rank_B,1:rank_B] = quad_form_sym(sigma_vu_join, L');
  sigma_z[rank_B+1:,rank_B+1:] = sigma_w;
  sigma_z[1:rank_B,rank_B+1:] = Lsigma_uv;
  sigma_z[rank_B+1:,1:rank_B] = Lsigma_uv';
  
  // Adjust the covariance matrix sigma_z to include transformed emulator and experimental error terms
  sigma_z_hat = rep_matrix(0,rank_B+m*p_eta, rank_B+m*p_eta);
  sigma_z_hat[1:rank_B,1:rank_B] = BTWyBinv/lambda_y;
  sigma_z_hat[rank_B+1:,rank_B+1:] = KTKinv/lambda_eta;
  sigma_z_hat = sigma_z_hat + sigma_z;
  
  // Perform Cholesky decomposition of sigma_z_hat for efficiency
  // either force symmetry somehow using cholesky_decompose function,
  L_z_hat = cholesky_decompose(sigma_z_hat);
  
  // Specify prior distribution of z_hat
  z_hat ~ multi_normal_cholesky(mu_z_hat, L_z_hat);
  
  // Specify priors here
  // Leave calibration parameters as uniform for now, though input these later using existing code
  rho_w[1:p_eta*q] ~ beta(1.0, 0.1);    // Correlation parameter for emulator
  // Higdon and Chong priors on lambda_w look fairly similar
  lambda_w[1:p_eta] ~ gamma(5.0, 5.0);      // Precision parameter for emulator, gamma (shape, rate) (Higdon)
  //lambda_w[1:p_eta] ~ gamma(10.0, 10.0);      // Precision parameter for emulator, gamma (shape, rate) (Chong)
  //lambda_v[1:p_delta] ~ gamma(1.0, 0.0001); // Precision parameter for discrepancy, gamma (shape, rate) (Higdon)
  lambda_v[1:p_delta] ~ gamma(10.0, 0.0001); // Precision parameter for discrepancy, gamma (shape, rate) (Chong)
  lambda_y ~ gamma(a_y_dash, b_y_dash);               // Precision parameter for observation error, gamma (shape, rate)
  lambda_eta ~ gamma(a_eta_dash, b_eta_dash);        // Precision parameter for PCA truncation error, gamma (shape, rate}.
  // define priors on calibration parameters
  for (i in 1:(q-1)){
    tf_gauss[i] ~ normal(tf_param_1[i], tf_param_2[i]);
  }
  // Remaining calibration parameters are uniformly distributed by default, between previously defined bounds
}
