# define a vector of calibration parameters and their bounds
row_vector<lower=-0.1,upper=1.1>[q] tf_gauss;

// row_vector<lower=0,upper=1>[q-1] tf_gauss; // Gaussian-distributed priors
// real<lower=tf_param_1[q],upper=tf_param_2[q]> tf_unif; // uniform-distributed prior

