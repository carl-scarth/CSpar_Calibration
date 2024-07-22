// define a vector of calibration parameters and their bounds
row_vector<lower=-0.3,upper=1.3>[7] tf_gauss; // Gaussian-distributed priors
real<lower=0.0,upper=1.0> tf_unif; // uniform-distributed prior
row_vector<lower=0.0>[2] tf_halfnorm; // half-normal-distributed prior
