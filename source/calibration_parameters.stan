// define a vector of calibration parameters and their bounds
row_vector<lower=-0.2,upper=1.2>[5] tf_gauss; // Gaussian-distributed priors
real<lower=0.0,upper=1.0> tf_unif; // uniform-distributed prior
row_vector<lower=-1.05, upper=1.05>[2] tf_halfnorm; // half-normal-distributed prior
