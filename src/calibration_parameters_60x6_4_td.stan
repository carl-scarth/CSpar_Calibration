// define a vector of calibration parameters and their bounds
row_vector<lower=-0.25,upper=1.25>[3] tf_gauss; // Gaussian-distributed priors
row_vector<lower=0,upper=1>[3] tf_unif; // uniform-distributed prior
