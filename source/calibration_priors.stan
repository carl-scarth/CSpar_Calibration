// Define specified prior distributions on the calibration parameters
for (i in 1:q){
  tf_gauss[i] ~ normal(tf_param_1[i], tf_param_2[i]);
}
// Remaining calibration parameters are uniformly distributed by default, 
// between previously defined bounds

