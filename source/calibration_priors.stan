// Define specified prior distributions on the calibration parameters
for (i in 1:3){
  tf_gauss[i] ~ normal(tf_param_1[i], tf_param_2[i]);
}
