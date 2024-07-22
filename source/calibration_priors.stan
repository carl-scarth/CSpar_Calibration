// Define specified prior distributions on the calibration parameters
for (i in 1:7){
  tf_gauss[i] ~ normal(tf_param_1[i], tf_param_2[i]);
}
for (i in 1:2){
  tf_halfnorm[i] ~ normal(0.0, tf_param_2[8+i]);
}

