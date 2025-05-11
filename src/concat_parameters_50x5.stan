// concatenate all calibration parameters into a single row vector 
// row_vector[q] tf; // Vector containing both uniform and Gaussian distributed priors

//tf = append_col(tf_gauss, tf_unif);
// tf = tf_gauss;
tf = append_col(tf_gauss[1:3],tf_unif);
tf = append_col(tf,tf_gauss[4]);
