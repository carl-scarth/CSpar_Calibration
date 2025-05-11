// concatenate all calibration parameters into a single row vector 
real tf_unif_trans;
real tf_halfnorm_trans;

//tf_unif_trans = rep_row_vector(0.0, 2);
tf_unif_trans = tf_param_1[4] + (tf_param_2[4]-tf_param_1[4])*tf_unif;
tf_halfnorm_trans = tf_param_1[5] + tf_halfnorm;

//tf_halfnorm_trans = tf_param_1[5:6] + tf_halfnorm;

tf = append_col(tf_gauss,tf_unif_trans);
tf = append_col(tf,tf_halfnorm_trans);
