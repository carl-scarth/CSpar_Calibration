// concatenate all calibration parameters into a single row vector 
row_vector[2] tf_unif_trans;
real tf_halfnorm_trans;

tf_unif_trans = rep_row_vector(0.0, 2);
tf_unif_trans[1] = tf_param_1[4] + (tf_param_2[4]-tf_param_1[4])*tf_unif[1];
tf_unif_trans[2] = tf_param_1[6] + (tf_param_2[6]-tf_param_1[6])*tf_unif[2];
tf_halfnorm_trans = tf_param_1[5] + tf_halfnorm;

tf = append_col(tf_gauss,tf_unif_trans[1]);
tf = append_col(tf,tf_halfnorm_trans);
tf = append_col(tf,tf_unif_trans[2]);
