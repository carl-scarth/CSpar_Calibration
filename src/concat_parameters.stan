// concatenate all calibration parameters into a single row vector 
real tf_unif_trans;
row_vector[2] tf_halfnorm_trans;

//tf_unif_trans = rep_row_vector(0.0, 2);
tf_unif_trans = tf_param_1[6] + (tf_param_2[6]-tf_param_1[6])*tf_unif;
tf_halfnorm_trans[1] = tf_param_1[7] + abs(tf_halfnorm[1]);
tf_halfnorm_trans[2] = tf_param_1[8] + abs(tf_halfnorm[2]);

//tf_halfnorm_trans = tf_param_1[5:6] + tf_halfnorm;

//tf = append_col(tf_gauss,tf_unif_trans);
//tf = append_col(tf,tf_halfnorm_trans);
tf[1:5] = tf_gauss;
tf[6] = tf_unif_trans;
tf[7:8] = tf_halfnorm_trans;

