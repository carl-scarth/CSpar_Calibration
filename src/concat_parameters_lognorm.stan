// concatenate all calibration parameters into a single row vector 
row_vector[2] tf_unif_trans;
//for (i in 1:2){
//  tf_halfn_trans[i] = tf_param_1[4+i] + tf_halfnorm[i];
//}
tf_unif_trans = tf_param_1[5:6] + (tf_param_2[5:6]-tf_param_1[5:6]) .* tf_unif;
tf = append_col(tf_gauss,tf_unif_trans);
