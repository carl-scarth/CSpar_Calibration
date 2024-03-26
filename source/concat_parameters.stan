// concatenate all calibration parameters into a single row vector 
row_vector[q] tf_unif_trans;
tf_unif_trans = tf_param_1[4:6] + (tf_param_2[4:6]-tf_param_1[4:6]) .* tf_unif;
tf = append_col(tf_gauss,tf_unif_trans);
