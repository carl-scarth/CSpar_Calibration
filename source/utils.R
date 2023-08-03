# Miscellaneous useful functions

estimate_mode <- function(x) {
  # estimates mode of distribution of samples in x
  d <- density(x)
  d$x[which.max(d$y)]
}

full_field_emulator_modes <- function(rho_w, lambda_w, lambda_eta) {
  # Extract number of inputs, and number of basis functions
  p_eta = ncol(lambda_w)
  q = ncol(rho_w)/p_eta
  # estimate modes of the posterior distribution
  modes = rep(0,(p_eta*(q+1))+1)
  for (i in 1:(p_eta*q)){
    modes[i] = estimate_mode(rho_w[,i])
  }
  for (i in 1:p_eta){
    modes[p_eta*q + i] = estimate_mode(lambda_w[,i])
  }
  modes[(p_eta*(q+1))+1] = estimate_mode(lambda_eta)
  # Create dataframe and write column names
  modes = as.data.frame(t(modes))
  for (i in 1:p_eta){
    colnames(modes)[((i-1)*q+1):(q*i)] <- sprintf("rho_w_%d_%d", i,1:q)
  }
  colnames(modes)[(p_eta*q+1):(p_eta*(q+1))] <- sprintf("lambda_w_%d", 1:p_eta)
  colnames(modes)[p_eta*(q+1)+1] <- "lambda_eta"
  return(modes)
}

write_output <- function(data, out_string, label_string){
  # Formats data and writes to a csv file with name, and column headings given 
  # in out_string, and label label_string
  data = as.data.frame(data)
  for (i in 1:ncol(data)){
    colnames(data)[i] = sprintf("%s_%d",out_string, i)
  }
  write.csv(data, sprintf("outputs/%s_%s.csv", out_string, label_string), row.names = FALSE)
}

write_output_samples <- function(data, out_string, label_string){
  # Formats data and writes multiple posterior samples of predictions to a csv
  # file with name, and column headings given in out_string, and label 
  # label_string. Assumes data is stored in a 3D array with dimensions 
  # [n_output, n_sam, n_pred]
  # Determine the number of predictions and posterior samples
  n_sam = dim(data)[2]
  n_pred = dim(data)[3]

  # Flatten the array into columns with format [sam1_pred1, sam2_pred1, ... sam_n_pred1, sam_1_pred_2 etc...]
  data = matrix(aperm(data, c(2, 3, 1)), nrow = nrow(data), byrow = TRUE)
  data = as.data.frame(data)
  # Name columns
  for (i in 1:n_pred){
    for (j in 1:n_sam){
      colnames(data)[(i-1)*n_sam+j] = sprintf("%s_%d_sam_%d",out_string,i,j)
    }
  }
  write.csv(data, sprintf("outputs/%s_%s.csv", out_string, label_string), row.names = FALSE)
}