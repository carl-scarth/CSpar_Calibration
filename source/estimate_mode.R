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