# Code for scaling input and output samples before and after fitting the model

# Function for normalising inputs onto the unit hypercube, using (scalar or 
# vector) minimum and maximum values x_min and x_max
# std indicates whether the x is a standard deviation, in which case it isn't
# necessary to subtract the minimum
normalise_inputs <- function(x, x_min, x_max, std = FALSE){
  # Put maximum and minimum values into correct format
  x_min = t(as.matrix(x_min))
  x_max = t(as.matrix(x_max))
  if (is.matrix(x)){
    N = nrow(x)
  } else {
    N = 1
  }
  if (!std){
    x = x - x_min[rep(1,N),]
  }
  x_norm = x/(x_max[rep(1,N),]-x_min[rep(1,N),])
  return(x_norm)
}

# Convert inputs which have been normalised to lie on the unit hypercube back
# onto their original scale using (scalar or vector) minimum and maximum values 
# std indicates whether the x is a standard deviation, in which case it isn't
# necessary to subtract the minimum
rescale_inputs <- function(x_norm, x_min, x_max, std = FALSE){
  # Put maximum and minimum values into the correct format
  x_min = t(as.matrix(x_min))
  x_max = t(as.matrix(x_max))
  if (is.matrix(x_norm)){
    N = nrow(x_norm)
  } else {
    N = 1
  }
  x = x_norm*(x_max[rep(1,N),]-x_min[rep(1,N),])
  if (!std){
    x = x + x_min[rep(1,N),]
  }
  return(x)
}

# Function for standardising vector output by the mean vector mu_y and scalar 
# standard deviation sigma_y to have unit mean and zero standard deviation
# std indicates whether y is a standard deviation, in which case it isn't
# necessary to centre before scaling
# There is an option to provide the mean vector and standard deviation, 
# otherwise these are calculated internally
standardise_vector_output <- function(y, mu_y = NULL, sigma_y = NULL, std = FALSE, q_y = 1, scale_by_component = F){
  if (!std) {
    if (is.null(mu_y)){
      mu_y = rowMeans(y)
    }
    y = sweep(y,1,mu_y,"-")
  }
  if (is.null(sigma_y)){
    # Note, it's better to calculate this after centring the data, as this 
    # affects the standard deviation (given the mean vector is used, rather than
    # the overall mean)
    if (q_y > 1 & scale_by_component) {
      sigma_y = c()
      n_y = nrow(y)/q_y
      for (i in 1:q_y) {
        y_i = y[((i-1)*n_y+1):(i*n_y),]
        sigma_y = c(sigma_y, rep(sd(y_i),n_y))
      }
    } else {
      sigma_y = sd(y)
    }
  }
  if (length(sigma_y) == 1) {
    y_scale = y/sigma_y
  } else {
    y_scale = sweep(y,1,sigma_y,"/") # Allows for input of vector y to scale the individual componenents individually
  }
  return(list(y_scale, mu_y, sigma_y))
}


# Function for the inverse standardisation of vector output by the mean vector
# mu_y, and scalar standard deviation sigma_y
# std indicates whether the y is a standard deviation, in which case it isn't
# necessary to add the mean
# I think this would actually also work for scalar valued mu_y etc, and vector 
# sigma_y
rescale_vector_output <- function(y_scale, mu_y, sigma_y, std = FALSE){
  if (length(sigma_y) == 1) {
    y = y_scale*sigma_y
  } else {
    y = sweep(y_scale,1,sigma_y,"*") # Allows for input of vector y to scale the individual componenents individually
  }
  if (!std) {
    #if (is.matrix)(y){
    #  y = y + t(replicate(ncol(y), mu_y))
      #y = y + matrix(replicate(ncol(y),mu_y),nrow=nrow(y))
    #} else {
    y = y + mu_y
    #}
  }
  return(y)
}