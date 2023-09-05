# Code for scaling input and output samples before and after fitting the model

# Function for normalising inputs onto the unit hypercube
normalise_inputs <- function(x, x_min, x_max, sd = FALSE){
  # Put maximum and minimum values into correct format
  x_min = t(as.matrix(x_min))
  x_max = t(as.matrix(x_max))
  if (is.matrix(x)){
    N = nrow(x)
  } else {
    N = 1
  }
  if (!sd){
    x = x - x_min[rep(1,N),]
  }
  x = x/(x_max[rep(1,N),]-x_min[rep(1,N),])
  return(x)
}

# Function for the inverse standardisation of vector output by the mean vector
# mu_y, and scalar standard deviation sigma_y
# sd indicates whether the y is a standard deviation, in which case it isn't
# necessary to add the mean
# I think this would actually also work for scalar valued mu_y etc, and vector 
# sigma_y
rescale_vector_output <- function(y_scale, mu_y, sigma_y, sd = FALSE){
  y = y_scale*sigma_y
  if (!sd) {
    y = y + mu_y
  }
  return(y)
}