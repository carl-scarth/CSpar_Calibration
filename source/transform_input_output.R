# Code for scaling input and output samples before and after fitting the model

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