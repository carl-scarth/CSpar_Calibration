import numpy as np
# Header for standardising inputs and outputs of Gaussian processes

# Function for normalising inputs onto the unit hypercube, using (scalar or 
# vector) minimum and maximum values x_min and x_max
# std indicates whether the x is a standard deviation, in which case it isn't
# necessary to subtract the minimum
def normalise_inputs(x, std = False, **kwargs):
  # Calcuate maximum and minimum values of data if required
  # Assumes data x is a N x d numpy array, with axis 0 being a 
  # sample of a d-dimensional input vector
  maxmin_out = False
  if "x_min" not in kwargs.keys():
    maxmin_out = True
    x_min = x.min(axis = 0)
    x_max = x.max(axis = 0)
  else:
    x_min = kwargs["x_min"]
    try:
        x_max = kwargs["x_max"]
    except:
        raise Exception("Please also specify x_max")

  if not std:
    x = x - x_min
  x_norm = x/(x_max-x_min)

  if maxmin_out:
    return((x_norm, x_min, x_max))
  else:
    return(x_norm)
  
# Function for standardising output using mean mu_y and standard deviation
# sigma_y to have unit mean and zero standard deviation.
# std indicates whether y is a standard deviation, in which case it isn't
# necessary to centre before scaling
# There is an option to provide the mean vector and standard deviation, 
# otherwise these are calculated internally
def standardise_output(y, mu_y = [], sigma_y = [], std = False):
  if not(mu_y and sigma_y):
    musd_out = True
  else:
    musd_out = False

  if not std:
    if not mu_y:
      mu_y = np.mean(y)
    
    y = y - mu_y

  if not sigma_y:
    sigma_y = np.std(y)

  y_scale = y/sigma_y
  if musd_out:
    return((y_scale, mu_y, sigma_y))
  else:
    return(y_scale)
  

# Function for standardising vector output by the mean vector mu_y and scalar 
# standard deviation sigma_y to have unit mean and zero standard deviation
# std indicates whether y is a standard deviation, in which case it isn't
# necessary to centre before scaling
# The the mean vector and standard deviation are calculated internally if 
# not specified
def standardise_vector_output(y, mu_y = [], sigma_y = [], std = False, q_y = 1):
  if not std:
    if not mu_y:
      mu_y = np.mean(y, axis= 1, keepdims=True)
    y = y - mu_y
  
  if not sigma_y:
    # Note, it's better to calculate this after centring the data, as this 
    # affects the standard deviation (given the mean vector is used, rather than
    # the overall mean)
    if q_y > 1:
      sigma_y = []
      n_y = y.shape[0]//q_y
      for i in range(q_y):
        sd_i = np.std(y[i*n_y:(i+1)*n_y])
        sigma_y.extend([sd_i for j in range(n_y)])
      sigma_y = np.array(sigma_y).reshape((-1,1))
    else:
      sigma_y = np.std(y)
  
  y_scale = y/sigma_y # Should work in for both scalar and vector sigma_y
  # Shape is an empty tuple for a scalar so this won't work if it's needed
  #if (sigma_y.shape[0] == 1):
  #  y_scale = y/sigma_y
  #} else {
  #  y_scale = sweep(y,1,sigma_y,"/") # Allows for input of vector y to scale the individual componenents individually
  #}
  return y_scale, mu_y, sigma_y


# Function for the inverse standardisation of vector by mean mu_y and 
# standard deviation sigma_y
# std indicates whether the y is a standard deviation, in which case it isn't
# necessary to add the mean
def rescale_output(y_scale, mu_y = 0.0, sigma_y = 1.0, std = False):
  y = y_scale*sigma_y
  if not std:
    y = y + mu_y

  return(y)


# Convert inputs which have been normalised to lie on the unit hypercube back
# onto their original scale using (scalar or vector) minimum and maximum values 
# std indicates whether the x is a standard deviation, in which case it isn't
# necessary to subtract the minimum
def rescale_inputs(x_norm, x_min, x_max, std = False):
  # Put maximum and minimum values into the correct format
  x = x_norm*(x_max-x_min)
  if not std:
    x = x + x_min
  
  return(x)