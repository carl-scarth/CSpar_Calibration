# Functions for plots of overlaid prior and posterior distributions

calibration_inp_hist <- function(tf, tf_param = NULL, inp_labels = NULL, new_window = FALSE){
  # Produce prior and posterior plots for samples of calibrated model inputs
  # tf = matrix of posterior samples
  # tf_param = data.frame where each row is a random variable, which has columns 
  # "distribution" indicating the prior distribution, and "param_1" and
  # "param_2" with parameters of the distribution
  
  if (new_window) {
    dev.new(noRStudioGD = TRUE) # Create a new window if needed
  }
  q = ncol(tf)
  
  # If labels aren't provided, get these from the parameter list if provided
  if (is.null(inp_labels)){
    if (!is.null(tf_param)){
      inp_labels = rownames(tf_param)
    }
  }

  # Produce a subplot for each input
  par(mfrow = c(1,q))
  for (i in 1:q){
    if (!is.null(inp_labels)){
      label = inp_labels[i]
    } else {
      label = paste("tf_",as.character(i),sep="")
    }
    # get limits for plot
    if (!is.null(tf_param)){
      if (tf_param$distribution[i] == "Gaussian" | tf_param$distribution[i] == "Lognormal"){
        # maximum and minimum of +/- 3 standard deviations
        mu = tf_param$param_1[i]
        std = tf_param$param_2[i]
        t_min = mu - 3*std
        t_max = mu + 3*std
      } else if (tf_param$distribution[i] == "Halfnormal") {
        mu = tf_param$param_1[i]
        std = tf_param$param_2[i]
        t_min = mu
        t_max = mu + 3*std
      } else if ((tf_param$distribution[i] == "Uniform") | (tf_param$distribution[i] == "Loguniform")){
        t_min = tf_param$param_1[i]
        t_max = tf_param$param_2[i]
      } else {
        stop("Error: Non-implemented distribution")
      }
    } else {
      t_min = min(tf[,i])
      t_max = max(tf[,i])
    }
    hist(tf[,i], 
         main = label,
         xlab = label,
         col = "brown1",
         breaks = 25,
         freq = FALSE,
         #xlim = c(t_min,t_max),
         cex.lab = 1.25,
         cex.axis = 1.25)
  
    # overlay plot of prior distribution if parameters are provided
    if (!is.null(tf_param)){
      # identify the appropriate distribution and properties
      t_plot = seq(t_min,t_max, length.out = 100)
      if (tf_param$distribution[i] == "Gaussian" | tf_param$distribution[i] == "Lognormal" ){
        # maximum and minimum of +/- 3 standard deviations
        prior_plot = dnorm(t_plot, mean = mu, sd = std)
      } else if (tf_param$distribution[i] == "Halfnormal") {
        prior_plot = 2*dnorm(t_plot, mean = mu, sd = std)
      } else if ((tf_param$distribution[i] == "Uniform") | (tf_param$distribution[i] == "Loguniform")){
        prior_plot = dunif(t_plot, min = t_min, max = t_max)
      }
      lines(t_plot,prior_plot,lwd=3,"col"="blue")
    }
  }
}

full_field_rho_hist <- function(rho_w, p_eta, inp_labels=NULL, prior_shape_1 = 1.0, prior_shape_2 = 0.1, new_window = FALSE){
  # Produce plots of posterior and prior distributions of correlation parameters
  # (rho) for full-field emulator. Arrange with each row corresponding a given
  # principal component, and each column to a given input.
  # inp_labels is an optional input with vector of input labels for plot titles
  # Option to pass shape parameters of the beta prior, with default values as 
  # used by Higdon et al.
  # new_window is a Boolean which dictates if the plot is produced in a new window
  
  # Create a new plot window if needed
  if (new_window){
    dev.new(noRStudioGD = TRUE)
  }
  q = ncol(rho_w)/p_eta # Identify number of inputs
  rho_plot = seq(0,1, length.out = 100) # Set of x values for line-plot of prior
  # Set-up subplots. Reduce margin size to prevent error from too many plots in 
  # RStudio plot window. There isn't an ideal way of automating this.
  par(mfrow = c(p_eta,q), mar = c(1,2,1,1), mgp = c(2,0.5,0))
  for (i in 1:p_eta){
    for (j in 1:q){
      if (is.null(inp_labels)){
        title_j = paste("rho_w_",as.character(i),"_",as.character(j)) 
      } else {
        # Only include title on first row
        if (i == 1){
          title_j = inp_labels[j]
        } else {
          title_j = NULL
        }
      }
      # Plot histogram of posterior distribution
      # Don't plot x axis as 
      hist(samples$rho_w[,(i-1)*q+j],
           main = title_j,
           xlab = paste("rho_w_",as.character(i),"_",as.character(j)),
           col = "firebrick1",
           breaks = 25,
           freq = FALSE,
           xlim = c(0,1), 
           cex.lab=1.5,
           cex.axis=1.5)
      # Overlay plot of prior distribution
      prior_plot = dbeta(rho_plot, shape1=prior_shape_1, shape2=prior_shape_2)
      lines(rho_plot, prior_plot, lwd=3, col="blue")
    }
  }
  # Reset plot window to default
  par(mar = c(5.1, 4.1, 4.1, 2.1), mgp = c(3,1,0))
}

full_field_lambda_hist <- function(lambda, p, prior_shape = 5.0, prior_rate = 5.0, new_window = FALSE){
  # Produce plots of posterior and prior distributions of GP precision
  # hyperparameters (lambda) for full-field output. Each plot corresponds to a
  # different basis function
  # Optional inputs prior_shape and prior_rate give the properties of the gamma
  # prior. Default values are the Higdon et al. values for emulator precision
  # new_window is a Boolean with dictates if the plot is produced in a new window
  
  # create a new window if needed
  if (new_window){
    dev.new(noRStudioGD = TRUE)
  }
  
  # If p is greater than 5, place subplots across multiple rows
  if (p <= 5){
    par(mfrow = c(1,p))
  } else {
    par(mfrow = c(p%/%5+1,5))
  }
  lambda_plot = seq(0, qgamma(0.9999, prior_shape, prior_rate), length.out = 100) # x values for prior line plot
  for (i in 1:p) {
    # plot posterior
    hist(lambda[,i],
      main = paste("lambda_w",as.character(i)),
      xlab = paste("lambda_w",as.character(i)),
      col = "firebrick1",
      breaks = 25,
      freq = FALSE,
      xlim = c(0, qgamma(0.9999, prior_shape, prior_rate)),
      cex.axis=1.5,
      cex.lab=1.5)
    # overlay plot of prior
    prior_plot = dgamma(lambda_plot,shape=prior_shape,rate=prior_rate)
  lines(lambda_plot,prior_plot,lwd=3,col="blue")
  }
  # reset plot window to default
  par(mfrow = c(1,1))
}

lambda_hist <- function(lambda, prior_shape = 5.0, prior_rate = 5.0, label = "lambda", adj_prior_shape = NULL, adj_prior_rate = NULL, new_window = FALSE){
  # Plots precision hyper-parameter lambda for scalar-valued quantities
  # Optional inputs prior_shape and prior_rate give the properties of the gamma
  # prior. Default values are the Higdon et al. values for emulator precision
  # If adj_prior_shape and adj_prior_rate are passed a second prior plot is 
  # overlaid with these parameter values.
  # new_window is a Boolean with dictates if the plot is produced in a new window
  
  # create a new window if needed
  if (new_window){
    dev.new(noRStudioGD = TRUE)
  }
  if (!is.matrix(lambda)){
    lambda = as.matrix(lambda)
  }
  n_lambda = ncol(lambda)
  par(mfrow = c(1,n_lambda))
  # Get maximum and minimum x-limit for plots
  max_x_prior <- qgamma(0.99, prior_shape, prior_rate)
  min_x_prior <- qgamma(0.01, prior_shape, prior_rate)
  max_x <- max_x_prior
  if (!is.null(adj_prior_shape)){
    max_x_adj_prior <- qgamma(0.99, adj_prior_shape, adj_prior_rate)
    min_x_adj_prior <- qgamma(0.01, adj_prior_shape, adj_prior_rate)
    if (max_x < max_x_adj_prior){
      max_x = max_x_adj_prior
    }
  }
  for (i in 1:n_lambda) {
    if (n_lambda > 1){
      header = paste(label, i, sep="_")
    } else {
      header = label
    }
    #plot posterior
    hist(lambda[,i],
      main = header,
      xlab = label,
      col = "firebrick1",
      breaks = 25,
      freq = FALSE,
      xlim = c(0,2^ceiling(log(max_x,2))), # round up upper limit
      cex.axis=1.5,
      cex.lab=1,5)
    lambda_plot <- seq(min_x_prior, max_x_prior, length.out = 1000) # set of x values for prior plot
    prior_plot = dgamma(lambda_plot, shape=prior_shape, rate=prior_rate)
    lines(lambda_plot,prior_plot,lwd=3,col="blue")
    # Produce additional prior plot if needed
    if (!is.null(adj_prior_shape)){
      adj_lambda_plot <- seq(min_x_adj_prior, max_x_adj_prior, length.out = 1000) # set of x values for prior plot
      adj_prior_plot = dgamma(adj_lambda_plot,shape=adj_prior_shape,rate=adj_prior_rate)
      lines(adj_lambda_plot,adj_prior_plot,lwd=3,col="green")
    }
  }
}

# If modelling discrepancy consider re-using the above code? It could be that the
# only changes would be axis labels

prior_posterior_pairs <- function(tf, tf_param, inp_labels = NULL, new_window = FALSE){
  # Produces scatter plots populated with samples from the prior and posterior 
  # distribution of the calibration parameters, where posterior samples are 
  # passed via matrix tf
  # tf_param = data.frame where each row is a random variable, which has columns 
  # "distribution" indicating the prior distribution, and "param_1" and
  # "param_2" with parameters of the distribution
  # inp_labels = optional labels for each random variable
  
  if (new_window) {
    dev.new(noRStudioGD = TRUE) # Create a new window if needed
  }
  q = ncol(tf)
  
  # If labels aren't provided, get these from the parameter list if provided
  if (is.null(inp_labels)){
    inp_labels = rownames(tf_param)
  }
  
  prior_sam = matrix(0, nrow = N_samples, ncol = q)
  for (i in 1:q){
    if (tf_param$distribution[i] == "Gaussian"){
      prior_sam[,i] = rnorm(N_samples, tf_param$param_1[i], tf_param$param_2[i])
    } else if ((tf_param$distribution[i] == "Uniform") | (tf_param$distribution[i] == "Loguniform")){
      prior_sam[,i] = runif(N_samples, tf_param$param_1[i], tf_param$param_2[i])
    } else {
      stop("Error: Non-implemented distribution")
    }
  }
  # Combine prior and posterior samples into a single data.frame
  prior_post_frame = as.data.frame(rbind(prior_sam, tf))
  colnames(prior_post_frame) = inp_labels
  # Create list of colours for the plot
  colours = c(rep("blue",N_samples),rep("red",N_samples))
  pairs(prior_post_frame,col=colours)
}


w_star_hist <- function(w_star, w_lim = NULL, new_window = FALSE){
  # Produces prior and posterior plots of the calibrated expansion coefficients 
  # in the reduced-dimensional representation of the model output
  # w_star = p x n_sam matrix of posterior samples of expansion coefficients,
  # where p is the number of coefficients and n_sam is the number of samples
  # w_lim = vector containing limits of the prior plot, if provided
  # Boolean indicating whether to plot in a new window
  if (new_window) {
    dev.new(noRStudioGD = TRUE) # Create a new window if needed
  }
  if (is.null(w_lim)){
    get_limits = TRUE
  } else {
    get_limits = FALSE
  }
  p = nrow(w_star)
  # If p is greater than 5, place subplots across multiple rows
  if (p <= 5){
    par(mfrow = c(1,p))
  } else {
    par(mfrow = c(p%/%5+1,5))
  }
  for (i in 1:p){
    if (get_limits){
      w_lim = c(min(w_star[i,]), max(w_star[i,]))
    }
    hist(out_list$w_star[i,], 
      main =  paste("Feature",as.character(i)),
      xlab = paste("w_star", as.character(i)),
      col = "brown1",
      breaks = 15,
      freq = FALSE,
      xlim = w_lim,
      cex.lab = 1.25,
      cex.axis = 1.25)
    # overlay plot of prior distribution
    w_plot = seq(w_lim[1], w_lim[2], length.out = 100)
    prior_plot = dnorm(w_plot,mean = 0, sd = 1)
    lines(w_plot,prior_plot,lwd=3,"col"="blue")
  }
}