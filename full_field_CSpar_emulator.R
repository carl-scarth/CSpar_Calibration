
# Fit an emulator to the full-field output from a C-Spar Finite Element model
# Input data are the Longitudinal Modulus E11, ply thickness, and the torsional
# stiffness of a spring representing an uncertain boundary condition
# Output data are the axial displacements of the nodes of the FE model at a fixed load.
# Follows emulator aspect of D. Higdon et al, "Computer Model Calibration Using 
# High-Dimensional Output",Journal of the American Statistical Association,2008.
# This code handles pre and post processing.

library(data.table)
library(rstan)
library(matrixStats)
library(MASS)

# Set current working directory. This should be modified to match the directory
# of the user
setwd("C:/Users/cs2361/Documents/CSpar_Calibration/")

# include functions which are called in this code
source("source/estimate_mode.R")
source("source/covariance_matrices.R")
source("source/dimension_reduction.R")
source("source/prior_posterior_plots.R")

#-------------------------------------------------------------------------------

# Set up parameters which govern the formulation

p_eta = 7 # Number of basis functions to be retained for the emulator from SVD
exp_tol = 1e-6 # Tolerance variance fraction used to assess SVD convergence
disp_str = "w" # String which identifies the displacement component of interest (u,v, or w)
# Define parameters of the gamma prior on the error associated with truncating
# the series expansion for the model output
a_eta = 1.0     # Shape parameter for the lambda_eta prior
b_eta = 0.0001  # Rate parameter for the lambda_eta prior 
iter = 4000 # Number of samples per chain
chains = 3 # Number of chains for simulation
print_svd_output = TRUE # Print diagnostic output of svd to the terminal?
export_modes = TRUE # Calculate modes of emulator hyperparameters and write to file?

#-------------------------------------------------------------------------------

# Set up simulation data

# Load in emulator training data input values from Design of Experiments. 
# in_file = "LHSDesign50x3" # File identifier string which is in both input and output csvs
in_file = "LHSDesign40x4" # File identifier string which is in both input and output csvs
XT_sim = fread(paste("inputs/",in_file,".csv", sep = ""))

# In this example I fit the emulator to the log of spring stiffness K, which is 
# a more natural choice of values across which outputs are expected for
# variations in this input
#XT_sim$K = log(XT_sim$K)
#colnames(XT_sim)[3] = "log_K"

# Determine useful quantities from model inputs and outputs. Variable names 
# match the notation of Higdon et al. 2008
q = ncol(XT_sim)          # number of calibration inputs
tc = as.matrix(XT_sim)    # Convert training data input points to a matrix for passing to stan
m = nrow(XT_sim)          # sample size of computer simulation data

# Load in  training data output displacement values from Abaqus. Here I've used
# a similar naming convention to the inputs to automate changes. The file name 
# can be changed manually if need be.
# Each row of XT_sim corresponds to a block of three columns of displacement 
# data, with a column for each component u,v,w 
# abaqus_displacements = fread(paste("inputs/",in_file,"_displacements.csv", sep=""))
abaqus_displacements = fread(paste("inputs/",in_file,"_fixed_100kN.csv", sep=""))
n_eta = nrow(abaqus_displacements) # number of output points per simulation

# Extract the displacement for the component of interest and store in a matrix
dt_simulation = matrix(NA,nrow = n_eta, ncol = m)
for (i in 1:m){
  dt_simulation[,i] = abaqus_displacements[[as.name(paste(disp_str,sprintf("_%d", i),sep=""))]]
}

#-------------------------------------------------------------------------------

# Load in test points at which predictions are required
# XT_pred = fread("inputs/LHSDesign50x3_1.csv")
XT_pred = fread("inputs/LHSDesign40x4_1.csv")
# Repeat the log transformation for the test data
# XT_pred$K = log(XT_pred$K)
# colnames(XT_pred)[3] = "log_K"
t_pred = as.matrix(XT_pred)
n_pred = nrow(t_pred) # number of predictions
  
#-------------------------------------------------------------------------------

# All training data points are normalised onto the unit hypercube [0,1]^q before
# being passed to stan. The same transformation is applied to the test points

# Determine the maximum and minimum value of each input within the training data
t_min = t(as.matrix(colMins(tc)))
t_max = t(as.matrix(colMaxs(tc)))
# Normalise the inputs using these maximum and minimum values
tc = (tc - t_min[rep(1,m),])/(t_max[rep(1,m),]-t_min[rep(1,m),])
t_pred = (t_pred - t_min[rep(1,n_pred),])/(t_max[rep(1,n_pred),]-t_min[rep(1,n_pred),])

# plot Design of Experiments and test points
pairs(tc,col="blue", pch=4, 
      main = "Normalised Training Data Input values",
      cex=1.5,
      cex.labels = 1.75,
      cex.axis=1.5,
      cex.lab=1.5)
pairs(t_pred,col = "blue", pch=4, 
      main = "Normalised Test Data Input Values",
      cex=1.5,
      cex.labels = 1.75,
      cex.axis=1.5,
      cex.lab=1.5)

#-------------------------------------------------------------------------------

# Performs SVD on the model outputs, and standardise so the coefficients of the 
# expansion have zero mean and unit variance

# Centre the simulation output for each node This guarantees that the
# coefficients have zero (sample) mean. 
mu_dt = rowMeans(dt_simulation)
# Write model mean to csv file
write.csv(mu_dt, paste("outputs/model_mean_",in_file,".csv"), row.names = FALSE) 
dt_all_cen = sweep(dt_simulation,1,mu_dt,"-")

# Divide by the standard deviation of the outputs
sd_dt = sd(as.matrix(dt_all_cen))
dt_all_cen = dt_all_cen/sd_dt
# Convert to matrix for stan
eta = as.matrix(dt_all_cen)
out_basis = svd_basis(eta, p_eta = p_eta, exp_tol = exp_tol, 
                      print_output = print_svd_output, csv_label = in_file)
K_eta = out_basis[[1]]
p_eta = out_basis[[2]] # Used to determine p_eta automatically if not provided as an argument, otherwise this is unchanged

#-------------------------------------------------------------------------------

# Reduce the dimension of the output data and calculate associated quantities 
# for input to Stan
processed_data = reduce_dimension(eta, K_eta, a_eta, b_eta)
# Adjusted prior parameters for the basis expansion truncation error
a_eta_dash = processed_data[[1]]
b_eta_dash = processed_data[[2]]
z_hat = processed_data[[3]] # Reduced-dimensional outputs
KTKinv = processed_data[[4]] # Inverse of the product of basis matrices

#-------------------------------------------------------------------------------

# Set up the environment for Stan, pass arguments, and run stan model
# Settings taken from: https://betanalpha.github.io/assets/case_studies/gaussian_processes.html#21_Simulating_From_A_Gaussian_Process
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
util = new.env()

# List of arguments to pass to stan
stan_data = list(m=m, q=q, n_eta=n_eta, p_eta=p_eta, a_eta_dash = a_eta_dash, 
                 b_eta_dash = b_eta_dash, z_hat = z_hat, tc = tc, KTKinv = KTKinv)
# Run stan
fit = stan(file = "source/full_field_emulator.stan",
           data = stan_data,
           iter = iter,
           chains = chains,
           model_name = "full_field_emulator")

#-------------------------------------------------------------------------------

# Post-process and plot the simulation data

# Produce trace plots
stan_trace(fit, pars = c("rho_w", "lambda_w","lambda_eta"))
# Print summary of results
print(fit, pars = c("rho_w", "lambda_w", "lambda_eta"))

# extract samples from stan output
samples <- extract(fit)
rho_w <- samples$rho_w           # Correlation lengths
lambda_w <- samples$lambda_w     # Emulator precisions
lambda_eta <- samples$lambda_eta # Expansion truncation error 
N_samples <- dim(rho_w)[1]       # Total number of samples post warm-up

# Extract label of inputs for plots
labels = colnames(XT_sim) # don't think I need this but keeping just in case

# If required, estimate the modes of emulator hyperparameters and write to a csv
if (export_modes){
  modes = full_field_emulator_modes(rho_w, lambda_w, lambda_eta)
  write.csv(modes, paste("outputs/emulator_modes_",in_file,".csv", sep=""), row.names = FALSE)
}

# Plot correlation parameter histograms
full_field_rho_hist(rho_w, p_eta, inp_labels = labels)
# Plot emulator precision parameters for emulator
full_field_lambda_hist(lambda_w, p_eta)
# Plot emulator trunction error precision
lambda_hist(lambda_eta, prior_shape = a_eta, prior_rate = b_eta, adj_prior_shape = a_eta_dash, adj_prior_rate = b_eta_dash)

#-------------------------------------------------------------------------------

# Code for making predictions

N_sam_plot = 10
# N_sam_plot = 1
# Pick N_sam_plot random samples without repetition
rand_ind = sample.int(N_samples,N_sam_plot)
#N_sam_plot = N_samples

w_star_mu = array(0, c(p_eta, N_sam_plot, n_pred))
w_star_sigma = array(0, c(p_eta, p_eta, N_sam_plot, n_pred))
w_star = array(0, c(p_eta, N_sam_plot, n_pred))
eta_sam = array(0, c(n_eta, N_sam_plot, n_pred))
eta_mu = array(0, c(n_eta, N_sam_plot, n_pred))
eta_sigma = array(0, c(n_eta, N_sam_plot, n_pred))

for (i in 1:N_sam_plot){
  print(i)
  sam_ind = rand_ind[i]
  # Construct covariance matrix for the current sample
  # I could probably package the below code into an external function
  
  # Extract the inferred parameters for the current sample
  beta_w_i = samples$beta_w[sam_ind,]
  lambda_w_i = samples$lambda_w[sam_ind,]
  lambda_eta_i = samples$lambda_eta[sam_ind]
  
  # Define the covariance matrix for the emulator weights, evaluated at the
  # training data points.
  sigma_z = matrix(0, m*p_eta, m*p_eta)
  for (j in 1:p_eta) {
    # Calculate the covariance matrix for the ith emulator weight
    sigma_z[((j-1)*m+1):(j*m), ((j-1)*m+1):(j*m)] = ARD_SE_cov(tc, lambda_w_i[j], beta_w_i[((j-1)*q+1):(j*q)], 0)
  }
  
  # Adjust the covariance matrix sigma_z to include transformed emulator and 
  # experimental error terms
  sigma_z_hat = matrix(0,m*p_eta, m*p_eta)
  sigma_z_hat = sigma_z + KTKinv/lambda_eta_i
  
  sigma_z_w_star = array(0, c(m*p_eta, p_eta, n_pred))
  # Determine covariance of training data with predictions
  for (j in 1:n_pred) {
    for (k in 1:p_eta) {
      sigma_z_w_star[((k-1)*m+1):(k*m),k,j] = ARD_SE_cov_non_sym(tc, t(as.matrix(t_pred[j,])), lambda_w_i[k], beta_w_i[((k-1)*q+1):(k*q)])
    }
  }
  
  # Define correlation of emulator predictions with themselves. In an ideal world
  # I'd also look at cross-correlations but I've chosen not to do that for the 
  # sake of efficiency. The same covariance matrix can therefore be used for all
  # predictions as the autocorrelation is 1 irrespective of the prediction 
  sigma_w_star = diag(1/lambda_w_i)
  
  # Two different methods for making predictions. The first method is quicker, but 
  # I think the second is more numerically stable, which seems to make a different
  # when the variance is small. Consider using this if the simulation is taking 
  # too long
  
  # Explicitly calculating inverse, then reusing for all predictions. 
  # Ainv = solve(sigma_z_hat)
  # w_mu = t(Lsigma_z_w_star) %*% Ainv %*% z_hat
  # w_sigma = sigma_w_star - (t(Lsigma_z_w_star) %*% Ainv %*% Lsigma_z_w_star)

  # Solving using solve
  Ainv_z_hat = solve(sigma_z_hat,z_hat)
  # Store mean and covariance matrices of discrepancy and adjusted prediction Gaussian processes
  for (j in 1:n_pred){
    w_star_mu[,i,j] = t(sigma_z_w_star[,,j]) %*% Ainv_z_hat
    w_star_sigma[,,i,j] = sigma_w_star - (t(sigma_z_w_star[,,j]) %*% solve(sigma_z_hat,sigma_z_w_star[,,j]))
    # Sample from the Gaussian process
    w_star[,i,j] = mvrnorm(n = 1, w_star_mu[,i,j], w_star_sigma[,,i,j])
    # Generate individual Samples
    eta_sam[,i,j] = (K_eta %*% w_star[,i,j])*sd_dt + mu_dt
    # Also look at samples of the mean for output
    eta_mu[,i,j] = (K_eta %*% w_star_mu[,i,j])*sd_dt + mu_dt
    # Calculate the diagonal terms of the covariance matrix. Take the standard 
    # deviation as this is more meaningful
    eta_sigma[,i,j] = sqrt(as.matrix(rowSums((K_eta %*% w_star_sigma[,,i,j]) * K_eta)))*sd_dt
  }
}

# Get average full-field displacement across all samples
#eta_sam_matrix = matrix(0, nrow = n_eta, ncol = n_pred)
#for (i in 1:n_pred){
#  eta_sam_matrix[,i] = rowMeans(eta_sam[,i,])
#}
eta_sam_matrix = apply(eta_sam, 3, rowMeans)
eta_mu_matrix = apply(eta_mu, 3, rowMeans)
eta_sigma_matrix = apply(eta_sigma, 3, rowMeans)
# Compare against known output

cross_val_disp = fread("inputs/LHSDesign40x4_1_fixed_100kN.csv")
cross_val = matrix(0,n_eta,n_pred)
for (i in 1:n_pred){
  cross_val[,i] = cross_val_disp[[as.name(sprintf('w_%d', i))]]
}

# Quick cross-validation exercise - the second node is the reference point at 
# which the load is applied. Look at this to check predictions are approximately
# the right magintude, which seems to be the case
View(rbind(eta_sam_matrix[2,],eta_mu_matrix[2,],cross_val[2,],eta_sigma_matrix[2,]))

# Plot histogram of reduced coefficients
dev.new(noRStudioGD = TRUE) # plot in new window
par(mfrow = c(1, p_eta))
pred_ind = 2
for (i in 1:p_eta){
  hist(w_star[i,,pred_ind], 
       main =  paste("Feature",as.character(i)),
       xlab = paste("w_star", as.character(i)),
       col = "brown1",
       breaks = 5,
       freq = FALSE,
       xlim = c(-3,3),
       cex.lab = 1.25,
       cex.axis = 1.25)
  # overlay plot of prior distribution
  w_plot = seq(-3,3, length.out = 100)
  prior_plot = dnorm(w_plot,mean = 0, sd = 1)
  lines(w_plot,prior_plot,lwd=3,"col"="blue")
}


# Write all output to text files for plotting
# Too big - consider just outputting samples for a single prediction
# write.csv(eta_sam, "outputs/eta_sam.csv", row.names = FALSE)
# write.csv(eta_sam_matrix, "outputs/eta_sam_mu.csv", row.names = FALSE)
write.csv(eta_mu_matrix, "outputs/eta_mu_mu.csv", row.names = FALSE)
write.csv(eta_sigma_matrix, "outputs/eta_sigma_mu.csv", row.names = FALSE)
