
# Fit an emulator to the full-field output from a C-Spar Finite Element model
# Follows emulator aspect of D. Higdon et al, "Computer Model Calibration Using 
# High-Dimensional Output",Journal of the American Statistical Association,2008.
# This code handles pre and post processing.

library(data.table)
library(rstan)
library(matrixStats)

# Set current working directory. This should be modified to match the directory
# of the user
setwd("C:/Users/cs2361/Documents/CSpar_Calibration/")

# include functions which are called in this code
source("source/transform_input_output.R")
source("source/utils.R")
source("source/dimension_reduction.R")
source("source/prior_posterior_plots.R")
source("source/gp_predictions.R")

#-------------------------------------------------------------------------------

# Set up parameters which govern the formulation

p_eta = 40 # Number of basis functions retained for the emulator from SVD
exp_tol = 1e-6 # Tolerance variance fraction used to assess SVD convergence
# disp_str = c("u","v","w") # String which identifies the displacement component of interest (u,v, or w)
disp_str = "w"
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
# in_file = "LHSDesign50x3" # File identifier for input and output csvs
in_file = "LHSDesign100x9" # File identifier string for input and output csvs
XT_sim = fread(paste("inputs/",in_file,".csv", sep = ""))

# In this example I fit the emulator to the log of spring stiffness K, which is 
# a more natural choice of values across which outputs are expected for
# variations in this input
#XT_sim$K = log(XT_sim$K)
#colnames(XT_sim)[3] = "log_K"

# Determine useful quantities from model inputs and outputs. Variable names 
# match the notation of Higdon et al. 2008
q = ncol(XT_sim)          # number of calibration inputs
tc = as.matrix(XT_sim)    # Convert to a matrix for passing to stan
m = nrow(XT_sim)          # sample size of computer simulation data

# Load in training data output displacement values from Abaqus. I've used
# a similar naming convention to the inputs to automate changes. 
# Each row of XT_sim corresponds to a block of three columns of displacement 
# data, with a column for each component u,v,w 
# abaqus_displacements = fread(paste("inputs/",in_file,"_displacements.csv", sep=""))
abaqus_displacements = fread(paste("inputs/",in_file,"_fixed_200kN.csv", sep=""))
n_eta = nrow(abaqus_displacements) # number of output points per simulation
# Extract the displacement for the component(s) of interest and store in a matrix
dt_simulation = matrix(NA,nrow = length(disp_str)*n_eta, ncol = m)
for (i in 1:m){
  for (j in 1:length(disp_str)) {
    dt_simulation[((j-1)*n_eta+1):(j*n_eta),i] = abaqus_displacements[[as.name(paste(disp_str[j],sprintf("_%d", i),sep=""))]]
  }
}
n_eta = n_eta*length(disp_str) # update n_eta definition for stan input

#-------------------------------------------------------------------------------

# Load in test points at which predictions are required
# XT_pred = fread("inputs/LHSDesign50x3_1.csv")
XT_pred = fread("inputs/LHSDesign100x9_1.csv")
# Repeat the log transformation for the test data
# XT_pred$K = log(XT_pred$K)
# colnames(XT_pred)[3] = "log_K"
t_pred = as.matrix(XT_pred)
n_pred = nrow(t_pred) # number of predictions
  
#-------------------------------------------------------------------------------

# All training data points are normalised onto the unit hypercube [0,1]^q before
# being passed to stan. The same transformation is applied to the test points

# Determine the maximum and minimum value of each input within the training data
t_min = colMins(tc)
t_max = colMaxs(tc)
# Normalise the inputs using these maximum and minimum values
tc = normalise_inputs(tc, t_min, t_max)
t_pred = normalise_inputs(t_pred, t_min, t_max)

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

# Outputs at the training data points are standardised to have zero mean vector
# (i.e. zero mean at each node), and overall standard deviation of 1, then 
# decomposed according to basis vectors K_eta

outlist = standardise_vector_output(dt_simulation)
eta = as.matrix(outlist[[1]]) # Convert to matrix for stan
mu_dt = outlist[[2]]
sd_dt = outlist[[3]]

# Write model mean to csv file
write_output(as.matrix(mu_dt,n_row = n_eta),"training_data_mean",in_file)

# Calculate reduced-dimensional basis of training data
#out_basis = svd_basis(eta, p_eta = p_eta, exp_tol = exp_tol, 
#                      print_output = print_svd_output, csv_label = in_file)
out_basis = svd_basis(eta, exp_tol = exp_tol, print_output = print_svd_output, csv_label = in_file)
K_eta = out_basis[[1]]
p_eta = out_basis[[2]] # Used to determine p_eta automatically if not provided as an argument, otherwise this is unchanged

#-------------------------------------------------------------------------------

# Reduce the dimension of the output data and calculate associated quantities 
# for input to Stan
processed_data = reduce_dimension_emulator(eta, K_eta, a_eta=a_eta, b_eta=b_eta)
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
stan_data = list(m=m, q=q, n_eta=n_eta, p_eta=p_eta, linear_mean = 0, a_eta_dash = a_eta_dash, 
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
beta_w <- samples$beta_w         # Transformed Correlation lengths
lambda_w <- samples$lambda_w     # Emulator precisions
lambda_eta <- samples$lambda_eta # Expansion truncation error 
N_samples <- dim(rho_w)[1]       # Total number of samples post warm-up

# Extract label of inputs for plots
labels = colnames(XT_sim)

# If required, estimate the modes of emulator hyperparameters and write to a csv
if (export_modes){
  modes = full_field_emulator_modes(rho_w, lambda_w, lambda_eta)
  #write.csv(modes, paste("outputs/emulator_modes_",in_file,".csv", sep=""), row.names = FALSE)
  write.csv(modes, paste("outputs/emulator_modes_",in_file,"_multi.csv", sep=""), row.names = FALSE)
}

# Plot correlation parameter histograms
full_field_rho_hist(rho_w, p_eta, inp_labels = labels)
# Plot emulator precision parameters for emulator
full_field_lambda_hist(lambda_w, p_eta)
# Plot emulator truncation error precision
lambda_hist(lambda_eta, prior_shape = a_eta, prior_rate = b_eta, label = "lambda_eta", adj_prior_shape = a_eta_dash, adj_prior_rate = b_eta_dash)

#-------------------------------------------------------------------------------

# Make predictions from fitted Gaussian process emulator
N_sam_pred = 500 # Required number of prediction samples
# Make predictions. Request only averages of the full-field across the posterior
# uncertainty
out_list = full_field_gp_pred(N_sam_pred, tc, z_hat, t_pred, beta_w, lambda_w, lambda_eta, K_eta, KTKinv, sam_gp = FALSE, output_coeff_sam = FALSE, output_ff_sam = FALSE, output_coeff_mean = FALSE, output_ff_mean = TRUE)

# Loop over each output,transform onto the correct scale, then write to csv for
# plotting outside of R
for (i in 1:length(out_list)){
  out_string = names(out_list[i])
  out_i = out_list[[out_string]]
  # If the quantity is in full-field, transform back onto it's original scale
  if (grepl("sigma", out_string, fixed=TRUE) & grepl("eta", out_string, fixed=TRUE)){
    out_i = rescale_vector_output(out_i, mu_dt, sd_dt, std=TRUE)
  } else if (grepl("eta", out_string, fixed=TRUE)) {
    out_i = rescale_vector_output(out_i, mu_dt, sd_dt)
  }
  # We don't want to output covariance matrices for w_star.
  if (!(grepl("sigma", out_string, fixed=TRUE) & grepl("w_star", out_string, fixed=TRUE))) {
    print(paste("writing outputs/", out_string, "_", in_file, ".csv", sep=""))
    write_output(out_i, out_string[1], in_file)
  } else {
    print(paste("Not outputting ", out_string, sep=""))
  }
}
