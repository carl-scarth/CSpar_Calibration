
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

p_eta = 7 # Number of basis functions retained for the emulator from SVD
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
# in_file = "LHSDesign50x3" # File identifier for input and output csvs
in_file = "LHSDesign40x4" # File identifier string for input and output csvs
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
beta_w <- samples$beta_w         # Transformed Correlation lengths
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

# Make predictions from fitted Gaussian process emulator
N_sam_pred = 1000 # Required number of prediction samples
# Make predictions. Request only averages of the full-field across the posterior
# uncertainty
out_list = full_field_gp_pred(N_sam_pred, tc, z_hat, t_pred, beta_w, lambda_w, lambda_eta, K_eta, KTKinv, sam_gp = TRUE, output_coeff_sam = FALSE, output_ff_sam = FALSE, output_coeff_mean = FALSE, output_ff_mean = TRUE)

# extract quantities of interst from output, transform back onto the original
# (un-standardised) scale, then write to csv
if ("eta_mu_mu" %in% names(out_list)){
  eta_mu_mu = out_list$eta_mu_mu
  eta_sigma_mu = out_list$eta_sigma_mu
  
  # Convert back on true scale
  eta_mu_mu = eta_mu_mu*sd_dt + mu_dt
  eta_sigma_mu = eta_sigma_mu*sd_dt
  
  write_output(eta_mu_mu, "eta_mu_mu", in_file)
  write_output(eta_sigma_mu, "eta_sigma_mu",in_file)
}
if ("eta_sam_mu" %in% names(out_list)){
  eta_sam_mu = out_list$eta_sam_mu
  eta_sam_mu = eta_sam_mu*sd_dt + mu_dt
  write_output(eta_sam_mu, "eta_sam_mu", in_file)
}
# Are there individual samples to be written to file?
if ("eta_mu" %in% names(out_list)){
  eta_mu = out_list$eta_mu
  eta_mu = eta_mu*sd_dt + mu_dt
  eta_sigma = out_list$eta_sigma
  eta_sigma = eta_sigma*sd_dt
  write_output_samples(eta_mu, "eta_mu", in_file)
  write_output_samples(eta_sigma, "eta_sigma", in_file)
}
if ("eta_sam" %in% names(out_list)){
  eta_sam = out_list$eta_sam
  eta_sam = eta_sam*sd_dt + mu_dt
  write_output_samples(eta_sam, "eta_sam", in_file)
}

#-------------------------------------------------------------------------------

# Example cross validation exercise, where comparison is made with another set 
# of samples with known output

# Compare against known output for separate set of samples
cross_val_disp = fread("inputs/LHSDesign40x4_1_fixed_100kN.csv")
cross_val = matrix(0,n_eta,n_pred)
for (i in 1:n_pred){
  cross_val[,i] = cross_val_disp[[as.name(sprintf('w_%d', i))]]
}
  
# Quick cross-validation exercise, which compares the displacement of a 
# reference point at the spar tip in the model with emulator predictions
View(rbind(eta_sam_mu[2,],eta_mu_mu[2,],cross_val[2,],eta_sigma_mu[2,]))
