
# Fit an emulator to the full-field output from a C-Spar Finite Element model
# Input data are the Longitudinal Modulus E11, ply thickness, and the torsional
# stiffness of a spring representing an uncertain boundary condition
# Output data are the axial displacements of the nodes of the FE model across
# a range of load steps of a nonlinear model
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
source("source/utils.R")
source("source/dimension_reduction.R")
source("source/prior_posterior_plots.R")
source("source/abaqus_json.R")
source("source/gp_predictions.R")
source("source/transform_input_output.R")

#-------------------------------------------------------------------------------

# Set up parameters which govern the formulation

p_eta = 11 # Number of basis functions retained for the emulator from SVD
exp_tol = 1e-6 # Tolerance variance fraction used to assess SVD convergence
disp_str = "w" # String which identifies the displacement component of interest (u,v, or w)
print_svd_output = TRUE # Print diagnostic output of svd to the terminal?
export_modes = TRUE # Calculate modes of emulator hyperparameters and write to file?
# Define parameters of the gamma prior on the error associated with truncating
# the series expansion for the model output
a_eta = 1.0     # Shape parameter for the lambda_eta prior
b_eta = 0.0001  # Rate parameter for the lambda_eta prior
iter = 4000 # Number of samples per chain
chains = 3 # Number of chains for simulation

#-------------------------------------------------------------------------------

# Set up simulation data

# Load in emulator training data input values from Design of Experiments. 
in_file = "LHSDesign40x4" # File identifier string for input and output files
XT_sim = fread(paste("inputs/",in_file,".csv", sep = ""))

# Determine useful quantities from model inputs and outputs. Variable names 
# match the notation of Higdon et al. 2008
q = ncol(XT_sim)          # number of calibration inputs
tc = as.matrix(XT_sim)    # Convert to a matrix for passing to stan
m = nrow(XT_sim)          # sample size of computer simulation data

# Load in training data output displacement values from Abaqus. I've used a
# similar naming convention to the inputs to automate changes. 
# Outputs are structured in a json across samples and load increments
abaqus_displacements <- fromJSON(file = paste("inputs/",in_file,"_output_struct.json", sep=""))
# Extract displacements from json, and store as matrix where each column is a 
# training sample, and the rows are the concatenation of displacement output
# across all output frames
sorted_data = extract_const_frame(abaqus_displacements,"w")
dt_simulation = sorted_data[[1]]
n_nodes = sorted_data[[2]] # Number of nodes
n_frames = sorted_data[[3]] # Number of frames
n_eta = nrow(dt_simulation) # total number of output points per simulation

#-------------------------------------------------------------------------------

# Load in test points at which predictions are required
XT_pred = fread("inputs/LHSDesign40x4_1.csv")
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

out_basis = svd_basis(eta, p_eta = p_eta, exp_tol = exp_tol, 
                      print_output = print_svd_output, csv_label = paste("nonlinear_",in_file,sep=""))

K_eta = out_basis[[1]]
p_eta = out_basis[[2]] # Used to determine p_eta automatically if not provided as an argument, otherwise this is unchanged

# Write basis to .json file
out_json = basis_mean_to_json(n_frames, n_nodes, K_eta, mu_dt)
write(out_json, paste("outputs/basis_nonlinear_",in_file,".json",sep=""))

#-------------------------------------------------------------------------------

# Reduce the dimension of the output data and determine associated quantities
# for input to Stan
processed_data = reduce_dimension_emulator(eta, K_eta, a_eta=a_eta, b_eta=b_eta)
# Adjusted prior parameters for the basis expansion truncation error
a_eta_dash = processed_data[[1]]
b_eta_dash = processed_data[[2]]
z_hat = processed_data[[3]] # Reduced-dimensional outputs
KTK_inv = processed_data[[4]] # Inverse of the product of basis matrices

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
  write.csv(modes, paste("outputs/nonlinear_emulator_modes_",in_file,".csv", sep=""), row.names = FALSE)
}

# Plot correlation parameter histograms
full_field_rho_hist(rho_w, p_eta, inp_labels = labels)
# Plot emulator precision parameters for emulator
full_field_lambda_hist(lambda_w, p_eta)
# Plot emulator trunction error precision
lambda_hist(lambda_eta, prior_shape = a_eta, prior_rate = b_eta, label = "lambda_eta", adj_prior_shape = a_eta_dash, adj_prior_rate = b_eta_dash)

#-------------------------------------------------------------------------------

# Make predictions from fitted Gaussian process emulator, and write these
# predictions to a .json file

N_sam_pred = 1000 # Required number of prediction samples
# Make predictions. Request only averages of the full-field across the posterior
# uncertainty
out_list = full_field_gp_pred(N_sam_pred, tc, z_hat, t_pred, beta_w, lambda_w, lambda_eta, K_eta, KTKinv, sam_gp = FALSE, output_coeff_sam = FALSE, output_ff_sam = FALSE, output_coeff_mean = FALSE, output_ff_mean = TRUE)

# extract quantities of interest from output, transform back onto the original
# (un-standardised) scale, then write to json
# List of full-field outputs which may be in out_list
out_strings = c("eta_mu_mu", "eta_sigma_mu", "eta_sam_mu", "eta_mu", "eta_sigma", "eta_sam")
# Loop over each output, check if it has been requested, then append to the list
# to be written to json
json_list <- list()
for (i in 1:length(out_strings)){
  if (out_strings[i] %in% names(out_list)){
    out_i = out_list[[out_strings[i]]]
    # Transform outputs back onto their individual scale
    # If the output is a standard deviation a different transformation is required
    if (grepl("sigma", out_strings[i], fixed=TRUE)){
      out_i = rescale_vector_output(out_i, mu_dt, sd_dt, std=TRUE)
    } else {
      out_i = rescale_vector_output(out_i, mu_dt, sd_dt)
    }
    # Append transformed values to list
    out_i = list(out_i)
    names(out_i) = out_strings[i]
    json_list = append(json_list, out_i)
  }
}

# Write json to file
out_json = gp_pred_to_json(json_list, n_frames, n_nodes, n_pred=n_pred, n_post_sam = N_sam_pred)
write(out_json, paste("outputs/gp_predictions_nonlinear_",in_file,".json",sep=""))
