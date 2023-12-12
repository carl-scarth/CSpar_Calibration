
# This code is an application of the multivariate calibration formulation 
# proposed in 'Computer Model Calibration using High-dimensional Output', by
# Higdon et al, JASA, to a CSpar finite element model with uncertain inputs 
# using DIC data from one experiment. A simplified version of Higdon et al. 
# with n = 1 is implemented here. Sampling is undertaken in stan. This code 
# handles pre and post processing.

library(data.table)
library(rstan)
library(matrixStats)

# Set current working directory. This should be modified to match the directory
# of the user
setwd("C:/Users/cs2361/Documents/CSpar_Calibration/")

# include functions which are called in this code
source("source/transform_input_output.R")
source("source/interpolate_data.R")
source("source/utils.R")
source("source/dimension_reduction.R")
source("source/prior_posterior_plots.R")
source("source/gp_predictions.R")

#-------------------------------------------------------------------------------

# Set up parameters which govern the formulation
# Might be able to delete some of these later
disp_str = "u" # String which identifies the displacement component of interest (u,v, or w)
DIC_coord_labels = c("x_proj","y_proj","z_proj") # Strings used to identify coordinates in DIC point_cloud
q_y = length(disp_str) # Number of observation error parameters to be inferred
# Define parameters for (gamma) prior distributions on the observation error
a_y = 5.0  # Shape parameter for the lambda_y prior
b_y = 0.05 # Rate parameter for the lambda_y prior
iter = 4000 # Nuber of samples per chain
chains = 3 # Number of chains for simulation
in_file = "LHSDesign40x4" # File identifier string for input and output csvs for model
exp_data_file = "Interpolated_DIC_200kN" # File identifier string for experimental data
# exp_data_file = "subsampled_data_200kN_8" # File identifier string for experimental data
surface_elements = "nominal_shell_mesh_outer_surface_elements" # File identifier string for surface mesh connectivity
interp_model = TRUE # If true, we need to interpolate model outputs to experimental coordinates. Otherwise it is assumed this has already been done

#-------------------------------------------------------------------------------

# Define prior distribution parameters for passing to stan

# Pre-processing for BC example (Mean and coefficients of variation for Gaussian
# inputs)
E11_mu = 115.6
t_ply_mu = 0.196
E11_cov = 6.0 
t_ply_cov = 5.0
# Bounds for log-uniform inputs
#K_lb = 100.0
#K_ub = 1.0e9

E22_mu = 9.24
E22_cov = 6.0
#nu12_mu = 0.335
#nu12_cov = 12.123
#nu23_mu = 0.487
#nu23_cov = 12.0
#G12_mu = 4.826
#G12_cov = 6.0


# Define data_frame of prior parameters (this could be done via csv?)
# tf_param <- data.frame(distribution = c("Gaussian","Gaussian","Loguniform"),
#                       param_1 = c(E11_mu, t_ply_mu, log(K_lb)),
#                       param_2 = c(E11_mu*E11_cov/100, t_ply_mu*t_ply_cov/100, log(K_ub)))
# row.names(tf_param) <- c("E11","t_ply","log_K")

#log_K_mu = 16.0
#log_K_sd = 1.0

# Additional pre-processing for flange-rotation example
flange_theta_mu = 0.0
flange_theta_sigma = 4.0
#flange_theta_sigma = 5.0/3.0
# Define data_frame of prior parameters
# tf_param <- data.frame(distribution = c("Gaussian","Gaussian","Gaussian","Gaussian"),
#                      param_1 = c(E11_mu, t_ply_mu, flange_theta_mu, flange_theta_mu ),
#                      param_2 = c(E11_mu*E11_cov/100, t_ply_mu*t_ply_cov/100, flange_theta_sigma, flange_theta_sigma))
# row.names(tf_param) <- c("E11","t_ply","LFlange_theta","RFlange_theta")
tf_param <- data.frame(distribution = c("Gaussian","Gaussian","Gaussian","Gaussian"),
                                             param_1 = c(E11_mu, t_ply_mu, flange_theta_mu, flange_theta_mu),
                                             param_2 = c(E11_mu*E11_cov/100, t_ply_mu*t_ply_cov/100, flange_theta_sigma, flange_theta_sigma))
row.names(tf_param) <- c("E11","t_ply","LFlange_theta","RFlange_theta")

#tf_param <- data.frame(distribution = c("Gaussian","Gaussian","Gaussian","Gaussian","Gaussian","Gaussian","Gaussian","Gaussian","Gaussian"),
#  param_1 = c(E11_mu, E22_mu, nu12_mu, nu23_mu, G12_mu, t_ply_mu, flange_theta_mu, flange_theta_mu, log_K_mu),
#  param_2 = c(E11_mu*E11_cov/100, E22_mu*E22_cov/100, nu12_mu*nu12_cov/100, nu23_mu*nu23_cov/100, G12_mu*G12_cov/100, t_ply_mu*t_ply_cov/100, flange_theta_sigma, flange_theta_sigma, log_K_sd)
#  )
#row.names(tf_param) <- c("E11","E22","nu12","nu23","G12","t_ply","LFlange_theta","RFlange_theta","log_K")
                       

                       
#-------------------------------------------------------------------------------

# Set up simulation data (need to retain this for normalisation of inputs)

# Load in emulator training data input values from Design of Experiments. 
# in_file = "LHSDesign50x3" # File identifier string for input and output csvs
XT_sim = fread(paste("inputs/",in_file,".csv", sep = ""))

# In this example I fit the emulator to the log of spring stiffness K, which is 
# a more natural choice of values across which outputs are expected for
# variations in this input
# XT_sim$K = log(XT_sim$K)
# colnames(XT_sim)[3] = "log_K"

# Determine useful quantities from model inputs and outputs. Variable names 
# match the notation of Higdon et al. 2008
q = ncol(XT_sim)          # number of calibration inputs
tc = as.matrix(XT_sim)    # Convert to a matrix for passing to stan
m = nrow(XT_sim)          # sample size of computer simulation data

# Each row of XT_sim corresponds to a block of three columns of displacement 
# data, with a column for each component u,v,w 
# abaqus_displacements = fread(paste("inputs/",in_file,"_displacements.csv", sep=""))
abaqus_displacements = fread(paste("inputs/",in_file,"_fixed_200kN.csv", sep=""))
n_eta = nrow(abaqus_displacements) # number of output points per simulation

# Extract the displacement for the component of interest and store in a matrix
#dt_simulation = matrix(NA,nrow = n_eta, ncol = m)
#for (i in 1:m){
#  dt_simulation[,i] = abaqus_displacements[[as.name(paste(disp_str,sprintf("_%d", i),sep=""))]]
#}

# Extract the displacement for the component(s) of interest and store in a matrix
dt_simulation = matrix(NA,nrow = length(disp_str)*n_eta, ncol = m)
for (i in 1:m){
  for (j in 1:length(disp_str)) {
    dt_simulation[((j-1)*n_eta+1):(j*n_eta),i] = abaqus_displacements[[as.name(paste(disp_str[j],sprintf("_%d", i),sep=""))]]
  }
}
n_eta = n_eta*length(disp_str) # update n_eta definition for stan input

#-------------------------------------------------------------------------------

# All training data points are normalised onto the unit hypercube [0,1]^q before
# being passed to stan. The same transformation is applied to the test points

# Determine the maximum and minimum value of each input within the training data
t_min = colMins(tc)
t_max = colMaxs(tc)
# Normalise the inputs using these maximum and minimum values
tc = normalise_inputs(tc, t_min, t_max)

# For all implemented distributions the transformation of the 1st parameter is the same
tf_param$p1_trans = normalise_inputs(tf_param$param_1, t_min, t_max)
# The transformation of the 2nd parameter is different if this is a standard deviation
tf_param$p2_trans = NA
for (i in 1:q){
  # Second parameter of a Gaussian is a standard deviation
  if (tf_param$distribution[i] == "Gaussian"){
    tf_param$p2_trans[i] = normalise_inputs(tf_param$param_2[i], t_min[i], t_max[i], std = TRUE)
  } else if ((tf_param$distribution[i] == "Uniform") | (tf_param$distribution[i] == "Loguniform")){
    tf_param$p2_trans[i] = normalise_inputs(tf_param$param_2[i], t_min[i], t_max[i])
  } else {
    stop("Error: Non-implemented distribution")
  }
}

# plot Design of Experiments and test points
pairs(tc,col="blue", pch=4, 
      main = "Normalised Training Data Input values",
      cex=1.5,
      cex.labels = 1.75,
      cex.axis=1.5,
      cex.lab=1.5)

#-------------------------------------------------------------------------------

# Outputs at the training data points are standardised to have zero mean vector
# (i.e. zero mean at each node), and overall standard deviation of 1, then 
# decomposed according to basis vectors K_eta

# Standardise the data to have zero mean and unit standard deviation
outlist = standardise_vector_output(dt_simulation)
eta = as.matrix(outlist[[1]]) # Convert to matrix for stan
mu_dt = outlist[[2]]
sd_dt = outlist[[3]]

# Load in basis vectors previously generated by emulator code
# K_eta = as.matrix(fread(paste("outputs/basis_",in_file,".csv", sep = "")))
K_eta = as.matrix(fread(paste("outputs/basis_",in_file,"_u.csv", sep = "")))
p_eta = ncol(K_eta) # Number of basis functions retained for the emulator from SVD

#-------------------------------------------------------------------------------

# Load in the experimental data, and standardise using the same method as used  
# for the model output. This requires interpolation of the mean model output to 
# the DIC point cloud coordinates
# Alongside displacements each row has entries for the element to which the 
# measurement has been matched, and its natural coordinates within the element
experimental_data = as.data.frame(fread(paste("inputs/", exp_data_file, ".csv", sep = "")))
n_y = nrow(experimental_data)# Number of observations
exp_str = paste(disp_str, "_rot", sep="")
exp_displacement = experimental_data[[as.name(exp_str)]] # Displacement component

# exp_displacement = experimental_data[[as.name(paste(disp_str, "_rot", sep=""))]] # Displacement component
if (interp_model) {
  y_element = py_to_R(experimental_data$Element) # Matched element indices. Must be converted from Python to R convention
  hr = as.matrix(experimental_data[c("h","r")]) # Matched natural coordinates
  # Interpolate the training data mean. Skip the first two nodes in the output as 
  # these are reference points which are not referenced by the connectivity file.
  # mu_y = as.vector(intp_nodes_to_cloud(y_element, hr, as.matrix(mu_dt), conn_file = paste("inputs/", surface_elements, ".csv", sep=""), skip_nodes=2))
  mu_y = rep(NA, q_y*n_y)
  for (i in 1:length(disp_str)){
    mu_dt_i = mu_dt[((i-1)*n_eta/q_y+1):(i*n_eta/q_y)]
    mu_y[((i-1)*n_y+1):(i*n_y)] = as.vector(intp_nodes_to_cloud(y_element, hr, as.matrix(mu_dt_i), conn_file = paste("inputs/", surface_elements, ".csv", sep=""), skip_nodes=2))
  }
} else {
  mu_y = mu_dt
}

# Calculate residuals of experimental error
residual = as.vector(exp_displacement) - mu_y
rel_error = (abs(residual)/abs(mu_y))*100
# Write to CSV
if (q_y == 1) {
  out_frame = cbind(experimental_data[DIC_coord_labels],training_data_mean = mu_y,residual,abs(residual),rel_error)
} else {
  out_frame = experimental_data[DIC_coord_labels]
  for (i in 1:q_y) {
    out_frame[paste("training_data_mean_", disp_str[i], sep="")] = mu_y[((i-1)*n_y+1):(i*n_y)]
    out_frame[paste("residual_", disp_str[i], sep="")] = residual[((i-1)*n_y+1):(i*n_y)]
    out_frame[paste("abs_residual_", disp_str[i], sep="")] = abs(residual[((i-1)*n_y+1):(i*n_y)])
    out_frame[paste("rel_error_", disp_str[i], sep="")] = rel_error[((i-1)*n_y+1):(i*n_y)]
  }
}
write.csv(out_frame, paste("outputs/",in_file,"_mean_error.csv", sep = ""), row.names = FALSE)

# Centre the experimental data and convert to vector to pass to stan
# as.vector reshapes from matrix to column vector
exp_displacement_cen = (as.vector(exp_displacement) - mu_y)/sd_dt
y = as.vector(exp_displacement_cen)

# We also need to interpolate the basis functions, K, (determined above using 
# SVD) to the DIC point cloud locations.
if (interp_model) {
  K_y = matrix(NA, nrow = n_y*q_y, ncol = p_eta)
  for (i in 1:q_y){
    K_y[((i-1)*n_y+1):(i*n_y),] = intp_nodes_to_cloud(y_element, hr, K_eta[((i-1)*n_eta/q_y+1):(i*n_eta/q_y),], conn_file = paste("inputs/", surface_elements, ".csv", sep=""), skip_nodes=2)
  }
  if (q_y == 1) {
    out_frame = as.data.frame(K_y)
    for (i in 1:p_eta) {colnames(out_frame)[i] = sprintf("K_y,%d",i)} 
  } else {
    out_frame = as.data.frame(matrix(NA, nrow = n_y, ncol = 0))
    for (i in 1:q_y) {
      for (j in 1:p_eta) {
        out_frame[[as.name(sprintf("K_y_%s_%d",disp_str[i],j))]] = K_y[((i-1)*n_y+1):(i*n_y),j]
      }
    }
  }
  out_frame = cbind(experimental_data[DIC_coord_labels], out_frame)
  # output interpolated bases for plotting
  write.csv(out_frame, paste("outputs/",in_file,"_interpolated_basis.csv", sep = ""), row.names = FALSE)
} else {
  K_y = as.matrix(K_eta)
}

#-------------------------------------------------------------------------------

# specify W_y, the (prior) precision of the observation error.

# The commented code below deals with specifying both an iid noise error, and 
# one due to an isotropic shift applied to every point
# For now I pass the identity matrix, and place a weaker prior on the observation
# error in Stan
# sigma_error = 0.01 # standard deviation associated with noise (a choice of 0.005 would also be reasonable if this is too large)
# sigma_shift = 0
# Sigma_y = diag(rep(sigma_error^2,n_y)) + matrix(sigma_shift^2, nrow = n_y, ncol = n_y)
# Standardise using the model output variance 
# Sigma_y = Sigma_y/(sd_dt^2)
# Convert covariance matrix to a precision 
# W_y = solve(Sigma_y)

# The below code doesn't work when dealing with large quantities of experimental
# data. I've implemented a method which instead works with a vector, assuming
# The precision matrix is diagonal. This will need to be adapted to work with 
# other types of precision matrix
# Directly pass the identity matrix
# W_y = diag(rep(1.0,n_y))
W_y = rep(1.0,n_y*q_y)

#-------------------------------------------------------------------------------

# Reduce the dimension of the output data and calculate associated quantities 
# for input to Stan
# First calculate quantities related to the emulator
processed_data = reduce_dimension_emulator(eta, K_eta)
# Adjusted prior parameters for the basis expansion truncation error
w_hat = processed_data[[1]] # Reduced-dimensional outputs
KTKinv = processed_data[[2]] # Inverse of the product of basis matrices

processed_data = reduce_dimension_calibration(y, K_y, W_y, q_y=q_y, a_y=a_y, b_y=b_y)
a_y_dash = processed_data[[1]]
b_y_dash = processed_data[[2]]
u_hat = processed_data[[3]]
# BTWyBinv = processed_data[[4]]
BTWyB = processed_data[[4]]

# Group together coefficient samples from experimental data and model into a 
# single vector
z_hat = c(u_hat,w_hat)


# Force the adjusted parameters of the observation error prior to specified 
# values to overcome over-constraint issues
a_y_dash = array(rep(a_y, q_y), dim=1) # Have to use array to prevent errors in stan with q_y = 1
b_y_dash = array(rep(b_y, q_y), dim=1)

#-------------------------------------------------------------------------------

# Load in fixed values of emulator parameters, determined separately, e.g. via
# MLE or MAP estimate
# emulator_parameters = c(as.matrix(read.table(paste("outputs/emulator_modes_",in_file,".csv", sep=""), sep = ",", header = TRUE)))
emulator_parameters = c(as.matrix(read.table(paste("outputs/emulator_modes_",in_file,"_u.csv", sep=""), sep = ",", header = TRUE)))
rho_w = emulator_parameters[1:(p_eta*q)]
lambda_w = emulator_parameters[(p_eta*q + 1):(p_eta*(q+1))]
lambda_eta = emulator_parameters[p_eta*(q+1)+1]

# Set up the environment for Stan, pass arguments, and run stan model
# Settings taken from: https://betanalpha.github.io/assets/case_studies/gaussian_processes.html#21_Simulating_From_A_Gaussian_Process
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
util = new.env()

# List of arguments to pass to stan
stan_data = list(m=m, q=q, n_eta=n_eta, n_y=n_y, p_eta=p_eta, q_y=q_y, a_y_dash = a_y_dash,
                 b_y_dash = b_y_dash, lambda_eta = lambda_eta, z_hat = z_hat,
                 tf_param_1=tf_param$p1_trans, tf_param_2=tf_param$p2_trans, 
                 rho_w = rho_w, lambda_w = lambda_w, tc = tc, KTKinv = KTKinv, 
                 BTWyB = BTWyB)
                 # BTWyBinv = BTWyBinv)

# Set via variable as in emulator case
fit = stan(file = "source/full_field_calibration_fixed_em.stan",
           data = stan_data,
           iter = iter,
           chains = chains,
           # chains = 1,
           # iter = 1,
           model_name = "full_field_calibration")

#-------------------------------------------------------------------------------

# Post-process and plot the simulation data

# plot trace plots
stan_trace(fit, pars = c("tf", "lambda_y"))
# Print summary of results
print(fit, pars = c("tf", "lambda_y")) # , "lambda_w", "lambda_v","lambda_y","lambda_eta"))

# extract samples from stan output
samples <- extract(fit)
tf <- samples$tf              # Calibrated model inputs
lambda_y <- samples$lambda_y  # Observation error magnitude
N_samples <- dim(tf)[1]       # Total number of samples post warm-up

# Extract label of inputs for plots
labels = colnames(XT_sim)

# Plot observation error precision
lambda_hist(lambda_y, prior_shape = a_y, prior_rate = b_y, label = "lambda_y") #, adj_prior_shape = a_y_dash, adj_prior_rate = b_y_dash)

# Transform calibrated inputs onto their original scale for plotting and output
tf_trans = rescale_inputs(tf, t_min, t_max)
# Plot prior and posterior distribution of the calibration parameters.
calibration_inp_hist(tf_trans, tf_param = tf_param)

# estimate means and modes of the posterior distribution and print to the screen
modes = rep(0,q)
for (i in 1:q){
  modes[i] = estimate_mode(tf_trans[,i])
}
print("calibration parameter modes = ")
print(modes)
print("calibration parameter means = ")
print(colMeans(tf_trans))

# Produce pairs plots of prior and posterior distributions
prior_posterior_pairs(tf_trans, tf_param)

#-------------------------------------------------------------------------------

# Make predictions from fitted Gaussian process emulator
N_sam_plot = N_samples
# Transform correlation lengths into appropriate format
beta_w = -4.0*log(rho_w)
# Make predictions using the calibrated Gaussian process
out_list = full_field_calibration_pred_fixed_em(N_sam_plot, tc, tf, z_hat, beta_w, lambda_w, lambda_eta, lambda_y, KTKinv, BTWyBinv, K = K_eta, K_y = K_y, nugget = F, sam_gp = T, output_coeff_sam = F, output_ff_sam = F, output_coeff_mean = F, output_ff_mean = T, output_ff_std = T)

# Loop over each output,transform onto the correct scale, then write to csv for
# plotting outside of R
for (i in 1:length(out_list)){
  out_string = names(out_list[i])
  out_i = out_list[[out_string]]
  # If the quantity is in full-field, transform back onto it's original scale
  # First, identify the correct mean vector to use for the transformation
  if (grepl("y", out_string, fixed = TRUE)){
    mu_out = mu_y
  } else {
    mu_out = mu_dt
  }
  if (grepl("sigma", out_string, fixed=TRUE) & grepl("eta", out_string, fixed=TRUE)){
    out_i = rescale_vector_output(out_i, mu_out, sd_dt, std=TRUE)
  } else if (grepl("eta", out_string, fixed=TRUE)) {
    # I encountered a bug here as mu_y was a matrix, not a vector. I've corrected 
    # this in the above code but haven't tested all the way through. I suggest 
    # that if this happens again it's worth checking this first
    out_i = rescale_vector_output(out_i, mu_out, sd_dt)
  }
  # We don't want to output covariance matrices. All other data sets have fewer
  # than two dimensions
  if (length(dim(out_i)) <= 2) {
    print(paste("writing outputs\\", out_string, "_", in_file, ".csv", sep=""))
    write_output(out_i, out_string, in_file)
  }
}

# Plot histogram of calibrated expansion coefficients
out_list = full_field_calibration_pred_fixed_em(N_sam_plot, tc, tf, z_hat, beta_w, lambda_w, lambda_eta, lambda_y, KTKinv, BTWyBinv, nugget = F, sam_gp = T, output_coeff_sam = T, output_ff_sam = F, output_coeff_mean = F, output_ff_mean = F, output_ff_std = F)
w_star_hist(out_list$w_star, w_lim = c(-3,3))
