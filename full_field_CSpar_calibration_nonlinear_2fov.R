
# Example for when there are two experimental fields of view

# This code is an application of the multivariate calibration formulation 
# proposed in 'Computer Model Calibration using High-dimensional Output', by
# Higdon et al, JASA, to a CSpar finite element model with uncertain inputs 
# using DIC data from one experiment. Output data are the axial displacements of
# the nodes of the FE model/DIC point cloud across a range of load steps. A
# simplified version of Higdon et al. with n = 1 is implemented here. Sampling 
# is undertaken in stan. This code handles pre and post processing.

library(data.table)
library(rstan)
library(matrixStats)

# Set current working directory. This should be modified to match the directory
# of the user
setwd("C:/Users/cs2361/Documents/CSpar_Calibration/")

# include functions which are called in this code
source("source/abaqus_json.R")
source("source/transform_input_output.R")
source("source/utils.R")
source("source/interpolate_data.R")
source("source/dimension_reduction.R")
source("source/prior_posterior_plots.R")
source("source/gp_predictions.R")

#-------------------------------------------------------------------------------

# Set up parameters which govern the formulation
# Might be able to delete some of these later
# p_eta = 7 # Number of basis functions retained for the emulator from SVD
# exp_tol = 1e-6 # Tolerance variance fraction used to assess SVD convergence
disp_str = c("u") # String which identifies the displacement component of interest (u,v, or w)
q_y = length(disp_str)
DIC_coord_labels = c("x_proj","y_proj","z_proj") # Strings used to identify coordinates in DIC point_cloud
a_y = 5.0  # Shape parameter for the lambda_y prior
b_y = 0.05 # Rate parameter for the lambda_y prior
# Define parameters of the gamma prior on the error associated with truncating
# the series expansion for the model output
# a_eta = 1.0     # Shape parameter for the lambda_eta prior
# b_eta = 0.0001  # Rate parameter for the lambda_eta prior 
iter = 5000 # Number of samples per chain
chains = 3 # Number of chains for simulation
in_file = "LHSDesign100x8" # File identifier string for input and output csvs
exp_data_file_1 = "Interpolated_DIC_multistep_160kN_LC" # Identifier of file containing DIC data
exp_data_file_2 = "Interpolated_DIC_multistep_160kN_RC" # Identifier of file containing DIC data
# surface_elements = "nominal_shell_mesh_outer_surface_elements" # File identifier string for surface mesh connectivity
surface_elements = "new_spar_mesh_outer_surface_elements" # File identifier string for surface mesh connectivity
p_sub = 2 # Use a subset of the basis functions
use_subset = F

#-------------------------------------------------------------------------------

# Define prior distribution parameters for passing to stan
# Note, input sd rather than cov
#tf_param = read.table(paste("inputs/", in_file, "_tf_param_training_data_relaxed_3.csv",sep=""), sep=",", header = TRUE, row.names = 1)
tf_param = read.table(paste("inputs/", in_file, "_tf_param_mod.csv",sep=""), sep=",", header = TRUE, row.names = 1)

# Pre-processing for BC example (Mean and coefficients of variation for Gaussian
# inputs)
# E11_mu = 115.6
# t_ply_mu = 0.196
# E11_cov = 6.0 
# t_ply_cov = 5.0
# Bounds for log-uniform inputs
# K_lb = 100.0
# K_ub = 1.0e9

# Define data_frame of prior parameters (this could be done via csv?)
# tf_param <- data.frame(distribution = c("Gaussian","Gaussian","Loguniform"),
#                       param_1 = c(E11_mu, t_ply_mu, log(K_lb)),
#                       param_2 = c(E11_mu*E11_cov/100, t_ply_mu*t_ply_cov/100, log(K_ub)))
# row.names(tf_param) <- c("E11","t_ply","log_K")

# Additional pre-processing for flange-rotation example
# flange_theta_mu = 0.0
# flange_theta_sigma = 4.0
# Define data_frame of prior parameters
# tf_param <- data.frame(distribution = c("Gaussian","Gaussian","Gaussian","Gaussian"),
#                       param_1 = c(E11_mu, t_ply_mu, flange_theta_mu, flange_theta_mu ),
#                       param_2 = c(E11_mu*E11_cov/100, t_ply_mu*t_ply_cov/100, flange_theta_sigma, flange_theta_sigma))
# row.names(tf_param) <- c("E11","t_ply","LFlange_theta","RFlange_theta")

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

# Load in training data output displacement values from Abaqus. I've used a
# similar naming convention to the inputs to automate changes. 
# Outputs are structured in a json across samples and load increments
abaqus_displacements <- fromJSON(file = paste("inputs/",in_file,"_output_struct_160kN_zeroed.json", sep=""))
# abaqus_displacements <- fromJSON(file = paste("inputs/",in_file,"_downsam_2.json", sep=""))
# Extract displacements from json, and store as matrix where each column is a 
# training sample, and the rows are the concatenation of displacement output
# across all output frames
sorted_data = extract_const_frame(abaqus_displacements,disp_str)
dt_simulation = sorted_data[[1]]
n_nodes = sorted_data[[2]] # Number of nodes
n_frames = sorted_data[[3]] # Number of frames
n_eta = nrow(dt_simulation) # total number of output points per simulation

#------------------------------------------------------------------------------

# All training data points are normalised onto the unit hypercube [0,1]^q before
# being passed to stan. The same transformation is applied to the test points

# Determine the maximum and minimum value of each input within the training data
t_min = colMins(tc)
t_max = colMaxs(tc)
# Normalise the inputs using these maximum and minimum values
tc = normalise_inputs(tc, t_min, t_max)

# For all implemented distributions the transformation of the 1st parameter is the same
tf_param$p1_trans = NA
# The transformation of the 2nd parameter is different if this is a standard deviation
tf_param$p2_trans = NA
for (i in 1:q){
  # Second parameter of a Gaussian is a standard deviation
  if (tf_param$distribution[i] == "Gaussian" | tf_param$distribution[i] == "Lognormal" | tf_param$distribution[i] == "Halfnormal"){
    tf_param$p1_trans[i] = normalise_inputs(tf_param$param_1[i], t_min[i], t_max[i])
    tf_param$p2_trans[i] = normalise_inputs(tf_param$param_2[i], t_min[i], t_max[i], std = TRUE)
  } else if ((tf_param$distribution[i] == "Uniform") | (tf_param$distribution[i] == "Loguniform")){
    tf_param$p1_trans[i] = normalise_inputs(tf_param$param_1[i], t_min[i], t_max[i])
    tf_param$p2_trans[i] = normalise_inputs(tf_param$param_2[i], t_min[i], t_max[i])
  } else if (tf_param$distribution[i] == "Gamma" | tf_param$distribution[i] == "Loggamma"){
    tf_param$p1_trans[i] = tf_param$param_1[i] # Shape parameter doesn't change under transformation
    tf_param$p2_trans[i] = tf_param$param_2[i]*(t_max[i] - t_min[i]) # I think this is correct...
    stop("Error: Scale implemented correctly, but not yet shifted. Consider passing shift as parameter and shifting back within stan")
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
out_list = standardise_vector_output(dt_simulation, q_y = q_y)
# out_list = standardise_vector_output(dt_simulation)
eta = as.matrix(out_list[[1]]) # Convert to matrix for stan
mu_dt = out_list[[2]]
sd_dt = out_list[[3]]

# Load in basis vectors previously generated by emulator code
# K_eta = as.matrix(fread(paste("outputs/basis_nonlinear_",in_file,".csv", sep = "")))
K_eta = as.matrix(fread(paste("outputs/basis_nonlinear_",in_file,"_u160z.csv", sep = "")))
p_eta = ncol(K_eta) # Number of basis functions retained for the emulator from SVD
if (use_subset) {
  K_eta = K_eta[,1:p_sub,drop=F]
}


#-------------------------------------------------------------------------------

# Load in the experimental data, and standardise using the same method as used  
# for the model output. This requires interpolation of the mean model output to 
# the DIC point cloud coordinates.
# Alongside displacements each row has entries for the element to which the 
# measurement has been matched, its natural coordinates within the element, and 
# applied load increment, which matches those of the model.

exp_str = paste(disp_str, "_rot", sep="")

# Extract relevant data from first field of view
experimental_data_1 = as.data.frame(fread(paste("inputs/", exp_data_file_1, ".csv", sep = "")))
n_y1 = nrow(experimental_data_1) # Total of observations
exp_displacement_1 = experimental_data_1[(exp_str)] # Displacement component
y_element_1 = py_to_R(experimental_data_1$Element) # Matched element indices. Must be converted from Python to R convention
# hr = as.matrix(experimental_data[c("h","r")]) # Matched natural coordinates
gh_1 = as.matrix(experimental_data_1[c("g","h")]) # Matched natural coordinates
# exp_displacement = experimental_data[[as.name(paste(disp_str, "_rot", sep=""))]] # Displacement component
frame_ind_1 = experimental_data_1$Increment # Index of which frame each point belongs to
exp_coords_1 = experimental_data_1[DIC_coord_labels]

# Extract relevant data from second field of view
experimental_data_2 = as.data.frame(fread(paste("inputs/", exp_data_file_2, ".csv", sep = "")))
n_y2 = nrow(experimental_data_2) # Total of observations
exp_displacement_2 = experimental_data_2[(exp_str)] # Displacement component
y_element_2 = py_to_R(experimental_data_2$Element) # Matched element indices. Must be converted from Python to R convention
# hr = as.matrix(experimental_data[c("h","r")]) # Matched natural coordinates
gh_2 = as.matrix(experimental_data_2[c("g","h")]) # Matched natural coordinates
# exp_displacement = experimental_data[[as.name(paste(disp_str, "_rot", sep=""))]] # Displacement component
frame_ind_2 = experimental_data_2$Increment # Index of which frame each point belongs to
exp_coords_2 = experimental_data_2[DIC_coord_labels]

# Combine datasets
# have to combine after interpolation
exp_displacement = rbind(exp_displacement_1,exp_displacement_2)
#y_element = c(y_element_1,y_element_2)
#gh = rbind(gh_1,gh_2)
n_y = n_y1 + n_y2
frame_ind = c(frame_ind_1,frame_ind_2)
exp_coords = rbind(exp_coords_1,exp_coords_2)

# Interpolate the training data mean. Skip the first two nodes in the output as 
# these are reference points which are not referenced by the connectivity file
mu_y1 = rep(NA, q_y*n_y1)
mu_y2 = rep(NA, q_y*n_y2)
for (i in 1:q_y){
  mu_dt_i = mu_dt[((i-1)*n_eta/q_y+1):(i*n_eta/q_y)]
  # mu_y[((i-1)*n_y+1):(i*n_y)] = as.vector(intp_nodes_to_cloud_inc(y_element, hr, as.matrix(mu_dt_i), frame_ind, n_frames, conn_file = paste("inputs/", surface_elements, ".csv", sep=""), skip_nodes=2))
  mu_y1[((i-1)*n_y1+1):(i*n_y1)] = as.vector(intp_nodes_to_cloud_inc(y_element_1, gh_1, as.matrix(mu_dt_i), frame_ind_1, n_frames, conn_file = paste("inputs/", surface_elements, ".csv", sep=""), skip_nodes=0))
  mu_y2[((i-1)*n_y2+1):(i*n_y2)] = as.vector(intp_nodes_to_cloud_inc(y_element_2, gh_2, as.matrix(mu_dt_i), frame_ind_2, n_frames, conn_file = paste("inputs/", surface_elements, ".csv", sep=""), skip_nodes=0))
}
mu_y = c(mu_y1,mu_y2)

# Calculate residuals of experimental error
# residual = exp_displacement - mu_y
residual = c(as.matrix(exp_displacement)) - mu_y
rel_error = (abs(residual)/abs(mu_y))*100

# Write to CSV
# write.csv(cbind(experimental_data[DIC_coord_labels],training_data_mean = mu_y,residual,abs(residual),rel_error, increment = experimental_data$Increment), paste("outputs/",in_file,"_mean_error_nonlinear.csv", sep = ""), row.names = FALSE)
if (q_y == 1) {
  out_frame = cbind(exp_coords,training_data_mean = mu_y,residual,abs(residual),rel_error, Increment = frame_ind)
} else {
  out_frame = exp_coords
  for (i in 1:q_y) {
    ind = c(((i-1)*n_y+1):(i*n_y),q_y*n_y+(((i-1)*n_y+1):(i*n_y)))
    out_frame[paste("training_data_mean_", disp_str[i], sep="")] = mu_y[ind]
    out_frame[paste("residual_", disp_str[i], sep="")] = residual[ind]
    out_frame[paste("abs_residual_", disp_str[i], sep="")] = abs(residual[ind])
    out_frame[paste("rel_error_", disp_str[i], sep="")] = rel_error[ind]
    #out_frame[paste("training_data_mean_", disp_str[i], sep="")] = c(mu_y[((i-1)*n_y+1):(i*n_y)],mu_y[n_y+(((i-1)*n_y+1):(i*n_y))])
    #out_frame[paste("residual_", disp_str[i], sep="")] = residual[((i-1)*n_y+1):(i*n_y)]
    #out_frame[paste("abs_residual_", disp_str[i], sep="")] = abs(residual[((i-1)*n_y+1):(i*n_y)])
    #out_frame[paste("rel_error_", disp_str[i], sep="")] = rel_error[((i-1)*n_y+1):(i*n_y)]
  }
  out_frame["Increment"] = frame_ind
}
write.csv(out_frame, paste("outputs/",in_file,"_mean_error_nonlinear.csv", sep = ""), row.names = FALSE)

# Centre the experimental data and convert to vector to pass to stan
if (length(sd_dt) == 1) {
  exp_displacement_cen = (c(as.matrix(exp_displacement)) - mu_y)/sd_dt
} else {
  sd_y = rep(NA, q_y*n_y)
  # Otherwise create vector of component-wise model standard deviations
  for (i in 1:2*q_y){
    sd_y[((i-1)*n_y+1):(i*n_y)] = rep(sd_dt[i*n_eta/q_y], n_y)
  }
  exp_displacement_cen = (c(as.matrix(exp_displacement)) - mu_y)/sd_y
}
y = exp_displacement_cen
#y = as.vector(exp_displacement_cen)

# We also need to interpolate the basis functions, K, (determined above using 
# SVD) to the DIC point cloud locations.
if (use_subset){
  K_y1 = matrix(NA, nrow = n_y1*q_y, ncol = p_sub)
  K_y2 = matrix(NA, nrow = n_y2*q_y, ncol = p_sub)
} else {
  K_y1 = matrix(NA, nrow = n_y1*q_y, ncol = p_eta)
  K_y2 = matrix(NA, nrow = n_y2*q_y, ncol = p_eta)
}

for (i in 1:q_y){
  K_y1[((i-1)*n_y1+1):(i*n_y1),] = intp_nodes_to_cloud_inc(y_element_1, gh_1, K_eta[((i-1)*n_eta/q_y+1):(i*n_eta/q_y),,drop=F], frame_ind_1, n_frames, conn_file = paste("inputs/", surface_elements, ".csv", sep=""), skip_nodes=0)
  K_y2[((i-1)*n_y2+1):(i*n_y2),] = intp_nodes_to_cloud_inc(y_element_2, gh_2, K_eta[((i-1)*n_eta/q_y+1):(i*n_eta/q_y),,drop=F], frame_ind_2, n_frames, conn_file = paste("inputs/", surface_elements, ".csv", sep=""), skip_nodes=0)
}
K_y = rbind(K_y1, K_y2)

if (q_y == 1) {
  out_frame = as.data.frame(K_y)
  for (i in 1:ncol(K_eta)) {colnames(out_frame)[i] = sprintf("K_y,%d",i)} 
} else {
  out_frame = as.data.frame(matrix(NA, nrow = n_y, ncol = 0))
  for (i in 1:q_y) {
    for (j in 1:ncol(K_eta)) {
      out_frame[[as.name(sprintf("K_y_%s_%d",disp_str[i],j))]] = K_y[((i-1)*n_y+1):(i*n_y),j]
    }
  }
}
out_frame = cbind(exp_coords, out_frame, Increment = frame_ind)

# output interpolated bases for plotting
write.csv(out_frame, paste("outputs/",in_file,"_interpolated_basis_nonlinear.csv", sep = ""), row.names = FALSE)

#-------------------------------------------------------------------------------

# specify W_y, the (prior) precision of the observation error.
# In the long run, consider varying the observation error with load magnitude.
# For now, keep constant.

sigma_noise = 1.0 # standard deviation associated with noise (a choice of 0.005 would also be reasonable if this is too large)
# sigma_shift = 0
# Sigma_y = diag(rep(sigma_error^2,n_y)) + matrix(sigma_shift^2, nrow = n_y, ncol = n_y)
# Standardise using the model output variance 
# Sigma_y = Sigma_y/(sd_dt^2)
# Convert covariance matrix to a precision 
# W_y = solve(Sigma_y)
# W_y = diag(rep(1/(sigma_noise^2), n_y))

# The below code doesn't work when dealing with large quantities of experimental
# data. I've implemented a method which instead works with a vector, assuming
# The precision matrix is diagonal. This will need to be adapted to work with 
# other types of precision matrix
# Directly pass the identity matrix
# W_y = diag(rep(1.0,n_y))
# W_y = rep(1.0,n_y)
W_y = rep(1/(sigma_noise^2),2*n_y*q_y)

#-------------------------------------------------------------------------------

# Reduce the dimension of the output data and calculate associated quantities 
# for input to Stan
# First calculate quantities related to the emulator
processed_data = reduce_dimension_emulator(eta, K_eta)
# Adjusted prior parameters for the basis expansion truncation error
w_hat = processed_data[[1]] # Reduced-dimensional outputs
KTKinv = processed_data[[2]] # Inverse of the product of basis matrices

processed_data = reduce_dimension_calibration_multi_fov(y, K_y, W_y, c(n_y1, n_y2))
BTB = processed_data[[1]]
BTy = processed_data[[2]]

# Force the adjusted parameters of the observation error prior to specified 
# values to overcome over-constraint issues
# a_y_dash = a_y
# b_y_dash = b_y
a_y_dash = array(rep(a_y, q_y), dim=2) # Have to use array to prevent errors in stan with q_y = 1
b_y_dash = array(rep(b_y, q_y), dim=2)

#-------------------------------------------------------------------------------

# Load in fixed values of emulator parameters, determined separately, e.g. via
# MLE or MAP estimate
# emulator_parameters = c(as.matrix(read.table(paste("outputs/nonlinear_emulator_modes_",in_file,".csv", sep=""), sep = ",", header = TRUE)))
emulator_parameters = c(as.matrix(read.table(paste("outputs/nonlinear_emulator_modes_",in_file,"_u160z_mod.csv", sep=""), sep = ",", header = TRUE)))
if (use_subset) {
  rho_w = emulator_parameters[1:(p_sub*q)]
  lambda_w = emulator_parameters[(p_eta*q + 1):(p_eta*q+p_sub)]
} else {
  rho_w = emulator_parameters[1:(p_eta*q)]
  lambda_w = emulator_parameters[(p_eta*q + 1):(p_eta*(q+1))]
}
lambda_eta = emulator_parameters[p_eta*(q+1)+1]


# Set up the environment for Stan, pass arguments, and run stan model
# Settings taken from: https://betanalpha.github.io/assets/case_studies/gaussian_processes.html#21_Simulating_From_A_Gaussian_Process
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
util = new.env()

if (use_subset) {
  p_eta = p_sub
}

# List of arguments to pass to stan

stan_data = list(m=m, q=q, n_eta=n_eta, n_y=n_y, p_eta=p_eta, q_y=2, a_y_dash = a_y_dash,
                 b_y_dash = b_y_dash, lambda_eta = lambda_eta, w_hat = w_hat,
                 tf_param_1=tf_param$p1_trans, tf_param_2=tf_param$p2_trans, 
                 rho_w = rho_w, lambda_w = lambda_w, tc = tc, KTKinv = KTKinv, 
                 BTB = BTB, BTy = BTy)

#stan_data = list(m=m, q=q, n_eta=n_eta, n_y=n_y, p_eta=p_eta, q_y=q_y, a_y_dash = a_y_dash,
#                 b_y_dash = b_y_dash, lambda_eta = lambda_eta, z_hat = z_hat,
#                 tf_param_1=tf_param$p1_trans, tf_param_2=tf_param$p2_trans, 
#                 rho_w = rho_w, lambda_w = lambda_w, tc = tc, KTKinv = KTKinv, 
#                 BTWyB = BTWyB)
                 #BTWyBinv = BTWyBinv)

# Set via variable as in emulator case
#fit = stan(file = "source/full_field_calibration_fixed_em.stan",
#           data = stan_data,
#           iter = iter,
#           chains = chains,
#           #chains = 1,
#           #iter = 1,
#           model_name = "full_field_calibration")
fit = stan(file = "source/full_field_calibration_fixed_em_multi.stan",
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
#labels[3] = "Eccentricity bias"
#labels[4] = "log K_truss"
#labels[5] = "Rad thinning"
#labels[6] = "log K_spring"

# Plot observation error precision
lambda_hist(lambda_y, prior_shape = a_y, prior_rate = b_y, label = "lambda_y") #, adj_prior_shape = a_y_dash, adj_prior_rate = b_y_dash)

# Transform calibrated inputs onto their original scale for plotting and output
tf_trans = rescale_inputs(tf, t_min, t_max)
# Plot prior and posterior distribution of the calibration parameters.
calibration_inp_hist(tf_trans, tf_param = tf_param, inp_labels = labels, nrows=2, new_window = T)

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

# update prediction code to match multi_dof

#-------------------------------------------------------------------------------

# Make predictions from fitted Gaussian process emulator
# N_sam_plot = N_samples
N_sam_plot = 1000
# Transform correlation lengths into appropriate format
beta_w = -4.0*log(rho_w)

# Take averages across posterior predictions of the calibrated Gaussian process
out_list = full_field_calibration_pred_fixed_em(N_sam_plot, tc, tf, z_hat, beta_w, lambda_w, lambda_eta, lambda_y, KTKinv, BTWyB, K = K_eta, K_y = K_y, nugget = F, sam_gp = T, output_coeff_sam = T, output_coeff_mean = F, output_ff_mean = T, output_ff_std = T)
# It is too expensive to store a large number of posterior samples. Instead 
# output a fewer number of samples and append to output list
N_sam_plot = 25
out_list = c(out_list, full_field_calibration_pred_fixed_em(N_sam_plot, tc, tf, z_hat, beta_w, lambda_w, lambda_eta, lambda_y, KTKinv, BTWyB, K = K_eta, K_y = K_y, nugget = F, sam_gp = T, output_ff_sam = T, output_coeff_mean = F, output_ff_mean = F))

# On next iteration, tidy up prediction code as follows:
# First package up all the w_star prediction into one function then do loop first
# automatically output all samples as there will always be enough memory for this,
# though could optionally output them later
# this would easily cut down on am ount of code because inverse transformation could
# also be packaged, thus reducing duplication of K_y, but wouldn't have to do
# the matrix solves twice so not more expensive. Do on next iteration
# After this I think it'll probably be possible to tidy further. Wait till
# I've progressed the calibration work to overcome some of the problems first

# Loop over each output then append to the list to be written to json
json_list <- list()
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
  # Predictions at experimental data points to be outputted separately as csvs
  if (!grepl("w_star", out_string, fixed = TRUE) & !grepl("y", out_string, fixed = TRUE)) {
    # Append transformed values to list
    out_i = list(out_i)
    names(out_i) = out_string
    json_list = append(json_list, out_i)
  } else if (out_string != "w_star_sigma") {
    print(out_string)
    # We don't want to output covariance matrices.
    # Output predictions at experimental coordinates to csv
    write_output(out_i, out_string, in_file)
  }
}

# Write json to file
out_json = gp_pred_to_json(json_list, n_frames, n_nodes, disp_str, n_post_sam = N_sam_plot)
# write(out_json, paste("outputs/gp_predictions_nonlinear_",in_file,".json",sep=""))
write(out_json, paste("E:Working_Folder/gp_predictions_nonlinear_",in_file,".json",sep=""))

# Now plot histogram of calibrated expansion coefficients
w_star_hist(out_list$w_star, w_lim = c(-3,3))
