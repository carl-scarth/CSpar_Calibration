
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
# library(MASS)
# library(colormap)

# Set current working directory. This should be modified to match the directory
# of the user
setwd("C:/Users/cs2361/Documents/CSpar_Calibration/")

# include functions which are called in this code
source("source/abaqus_json.R")
source("source/transform_input_output.R")
source("source/utils.R")
source("source/interpolate_data.R")
# source("source/dimension_reduction.R")
# source("source/prior_posterior_plots.R")
# source("source/gp_predictions.R")

#-------------------------------------------------------------------------------

# Set up parameters which govern the formulation
# Might be able to delete some of these later
# p_eta = 7 # Number of basis functions retained for the emulator from SVD
# exp_tol = 1e-6 # Tolerance variance fraction used to assess SVD convergence
disp_str = "w" # String which identifies the displacement component of interest (u,v, or w)
DIC_coord_labels = c("x_proj","y_proj","z_proj") # Strings used to identify coordinates in DIC point_cloud

# Define parameters of the gamma prior on the error associated with truncating
# the series expansion for the model output
# a_eta = 1.0     # Shape parameter for the lambda_eta prior
# b_eta = 0.0001  # Rate parameter for the lambda_eta prior 
iter = 4000 # Number of samples per chain
chains = 3 # Number of chains for simulation
exp_data_file = "interpolated_DIC" # Identifier of file containing DIC data
surface_elements = "nominal_shell_mesh_outer_surface_elements" # File identifier string for surface mesh connectivity

#-------------------------------------------------------------------------------

# Define prior distribution parameters for passing to stan

# Pre-processing for BC example (Mean and coefficients of variation for Gaussian
# inputs)
E11_mu = 115.6
t_ply_mu = 0.196
E11_cov = 6.0 
t_ply_cov = 5.0
# Bounds for log-uniform inputs
# K_lb = 100.0
# K_ub = 1.0e9

# Define data_frame of prior parameters (this could be done via csv?)
# tf_param <- data.frame(distribution = c("Gaussian","Gaussian","Loguniform"),
#                       param_1 = c(E11_mu, t_ply_mu, log(K_lb)),
#                       param_2 = c(E11_mu*E11_cov/100, t_ply_mu*t_ply_cov/100, log(K_ub)))
# row.names(tf_param) <- c("E11","t_ply","log_K")

# Additional pre-processing for flange-rotation example
flange_theta_mu = 0.0
flange_theta_sigma = 4.0
# Define data_frame of prior parameters
tf_param <- data.frame(distribution = c("Gaussian","Gaussian","Gaussian","Gaussian"),
                       param_1 = c(E11_mu, t_ply_mu, flange_theta_mu, flange_theta_mu ),
                       param_2 = c(E11_mu*E11_cov/100, t_ply_mu*t_ply_cov/100, flange_theta_sigma, flange_theta_sigma))
row.names(tf_param) <- c("E11","t_ply","LFlange_theta","RFlange_theta")

#-------------------------------------------------------------------------------
# Set up simulation data (need to retain this for normalisation of inputs)

# Load in emulator training data input values from Design of Experiments. 
# in_file = "LHSDesign50x3" # File identifier string for input and output csvs
in_file = "LHSDesign40x4" # File identifier string for input and output csvs
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
abaqus_displacements <- fromJSON(file = paste("inputs/",in_file,"_output_struct.json", sep=""))
# Extract displacements from json, and store as matrix where each column is a 
# training sample, and the rows are the concatenation of displacement output
# across all output frames
sorted_data = extract_const_frame(abaqus_displacements,"w")
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
K_eta = as.matrix(fread(paste("outputs/basis_nonlinear_",in_file,".csv", sep = "")))
p_eta = ncol(K_eta) # Number of basis functions retained for the emulator from SVD

#-------------------------------------------------------------------------------

# Load in the experimental data, and standardise using the same method as used  
# for the model output. This requires interpolation of the mean model output to 
# the DIC point cloud coordinates.
# Alongside displacements each row has entries for the element to which the 
# measurement has been matched, its natural coordinates within the element, and 
# applied load increment, which matches those of the model.
experimental_data = as.data.frame(fread(paste("inputs/", exp_data_file, ".csv", sep = "")))
n_y = nrow(experimental_data) # Total of observations
y_element = py_to_R(experimental_data$Element) # Matched element indices. Must be converted from Python to R convention
hr = as.matrix(experimental_data[c("h","r")]) # Matched natural coordinates
exp_displacement = experimental_data[[as.name(paste(disp_str, "_rot", sep=""))]] # Displacement component
frame_ind = experimental_data$Increment # Index of which frame each point belongs to
# Interpolate the model output across all increments
mu_y = intp_nodes_to_cloud_inc(y_element, hr, as.matrix(mu_dt), frame_ind, n_frames, conn_file = paste("inputs/", surface_elements, ".csv", sep=""), skip_nodes=2)

# Calculate residuals of experimental error
residual = exp_displacement - mu_y
rel_error = (abs(residual)/abs(mu_y))*100

# Write to CSV
write.csv(cbind(experimental_data[DIC_coord_labels],training_data_mean = mu_y,residual,abs(residual),rel_error, increment = experimental_data$Increment), "outputs/mean_error.csv", row.names = FALSE)

# Centre the experimental data and convert to vector to pass to stan
exp_displacement_cen = (exp_displacement - mu_y)/sd_dt
y = as.vector(exp_displacement_cen)

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------

# specify W_y, the (prior) precision of the observation.

# The commented code below deals with specifying both an iid noise error, and 
# one due to an isotropic shift applied to every point
# For now I pass the identity matrix, and place a weaker prior on the observation
# error in Stan
# In the long run it might be worth considering an error term which can change
# in magnitude across increments. Consider adding this later.
# sigma_error = 0.01 # standard deviation associated with noise (a choice of 0.005 would also be reasonable if this is too large)
# sigma_shift = 0
# Sigma_y = diag(rep(sigma_error^2,n_y)) + matrix(sigma_shift^2, nrow = n_y, ncol = n_y)
# Standardise using the model output variance 
# Sigma_y = Sigma_y/(sd_dt^2)
# Convert covariance matrix to a precision 
# W_y = solve(Sigma_y)

# Directly pass the identity matrix
W_y = diag(rep(1.0,n_y))