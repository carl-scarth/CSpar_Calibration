
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
# source("source/utils.R")
# source("source/dimension_reduction.R")
# source("source/prior_posterior_plots.R")
# source("source/gp_predictions.R")

#-------------------------------------------------------------------------------

# Set up parameters which govern the formulation
# Might be able to delete some of these later
p_eta = 7 # Number of basis functions retained for the emulator from SVD
# exp_tol = 1e-6 # Tolerance variance fraction used to assess SVD convergence
disp_str = "w" # String which identifies the displacement component of interest (u,v, or w)
# Define parameters of the gamma prior on the error associated with truncating
# the series expansion for the model output
# a_eta = 1.0     # Shape parameter for the lambda_eta prior
# b_eta = 0.0001  # Rate parameter for the lambda_eta prior 
iter = 4000 # Number of samples per chain
chains = 3 # Number of chains for simulation
# print_svd_output = TRUE # Print diagnostic output of svd to the terminal?

p_eta = 5 # Number of basis functions to be retained for the emulator from SVD

#-------------------------------------------------------------------------------

# Define prior distribution parameters for passing to stan

# Pre-processing for BC example (Mean and coefficients of variation for Gaussian
# inputs)
E11_mu = 115.6
t_ply_mu = 0.196
E11_cov = 6.0 
t_ply_cov = 5.0
# Bounds for log-uniform inputs
K_lb = 100.0
K_ub = 1.0e9

# Define data_frame of prior parameters (this could be done via csv?)
tf_param = data.frame(E11 = c("Gaussian",E11_mu, E11_mu*E11_cov/100),
                      t_ply = c("Gaussian",t_ply_mu,t_ply_mu*t_ply_cov/100),
                      K = c("Loguniform",log(K_lb),log(K_ub)))

# Additional pre-processing for flange-rotation example
flange_theta_mu = 0.0
flange_theta_sigma = 4.0
# Define data_frame of prior parameters
tf_param = data.frame(E11 = c("Gaussian",E11_mu, E11_mu*E11_cov/100),
                      t_ply = c("Gaussian",t_ply_mu,t_ply_mu*t_ply_cov/100),
                      LFlange_theta = c("Gaussian",flange_theta_mu,flange_theta_sigma),
                      RFlange_theta = c("Gaussian",flange_theta_mu,flange_theta_sigma))

#-------------------------------------------------------------------------------
# Set up simulation data

# Load in emulator training data (input values) from Design of Experiments. 
# This input includes different values across uncontrolled calibration inputs.
XT_sim = fread("inputs/LHSDesign30x3.csv")
# take natural logarithm of spring stiffness
XT_sim$K = log(XT_sim$K)
colnames(XT_sim)[3] = "log_K"

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
