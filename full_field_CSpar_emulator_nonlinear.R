
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
library(rjson)

# Set current working directory. This should be modified to match the directory
# of the user
setwd("C:/Users/cs2361/Documents/CSpar_Calibration/")

# include functions which are called in this code
source("source/estimate_mode.R")
source("source/covariance_matrices.R")
source("source/dimension_reduction.R")
source("source/abaqus_json.R")

#-------------------------------------------------------------------------------

# Set up parameters which govern the formulation

p_eta = 11 # Number of basis functions to be retained for the emulator from SVD
exp_tol = 1e-6 # Tolerance variance fraction used to assess SVD convergence
disp_str = "w" # String which identifies the displacement component of interest (u,v, or w)
print_output = TRUE # Print diagnostic output to the terminal?
# Define parameters of the gamma prior on the error associated with truncating
# the series expansion for the model output
a_eta = 1.0     # Shape parameter for the lambda_eta prior
b_eta = 0.0001  # Rate parameter for the lambda_eta prior

#-------------------------------------------------------------------------------

# Set up simulation data

# Load in emulator training data input values from Design of Experiments. 
in_file = "LHSDesign40x4" # File identifier string which is in both input and output files
XT_sim = fread(paste("inputs/",in_file,".csv", sep = ""))

# Determine useful quantities from model inputs and outputs. Variable names 
# match the notation of Higdon et al. 2008
q = ncol(XT_sim)          # number of calibration inputs
tc = as.matrix(XT_sim)    # Convert training data input points to a matrix for passing to stan
m = nrow(XT_sim)          # sample size of computer simulation data

# Load in  training data output displacement values from Abaqus. Here I've used
# a similar naming convention to the inputs to automate changes. The file name 
# can be changed manually if need be.
# Here I've stored outputs structured in a json across samples and load increments
abaqus_displacements <- fromJSON(file = paste("inputs/",in_file,"_output_struct.json", sep=""))
# DataFrame might be a more intuitive format for R...
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
dt_all_cen = sweep(dt_simulation,1,mu_dt,"-")

# Divide by the standard deviation of the outputs
sd_dt = sd(as.matrix(dt_all_cen))
dt_all_cen = dt_all_cen/sd_dt
# Convert to matrix for stan
eta = as.matrix(dt_all_cen)

# Perform SVD on the centred data
dt_svd = svd(dt_all_cen)

# Extract the first p_eta basis functions from the svd, and standardise so the
# coefficients (columns of sqrt(m-1)*v) have unit variance
K_eta = dt_svd$u[,1:p_eta]*matrix(dt_svd$d[1:p_eta],nrow = n_eta, ncol = p_eta, byrow = TRUE)/sqrt(m-1)

if (print_output){
  # Sanity check that weights of SVD have zero mean and unit variance
  print("mean of reduced dimension output w = ")
  print(colMeans(dt_svd$v))
  print("standard deviations of reduced dimension output w = ")
  print(colSds(dt_svd$v*sqrt(m-1)))
}


# Write a JSON file with the output (still tempted to use DataFrame as the code
# is a little convoluted) 
out_json = basis_mean_to_json(n_frames, n_nodes, K_eta, mu_dt)
write(out_json, paste("outputs/basis_nonlinear_",in_file,".json", sep=""))


# Plot magnitude of singular value d with increasing number of basis functions
# to indicate how much each base contributes to the output variance.
dr_sqr = sum(dt_svd$d^2)
d_r = dt_svd$d[1:p_eta]
d_r_norm = d_r^2/dr_sqr
par(mfrow = c(1,2))
plot(1:p_eta,d_r_norm,"type"="p","col"="red","pch"=4,"lwd"=3,cex=1.5,
     'xlab' = "Feature",'ylab'="normalised d_i",cex.axis=1.75,cex.lab=1.75)
# Omit first point for greater clarity on convergence
plot(2:p_eta,d_r_norm[-1],"type"="p","col"="red","pch"=4,"lwd"=3,cex=1.5,
     'xlab' = "Feature",'ylab'="normalised d_i",cex.axis=1.75,cex.lab=1.75)
title("Convergence with Number of Basis Functions", line = -2, outer = TRUE, cex.main=1.75)

# How many functions are needed achieve a tolerance fraction of the variance?
basis_tol = 1:p_eta
basis_tol = basis_tol[cumsum(d_r_norm) > (1-exp_tol)]
basis_tol = basis_tol[1]
if (print_output){
  print("Number of basis functions required to represent output within tolerance = ")
  print(basis_tol)
}

#-------------------------------------------------------------------------------

# Reduce the dimension of the output data and determine associated quantities
# for input to Stan
processed_data = reduce_dimension(eta, K_eta, a_eta, b_eta)
# Adjusted prior parameters for the basis expansion truncation error
a_eta_dash = processed_data[[1]]
b_eta_dash = processed_data[[2]]
z_hat = processed_data[[3]] # Reduced-dimensional outputs
KTK_inv = processed_data[[4]] # Inverse of the product of basis matrices

#-------------------------------------------------------------------------------