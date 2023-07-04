
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
source("source/estimate_mode.R")
source("source/covariance_matrices.R")

#-------------------------------------------------------------------------------

# Set up simulation data

# Load in emulator training data input values from Design of Experiments. 
in_file = "LHSDesign50x3" # File identifier string which is in both input and output csvs
XT_sim = fread(paste("inputs/",in_file,".csv", sep = ""))

# In this example I fit the emulator to the log of spring stiffness K, which is 
# a more natural choice of values across which outputs are expected for
# variations in this input
XT_sim$K = log(XT_sim$K)
colnames(XT_sim)[3] = "log_K"

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
abaqus_displacements = fread(paste("inputs/",in_file,"_displacements.csv", sep=""))
n_eta = nrow(abaqus_displacements) # number of output points per simulation

# Extract the displacement for the component of interest and store in a matrix
dt_simulation = matrix(NA,nrow = n_eta, ncol = m)
for (i in 1:m){
  dt_simulation[,i] = abaqus_displacements[[as.name(paste(disp_str,sprintf("_%d", i),sep=""))]]
}