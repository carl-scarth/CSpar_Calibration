
# Full-field calibration method applied to C-spar with no calibration

# Set up R

library(data.table)
library(rstan)
#library(maximin)
library(matrixStats)
#library(colormap)
library(MASS)

# Set current working directory. This should be modified to match the directory
# at which the stan code and any data is stored
setwd("C:/Users/cs2361/Documents/CSpar_Calibration/")

# include functions which are called in this code
source("source/estimate_mode.R")
source("source/covariance_matrices.R")

#-------------------------------------------------------------------------------

# Set up parameters which govern the formulation

p_eta = 1 # Number of basis functions to be retained for the emulator from SVD

# define vector of prior means for the calibration inputs
E11_mu = 115.6
# E22_mu = 9.24
# nu12_mu = 0.335
# nu23_mu = 0.487
# G12_mu = 4.826
t_ply_mu = 0.196

# define prior coefficients of variation for the calibration parameters
E11_cov = 6.0 
# E22_cov = 6.0 
# nu12_cov = 12.123 
# nu23_cov = 12.0
# G12_cov = 6.0
t_ply_cov = 5.0

# Rotational spring stiffness is defined by bounds
K_lb = 100.0
K_ub = 1.0e9
  
# Combine prior means and coefficients of variation into a single vector.
# tf_mu = c(E11_mu,E22_mu,nu12_mu,nu23_mu,G12_mu,t_ply_mu)
# tf_cov = c(E11_cov,E22_cov,nu12_cov,nu23_cov,G12_cov,t_ply_cov)
# Calculate prior standard deviations from mean and COV
#tf_sigma = tf_mu*tf_cov/100

# Combine prior distribution parameters into a single vector
# take natural logarithm of spring stiffness, as this is unformly-distributed
tf_param_1 = c(E11_mu, t_ply_mu, log(K_lb))
tf_param_2 = c(E11_mu*E11_cov/100, t_ply_mu*t_ply_cov/100, log(K_ub))

#-------------------------------------------------------------------------------

# Set up simulation data

# Load in emulator training data (input values) from Design of Experiments. 
# This input includes different values across uncontrolled calibration inputs.
XT_sim = fread("inputs/LHSDesign30x3.csv")
# take natural logarithm of spring stiffness
XT_sim$K = log(XT_sim$K)
colnames(XT_sim)[3] = "log_K"

# Determine useful quantities from model inputs and outputs. Variable names 
# match the notation of Higdon et al.
q = ncol(XT_sim)          # number of calibration inputs
tc = as.matrix(XT_sim)    # Convert training data input points to a matrix for passing to stan
m = nrow(XT_sim)          # sample size of computer simulation data

# Load in emulator training data (outputs - Abaqus nodal displacement data).
# Each row matches inputs for the corresponding row in XT_sim
displacement_data = fread("inputs/LHSDesign30x3_displacements.csv")#, header = FALSE, sep = ",")
n_eta = nrow(displacement_data) # number of output points per simulation

dt_all_simulation = matrix(NA,nrow = n_eta, ncol = m)
# Extract the axial displacement, w
for (i in 1:m){
  dt_all_simulation[,i] = displacement_data[[as.name(sprintf('w_%d', i))]]
}

#mins = matrix(colMins(as.matrix(dt_all_simulation)),nrow = 3,ncol = 60)
# maxs = matrix(colMaxs(as.matrix(dt_all_simulation)),nrow = 3,ncol = 60)
#dt_all_simulation = matrix(as.matrix(dt_all_simulation),nrow = n_eta*3,ncol = m)


# For the time being just work with the axial (z) displacement
# think a bit about how to process full displacement vector (might require
# alterations to the statistical model)
#coord_ind = 3 # Coordinate of interest
#dt_all_simulation = dt_all_simulation[((coord_ind-1)*n_eta+1):(coord_ind*n_eta),]

#-------------------------------------------------------------------------------

# All training data points are normalised onto the unit hypercube [0,1]^q before
# being passed to stan. This code transforms all of these points, also 
# transforming the prior mean and standard deviation of the calibration 
# parameters for the sake of consistency

# Determine the maximum and minimum value of each input within the training data
t_min = t(as.matrix(colMins(tc)))
t_max = t(as.matrix(colMaxs(tc)))
# Normalise the training data points using these maximum and minimum values
tc = (tc - t_min[rep(1,m),])/(t_max[rep(1,m),]-t_min[rep(1,m),])
# also transform the prior mean and standard deviations using the same values for consistency
tf_param_1 = as.vector((tf_param_1 - t_min)/(t_max-t_min))
tf_param_2[3] = tf_param_2[3] - t_min[3]
tf_param_2 = as.vector(tf_param_2/(t_max-t_min))

#-------------------------------------------------------------------------------

# plot Design of Experiments used to generate training data
# pairs(~E11 + E22 + nu12 + nu23 + G12 + t_ply, data = tc,col = "blue",pch=4)
pairs(~E11 + t_ply + log_K, data = tc,col = "blue",pch=4)

#-------------------------------------------------------------------------------

# This portion of code performs SVD on the model data, which is standardised 
# such that the weight coefficients of the resulting expansion will have zero 
# mean and unit variance

# Centre the simulation output for each element. This guarantees that the weights (w)
# will have zero (sample) mean. Note that this same transformation must also be 
# applied to experimental data y, likely requiring some interpolation. 
# An alternative approach would be to centre the data using the overall mean, 
# calculating the sample mean of the transformed output, then specify this 
# (constant) mean value to each Gaussian process. This is arguably easier and 
# might reduce this possibility of errors being introduced.
mu_dt = rowMeans(dt_all_simulation)
dt_all_cen = sweep(dt_all_simulation,1,mu_dt,"-")

# Divide by the overall standard deviation of the simulation data to standardise
# to unit variance
sd_dt = sd(as.matrix(dt_all_cen))
dt_all_cen = dt_all_cen/sd_dt
# Store standardised simulation output as eta and convert to matrix to pass to stan
eta = as.matrix(dt_all_cen)

# Perform SVD on centred data
dt_svd = svd(dt_all_cen)

#------------------------------------------------------------------------------------------------------------------

# Extract the first p_eta basis functions from the svd. 
# standardised such that the weights (columns of sqrt(m-1)*v) have unit variance
# K_eta = dt_svd$u[,1:p_eta]*matrix(rep(dt_svd$d[1:p_eta],2),nrow = n_eta, ncol = p_eta, byrow = TRUE)/sqrt(m-1)
K_eta = dt_svd$u[,1:p_eta]*matrix(dt_svd$d[1:p_eta],nrow = n_eta, ncol = p_eta, byrow = TRUE)/sqrt(m-1)
# K_eta = as.matrix(dt_svd$u[,1:p_eta]/sqrt(m-1)) # This one has worked the best so far
# K_eta = as.matrix(dt_svd$u[,1:p_eta]) # Also quite good


# Sanity check that weights of SVD have zero mean and unit variance
print("mean of reduced dimension output w = ")
print(colMeans(dt_svd$v))
print("standard deviations of reduced dimension output w = ")
print(colSds(dt_svd$v*sqrt(m-1)))

#------------------------------------------------------------------------------------------------------------------------------------------

# Alternative normalisation - just to vhevk

## Determine samples of unknown coefficients of expansion for the training samples.This exploits orthogonality
## of the basis vectors, as multiplying the output by one of the vectors removes all contributions from the 
## other basis vectors, noting that t(u_r[,i]) %*% u_r[,j] = `1 if i = j, and 0 otherwise. Identify training
## samples of the coefficients for each basis vector (row), and each training sample (column)
#a_sim = t(as.matrix(K_eta)) %*% as.matrix(dt_all_cen)

#rowMeans(a_sim)
#a_sd = rowSds(a_sim)
# Scale the basis functions using the standard deviations, such that the weights
# have unit (sample) standard deviation
#K_eta = sweep(K_eta,2,a_sd,"/")
# Code for if only using p_eta = 1, delete/don't run the below line if that is the case
#K_eta = as.matrix(K_eta/a_sd)
# Sanity check
#a_sim_t = t(as.matrix(K_eta)) %*% as.matrix(dt_all_cen)
#rowSds(a_sim_t)

#-------------------------------------------------------------------------------

# Load in the experimental data, and centre using the same method as used to 
# standardise the model output. This requires some interpolation to determine the
# mean model output at each measurement location

# Read in DIC data. Each row has entries for which element the measurement has 
# been mapped to, and the natural coordinates for that point within the element
# experimental_data = read.table("inputs/Ext_LCorner_Image_0115_0.tiff_aligned_adjusted.csv", sep = ',', header = TRUE)
experimental_data = read.table("inputs/Ext_LCorner_Image_0115_0.tiff_nat_coord_rad_trim.csv", sep = ',', header = TRUE)
n_y = nrow(experimental_data)# Number of observations
y_element = experimental_data$Element
hr = as.matrix(experimental_data[,c("h","r")])
exp_displacement = experimental_data$W

# I also need to load mesh connectivity information about the outer surface of
# the spar. This must be consistent across all of the training data for this 
# method to work.
connectivity = read.table("inputs/outer_surface_elements.csv", sep = ',', header = TRUE)

# Define the position of the nodes in natural coordinates (this is 4,1,5,8) as
# this is what I used to determine the coordinates.
# It would probably make more sense to use standard ordering for a quad element,
# but I would need to re-do the mapping
HR = matrix(c(1.0, -1.0, -1.0, 1.0, -1.0, -1.0, 1.0, 1.0),nrow=2,ncol=4,byrow = TRUE)

# For consistency with simulation output, the DIC data must be centred
# using with the mean model output at their locations. As this is  unknown 
# precisely, this must be approximated through interpolation. Here I use the
# same linear interpolation functions of the isoparametric elements used by
# Abaqus
mu_y = rep(NA,n_y)
for (i in 1:n_y) {
  # Might be clearer just to write out basis equations in full... 
  bases = 1.0 + matrix(hr[i,],2,4)*HR
  # The complete basis functions for the quad element are given by the product
  # of those in h and r, contained in the rows of "bases"
  # Note that the element index is in Python indexing convention, but connectivities are in the Abaqus convention
  # plus 2 as the first two nodes are the reference points, which are numbered
  # according to a different system as they are defined directly on the assembly
  mu_y[i] = sum(bases[1,]*bases[2,]*mu_dt[as.numeric(connectivity[y_element[i]+1,])+2])/4.0
}

# Centre the experimental data using the interpolated mean model output
# output mean at data point, residual and relative error (with mean) across data points (consider other full-field metrics)
residual = abs(exp_displacement - mu_y)
rel_error = (residual/abs(mu_y))*100
write.csv(cbind(experimental_data[c("X","Y","Z")],mu_y,residual,rel_error), "outputs/mean_error.csv", row.names = FALSE)
exp_displacement_cen = (exp_displacement - mu_y)/sd_dt
# Store experimental data as y and convert to vector to pass to stan
y = as.vector(exp_displacement_cen)
 
#-------------------------------------------------------------------------------

# In this Section, W_y, the (prior) precision of the observation error is specified
# Here we give a value to the expected standard deviation of the observation error,
# which will later be weighted by parameter lambda_y, which will be given a 
# strong prior centred roughly around a value of 1. Here I assumed that error at
# each of the DIC observations is iid with standard deviation sigma_error. A 
# value of 0.01 has been chosen considering samples loaded at 0kN.

# In practice there will also be an error correlated across all observations due 
# to the fact that I had to define a zero datum for the DIC data. This can be
# added to every entry of the covariance matrix as sigma_shift, but if not
# careful this will dominate. I'm unsure of a precise value, but I suggest
# investigating differing values of 0, 0.0033, 0.01, and 0.025

# From the data, the DIC error appears correlated, with regions of higher noise
# appearing in clusters rather than as white noise. Ultimately, it would be good
# to represent these as random fields about which we infer the amplitude and
# correlation length parameters, but this requires a different statistical model

# Finally, there is an error which I do not account for, due to misalignment of
# the DIC data and model. This error would obviously be correlated, applying to
# all data points, but would not take the form of a uniform value. In some places
# such a shift would result in higher values than the "truth", in some cases
# lower. I'm unsure of how to account for this in practice.

# There will also be a "time" correlated error in-line with Higdon et al, but 
# there is no need to account for this here for one load value

# sigma_error = 0.01 # standard deviation associated with noise (a choice of 0.005 would also be reasonable if this is too large)
# sigma_shift = 0.0033 # standard deviation due to a shift in the zero value of the DIC data
# sigma_shift = 0
# Sigma_y = diag(rep(sigma_error^2,n_y)) + matrix(sigma_shift^2, nrow = n_y, ncol = n_y)

# Because the experimental data has been standardised by dividing through by the 
# standard deviation of the model output, the covariance of the measurement error
# must also be divided through by the variance of the model output for consistency
# Sigma_y = Sigma_y/(sd_dt^2)
# Convert covariance matrix to a precision by taking the inverse, as this is what
# is specified to stan as in Higdon et al.
# W_y = solve(Sigma_y) # There are almost certainly more efficient ways of implementing this...
W_y = diag(rep(1,n_y))

#-------------------------------------------------------------------------------

# The statistical model for the experimental data is y = Kw(theta) + Dv + e. For this 
# model it is necessary to evaluate the basis fuctions, K (as determined above 
# using SVD) at the physical locations of the experimental measurements in y. 
# To do this is is necessary do interpolate the basis functions, as was the case
# with the mean model output.

# Firstly interpolate emulator basis vectors, K.
K_y = matrix(NA,n_y,p_eta)
for (i in 1:n_y) {
  # Might be clearer just to write out basis equations in full... 
  bases = 1.0 + matrix(hr[i,],2,4)*HR
  # The complete basis functions for the quad element are given by the product
  # of those in h and r, contained in the rows of "bases"
  # Note that the element index is in Python indexing convention, but connectivities are in the Abaqus convention
  K_y[i,] = colSums(matrix(bases[1,]*bases[2,],4,p_eta)*K_eta[as.numeric(connectivity[y_element[i]+1,])+2,])/4.0
}

output_frame = cbind(experimental_data[c("X","Y","Z")],K_y)
for (i in 1:p_eta) {
  colnames(output_frame)[3+i] = sprintf("K_y,%d",i)
}
# output interpolated bases for plotting
write.csv(output_frame, "outputs/interpolated_basis.csv", row.names = FALSE)

#-------------------------------------------------------------------------------

# This portion of the code deals with the matrix algebra from Section 2.2.4 of 
# Higdon et al., calculating the necessary quantities for passing to stan, where
# the sampling is undertaken

# Assemble B matrix, which maps between the experimental data points,y, and the
# reduced-dimension coefficients, v (discrepancy) and u (emulator) corresponding 
# to these points. For n = 1 experiment this is simply the concatenation of D_y 
# and K_y
B = K_y

# Calculate B'*W_y*B
BTWyB = t(K_y)%*%W_y%*%K_y
BTWyBinv = solve(BTWyB)

# Calculate inverse of K'*K. This is used for multiple calculations and so it is
# more efficient to store the inverse than to solve the equations via other means
# Here K is the matrix of emulator basis functions evaluated for full model 
# output, arranged as specified at the end of Section 2.2.2 of Higdon et al.
# Note that this is a very large matrix (m x n_eta) by (m x p_eta), and as such
# calculating K'*K is computationally expensive, and requires a prohibitively 
# large amount of memory. To make this code possible I've calculated closed-form
# expressions for this matrix product, and inputted K'*K directly.

# Populate dense matrix of products of basis vectors, 
# with KTK_dense[i,j] = k_i'*k_j
# These product of individual basis vectors make up the entries of K'*K
KTK_dense = t(K_eta)%*%K_eta
# Note that the below code is the general expression for this matrix product if
# basis vectors are not orthogonal. As the basis vectors obtained by SVD are 
# orthogonal such that k_i'*k_j = 0 if i not equal to j, the code for populating
# off-diagonal terms has been commented out. This may be adapted to non-
# orthogonal bases by un-commenting the nested for loop.
KTK = matrix(0,m*p_eta,m*p_eta)
for (i in 1:p_eta){
  # Populate diagonal terms
  KTK[((i-1)*m+1):(i*m),((i-1)*m+1):(i*m)] = diag(m)*KTK_dense[i,i]
  # Note, Need to use seq_len as colon operator in R does not assume the list of 
  # indices must always increase, and j in 1:(i+1) throws an error when i = p_eta
  # for (j in (i + seq_len(p_eta-i))){
  #  # Populate off-diagonal terms. Note that these will be zero if the ks are orthogonal
  #  KTK[((i-1)*m+1):(i*m),((j-1)*m+1):(j*m)] = diag(m)*KTK_dense[i,j]
  #  KTK[((j-1)*m+1):(j*m),((i-1)*m+1):(i*m)] = diag(m)*KTK_dense[i,j]
  #}
}
# Determine inverse of K'*K. Note that as this matrix is diagonal for
# orthogonal basis functions in such cases it would be more efficient to 
# directly input the inverse, by inputting the reciprocal of the entries in the 
# above loop.
KTKinv = solve(KTK)

# Calculate K'*eta, which is used to calculate OLS solution w_hat. This 
# calculation has also been implemented by directly inputting the closed-form
# expression for this product, due to the large size of K.
# Note that in Higdon et al, eta is reshaped from a n_eta x m matrix into a 
# (m x n_eta) x 1 vector, as eta = [eta_1',eta_2',...,eta_m']'.
# K'*eta therefore results in a stack of
# [k_1'*eta_1;...;k_1'*eta_m;k_2'*eta_1;k_2'*eta_m;...;k_p_eta'*eta_m]
# It is more efficient to do this product by keeping K_eta stored as a matrix
# and likewise retaining eta as a matrix (for the purposes of this comment let's
# call this matrix eta_mat), then calculating K'*eta by obtaining eta_mat'*K_eta
# reshaped column-wise into a vector, as below
KTeta = t(eta)%*%(K_eta)
KTeta = as.vector(KTeta)

# Determine pseudo-inverses of basis matrices to get reduced-dimensional 
# representation of model output and experimental data
w_hat = KTKinv%*%KTeta # reduced-dimensional representation of emulator training data
# Calculate reduced dimensional representation of discrepancy and emulator 
# corresponding to experimental data.
# Note. I have not implemented an efficient version of the below matrix algebra
# as the quantity of experimental data is small in this application. When using
# DIC data, I will likely have to implement more closed-form expressions for the
# matrix products as in the case of K.I've implemented the more efficient version
# for b_y_dash - just not the pseudo-inverse
BTWyy = t(B)%*%W_y%*%y
# vu_hat = BTWyBinv%*%t(B)%*%W_y%*%y 
u_hat = BTWyBinv%*%BTWyy

# Group together all samples from experimental data and model output into a 
# single vector
z_hat = as.vector(rbind(u_hat,w_hat))

# Define parameters for (gamma) prior distributions on the emulator error and 
# observation error terms
a_eta = 1.0     # Shape parameter for the lambda_eta prior
b_eta = 0.0001  # Rate parameter for the lambda_eta prior 
a_y = 5.0       # Shape parameter for the lambda_y prior
# b_y = 5.0       # Rate parameter for the lambda_y prior
b_y = 0.05

a_eta_dash = a_eta+(0.5*(m*(n_eta-p_eta))) # Adjusted shape parameter for the lambda_eta prior (Eq. 11 Higdon et al.)
a_y_dash = a_y+(0.5*(n_y-p_eta)) # Adjusted shape parameter for the lambda_y prior.

# Stack reshaped model output eta into a single vector. This is the format used 
# in Higdon et al., although it wasn't necessary for the above calculations it
# is used to adjust the prior parameters.
eta_vec = as.vector(eta)
# Re-arraged version of b_eta_dash from that in Eq. (11) for the sake of 
# computational efficiency, using the fact that:
# eta'*(I - K*(K'*K)^-1*K')*eta = eta'*eta - (K'*eta)'*w_hat
b_eta_dash = as.numeric(b_eta + (0.5*(t(eta_vec)%*%eta_vec - t(KTeta) %*% w_hat)))
#b_y_dash = as.numeric(b_y + (0.5*(t(y)%*%(W_y - W_y%*%B%*%BTWyBinv%*%t(B)%*%W_y)%*%y)))
b_y_dash = as.numeric(b_y + (0.5*(t(y)%*%W_y%*%y - t(BTWyy) %*% u_hat)))

#-------------------------------------------------------------------------------


# This segment of code deals with setting up the environment for stan, passing 
# data to stan, and running the correct stan code depending upon which method is
# required for sampling

# Set up the environment for the Stan model to run in parallel. Taken from:
# https://betanalpha.github.io/assets/case_studies/gaussian_processes.html#21_Simulating_From_A_Gaussian_Process
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
parallel:::setDefaultClusterOptions(setup_strategy = "sequential")
util = new.env()
par(family="CMU Serif", las=1, bty="l", cex.axis=1, cex.lab=1, cex.main=1,
    xaxs="i", yaxs="i", mar = c(5, 5, 3, 5))

# Set up stan code
stan_data = list(m=m, q=q, n_eta=n_eta, n_y=n_y, p_eta=p_eta, a_eta_dash = a_eta_dash, a_y_dash = a_y_dash, b_eta_dash = b_eta_dash, b_y_dash = b_y_dash, z_hat = z_hat, tf_param_1=tf_param_1, tf_param_2=tf_param_2, tc = tc, KTKinv = KTKinv, BTWyBinv = BTWyBinv)
fit = stan(file = "source/Full_Field_Calibration_Higdon_Prior_Spec_No_Disc.stan",
           data = stan_data,
           iter = 4000,
           chains = 3)
           # iter = 1,
           # chains = 1)


# Good up to here ^^^ Run just one to check

#-------------------------------------------------------------------------------

# This section of code deals with post-processing of the data coming out of the 
# simulations

# plot trace plots of simulation samples
# dev.new(noRStudioGD = TRUE)  # generate plots in separate window
stan_trace(fit, pars = c("tf", "rho_w", "lambda_w","lambda_y","lambda_eta"))

# Summarise results to check convergence
print(fit, pars = c("tf", "rho_w", "lambda_w","lambda_y","lambda_eta"))

# extract samples from stan output

samples <- extract(fit)
N_samples = dim(samples$rho_w)[1] # get total number of samples

# Extract label of inputs for plots
labels = colnames(XT_sim)

# Produce plots of posterior and prior distributions of correlation parameters
# (rho) for emulator. Here the rows corresponding to the different principal
# components, whereas the columns correspond to the different calibration inputs
# A values of rho close to 1 implies that an input does strongly affect the 
# model output
dev.new(noRStudioGD = TRUE)
par(mfrow = c(p_eta,q))
rho_plot = seq(0,1, length.out = 100)
for (i in 1:p_eta){
  for (j in 1:q){
    # Plot histogram of posterior distribution
    hist(samples$rho_w[,(i-1)*q+j],
         main = labels[j],
         #main = paste("rho_w,",as.character(i),",",as.character(j)),
         xlab = paste("rho_w,",as.character(i),",",as.character(j)),
         col = "firebrick1",
         breaks = 25,
         freq = FALSE,
         xlim = c(0,1), 
         cex.lab=1.5,
         cex.axis=1.5)
    # Plot prior distribution
    prior_plot = dbeta(rho_plot,shape1=1,shape2=0.1)
    lines(rho_plot,prior_plot,lwd=3,col="blue")
  }
}

# Plot precision parameters for emulator. Each plot corresponds to a different
# principal component
dev.new(noRStudioGD = TRUE)
par(mfrow = c(1,p_eta))
lambda_plot = seq(0,2.5, length.out = 100)
for (i in 1:p_eta) {
  # plot posterior
  hist(samples$lambda_w[,i],
       main = paste("lambda_w",as.character(i)),
       xlab = paste("lambda_w",as.character(i)),
       col = "firebrick1",
       breaks = 25,
       freq = FALSE,
       # xlim = c(0,2.5),
       cex.axis=1.5,
       cex.lab=1.5)
  # plot prior
  prior_plot = dgamma(lambda_plot,shape=5,rate=5)
  lines(lambda_plot,prior_plot,lwd=3,col="blue")
}

# Plot precision parameters for the discrepancy. Each plot corresponds to a 
# different basis function
# dev.new(noRStudioGD = TRUE)
# par(mfrow = c(1,p_delta))
# lambda_plot = seq(0,1E5, length.out = 100)
# for (i in 1:p_delta){
#   # plot posterior
#  hist(samples$lambda_v[,i],
#       main = paste("lambda_v",as.character(i)),
#       xlab = paste("lambda_v",as.character(i)),
#       col = "firebrick1",
#       breaks = 25,
#       freq = FALSE,
#       #xlim = c(0,75),
#       cex.axis=1.5,
#       cex.lab=1.5)
#  # plot prior
#  prior_plot = dgamma(lambda_plot,shape=1.0,rate=0.0001)
#  lines(lambda_plot,prior_plot,lwd=3,col="blue")
#}

# Plot precision of the observation error and PCA truncation error
lambda_plot = seq(0,10.0, length.out = 100)
dev.new(noRStudioGD = TRUE)
par(mfrow = c(1,2))
#plot posterior
hist(samples$lambda_y,
     main = "lambda_y",
     xlab = "lambda_y",
     col = "firebrick1",
     breaks = 25,
     freq = FALSE,
#     xlim = c(0,10.0),
     cex.axis=1.5,
     cex.lab=1.5)
# plot prior
prior_plot = dgamma(lambda_plot,shape=a_y,rate=b_y)
lines(lambda_plot,prior_plot,lwd=3,col="blue")
# plot of prior adjusted for equivalent, reduced-dimension normal-gamma model
adj_prior_plot = dgamma(lambda_plot,shape=a_y_dash,rate=b_y_dash)
lines(lambda_plot,adj_prior_plot,lwd=3,col="green")

hist(samples$lambda_eta,
     main = "lambda_eta",
     xlab = "lambda_eta",
     col = "firebrick1",
     breaks = 25,
     freq = FALSE,
     xlim = c(0,5e6),
     cex.axis=1.5,
     cex.lab=1,5)
# plot prior
lambda_plot = seq(0,5E6, length.out = 1000)
prior_plot = dgamma(lambda_plot,shape=a_eta,rate=b_eta)
lines(lambda_plot,prior_plot,lwd=3,col="blue")
adj_prior_plot = dgamma(lambda_plot,shape=a_eta_dash,rate=b_eta_dash)
lines(lambda_plot,adj_prior_plot,lwd=3,col="green")

# Before plotting their distributions, the calibration parameters must be 
# transformed back onto their original scale, recalling that they were initially
# transformed onto [0,1]
tf_trans = samples$tf
tf_trans = tf_trans*matrix(rep(t_max-t_min,N_samples),ncol=q,byrow=TRUE) + matrix(rep(t_min,N_samples),ncol=q,byrow=TRUE)

# Plot prior and posterior distribution of the calibration parameters.
# Note that the calibration parameters were initially defined as deviations from
# their nominal values in this application. Here I've added these nominal values
# back onto the prior and posterior plots, such that they are more readable by 
# engineers
dev.new(noRStudioGD = TRUE) # plot in new window
par(mfrow = c(1,q))
# plot histogram of posterior
for (i in 1:q){
  hist(tf_trans[,i], 
       main = labels[i],
       xlab = labels[i],
       col = "brown1",
       breaks = 25,
       freq = FALSE,
       xlim = c(t_min[i],t_max[i]),
       cex.lab = 1.25,
       cex.axis = 1.25)
  # overlay plot of prior distribution
  t_plot = seq(t_min[i],t_max[i], length.out = 100)
  if (i < q){
    print(i)
    prior_plot = dnorm(t_plot,mean = tf_param_1[i]*(t_max[i]-t_min[i]) + t_min[i], sd = tf_param_2[i]*(t_max[i]-t_min[i]))
    } else {
    prior_plot = dunif(t_plot, min = (tf_param_1[i]*(t_max[i]-t_min[i]) + t_min[i]), max = (tf_param_2[i]*(t_max[i]-t_min[i]) + t_min[i]))
  }
  lines(t_plot,prior_plot,lwd=3,"col"="blue")
}


# estimate mode of the posterior distribution
modes = rep(0,q)
for (i in 1:q){
  modes[i] = estimate_mode(tf_trans[,i])
}
print("calibration parameter modes = ")
print(modes)

print("calibration parameter means = ")
print(colMeans(tf_trans))

#-------------------------------------------------------------------------------


# Code for making predictions - here I'm assuming we're using the transformed 
# version for matrices which are not full rank

# N_sam_plot = 1000
N_sam_plot = 1
# Pick N_sam_plot random samples without repetition
rand_ind = sample.int(N_samples,N_sam_plot)
#N_sam_plot = N_samples
#rand_ind = 1:N_samples


v_star_mu = matrix(0, p_delta, N_sam_plot)
w_star_mu = matrix(0, p_eta, N_sam_plot)
v_star_sigma = array(0, c(p_delta, p_delta, N_sam_plot))
w_star_sigma = array(0, c(p_eta, p_eta, N_sam_plot))

v_star = matrix(0, p_delta, N_sam_plot)
w_star = matrix(0, p_eta, N_sam_plot)

delta_sam = matrix(0, n_eta, N_sam_plot)
delta_y = matrix(0, n_y, N_sam_plot)
eta_sam = matrix(0, n_eta, N_sam_plot)
eta_y = matrix(0, n_y, N_sam_plot)


sam_ind = 1
# Sort out predictions for full-rank code from here!!
for (i in 1:N_sam_plot){
  sam_ind = rand_ind[i]
  # Construct covariance matrix for the current sample
  
  # I could probably package the below code into an external function
  
  # Extract the inferred parameters for the current sample
  beta_w_i = samples$beta_w[sam_ind,]
  #lambda_v_i = samples$lambda_v[sam_ind,]
  # Code for if there are only one set of discrepancy parameters to learn about
  # Probably need some kind of conditional statement
  lambda_v_i = samples$lambda_v[sam_ind]
  
  lambda_w_i = samples$lambda_w[sam_ind,]
  tf_i = samples$tf[sam_ind, ]
  lambda_y_i = samples$lambda_y[i]
  lambda_eta_i = samples$lambda_eta[i]
  
  # THE ABOVE SHOULD WORK FOR ALL PREDICTION METHODS CHECK
  
  # Firstly we need to reconstruct the covariance matrix for the training data points,
  # Sigma_w_hat
  
  # Covariance matrix discrepancy weights.Here I assume that each weight is an
  # independent Gaussian process (i.e. F = p_delta, |G_i| = 1 for all i = 1,...p_delta)
  # Code for set of discrepancy parameters per basis functio
  # sigma_v = diag(1/lambda_v_i)
  # Code for if there is only one set of discrepancy parameters
  sigma_v = diag(rep(1/lambda_v_i,p_delta))
  
  # Define the covariance matrix for the emulator weights, evaluated at the
  # experimental data points. As we only have one experiment, and
  # R(theta,theta) = 1, this also reduces to a diagonal matrix.
  sigma_u = diag(1/lambda_w_i)

  # Combine contributions from sigma_u and sigma_v into a single block-diagonal
  # matrix for subsequent transformation onto their full-rank equivalent combining
  # terms from both
  sigma_vu_join = matrix(0, p_eta+p_delta, p_eta+p_delta)
  sigma_vu_join[1:p_delta, 1:p_delta] = sigma_v
  sigma_vu_join[-(1:p_delta), -(1:p_delta)] = sigma_u
  
  # Define the covariance matrix for the emulator weights, evaluated at the
  # training data points.
  sigma_w = matrix(0, m*p_eta, m*p_eta)
  
  for (j in 1:p_eta) {
    # Calculate the covariance matrix for the ith emulator weight
    sigma_w[((j-1)*m+1):(j*m), ((j-1)*m+1):(j*m)] = ARD_SE_cov(tc, lambda_w_i[j], beta_w_i[((j-1)*q+1):(j*q)], 0)
  }
  
  # Evaluate cross-covariance of emulator weights between training data and experimental data points
  sigma_uw = matrix(0, p_eta, m*p_eta)
  for (j in 1:p_eta) {
    sigma_uw[j,((j-1)*m+1):(j*m)] = ARD_SE_cov_non_sym(t(as.matrix(tf_i)), tc, lambda_w_i[j], beta_w_i[((j-1)*q+1):(j*q)])
  }
  
  # The covariance matrices used in predictions will differ depending on whether
  # or not B is full rank
  if ((rank_B < (p_eta+p_delta)) && exists("B_tilde")){
    # Transform contributions from sigma uv to into their full rank equivalent
    Lsigma_uw = L %*% rbind(matrix(0, p_delta, m*p_eta), sigma_uw)
    # Assemble covariance matrix for joint experimental and model data
    sigma_z = matrix(0,rank_B+m*p_eta, rank_B+m*p_eta)
    sigma_z[1:rank_B,1:rank_B] = L %*% sigma_vu_join %*% t(L)
    sigma_z[-(1:rank_B),-(1:rank_B)] = sigma_w
    sigma_z[1:rank_B,-(1:rank_B)] = Lsigma_uw
    sigma_z[-(1:rank_B),1:rank_B] = t(Lsigma_uw)
  } else {
    # Assemble covariance matrix for joint experimental and model data
    sigma_z = matrix(0,p_delta+(m+1)*p_eta, p_delta+(m+1)*p_eta)
    sigma_z[1:p_delta,1:p_delta] = sigma_v
    sigma_z[(p_delta+1):(p_delta+p_eta),(p_delta+1):(p_delta+p_eta)] = sigma_u
    sigma_z[-(1:(p_delta+p_eta)),-(1:(p_delta+p_eta))] = sigma_w
    sigma_z[(p_delta+1):(p_delta+p_eta),-(1:(p_delta+p_eta))] = sigma_uw
    sigma_z[-(1:(p_delta+p_eta)),(p_delta+1):(p_delta+p_eta)] = t(sigma_uw)
  }
  # Adjust the covariance matrix sigma_z to include transformed emulator and 
  # experimental error terms
  sigma_z_hat = matrix(0,rank_B+m*p_eta, rank_B+m*p_eta)
  sigma_z_hat[1:rank_B,1:rank_B] = BTWyBinv/lambda_y_i
  sigma_z_hat[-(1:rank_B),-(1:rank_B)] = KTKinv/lambda_eta_i
  sigma_z_hat = sigma_z_hat + sigma_z

  # Likewise covariance matrix with predictions will depend upon whether full-
  # rank format is used
  if ((rank_B < (p_eta+p_delta)) && exists("B_tilde")){
    # Define cross-correlation terms for z with discrepancy prediction, v_star. Note
    # That this does not need to be adjusted to as sigma_z_hat was, as the error 
    # terms do not apply to predictions
    sigma_z_v_star = rbind(sigma_v,matrix(0,p_eta,p_delta))
    # Need to multiply by L to obtain full-rank equivalent
    Lsigma_z_v_star = L %*% sigma_z_v_star
    # The resulting matrix is then  padded with zeros for correlation with w
    Lsigma_z_v_star = rbind(Lsigma_z_v_star,matrix(0,m*p_eta,p_delta))
  
    # Define cross-correlation of emulator predictions w_star with emulator training
    # data. Note that in this application, because we are making a prediction at the
    # same value of controlled input x, and uncontrolled inputs are set to their
    # calibrated values as was the case for the experimental data, then this is 
    # given by sigma_u_w', Note that this won't be the case if we want to make a
    # prediction at a different x, or want to evaluate just the emulator across a 
    # range of t values
    sigma_w_w_star = t(sigma_uw)

    # Define cross-correlation terms for z with emulator prediction w_star
    sigma_z_w_star = rbind(matrix(0,p_delta,p_eta), sigma_u)
    Lsigma_z_w_star = L %*% sigma_z_w_star
    # Combine with the correlation matrix for w to get the full correlation with z
    Lsigma_z_w_star = rbind(Lsigma_z_w_star,sigma_w_w_star)
  } else {
    # Define cross-correlation terms for z with discrepancy prediction, v_star. 
    # Note that this does not need to be adjusted to as sigma_z_hat was, as the 
    # error terms do not apply to predictions
    sigma_z_v_star = rbind(sigma_v,matrix(0,(m+1)*p_eta,p_delta))
    
    # Define cross-correlation of emulator predictions w_star with emulator training
    # data. Note that in this application, because we are making a prediction at the
    # same value of controlled input x, and uncontrolled inputs are set to their
    # calibrated values as was the case for the experimental data, then this is 
    # given by sigma_u_w', Note that this won't be the case if we want to make a
    # prediction at a different x, or want to evaluate just the emulator across a 
    # range of t values
    sigma_w_w_star = t(sigma_uw)
    # Define cross-correlation terms for z with emulator prediction w_star
    sigma_z_w_star = rbind(matrix(0,p_delta,p_eta), sigma_u, sigma_w_w_star)
  }
  
  # Define correlation of discrepancy prediction v_star with itself
  sigma_v_star = sigma_v
  
  # Define correlation of emulator prediction with itself
  sigma_w_star = sigma_u
  
  # Consider adding a nugget...
  
  # Two different methods for making predictions. The first method is quicker, but 
  # I think the second is more numerically stable, which seems to make a different
  # when the variance is small. Consider using this if the simulation is taking 
  # too long
  
  # Explicitly calculating inverse, then reusing for all predictions. 
  # Ainv = solve(sigma_z_hat)
  # v_mu = t(Lsigma_z_v_star) %*% Ainv %*% z_hat
  # w_mu = t(Lsigma_z_w_star) %*% Ainv %*% z_hat
  # v_sigma = sigma_v_star - (t(Lsigma_z_v_star) %*% Ainv %*% Lsigma_z_v_star)
  # w_sigma = sigma_w_star - (t(Lsigma_z_w_star) %*% Ainv %*% Lsigma_z_w_star)

  # Solving using solve
  Ainv_z_hat = solve(sigma_z_hat,z_hat)
  # Store mean and covariance matrices of discrepancy and adjusted prediction Gaussian processes
  if ((rank_B < (p_eta+p_delta)) && exists("B_tilde")){
    v_star_mu[,i] = t(Lsigma_z_v_star) %*% Ainv_z_hat
    w_star_mu[,i] = t(Lsigma_z_w_star) %*% Ainv_z_hat
    v_star_sigma[,,i] = sigma_v_star - (t(Lsigma_z_v_star) %*% solve(sigma_z_hat,Lsigma_z_v_star))
    w_star_sigma[,,i] = sigma_w_star - (t(Lsigma_z_w_star) %*% solve(sigma_z_hat,Lsigma_z_w_star))
  } else {
    v_star_mu[,i] = t(sigma_z_v_star) %*% Ainv_z_hat
    w_star_mu[,i] = t(sigma_z_w_star) %*% Ainv_z_hat
    v_star_sigma[,,i] = sigma_v_star - (t(sigma_z_v_star) %*% solve(sigma_z_hat,sigma_z_v_star))
    w_star_sigma[,,i] = sigma_w_star - (t(sigma_z_w_star) %*% solve(sigma_z_hat,sigma_z_w_star))
  }

  # Sample from Gaussian process (might not be necessary - could just cheat and 
  # look at the mean
  v_star[,i] = mvrnorm(n = 1, v_star_mu[,i], v_star_sigma[,,i])
  w_star[,i] = mvrnorm(n = 1, w_star_mu[,i], w_star_sigma[,,i])
  
  # Generate individual Samples
  delta_sam[,i] = (D_eta %*% v_star[,i])*sd_dt
  delta_y[,i] = (D_y %*% v_star[,i])*sd_dt
  eta_sam[,i] = (K_eta %*% w_star[,i])*sd_dt + mu_dt
  eta_y[,i] = (K_y %*% w_star[,i])*sd_dt + mu_y
}

# Good to here!!!

# Plot histogram of reduced coefficients
dev.new(noRStudioGD = TRUE) # plot in new window
par(mfrow = c(1, p_eta))
for (i in 1:p_eta){
  hist(w_star[i,], 
       main =  paste("Feature",as.character(i)),
       xlab = paste("w_star", as.character(i)),
       col = "brown1",
       breaks = 25,
       freq = FALSE,
       xlim = c(-4,4),
       cex.lab = 1.25,
       cex.axis = 1.25)
  # overlay plot of prior distribution
  w_plot = seq(-4,4, length.out = 100)
  prior_plot = dnorm(w_plot,mean = 0, sd = 1)
  lines(w_plot,prior_plot,lwd=3,"col"="blue")
}

dev.new(noRStudioGD = TRUE) # plot in new window
par(mfrow = c(1, p_delta))
for (i in 1:p_delta){
  hist(w_star[i,], 
       main =  paste("Feature",as.character(i)),
       xlab = paste("v_star", as.character(i)),
       col = "brown1",
       breaks = 25,
       freq = FALSE,
       xlim = c(-4,4),
       cex.lab = 1.25,
       cex.axis = 1.25)
  # overlay plot of prior distribution
  w_plot = seq(-4,4, length.out = 100)
  prior_plot = dnorm(w_plot,mean = 0, sd = 1)
  lines(w_plot,prior_plot,lwd=3,"col"="blue")
}

# Also look at samples of the mean for output
# Multiply stored mean values by the appropriate basis matrix and transform back
# onto their correct scales
delta_mu = (D_eta %*% v_star_mu)*sd_dt
delta_y_mu = (D_y %*% v_star_mu)*sd_dt
eta_mu = (K_eta %*% w_star_mu)*sd_dt + mu_dt
eta_y_mu = (K_y %*% w_star_mu)*sd_dt + mu_y

# Calculate the diagonal terms of the covariance matrix. Just do this as the
# full covariance is too big. Take the square root (standard deviation) as this
# is more meaningful
delta_sigma = matrix(0, n_eta, N_sam_plot)
delta_y_sigma = matrix(0, n_y, N_sam_plot)
eta_sigma = matrix(0, n_eta, N_sam_plot)
eta_y_sigma = matrix(0, n_y, N_sam_plot)
for (i in 1:N_sam_plot){
  delta_sigma[,i] = sqrt(as.matrix(rowSums((D_eta %*% v_star_sigma[,,i]) * D_eta)))*sd_dt
  delta_y_sigma[,i] = sqrt(as.matrix(rowSums((D_y %*% v_star_sigma[,,i]) * D_y)))*sd_dt
  eta_sigma[,i] = sqrt(as.matrix(rowSums((K_eta %*% w_star_sigma[,,i]) * K_eta)))*sd_dt
  eta_y_sigma[,i] = sqrt(as.matrix(rowSums((K_y %*% w_star_sigma[,,i]) * K_y)))*sd_dt
}

# Integrate uncertainty out of mean and standard deviation by taking average
eta_mu_mu = rowMeans(eta_mu)
delta_mu_mu = rowMeans(delta_mu)
eta_y_mu_mu = rowMeans(eta_y_mu)
delta_y_mu_mu = rowMeans(delta_y_mu)

eta_sigma_mu = rowMeans(eta_sigma)
delta_sigma_mu = rowMeans(delta_sigma)
eta_y_sigma_mu = rowMeans(eta_y_sigma)
delta_y_sigma_mu = rowMeans(delta_y_sigma)



# Also separate emulator

# Here it would be good to take some averages!! We haven't averaged out the calibration parameters!!



# Write all output to text files for plotting
write.csv(eta_mu, "outputs/eta_mu.csv", row.names = FALSE)
write.csv(delta_mu, "outputs/delta_mu.csv", row.names = FALSE)
write.csv(eta_y_mu, "outputs/eta_y_mu.csv", row.names = FALSE)
write.csv(delta_y_mu, "outputs/delta_y_mu.csv", row.names = FALSE)
write.csv(delta_sam, "outputs/delta_sam.csv", row.names = FALSE)
write.csv(eta_sam, "outputs/eta_sam.csv", row.names = FALSE)
write.csv(delta_y, "outputs/delta_y.csv", row.names = FALSE)
write.csv(eta_y, "outputs/eta_y.csv", row.names = FALSE)
write.csv(delta_sigma, "outputs/delta_sigma.csv", row.names = FALSE)
write.csv(delta_y_sigma, "outputs/delta_y_sigma.csv", row.names = FALSE)
write.csv(eta_sigma, "outputs/eta_sigma.csv", row.names = FALSE)
write.csv(eta_y_sigma, "outputs/eta_y_sigma.csv", row.names = FALSE)

write.csv(eta_mu_mu, "outputs/eta_mu_mu.csv", row.names = FALSE)
write.csv(delta_mu_mu, "outputs/delta_mu_mu.csv", row.names = FALSE)
write.csv(eta_y_mu_mu, "outputs/eta_y_mu_mu.csv", row.names = FALSE)
write.csv(delta_y_mu_mu, "outputs/delta_y_mu_mu.csv", row.names = FALSE)

write.csv(eta_sigma_mu, "outputs/eta_sigma_mu.csv", row.names = FALSE)
write.csv(delta_sigma_mu, "outputs/delta_sigma_mu.csv", row.names = FALSE)
write.csv(eta_y_sigma_mu, "outputs/eta_y_sigma_mu.csv", row.names = FALSE)
write.csv(delta_y_sigma_mu, "outputs/delta_y_sigma_mu.csv", row.names = FALSE)
