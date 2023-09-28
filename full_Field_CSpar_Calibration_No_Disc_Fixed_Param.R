
# Full-field calibration method applied to C-spar with no calibration

# Set up R

library(data.table)
library(rstan)
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

p_eta = 4 # Number of basis functions to be retained for the emulator from SVD

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

# Need to retain this for normalisation of parameter values

## Load in emulator training data (input values) from Design of Experiments. 
## This input includes different values across uncontrolled calibration inputs.
# XT_sim = fread("inputs/LHSDesign30x3.csv")
XT_sim = fread("inputs/LHSDesign50x3.csv")
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
# Retain for predictions
# displacement_data = fread("inputs/LHSDesign30x3_displacements.csv")#, header = FALSE, sep = ",")
displacement_data = fread("inputs/LHSDesign50x3_displacements.csv")#, header = FALSE, sep = ",")
n_eta = nrow(displacement_data) # number of output points per simulation

dt_all_simulation = matrix(NA,nrow = n_eta, ncol = m)
# Extract the axial displacement, w
for (i in 1:m){
  dt_all_simulation[,i] = displacement_data[[as.name(sprintf('w_%d', i))]]
  # dt_all_simulation[,i] = displacement_data[[as.name(sprintf('u_%d', i))]]
}


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

#------------------------------------------------------------------------------------------------------------------

# Extract Load in the basis vectors from file. 

# K_eta = as.matrix(read.table("outputs/basis_unit.csv", sep = ',', header = TRUE))
K_eta = as.matrix(read.table("outputs/basis_50x3.csv", sep = ',', header = TRUE))
# K_eta = as.matrix(read.table("outputs/basis_50x3_u.csv", sep = ',', header = TRUE))
K_eta = as.matrix(K_eta[,1:p_eta])

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

sigma_error = 0.01 # standard deviation associated with noise (a choice of 0.005 would also be reasonable if this is too large)
# sigma_shift = 0.0033 # standard deviation due to a shift in the zero value of the DIC data
sigma_shift = 0
Sigma_y = diag(rep(sigma_error^2,n_y)) + matrix(sigma_shift^2, nrow = n_y, ncol = n_y)

# Because the experimental data has been standardised by dividing through by the 
# standard deviation of the model output, the covariance of the measurement error
# must also be divided through by the variance of the model output for consistency
Sigma_y = Sigma_y/(sd_dt^2)
# Convert covariance matrix to a precision by taking the inverse, as this is what
# is specified to stan as in Higdon et al.
# W_y = solve(Sigma_y) # There are almost certainly more efficient ways of implementing this...
# W_y = diag(rep(0.01,n_y))
W_y = diag(rep(1.0,n_y))


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
  matrix(bases[1,]*bases[2,],4,p_eta)*K_eta[as.numeric(connectivity[y_element[i]+1,])+2,]
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
# a_eta = 1.0     # Shape parameter for the lambda_eta prior
# b_eta = 0.0001  # Rate parameter for the lambda_eta prior 
a_y = 5.0       # Shape parameter for the lambda_y prior
b_y = 5.0       # Rate parameter for the lambda_y prior
# b_y = 0.05
#b_y = 0.000000000000000005

# a_eta_dash = a_eta+(0.5*(m*(n_eta-p_eta))) # Adjusted shape parameter for the lambda_eta prior (Eq. 11 Higdon et al.)
a_y_dash = a_y+(0.5*(n_y-p_eta)) # Adjusted shape parameter for the lambda_y prior.
#
# Stack reshaped model output eta into a single vector. This is the format used 
# in Higdon et al., although it wasn't necessary for the above calculations it
# is used to adjust the prior parameters.
# eta_vec = as.vector(eta)
# Re-arraged version of b_eta_dash from that in Eq. (11) for the sake of 
# computational efficiency, using the fact that:
# eta'*(I - K*(K'*K)^-1*K')*eta = eta'*eta - (K'*eta)'*w_hat
# b_eta_dash = as.numeric(b_eta + (0.5*(t(eta_vec)%*%eta_vec - t(KTeta) %*% w_hat)))
#b_y_dash = as.numeric(b_y + (0.5*(t(y)%*%(W_y - W_y%*%B%*%BTWyBinv%*%t(B)%*%W_y)%*%y)))
b_y_dash = as.numeric(b_y + (0.5*(t(y)%*%W_y%*%y - t(BTWyy) %*% u_hat)))
a_y_dash = 5.0
b_y_dash = 0.05

#-------------------------------------------------------------------------------

# Load in emulator hyperparameters from previous fit

p_eta_em = 5
# emulator_parameters = c(as.matrix(read.table("outputs/emulator_modes_unit.csv", sep = ',', header = TRUE)))
emulator_parameters = c(as.matrix(read.table("outputs/emulator_modes_50x3.csv", sep = ',', header = TRUE)))
# emulator_parameters = c(as.matrix(read.table("outputs/emulator_modes_50x3_u.csv", sep = ',', header = TRUE)))
rho_w = emulator_parameters[1:(p_eta*q)]
lambda_w = emulator_parameters[(p_eta_em*q + 1):(p_eta_em*(q+1))]
lambda_w = as.vector(lambda_w[1:p_eta])
lambda_eta = emulator_parameters[p_eta_em*(q+1)+1]

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
stan_data = list(m=m, q=q, n_eta=n_eta, n_y=n_y, p_eta=p_eta, a_y_dash = a_y_dash, b_y_dash = b_y_dash, lambda_eta = lambda_eta, z_hat = z_hat, tf_param_1=tf_param_1, tf_param_2=tf_param_2, rho_w = rho_w, lambda_w = lambda_w, tc = tc, KTKinv = KTKinv, BTWyBinv = BTWyBinv)
fit = stan(file = "source/Full_Field_Calibration_Higdon_Prior_Spec_No_Disc_fixed_param.stan",
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
stan_trace(fit, pars = c("tf", "lambda_y"))

# Summarise results to check convergence
print(fit, pars = c("tf", "lambda_y"))

# extract samples from stan output

samples <- extract(fit)
N_samples = dim(samples$tf)[1] # get total number of samples

# Extract label of inputs for plots
labels = colnames(XT_sim)

# Plot precision of the observation error and PCA truncation error
lambda_plot = seq(0,10.0, length.out = 1000)
dev.new(noRStudioGD = TRUE)
#plot posterior
hist(samples$lambda_y,
     main = "lambda_y",
     xlab = "lambda_y",
     col = "firebrick1",
     breaks = 25,
     freq = FALSE,
     # xlim = c(0,10.0),
     cex.axis=1.5,
     cex.lab=1.5)
# plot prior
prior_plot = dgamma(lambda_plot,shape=a_y,rate=b_y)
lines(lambda_plot,prior_plot,lwd=3,col="blue")
# plot of prior adjusted for equivalent, reduced-dimension normal-gamma model
adj_prior_plot = dgamma(lambda_plot,shape=a_y_dash,rate=b_y_dash)
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
colnames(tf_trans) = c("E11","t_ply","logK")
dev.new(noRStudioGD = TRUE) # plot in new window
#pairs(~E11 + t_ply + logK, data = tf_trans, 
#      col = "blue",
#      pch = 18)

prior_sam = cbind(rnorm(N_samples, mean = E11_mu, sd = E11_mu*E11_cov/100), rnorm(N_samples, mean = t_ply_mu, sd = t_ply_cov*t_ply_mu/100), runif(N_samples,log(K_lb),log(K_ub)))

dev.new(noRStudioGD = TRUE, width=5, height=4) # plot in new window
par(mfrow = c(1,3))
plot(prior_sam[,1],prior_sam[,2],
     col = "blue",
     xlab = expression("E"[11]*" (GPa)"),
     ylab = expression("t"["ply"]*" (mm)"),#'_ply (mm)",
     cex.lab = 2,
     cex.axis = 2,
     pch = 18)
points(tf_trans[,1],tf_trans[,2],
       col = "red",
       pch = 18)


plot(prior_sam[,1],prior_sam[,3],
     col = "blue",
     xlab = expression("E"[11]*" (GPa)"),
     ylab = "log(K)",
     cex.lab = 2,
     cex.axis = 2,
     pch = 18)
points(tf_trans[,1],tf_trans[,3],
       col = "red",
       pch = 18)


plot(prior_sam[,2],prior_sam[,3],
     col = "blue",
     xlab = expression("t"["ply"]*" (mm)"),#'_ply (mm)",
     ylab = "log(K)",
     cex.lab = 2,
     cex.axis = 2,
     pch = 18)
points(tf_trans[,2],tf_trans[,3],
       col = "red",
       pch = 18)


#--------------------------------------------------------------------------------

# Code for making predictions - here I'm assuming we're using the transformed 
# version for matrices which are not full rank

N_sam_plot = N_samples
#N_sam_plot = 100
# Pick N_sam_plot random samples without repetition
rand_ind = sample.int(N_samples,N_sam_plot)

w_star_mu = matrix(0, p_eta, N_sam_plot)
w_star_sigma = array(0, c(p_eta, p_eta, N_sam_plot))
w_star = matrix(0, p_eta, N_sam_plot)
eta_sam = matrix(0, n_eta, N_sam_plot)
eta_y = matrix(0, n_y, N_sam_plot)

beta_w_i = -4.0 * log(rho_w)
lambda_w_i = lambda_w
lambda_eta_i = lambda_eta
# Define the covariance matrix for the emulator weights, evaluated at the
# experimental data points. As we only have one experiment, and
# R(theta,theta) = 1, this also reduces to a diagonal matrix.
sigma_u = diag(1/lambda_w_i, nrow=p_eta, ncol=p_eta)
# Define the covariance matrix for the emulator weights, evaluated at the
# training data points.
sigma_w = matrix(0, m*p_eta, m*p_eta)
for (j in 1:p_eta) {
  # Calculate the covariance matrix for the ith emulator weight
  sigma_w[((j-1)*m+1):(j*m), ((j-1)*m+1):(j*m)] = ARD_SE_cov(tc, lambda_w_i[j], beta_w_i[((j-1)*q+1):(j*q)], 0)
}

# Sort out predictions for full-rank code from here!!
for (i in 1:N_sam_plot){
  print(i)
  sam_ind = rand_ind[i]
  # Construct covariance matrix for the current sample
  
  # Code for if there are only one set of discrepancy parameters to learn about

  tf_i = samples$tf[sam_ind, ]
  lambda_y_i = samples$lambda_y[i]
  
  # Evaluate cross-covariance of emulator weights between training data and experimental data points
  sigma_uw = matrix(0, p_eta, m*p_eta)
  for (j in 1:p_eta) {
    sigma_uw[j,((j-1)*m+1):(j*m)] = ARD_SE_cov_non_sym(t(as.matrix(tf_i)), tc, lambda_w_i[j], beta_w_i[((j-1)*q+1):(j*q)])
  }
  
  # The covariance matrices used in predictions will differ depending on whether
  # or not B is full rank
  sigma_z = matrix(0,(m+1)*p_eta, (m+1)*p_eta)
  sigma_z[1:p_eta,1:p_eta] = sigma_u
  sigma_z[-(1:p_eta),-(1:p_eta)] = sigma_w
  sigma_z[1:p_eta,-(1:p_eta)] = sigma_uw
  sigma_z[-(1:p_eta),1:p_eta] = t(sigma_uw)
    
  # Adjust the covariance matrix sigma_z to include transformed emulator and 
  # experimental error terms
  sigma_z_hat = matrix(0,(m+1)*p_eta, (m+1)*p_eta)
  sigma_z_hat[1:p_eta,1:p_eta] = BTWyBinv/lambda_y_i
  sigma_z_hat[-(1:p_eta),-(1:p_eta)] = KTKinv/lambda_eta_i
  sigma_z_hat = sigma_z_hat + sigma_z

  # Define cross-correlation of emulator predictions w_star with emulator training
  # data. Note that in this application, because we are making a prediction at the
  # same value of controlled input x, and uncontrolled inputs are set to their
  # calibrated values as was the case for the experimental data, then this is 
  # given by sigma_u_w', Note that this won't be the case if we want to make a
  # prediction at a different x, or want to evaluate just the emulator across a 
  # range of t values
  sigma_w_w_star = t(sigma_uw)
  # Define cross-correlation terms for z with emulator prediction w_star
  sigma_z_w_star = rbind(sigma_u, sigma_w_w_star)
  
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
  w_star_mu[,i] = t(sigma_z_w_star) %*% Ainv_z_hat
  w_star_sigma[,,i] = sigma_w_star - (t(sigma_z_w_star) %*% solve(sigma_z_hat,sigma_z_w_star))

  # Sample from Gaussian process (might not be necessary - could just cheat and 
  # look at the mean
  w_star[,i] = mvrnorm(n = 1, w_star_mu[,i], w_star_sigma[,,i])
  
  # Generate individual Samples
  eta_sam[,i] = (K_eta %*% w_star[,i])*sd_dt + mu_dt
  eta_y[,i] = (K_y %*% w_star[,i])*sd_dt + mu_y
}

# Good to here. No need to check. Fix Below for plot details comparison of residuals etc
sadsadsadasdasdsadsadsdsadjsakhda

# Plot histogram of reduced coefficients
dev.new(noRStudioGD = TRUE) # plot in new window
par(mfrow = c(1, p_eta))
for (i in 1:p_eta){
  hist(w_star[i,], 
       main =  paste("Feature",as.character(i)),
       xlab = paste("w_star", as.character(i)),
       col = "brown1",
       breaks = 15,
       freq = FALSE,
       # xlim = c(-4,4),
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
eta_mu = (K_eta %*% w_star_mu)*sd_dt + mu_dt
eta_y_mu = (K_y %*% w_star_mu)*sd_dt + mu_y

# Calculate the diagonal terms of the covariance matrix. Just do this as the
# full covariance is too big. Take the square root (standard deviation) as this
# is more meaningful
eta_sigma = matrix(0, n_eta, N_sam_plot)
eta_y_sigma = matrix(0, n_y, N_sam_plot)
for (i in 1:N_sam_plot){
  eta_sigma[,i] = sqrt(as.matrix(rowSums((K_eta %*% w_star_sigma[,,i]) * K_eta)))*sd_dt
  eta_y_sigma[,i] = sqrt(as.matrix(rowSums((K_y %*% w_star_sigma[,,i]) * K_y)))*sd_dt
}

# Integrate uncertainty out of mean and standard deviation by taking average
eta_mu_mu = rowMeans(eta_mu)
eta_y_mu_mu = rowMeans(eta_y_mu)
eta_sigma_mu = rowMeans(eta_sigma)
eta_y_sigma_mu = rowMeans(eta_y_sigma)
eta_sam_mu = rowMeans(eta_sam)
eta_sam_y_mu = rowMeans(eta_y)
eta_sam_sd = rowSds(eta_sam)
  
# Write all output to text files for plotting
write.csv(eta_mu, "outputs/eta_mu.csv", row.names = FALSE)
write.csv(eta_y_mu, "outputs/eta_y_mu.csv", row.names = FALSE)
write.csv(eta_sam, "outputs/eta_sam.csv", row.names = FALSE)
write.csv(eta_y, "outputs/eta_y.csv", row.names = FALSE)
write.csv(eta_sigma, "outputs/eta_sigma.csv", row.names = FALSE)
write.csv(eta_y_sigma, "outputs/eta_y_sigma.csv", row.names = FALSE)

write.csv(eta_mu_mu, "outputs/eta_mu_mu.csv", row.names = FALSE)
write.csv(eta_y_mu_mu, "outputs/eta_y_mu_mu.csv", row.names = FALSE)
write.csv(eta_sigma_mu, "outputs/eta_sigma_mu.csv", row.names = FALSE)
write.csv(eta_y_sigma_mu, "outputs/eta_y_sigma_mu.csv", row.names = FALSE)
write.csv(eta_sam_mu, "outputs/eta_sam_mu.csv", row.names = FALSE)

write.csv(eta_sam_sd, "outputs/eta_sam_sd.csv", row.names = FALSE)
write.csv(eta_sam_y_mu, "outputs/eta_y_sam_mu.csv", row.names = FALSE)
