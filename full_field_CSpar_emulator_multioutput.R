# Fit an emulator to the full displacement vector output from a C-Spar Finite 
# Element model
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

p_eta = 15 # Number of basis functions retained for the emulator from SVD
exp_tol = 1e-6 # Tolerance variance fraction used to assess SVD convergence
# Delete or replace with something else. Or do disp_str = "u, "v", "w" to enable generic code
# disp_str = "w" # String which identifies the displacement component of interest (u,v, or w)
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
abaqus_displacements = fread(paste("inputs/",in_file,"_fixed_200kN.csv", sep=""))

disp_str = c("u","v","w")
print(disp_str[1])
disp_str = "w"
n_eta = nrow(abaqus_displacements) # number of output points per simulation
# We can treat disp_str as a vector - then loop over it's length within the 
# below loop - this should allow the same code to be used for this, so I can just
# delete this version from github
# Would have to multiply n_eta by the length of disp_str
# Implement then check below code
# Extract the displacement for the component of interest and store in a matrix
dt_simulation = matrix(NA,nrow = length(disp_str)*n_eta, ncol = m)
for (i in 1:m){
  for (j in 1:length(disp_str)) {
    dt_simulation[((j-1)*n_eta+1):(j*n_eta),i] = abaqus_displacements[[as.name(paste(disp_str[j],sprintf("_%d", i),sep=""))]]
  }
}
n_eta = n_eta*length(disp_str) # update n_eta definition for stan input

# The above code correctly formats the output regardless of the number of elements
# of disp_str
# I think the code is identical from here on out. check if this is the case
# If so, run the remaining portion of the emulator code - see if it works. 
# The just modify the appropriate lines of the main emulator code to match the
# above, and delete this file.


#-------------------------------------------------------------------------------

# Load in test points at which predictions are required
t_pred = fread("inputs/LHSDesign50x3_1.csv")
t_pred$K = log(t_pred$K)
colnames(t_pred)[3] = "log_K"
t_pred = as.matrix(t_pred)
n_pred = nrow(t_pred)
  
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
# also transform the test points for consistencu
t_pred = (t_pred - t_min[rep(1,n_pred),])/(t_max[rep(1,n_pred),]-t_min[rep(1,n_pred),])

#-------------------------------------------------------------------------------

# plot Design of Experiments used to generate training data
# pairs(~E11 + E22 + nu12 + nu23 + G12 + t_ply, data = tc,col = "blue",pch=4)
pairs(~E11 + t_ply + log_K, data = tc,col = "blue",pch=4)
pairs(~E11 + t_ply + log_K, data = t_pred,col = "blue",pch=4)

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

# Extract the first p_eta basis functions from the svd. 
# standardised such that the weights (columns of sqrt(m-1)*v) have unit variance
# K_eta = dt_svd$u[,1:p_eta]*matrix(rep(dt_svd$d[1:p_eta],2),nrow = n_eta, ncol = p_eta, byrow = TRUE)/sqrt(m-1)
K_eta = dt_svd$u[,1:p_eta]*matrix(dt_svd$d[1:p_eta],nrow = n_eta, ncol = p_eta, byrow = TRUE)/sqrt(m-1)
# K_eta = dt_svd$u[,1:p_eta]

# Sanity check that weights of SVD have zero mean and unit variance
print("mean of reduced dimension output w = ")
print(colMeans(dt_svd$v))
print("standard deviations of reduced dimension output w = ")
print(colSds(dt_svd$v*sqrt(m-1)))

# write basis functions to file for external plotting
write.csv(K_eta, "outputs/basis_50x3_multi.csv", row.names = FALSE)
# Also write the mean vector
write.csv(mu_dt, "outputs/model_mean.csv", row.names = FALSE)

# Plot magnitude of singular value d with increasing number of basis functions. 
# This gives an indication of how much each base contributes to the output variance.
# I've used this plot to select which value of p_eta to choose with this method
d_r = dt_svd$d[1:p_eta]
d_r_norm = d_r^2/sum(d_r^2)

par(mfrow = c(1,2))
plot(1:p_eta,d_r_norm,"type"="p","col"="red","pch"=4,"lwd"=3,cex=1.5,'xlab' = "Feature",'ylab'="normalised d_i",cex.axis=1.75,cex.lab=1.75)
plot(2:p_eta,d_r_norm[-1],"type"="p","col"="red","pch"=4,"lwd"=3,cex=1.5,'xlab' = "Feature",'ylab'="normalised d_i",cex.axis=1.75,cex.lab=1.75)
# How many basis functions are required to get within a tolerance fraction of the
# variance
tol = 1e-6
cumsum(d_r_norm) > (1-tol)
basis_tol = 1:p_eta
basis_tol = basis_tol[cumsum(d_r_norm) > (1-tol)]
basis_tol = basis_tol[1]
print("Number of basis functions required to represent output within tolerance = ")
print(basis_tol)

#-------------------------------------------------------------------------------

# This portion of the code deals with the matrix algebra from Section 2.2.4 of 
# Higdon et al., calculating the necessary quantities for passing to stan, where
# the sampling is undertaken

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
z_hat = as.vector(KTKinv%*%KTeta) # reduced-dimensional representation of emulator training data

# Define parameters for (gamma) prior distributions on the emulator error and 
# observation error terms
a_eta = 1.0     # Shape parameter for the lambda_eta prior
b_eta = 0.0001  # Rate parameter for the lambda_eta prior 

a_eta_dash = a_eta+(0.5*(m*(n_eta-p_eta))) # Adjusted shape parameter for the lambda_eta prior (Eq. 11 Higdon et al.)

# Stack reshaped model output eta into a single vector. This is the format used 
# in Higdon et al., although it wasn't necessary for the above calculations it
# is used to adjust the prior parameters.
eta_vec = as.vector(eta)
# Re-arraged version of b_eta_dash from that in Eq. (11) for the sake of 
# computational efficiency, using the fact that:
# eta'*(I - K*(K'*K)^-1*K')*eta = eta'*eta - (K'*eta)'*w_hat
b_eta_dash = as.numeric(b_eta + (0.5*(t(eta_vec)%*%eta_vec - t(KTeta) %*% z_hat)))

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
stan_data = list(m=m, q=q, n_eta=n_eta, p_eta=p_eta, a_eta_dash = a_eta_dash, b_eta_dash = b_eta_dash, z_hat = z_hat, tc = tc, KTKinv = KTKinv)
fit = stan(file = "source/Full_Field_Emulator_Higdon.stan",
           data = stan_data,
           iter = 4000,
           chains = 3,
           model_name = "full_field_emulator")

#-------------------------------------------------------------------------------

# This section of code deals with post-processing of the data coming out of the 
# simulations

# plot trace plots of simulation samples
stan_trace(fit, pars = c("rho_w", "lambda_w","lambda_eta"))

# Summarise results to check convergence
print(fit, pars = c("rho_w", "lambda_w", "lambda_eta"))

# extract samples from stan output

samples <- extract(fit)
N_samples = dim(samples$rho_w)[1] # get total number of samples

# Extract label of inputs for plots
labels = colnames(XT_sim) # don't think I need this but keeping just in case

# estimate mode of the posterior distribution
modes = rep(0,(p_eta*(q+1))+1)
for (i in 1:(p_eta*q)){
  modes[i] = estimate_mode(samples$rho_w[,i])
}
for (i in 1:p_eta){
  modes[p_eta*q + i] = estimate_mode(samples$lambda_w[,i])
}
modes[(p_eta*(q+1))+1] = estimate_mode(samples$lambda_eta)

# write modes to csv for use in subsequent modelling
write.csv(modes, "outputs/emulator_modes_50x3_multi.csv", row.names = FALSE)

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
       xlim = c(0,2.5),
       cex.axis=1.5,
       cex.lab=1.5)
  # plot prior
  prior_plot = dgamma(lambda_plot,shape=5,rate=5)
  lines(lambda_plot,prior_plot,lwd=3,col="blue")
}

# Plot precision of the PCA truncation error
dev.new(noRStudioGD = TRUE)
#plot posterior
hist(samples$lambda_eta,
     main = "lambda_eta",
     xlab = "lambda_eta",
     col = "firebrick1",
     breaks = 25,
     freq = FALSE,
     xlim = c(0,1e8),
     cex.axis=1.5,
     cex.lab=1,5)
# plot prior
lambda_plot = seq(0,1E8, length.out = 1000)
prior_plot = dgamma(lambda_plot,shape=a_eta,rate=b_eta)
lines(lambda_plot,prior_plot,lwd=3,col="blue")
adj_prior_plot = dgamma(lambda_plot,shape=a_eta_dash,rate=b_eta_dash)
lines(lambda_plot,adj_prior_plot,lwd=3,col="green")

#-------------------------------------------------------------------------------

# Code for making predictions - here I'm assuming we're using the transformed 
# version for matrices which are not full rank

# N_sam_plot = 1
# Pick N_sam_plot random samples without repetition
N_sam_plot = 10
rand_ind = sample.int(N_samples,N_sam_plot)

w_star_mu = array(0, c(p_eta, N_sam_plot, n_pred))
w_star_sigma = array(0, c(p_eta, p_eta, N_sam_plot, n_pred))
w_star = array(0, c(p_eta, N_sam_plot, n_pred))
eta_sam = array(0, c(n_eta, N_sam_plot, n_pred))
eta_mu = array(0, c(n_eta, N_sam_plot, n_pred))
eta_sigma = array(0, c(n_eta, N_sam_plot, n_pred))

# Sort out predictions for full-rank code from here!!
for (i in 1:N_sam_plot){
  print(i)
  sam_ind = rand_ind[i]
  # Construct covariance matrix for the current sample
  # I could probably package the below code into an external function
  
  # Extract the inferred parameters for the current sample
  beta_w_i = samples$beta_w[sam_ind,]
  lambda_w_i = samples$lambda_w[sam_ind,]
  lambda_eta_i = samples$lambda_eta[sam_ind]
  
  # Define the covariance matrix for the emulator weights, evaluated at the
  # training data points.
  sigma_z = matrix(0, m*p_eta, m*p_eta)
  for (j in 1:p_eta) {
    # Calculate the covariance matrix for the ith emulator weight
    sigma_z[((j-1)*m+1):(j*m), ((j-1)*m+1):(j*m)] = ARD_SE_cov(tc, lambda_w_i[j], beta_w_i[((j-1)*q+1):(j*q)], 0)
  }
  
  # Adjust the covariance matrix sigma_z to include transformed emulator and 
  # experimental error terms
  sigma_z_hat = matrix(0,m*p_eta, m*p_eta)
  sigma_z_hat = sigma_z + KTKinv/lambda_eta_i
  
  sigma_z_w_star = array(0, c(m*p_eta, p_eta, n_pred))
  # Determine covariance of training data with predictions
  for (j in 1:n_pred) {
    for (k in 1:p_eta) {
      sigma_z_w_star[((k-1)*m+1):(k*m),k,j] = ARD_SE_cov_non_sym(tc, t(as.matrix(t_pred[j,])), lambda_w_i[k], beta_w_i[((k-1)*q+1):(k*q)])
    }
  }
  
  # Define correlation of emulator predictions with themselves. In an ideal world
  # I'd also look at cross-correlations but I've chosen not to do that for the 
  # sake of efficiency. The same covariance matrix can therefore be used for all
  # predictions as the autocorrelation is 1 irrespective of the prediction 
  sigma_w_star = diag(1/lambda_w_i)
  
  # Two different methods for making predictions. The first method is quicker, but 
  # I think the second is more numerically stable, which seems to make a different
  # when the variance is small. Consider using this if the simulation is taking 
  # too long
  
  # Explicitly calculating inverse, then reusing for all predictions. 
  # Ainv = solve(sigma_z_hat)
  # w_mu = t(Lsigma_z_w_star) %*% Ainv %*% z_hat
  # w_sigma = sigma_w_star - (t(Lsigma_z_w_star) %*% Ainv %*% Lsigma_z_w_star)

  # Solving using solve
  Ainv_z_hat = solve(sigma_z_hat,z_hat)
  # Store mean and covariance matrices of discrepancy and adjusted prediction Gaussian processes
  for (j in 1:n_pred){
    w_star_mu[,i,j] = t(sigma_z_w_star[,,j]) %*% Ainv_z_hat
    w_star_sigma[,,i,j] = sigma_w_star - (t(sigma_z_w_star[,,j]) %*% solve(sigma_z_hat,sigma_z_w_star[,,j]))
    # Sample from the Gaussian process
    w_star[,i,j] = mvrnorm(n = 1, w_star_mu[,i,j], w_star_sigma[,,i,j])
    # Generate individual Samples
    eta_sam[,i,j] = (K_eta %*% w_star[,i,j])*sd_dt + mu_dt
    # Also look at samples of the mean for output
    eta_mu[,i,j] = (K_eta %*% w_star_mu[,i,j])*sd_dt + mu_dt
    # Calculate the diagonal terms of the covariance matrix. Take the standard 
    # deviation as this is more meaningful
    eta_sigma[,i,j] = sqrt(as.matrix(rowSums((K_eta %*% w_star_sigma[,,i,j]) * K_eta)))*sd_dt
  }
}

# Get average full-field displacement across all samples
#eta_sam_matrix = matrix(0, nrow = n_eta, ncol = n_pred)
#for (i in 1:n_pred){
#  eta_sam_matrix[,i] = rowMeans(eta_sam[,i,])
#}
eta_sam_matrix = apply(eta_sam, 3, rowMeans)
eta_mu_matrix = apply(eta_mu, 3, rowMeans)
eta_sigma_matrix = apply(eta_sigma, 3, rowMeans)
# Compare against known output
cross_val_disp = fread("inputs/LHSDesign50x3_1_displacements.csv")
cross_val = matrix(0,1,2*n_pred)
for (i in 1:n_pred){
  cross_val[,i] = cross_val_disp[[as.name(sprintf('w_%d', i))]][2]
  cross_val[,(i+n_pred)] = max(cross_val_disp[[as.name(sprintf('u_%d', i))]])
}

# Quick cross-validation exercise - the second node is the reference point at 
# which the load is applied. Look at this to check predictions are approximately
# the right magintude, which seems to be the case
View(rbind(cbind(matrix(eta_sam_matrix[2,],nrow=1),matrix(colMaxs(eta_sam_matrix[-(1:(n_eta/2)),]),nrow=1)),cbind(matrix(eta_mu_matrix[2,],nrow=1),matrix(colMaxs(eta_mu_matrix[-(1:(n_eta/2)),]),nrow=1)),cross_val))

# Plot histogram of reduced coefficients
dev.new(noRStudioGD = TRUE) # plot in new window
par(mfrow = c(1, p_eta))
pred_ind = 1
for (i in 1:p_eta){
  hist(w_star[i,,pred_ind], 
       main =  paste("Feature",as.character(i)),
       xlab = paste("w_star", as.character(i)),
       col = "brown1",
       breaks = 5,
       freq = FALSE,
       xlim = c(-4,4),
       cex.lab = 1.25,
       cex.axis = 1.25)
  # overlay plot of prior distribution
  w_plot = seq(-4,4, length.out = 100)
  prior_plot = dnorm(w_plot,mean = 0, sd = 1)
  lines(w_plot,prior_plot,lwd=3,"col"="blue")
}


# Write all output to text files for plotting
# Too big - consider just outputting samples for a single prediction
# write.csv(eta_sam, "outputs/eta_sam.csv", row.names = FALSE)
write.csv(eta_sam_matrix, "outputs/eta_sam_mu.csv", row.names = FALSE)
write.csv(eta_mu_matrix, "outputs/eta_mu_mu.csv", row.names = FALSE)
write.csv(eta_sigma_matrix, "outputs/eta_sigma_mu.csv", row.names = FALSE)
