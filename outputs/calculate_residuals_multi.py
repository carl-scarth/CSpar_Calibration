import numpy as np

# Load in sample mean and standard deviation from gaussian process
output_cloud = np.loadtxt('eta_y_sam_mu.csv', delimiter = ',', skiprows = 1, ndmin=2)
output_sigma = np.loadtxt('eta_y_sigma_mu.csv', delimiter = ',', skiprows = 1)
output_sigma_2 = np.loadtxt('eta_y.csv',delimiter = ',', skiprows = 1)
output_sigma_2 = np.std(output_sigma_2, axis = 1)

# Load in experimental data
y = np.loadtxt('Ext_LCorner_Image_0115_0.tiff_nat_coord_rad_trim.csv', delimiter = ',', skiprows = 1)
xyz = y[:,0:3]
uvw = y[:,7:10]
N = len(output_cloud)

# Calculate residuals
res_w = np.absolute(output_cloud[0:(N//2),:]-uvw[:,2].reshape(N//2,1))
res_u = np.absolute(output_cloud[(N//2):,:]-uvw[:,0].reshape(N//2,1))

# Write output alongside coordinates
np.savetxt('GP_mean_residuals.csv',np.concatenate((xyz,res_w,res_u,output_sigma.reshape(2,N//2).T,output_sigma_2.reshape(2,N//2).T),1), delimiter =',', header = 'X, Y, Z, residual_w, residual_u, sigma_w, sigma_u, sigma_2_w, sigma_2_u', comments = '')