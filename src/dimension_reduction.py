import numpy as np
from numpy.linalg import svd, inv
import pandas as pd
import warnings
import matplotlib.pyplot as plt

def svd_basis(eta, p_eta = None, exp_tol = None, print_output = False, export_basis = True, csv_label = None):
  
    # Perform an SVD upon eta to generate a basis of principal components and 
    # return the first p_eta basis functions. If p_eta is not specified this 
    # is computed by determining how many functions are required to capture 
    # exp_tol fraction of the variance in the training data.
  
    # eta = n_eta x m matrix of simulation outputs, where n_eta is the number of
    #     output values per simulation, and m is the number of training data runs
    # p_eta = integer value of number of basis functions to be returned
    # exp_tol = desired fraction of the variance which is to be returned
    # print_output = Boolean indicating whether or not to print output to the 
    #     console, and produce convergence plots
    # export_basis = Boolean indicating whether or not to export basis functions
    #     to a csv for plotting externally
    # csv_label = String for labelling exported csv file 
  
    # Get output dimensions
    m = eta.shape[1]
  
    # Perform SVD on the centred data
    U, S, V = svd(eta, full_matrices = False)

    # Calculate the fraction of total variance captured by each singular value
    if exp_tol or print_output:
        dr_sqr = np.sum(S**2)
        # Fraction of the total variance captured by each singular value
        d_r_norm = S**2/dr_sqr

        # Determine the number of basis functions required to capture a specified
        # variance fraction
        if exp_tol:
            # How many functions are needed achieve a tolerance fraction of the variance?
            basis_tol = np.arange(m)
            basis_tol = basis_tol[np.cumsum(d_r_norm) > (1-exp_tol)]
            basis_tol = basis_tol[1]
            if not p_eta:
                p_eta = basis_tol
      
            if print_output:
                print("\nNumber of basis functions required to represent output within tolerance = " + str(basis_tol))
    
        if print_output:
            # Sanity check that weights of SVD have zero mean and unit variance
            print("\nmean of reduced dimension output w = ")
            print(np.mean(V, axis = 1))
            print("\nstandard deviations of reduced dimension output w = ")
            print(np.std(V*np.sqrt(m), axis = 1))
  
    # If neither p_eta nor exp_tol has been provided, retain all basis functions
    if not p_eta:
        warnings.warn("Neither p_eta nor exp_tol has been specified, and so all bases are retained")
        p_eta = m

    # If needed plot magnitude of singular value with increasing number of bases
    if print_output:
        fig, axes = plt.subplots(1,2)
        axes[0].plot(np.arange(p_eta)+1,d_r_norm[0:p_eta],"rx", markersize = 10)
        axes[0].set_xlabel("Feature")
        axes[0].set_ylabel("normalised $d_i$")
        
        # Omit first point for greater clarity on convergence
        axes[1].plot(np.arange(1,p_eta)+1,d_r_norm[1:p_eta],"rx", markersize = 10)
        axes[1].set_xlabel("Feature")
        axes[1].set_ylabel("normalised $d_i$")
        fig.suptitle("Convergence with Number of Basis Functions")
    
    # Extract the first p_eta basis functions from the svd, and standardise so the
    # coefficients (columns of sqrt(m)*v) have unit variance
    K = U[:,0:p_eta]*S[0:p_eta]/np.sqrt(m)
    
    # write basis functions and mean vector to file for external plotting if needed
    if export_basis:
        col_labels = ["K_basis_"+str(i+1) for i in range(p_eta)]
        out_frame = pd.DataFrame(K,columns = col_labels)
        if csv_label:
            out_frame.to_csv("outputs/basis_"+csv_label+".csv", index=False)
        else:
            out_frame.to_csv("outputs/basis.csv", index=False)    
  
    return K, p_eta

def reduce_dim_emulator(eta, K, a_eta = None, b_eta = None, orthog_K = True):
    # Perform the matrix algebra from Section 2.2.4 of Higdon et al., 
    # Focus on quantities relating to the emulator

    # eta = n_eta x m numpy array of simulation outputs, where n_eta is the number of
    # output values per simulation, and m is the number of training data runs
    # K = n_eta x p_eta numpy array of basis functions used to decompose eta, where
    # p_eta is the number of basis functions
    # a_eta = Shape parameter of the gamma prior
    # b_eta = Rate parameter of the gamma prior
    # orthog_K = Boolean stating if the basis functions K are orthogonal
    
    p_eta = K.shape[1]   # Number of basis functions
    n_eta, m = eta.shape # Number of training samples

    # Populate dense matrix of products of products with KTK_dense[i,j] = k_i'*k_j
    KTK_dense = np.matmul(K.T,K)
    print(KTK_dense)
    print("")
    # If k_i are orthogonal then K^T*K is diagonal and the inverse is specified 
    # directly via the closed-form expression. 
    if orthog_K:
        KTKinv = np.zeros((m*p_eta, m*p_eta))
        for i in range(p_eta):
            KTKinv[i*m:(i+1)*m, i*m:(i+1)*m] = np.eye(m)*(1/KTK_dense[i,i])

    else:
        # Otherwise populate matrix and calculate inverse directly
        KTK = np.zeros((m*p_eta, m*p_eta))
        for i in range(p_eta):
            # Populate diagonal terms
            KTK[i*m:(i+1)*m, i*m:(i+1)*m] = np.eye(m)*KTK_dense[i,i]
            # Compute off-diagonal terms. 
            for j in range(i):
                KTK[i*m:(i+1)*m, j*m:(j+1)*m] = np.eye(m)*KTK_dense[i,j]
                KTK[j*m:(j+1)*m, i*m:(i+1)*m] = np.eye(m)*KTK_dense[i,j]

        # Determine inverse of K'*K.
        KTKinv = inv(KTK)
    
    # Calculate K^T*eta and arrange into the required format
    KTeta = np.matmul(K.T, eta).reshape((-1))
   
    # Do pseudo-inverse of basis matrix to reduced dimension of model output
    z_hat = np.matmul(KTKinv, KTeta)
    
    # If required, adjust prior parameters of the expansion truncation error to 
    # account for the dimension reduction (Eq. 11 Higdon et al.)
    if a_eta is not None and b_eta is not None:
        a_eta_dash = a_eta + (0.5*m*(n_eta-p_eta))
        # Stack model output eta into a single vector
        eta_vec = eta.reshape((-1))
        # Re-arranged b_eta_dash from Higdon et al. for computational efficiency using
        # eta'*(I - K*(K'*K)^-1*K')*eta = eta'*eta - (K'*eta)'*w_hat
        # More efficient if eta_vec kept as one d using dot product?
        b_eta_dash = b_eta + (0.5*(np.dot(eta_vec.T, eta_vec) - np.matmul(KTeta.T, z_hat)))
        return z_hat, KTKinv, a_eta_dash, b_eta_dash
    else:
        return z_hat, KTKinv
    