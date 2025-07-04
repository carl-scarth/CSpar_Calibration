import os
import pandas as pd
import numpy as np
import json
import sys
import matplotlib.pyplot as plt
import pymc as pm
import arviz as az
import pytensor.tensor as tt

src_path = "src"
sys.path.insert(0, src_path)

from abaqus_json import extract_const_frame, basis_mean_to_json
from transform_input_output import normalise_inputs, standardise_vector_output, rescale_output, rescale_inputs
from dimension_reduction import svd_basis, reduce_dim_emulator
from vtk_output import plot_basis_nonlinear_vtk, plot_basis_nonlinear_point

class ff_mean_Zero(pm.gp.mean.Mean):
    def __init__(self, p_eta):
        super().__init__()
        self.p_eta = p_eta
        print(self.p_eta)
        
    def __call__(self, X):
        return(tt.zeros(X.shape[0]*self.p_eta))

# Create a custom covariance function for full-field emulator
class ff_cov(pm.gp.cov.Covariance):
    def __init__(self, 
                 input_dim: int,  # Number of inputs
                 m: int,                # Number of training samples
                 p_eta: int,            # Number of basis functions
                 ls: tt.tensor,         # Correlation lengths
                 sigma_em: tt.tensor,   # Scale parameter
                 lambda_eta: tt.tensor = None, # Precision of SVD truncation error
                 KTK_inv = None,        # Inverse of inner product of basis matrix 
                 untile_inp = False):   # Does the training data need to be untiled before evaluating       
        
        super().__init__(input_dim) # Inherit characteristics of the pymc covariance class (declared in class definition)
        # Note: X and Xs aren't assigned until calling a method
        self.q = input_dim
        self.m = m
        self.p_eta = p_eta
        self.ls = ls
        self.sigma_em = sigma_em
        if KTK_inv is None:
            self.KTK_inv = np.zeros((p_eta*m, p_eta*m))
        else:
            self.KTK_inv = KTK_inv
        if lambda_eta is None:
            self.lambda_eta = 1.0
        else:
            self.lambda_eta = lambda_eta
        self.untile_inp = untile_inp

    def diag(self):
        # Diagonal terms
        #diag_em = tt.zeros(m*p_eta)
        #for i in range(self.p_eta):
        #    diag_em[i*self.m:(i+1)*self.m] = tt.square(self.sigma_em[i])
        # More concise code - dunno if it'll work
        diag_em = [tt.alloc(tt.square(self.sigma_em[i]), self.m) for i in range(self.p_eta)]
        diag_em = tt.stack(diag_em)

        return diag_em
    
    def full(self, X, Xs = None):
        # Xs is passed as a second input to the covariance to return the cross-covariance between X and Xs
        # ls does not stand for l^2, so actually covariance is exp(-(X-X')/(2*LS^2))
        if Xs is None:
            # For compatibility with the output data it may be necessary to tile the training data. If so, untile 
            # before evaluating the covariance
            if self.untile_inp:
                X = X[0:m,:]
            
            cov_full = tt.zeros([m*p_eta, m*p_eta])
            # Populate entries for each principal component
            for i in range(self.p_eta):
                cov_em_i = pm.gp.cov.Constant(self.sigma_em[i]**2) * pm.gp.cov.ExpQuad(self.q,ls=self.ls[i*self.q:(i+1)*self.q])
                # Can't use indexing with tensors, instead use subtensor.set_subtensor
                cov_full = tt.subtensor.set_subtensor(cov_full[i*self.m:(i+1)*self.m, i*self.m:(i+1)*self.m],cov_em_i(X))

        else:
            # Assumes Xs is to be used for prediction, and so there's no need to account for the observation error
            # Write when I need this
            pass
            #Xs = tt.as_tensor_variable(Xs)
            #cov_full = cov_em(X,Xs)
        
        return(cov_full)

if __name__ == "__main__":

    # Set up parameters which govern the formulation

    in_file = "LHSDesign100x8" # File identifier string for input and output files
    pred_file = "LHSDesign10x8" # File identifier for input values where predictions are required
    in_suffix = "150kN_zeroed" # Other text added to the file identifier

    p_eta = 20 # Number of basis functions retained for the emulator from SVD
    exp_tol = 1e-4 # Tolerance variance fraction used to assess SVD convergence
    disp_str = ["u", "w"] # String which identifies the displacement component of interest (u,v, or w) or list of multiple
    q_y = len(disp_str)
    print_svd_output = True # Print diagnostic output of svd to the terminal?
    reduce_gamma_param = True # Alter precision prior gamma parameters to account for the dimension reduction
    export_modes = True # Calculate modes and means of emulator hyperparameters and write to file?
    output_json = True  # Output structured data via json
    output_vtk = True   # Output plottable data via vtk
    debug = False       # Print output for debugging

    # Define parameters of the gamma prior on the error associated with truncating
    # the series expansion for the model output
    a_eta = 1.0     # Shape parameter for the lambda_eta prior
    b_eta = 0.0001  # Rate parameter for the lambda_eta prior
    N_iter = 2500 # Number of samples per chain

    #-------------------------------------------------------------------------------

    # Set up simulation data

    # Load in emulator training data input values from Design of Experiments. 
    XT_sim = pd.read_csv(os.path.join("inputs",in_file+".csv"))

    # Load in test points at which predictions are required
    XT_pred = pd.read_csv(os.path.join("inputs",pred_file+".csv"))
    n_pred = XT_pred.shape[0] # number of predictions

    # Determine useful quantities from model inputs and outputs. Variable names 
    # match the notation of Higdon et al. 2008
    m, q = XT_sim.shape # sample size of simulation data, and number of calibration inputs

    # Load in training data output displacement values from Abaqus.
    # Outputs are structured in a json across samples and load increments
    with open(os.path.join("inputs", "_".join((in_file,"output_struct",in_suffix))+'.json'),'r') as f:
        abaqus_dict = json.loads(f.readline())

    dt_simulation, n_nodes, n_frames = extract_const_frame(abaqus_dict, disp_str)
    n_eta = dt_simulation.shape[0] # total number of output points per simulation

    #-------------------------------------------------------------------------------

    # Standardise the data

    # Normalise inputs such that training data is on the unit hypercube
    tc, t_min, t_max = normalise_inputs(XT_sim)
    tc = tc.to_numpy() # Covariance function doesn't like pandas

    # Normalise test data in the same way as the training data for consistency
    t_pred = normalise_inputs(XT_pred, x_min=t_min, x_max=t_max)

    # Standardise the outputs to have zero mean (for each row) and unit standard
    # deviation (for each displacement)
    eta, mu_dt, sd_dt = standardise_vector_output(dt_simulation, q_y = q_y)

    #-------------------------------------------------------------------------------

    # Perform dimension reduction on data via SVD
    K_eta, p_eta = svd_basis(eta, exp_tol = exp_tol, print_output = print_svd_output, csv_label = "nonlinear_"+in_file)

    # Output training data mean and basis vectors to json if required
    if output_json or output_vtk:
        # Structure output as dictionary
        basis_dict = basis_mean_to_json(n_frames, n_nodes, K_eta, mu_dt, disp_str)

    # Write output to json
    if output_json:
        with open(os.path.join("outputs","basis_nonlinear_"+in_file+".json"),'w') as f:
            f.write(json.dumps(basis_dict))

    # Produce 3D vtk plot of output
    if output_vtk:
        # Last four nodes in basis outputs are reference points which aren't listed in the node file. Skip these.
        plot_basis_nonlinear_vtk(in_dict = basis_dict, file_str = in_file, mesh_str = "inputs\\new_spar_mesh", out_folder = "outputs", skip_nodes = [-4,-3,-2,-1])

    if print_svd_output:
        # Produce line plots basis functions at location of maximum (absolute) displacement
        max_ind = {label : np.argmax(np.abs(value)) for label, value in basis_dict["Frame"][-1]["Training_Data_Mean"].items()}
        plot_basis_nonlinear_point(max_ind, in_dict = basis_dict, location_title = "point of maximum displacement")
        # If longittudinal displacement w is in the output, also extract the second row from
        # each basis, gives corresponds to displacement at the reference point. 
        if "w" in disp_str:
            plot_basis_nonlinear_point({"w" : -1}, in_dict = basis_dict, location_title = "Reference Point")

    # Reduce the dimension of the output data and determine pseudo-inverse matrix
    # Adjust the prior parameters if needed
    if reduce_gamma_param:
        z_hat, KTK_inv, a_eta_dash, b_eta_dash = reduce_dim_emulator(eta, K_eta, a_eta=a_eta, b_eta=b_eta)
    else:
        z_hat, KTK_inv = reduce_dim_emulator(eta, K_eta, a_eta=a_eta, b_eta=b_eta)
        a_eta_dash = a_eta
        b_eta_dash = b_eta

    # Replicate training data for compatibility with structure of reduced dimensional data
    tc_rep = np.tile(tc, (p_eta, 1))
    with pm.Model() as full_field_em:
        # Priors on emulator hyperparameters and noise parameter
        # Updata to improve convergence. Stan guidance on GP priors:
        #'https://github.com/stan-dev/stan/wiki/Prior-Choice-Recommendations#priors-for-gaussian-processes'
        # Example I'm following:
        # https://www.pymc.io/projects/examples/en/latest/gaussian_processes/GP-Latent.html
            
        # Scale parameter
        sigma_em = pm.HalfNormal("sigma_em", sigma = 1.0, shape = p_eta)
        # Precision of error due to SVD truncation
        lambda_eta = pm.Gamma("lambda_eta", alpha = a_eta_dash, beta = b_eta_dash, shape = 1)                    
        # Define correlation length parameters
        ls = pm.InverseGamma("ls", alpha = 4.0, beta = 4.0, shape = q*p_eta)
        # Define mean and covariance functions
        #mean_func = ff_mean_Zero(p_eta)
        mean_func = pm.gp.mean.Zero()
        cov_func = ff_cov(q, m, p_eta, ls, sigma_em, lambda_eta=lambda_eta, KTK_inv=KTK_inv, untile_inp = True)
        if debug:
            print(f"sigma_em = {sigma_em.eval()}\n")
            print(f"ls = {ls.eval()}\n")
            print(f"Mean vector = \n{mean_func(tc_rep).eval()}\n")
            print(f"Covariance matrix = \n{cov_func(tc_rep).eval()}\n")
            print((cov_func(tc_rep).shape.eval()))
            print((mean_func(tc_rep).shape.eval()))

        # Implementation using multivariate normal: (I think this vvvv isn't working because eval returns the same sample each iteration)
        # y_ = pm.MvNormal("y", mu=mean_func(tc_rep).eval(), cov=cov_func(tc).eval(), observed = z_hat)
        
        # Alternative example using marginal likelihood GP with nugget (more stable)
        gp = pm.gp.Marginal(mean_func = mean_func, cov_func = cov_func)
        y_ = gp.marginal_likelihood("y", X = tc_rep, y = z_hat, sigma = 1e-8)

        # Draw samples from posterior
        idata = pm.sample(N_iter, target_accept=0.9)

    # Produce trace plots    
    az.plot_trace(idata, combined=True, figsize=(10, 8))
    plt.show()
    # Some issues with rhat - low - try seeing try different priors, fewer basis. Commit for now as code at least runs.
    asdsad


# think about having option for individual principal components vs all at once
# Might need to re-think priors - investigate other options to gamma

## plot Design of Experiments and test points
#pairs(tc,col="blue", pch=4, 
#      main = "Normalised Training Data Input values",
##      cex=1.5,
#      cex.labels = 1.75,
#      cex.axis=1.5,
#      cex.lab=1.5)
#pairs(t_pred,col = "blue", pch=4, 
#      main = "Normalised Test Data Input Values",
#      cex=1.5,
#      cex.labels = 1.75,
#      cex.axis=1.5,
#      cex.lab=1.5)

# Future Versions
# - FastMAP
# Use Higdon prior parameters

# NOTE - IT MIGHT MAKE MORE SENSE TO DEFINE A SEPARATE COVARIANCE FUNCTION FOR THE HIGDON FORMAT VS ONE WITH DIFFERENT PRIORS
# Consider getting rid of type hints in covariance function? - potentially inconsistent with default None argument