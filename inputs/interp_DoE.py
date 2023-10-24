# Interpolates model outputs for Design of Experiments to DIC coordinates
from math import prod
import numpy as np
import pandas as pd

def intp_nodes_to_cloud(el_ind, gh, f_node, conn, GH = [], skip_nodes = 0):
    # Interpolate nodal quantities f_node from nodes of a quad element mesh to 
    # n_cloud point cloud coordinates. 
    # el_ind = n_cloud vector of integer element indices, matching an element in
    # the connectivity matrix
    # gh = n_cloud x 2 matrix of element natural coordinates, with columns g and h
    # f_node = n_nodes x n_f matrix of nodal values of n_f quantities of interest
    # to be interpolated to point cloud coordinates
    # conn = n_nodes x 4 matrix containing element connectivities
    # GH = 2 x 4 matrix with natural coordinates of the 4 nodes connected to each
    # element with rows g and h. If not specified default Abaqus values are used.
    # skip_nodes = number of nodes to skip at the beginning of f_node This is 
    # useful when the first set of points in the output are not referenced by the 
    # connectivity file

    if (skip_nodes > 0):
        f_node = f_node[skip_nodes:,:]
    
    n_cloud = gh.shape[0] # Number of experimental data points
    n_f = f_node.shape[1] # Number of output quantities to interpolate

    # Define default Abaqus values for nodal natural coordinates, if not inputted
    if len(GH) < 1:
        #GH = np.array([[-1.0,1.0,1.0,-1.0],[-1.0,-1.0,1.0,1.0]])
        GH = np.array([[-1.0,-1.0],[1.0,-1.0],[1.0,1.0],[-1.0,1.0]])
    
    # Interpolate using the functions for Abaqus isoparametric quad elements
    f_cloud = np.empty((n_cloud,n_f))
    for i, el in enumerate(el_ind):
        # bases = 1.0 + gh[i,:]*GH.T
        bases = 1.0 + gh[i,:]*GH
        prod_bases = np.prod(bases, axis = 1, keepdims = True)*f_node[conn[el,:],:]
        # bases = 1.0 + gh[i,:].T*GH # 1-dimensional basis functions
        # prod_bases = np.prod(bases, axis=0)*f_node[conn[el,:],:]
        # f_cloud[i,:] = np.sum(prod_bases, axis=1)/4.0
        f_cloud[i,:] = np.sum(prod_bases, axis=0)/4.0

    return(f_cloud)

if __name__ == "__main__":
    # Load inputs from file
    conn_file = "nominal_shell_mesh_outer_surface_elements" # Connectivity
    DoE_file = "LHSDesign40x4_fixed_200kN"                  # Model output
    exp_data_file = "Interpolated_DIC_200kN"                # Experimental data
    conn = pd.read_csv(conn_file+".csv").to_numpy(dtype=int)
    conn = conn - 1 # Convert to Python indexing from Abaqus
    model_out = pd.read_csv(DoE_file+".csv")
    exp_data = pd.read_csv(exp_data_file+".csv")

    # Extract element index and natural coordinates
    el_ind = exp_data["Element"].to_numpy(dtype=int)
    gh = exp_data[["h","r"]].to_numpy()

    interp_outputs = intp_nodes_to_cloud(el_ind, gh, model_out.to_numpy(), conn, GH = [], skip_nodes = 2)
    interp_outputs = pd.DataFrame(interp_outputs, columns=model_out.columns.values)
    interp_outputs.to_csv(DoE_file+"_interp.csv", sep=",", index=False)

    # concatenate with coordinates for plotting
    interp_outputs = pd.concat((exp_data[["x_0_rot","y_0_rot","z_0_rot"]],interp_outputs),axis=1)
    interp_outputs.to_csv(DoE_file+"_interp_plot.csv", sep=",", index=False)