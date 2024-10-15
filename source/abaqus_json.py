import numpy as np

def extract_const_frame(in_dict, disp_str):
    # Assume number of frames is the same for each sample, and the number of nodes 
    # the same for each frame.
    m = len(in_dict["Sample"]) # Number of samples
    n_frames = len(in_dict["Sample"][0]["Frame"]) # Number of frames
    n_nodes = len(in_dict["Sample"][0]["Frame"][0]["Displacements"]) # Number of nodes  
    out_mat = np.empty((n_nodes*n_frames*len(disp_str),m))
    
    # Assign a column index depending on which displacement component is of interest
    col_dict = {"u" : 0, "v" : 1, "w" : 2}
    col_ind = [col_dict[disp] for disp in disp_str]
  
    # Loop over output for all samples
    for i, sample in enumerate(in_dict["Sample"]):
        if (len(sample["Frame"]) != n_frames):
            raise Exception("There are a different number of frames across samples")
        disp_i = np.empty((0, 3))
        # Loop over each frame and extract the displacement
        for frame in sample["Frame"]:
            frame_j = np.array(frame["Displacements"])
            if frame_j.shape[0] != n_nodes:
                raise Exception("There are a different number of nodes across samples")

            # Concatenate the outputs of each frame
            disp_i = np.concatenate((disp_i, frame_j), axis = 0)

        # Store the component(s) of interest to the matrix of model outputs
        disp_i = disp_i[:,col_ind].flatten(order = 'F') # Reshape in column order
        out_mat[:,i] = disp_i
  
    return(out_mat,n_nodes,n_frames)