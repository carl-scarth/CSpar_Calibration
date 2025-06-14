from odbAccess import *
import numpy as np
import json
import sys

def update_outdict(frame, index, step, fix_to_ground):
    # Extract displacements from a given frame
    displacements = [value.data.tolist() for value in frame.fieldOutputs['U'].values]

    # Write dictionary for output in json format
    out_frame = {   "Index"         : index,
                    "Step"          : step,
                    "Increment"     : frame.frameValue,
                    "RFs"           : get_RFs(frame, fix_to_ground),
                    "Displacements" : displacements}

    return(out_frame)

def get_RFs(frame, fix_to_ground):
    # Extract the reaction force. The correct nodeset will depend on whether the model is fixed at the ground or not
    if fix_to_ground:
        RFs = frame.fieldOutputs['RF'].getSubset(region = o1.rootAssembly.instances['M'].nodeSets['FIXED_SUP_RP'])
    else:
        RFs = frame.fieldOutputs['RF'].getSubset(region = o1.rootAssembly.instances['M'].nodeSets['GROUND'])
    return(RFs.values[0].data.tolist())

# output to command prompt
sys.stdout = sys.__stdout__

# Extract inputs from the commandline string.
analysis = sys.argv[-1]
fix_to_ground = bool(sys.argv[-2])

# Open odb file. Assumes the scrip is run from the parent directory, and not from the Abaqus subfolder
o1 = session.openOdb(name = analysis + '.odb')
# Get number of nodes and frames
n_nodes = len(o1.steps['Step-1'].frames[0].fieldOutputs['U'].values) # Number of nodes
n_frames = len(o1.steps['Step-1'].frames) + 2 # Number of Frames (add two for first two steps)

# Initialise output dictionary
out_dict = {"Frame" : [], "Nodes" : []}
# Get nodal coordinates
out_dict["Nodes"] = [value.instance.getNodeFromLabel(value.nodeLabel).coordinates.tolist() for value in o1.steps['Step-0'].frames[0].fieldOutputs['U'].values]
# Extract output from initial steps (don't need anything from step 0.5 as this is 0 for step 1)
out_dict["Frame"].append(update_outdict(o1.steps['Step-0'].frames[0], 0, 0, fix_to_ground))
out_dict["Frame"].append(update_outdict(o1.steps['Step-0'].frames[-1], 1, 0, fix_to_ground))
# Loop over all frames in Step 1 and extract output
for i, frame in enumerate(o1.steps['Step-1'].frames):
    out_dict["Frame"].append(update_outdict(frame, i+2, 1, fix_to_ground))

# Output dictionary as json
with open(analysis + "_output.json",'w') as f:
    f.write(json.dumps(out_dict))