import numpy as np
import pandas as pd



file_str = "LHSDesign60x6"
# file_str = "Problem_run"
max_load = -250.0
max_inc = 0.05
displacements = np.loadtxt("inputs\\" + file_str + "_displacements_load=" + str(max_load) + "_max_inc=" + str(max_inc) +".csv",delimiter=",",skiprows=1)
print(displacements.shape)
sadsad
#displacements = pd.read_csv("inputs\\LHSDesign60x6_displacements_load=-250.0_max_inc=0.1.csv",sep=",")
RFs = np.loadtxt("inputs\\" + file_str + "_RFs_load=" + str(max_load) + "_max_inc=" + str(max_inc) + ".csv", delimiter=',', skiprows=1)
with open ("inputs\\" + file_str + "_incs_load=" + str(max_load) + "_max_inc=" + str(max_inc) + ".txt",'r') as f:
    increments = [[float(increment) for increment in line.strip().split(',')] for line in f.readlines()]

# Play around with different data structures: json or dataframe
displacements_struct = {"Sample" : []}
# Creates an empty dataframe with all the correct column names. Not sure how to add to these in the loop
displacements_df = pd.DataFrame(columns=["Sample","Frame","Increment","RF1","RF2","RF3","Node","U","V","W"])
frame_cnt = 0
for i, sample in enumerate(increments):
    displacements_struct["Sample"].append({"Frame" : []})
    # displacements_df["Sample"].append(i) # Doesn't work. Can append a series though. Look at other ways of constructing pandas dataframes on the fly
    for frame in sample:
        displacements_struct["Sample"][i]["Frame"].append({"Increment" : frame,
                                                           "RFs" : RFs[3*frame_cnt:3*(frame_cnt+1)].tolist(),
                                                           "Displacements" : displacements[:,3*frame_cnt:3*(frame_cnt+1)]})
        frame_cnt += 1
        
    #print(sample)

# displacements_struct["sample"][4]["frame"][-1]

displacements_struct["sample"][55] # This is the problem sample


displacements_struct_all = displacements_struct
disp_i = displacements_struct_all["sample"][56]["frame"][-1]#["Displacements"]
disp_RP = []
for sample in displacements_struct_all["sample"]:
    disp_RP.append(sample["frame"][-1]["Displacements"][1,2])


unique_disp_RP = set(disp_RP)
len(unique_disp_RP)
# An idea... Figure out exactly what to do with this later....

# Provided that I have reasonable displacement data even for the cases in which
# the analysis didn't complete, I can use this as emulator training data. Check
# out force displacement curves to see if this is the case
# Something funny happened with sample #56 (python index 55). Have a look at the
# output I re=ran tp see what happened. IS it possible that python just re-loaded
# saved outputs from previous sample. Check for this.
# Plot load displacement curve at RP to get indication of the sort of behaviour
# Also, research a bit more about dataframes
# Fixed value 200kN, and nonlinear emulator codes can now be written around this.
# For nonlinear treat load both as an input, but also do in the svd. See what
# works better

asdsdsad

# Plot force-displacement of reference pointe eventually?

infile = "inputs\\LHSDesign60x6"
x_DoE = np.loadtxt(infile + ".csv", delimiter = ",", skiprows = 1, ndmin = 2)
print(x_DoE[55,:])
np.savetxt("Problem_run.csv", x_DoE[55,:].reshape((1,-1)), delimiter=",", header="Blank Line", )

N, d = x_DoE.shape # Number of samples and number of inputs
N_complete = 56
x_DoE = x_DoE[N_complete:,:]
iterable = enumerate(x_DoE, start=N_complete)
for i, x_i in iterable:
    print(i)
    print(x_i)

with open ('RF_buffer.csv','r') as f:
    RFs = [[float(RF) for RF in line.strip().split()] for line in f.readlines()]    


n_incs = 0
for run in increments:
    n_incs += len(run)