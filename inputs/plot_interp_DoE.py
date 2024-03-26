import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("C:\\Users\\cs2361\\Documents\\Bayesian_Model_Calibration\\source")
from utils import * 

set_plot_params()

infile = "LHSDesign60x6_4"
exp_file = "CS02P_mean_disp_412.5mm"
model_fd = pd.read_csv(infile+"_interp_mid.csv")

exp_fd = pd.read_csv(exp_file+".csv")
sort_data = True # Option to sort data in ascending input order (only relevant for 1D study, though will also work by sorting according to first column)
if sort_data:
    in_data = pd.read_csv(infile+".csv")
    # in_data.columns.values[0] = "alpha"
    col_label = in_data.columns.values[0]
    in_data = in_data.sort_values(col_label)
    sort_ind = in_data.index.values.tolist()
    model_force = model_fd[["force_" + str(ind) for ind in sort_ind]]
    model_disp = model_fd[["disp_" + str(ind) for ind in sort_ind]]
    model_crosshead = model_fd[["crosshead_" + str(ind) for ind in sort_ind]]
    colors = plt.cm.jet(np.linspace(0,1,model_disp.shape[1])) # get colourmap
else:
    model_force = model_fd[["force_" + str(ind) for ind in range(model_fd.shape[1]//3)]]
    model_disp = model_fd[["disp_" + str(ind) for ind in range(model_fd.shape[1]//3)]]
    model_crosshead = model_fd[["crosshead_" + str(ind) for ind in range(model_fd.shape[1]//3)]]

fig, ax = plt.subplots()
for i, (force,disp) in enumerate(zip(model_force, model_disp)):
    if sort_data:
        ax.plot(-model_disp[disp].to_numpy(), model_force[force].to_numpy(), color = colors[i], linewidth=1.0, label = "{} = {:.2f}".format(col_label, in_data.loc[sort_ind[i], col_label]))
    elif i == 0:
        ax.plot(-model_disp[disp].to_numpy(), model_force[force].to_numpy(), "-b", linewidth=1.0, label = "Prior sample")
    else:
        ax.plot(-model_disp[disp].to_numpy(), model_force[force].to_numpy(), "-b", linewidth=1.0)

ax.plot(-exp_fd["Mean DIC Displacement"].to_numpy(), -exp_fd["Load"].to_numpy(), "rx", label = "DIC")
ax.set_title("Prior Force-displacement vs DIC")
ax.set_xlabel("Displacement (mm)")
ax.set_ylabel("Force (kN)")
ax.legend()

fig2, ax2 = plt.subplots(1,2)
for i, (crosshead,disp) in enumerate(zip(model_crosshead, model_disp)):
    if sort_data:
        ax2[0].plot(-model_crosshead[crosshead].to_numpy(), -model_disp[disp].to_numpy(), color = colors[i], linewidth=1.0, label = "{} = {:.2f}".format(col_label, in_data.loc[sort_ind[i], col_label]))
    elif i == 0:
        ax2[0].plot(-model_crosshead[crosshead].to_numpy(), -model_disp[disp].to_numpy(), "-b", linewidth=1.0, label = "Prior sample")
    else:
        ax2[0].plot(-model_crosshead[crosshead].to_numpy(), -model_disp[disp].to_numpy(), "-b", linewidth=1.0)

print(exp_fd)

ax2[0].plot(-exp_fd["Crosshead_disp_estimate"].to_numpy(), -exp_fd["Mean DIC Displacement"].to_numpy(), "rx", label = "DIC")
ax2[0].set_xlabel("Crosshead displacement (mm)")
ax2[0].set_ylabel("Web tip displacment (mm)")
ax2[0].set_ylim((0,5))
ax2[0].set_xlim((0,4))
ax2[0].legend()

for i, (crosshead,force) in enumerate(zip(model_crosshead, model_force)):
    if sort_data:
        ax2[1].plot(-model_crosshead[crosshead].to_numpy(), model_force[force].to_numpy(), color = colors[i], linewidth=1.0, label = "{} = {:.2f}".format(col_label, in_data.loc[sort_ind[i], col_label]))
    elif i == 0:
        ax2[1].plot(-model_crosshead[crosshead].to_numpy(), model_force[force].to_numpy(), "-b", linewidth=1.0, label = "Prior sample")
    else:
        ax2[1].plot(-model_crosshead[crosshead].to_numpy(), model_force[force].to_numpy(), "-b", linewidth=1.0)

# ax2[1].plot(-exp_fd["Crosshead_disp_estimate"].to_numpy(), -exp_fd["Load"].to_numpy(), "rx", label = "DIC")
ax2[1].set_xlabel("Crosshead displacement (mm)")
ax2[1].set_ylabel("Force (kN)")
ax2[1].legend()
fig.suptitle("Displacement on web and applied force against crosshead displacement")
plt.show()
# Sample 2 (python)