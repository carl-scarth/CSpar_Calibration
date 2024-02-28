import pandas as pd
import matplotlib.pyplot as plt
import sys

sys.path.append("C:\\Users\\cs2361\\Documents\\Bayesian_Model_Calibration\\source")
from utils import * 

set_plot_params()

infile = "LHSDesign60x6_2"
exp_file = "CS02P_mean_disp_412.5mm"
model_fd = pd.read_csv(infile+"_interp.csv")
exp_fd = pd.read_csv(exp_file+".csv")


fig, ax = plt.subplots()
for i in range(model_fd.shape[1]//2):
    if i == 0:
        plt.plot(-model_fd["disp_"+str(i)].to_numpy(), model_fd["force_"+str(i)].to_numpy(), "-b", linewidth=1.0, label = "Prior sample")
    else:
        plt.plot(-model_fd["disp_"+str(i)].to_numpy(), model_fd["force_"+str(i)].to_numpy(), "-b", linewidth=1.0)
ax.plot(-exp_fd["Mean DIC Displacement"].to_numpy(), -exp_fd["Load"].to_numpy(), "rx", label = "DIC")
ax.set_title("Prior Force-displacement vs DIC")
ax.set_xlabel("Displacement (mm)")
ax.set_ylabel("Force (kN)")
ax.legend()
plt.show()