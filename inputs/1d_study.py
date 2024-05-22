import numpy as np

# Parametric study investigating rotational spring stiffness behaviour in Abaqus
E11 = 140.9
t_ply = 0.125

N = 40
K = np.linspace(0.0,16.0,N).reshape([N,1])
# K = np.exp(K)
x_out = np.concatenate((np.full([N,1],E11),np.full([N,1],t_ply),K),1)
np.savetxt('translator_study.csv', x_out, delimiter=",", header = 'E_11,t_ply,log_K_rig', comments = "")