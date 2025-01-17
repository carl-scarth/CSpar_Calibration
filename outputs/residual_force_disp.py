import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import sys

sys.path.append("C:\\Users\\cs2361\\Documents\\Bayesian_Model_Calibration\\source")
from utils import * 

def residuals_plot(exp_data_file, all_exp_data_file, selected_points, disp_str, load, cam_str = "", gp_pred_file=None, gp_mu_file=None, gp_sd_file=None, downsam_rate = 2, obs_err_dt = None):
    # Load in experimental data
    exp_data = pd.read_csv(exp_data_file)
    # Sort into correct format
    # exp_data.columns.values[2] = "point_ind"
    # exp_data["point_ind"] = exp_data["point_ind"].astype("int")
    
    #all_exp_data = pd.read_csv(all_exp_data_file)
    #all_exp_data["Load"] = -all_exp_data["Load"]
    #all_exp_data["Mean DIC Displacement"] = -all_exp_data["Mean DIC Displacement"]
    #all_exp_data = all_exp_data[(all_exp_data["Load"].abs()<=(max_force*1.025)) & (all_exp_data["Mean DIC Displacement"].abs()<=3.5)]
    #all_exp_data = all_exp_data.iloc[[i for i in range(all_exp_data.shape[0]) if i%downsam_rate == 0]]
    
    # Also load GP samples
    #gp_sam = pd.read_csv(gp_pred_file, header=None)
    gp_mu = pd.read_csv(gp_mu_file, header=None)
    gp_sigma = pd.read_csv(gp_sd_file, header=None)

    # Select relevant points and displacement from the experimental data
    if len(cam_str) > 0:
        exp_data = exp_data[(exp_data["point_ind"].isin(selected_points)) & (exp_data["Camera_pair"] == cam_str)].sort_values("Increment")
    else:
        exp_data = exp_data[exp_data["point_ind"].isin(selected_points)].sort_values("Increment")
    increments = exp_data["Increment"].unique()
    
    # Get mean displacement
    disp_exp = np.array([(exp_data[exp_data["Increment"] == increment][disp_str].mean()) for increment in increments])
    disp_exp = -disp_exp

    # Calculate residual
    #gp_sam = - gp_sam
    gp_mu = - gp_mu
    #res = disp_exp.reshape([-1,1]) - gp_sam.to_numpy()
    #res_mu = np.mean(res,axis=1)
    #res_sd = np.std(res,axis=1)

    #interaction = res*gp_sam
    #int_mu = interaction.mean(axis=1).to_numpy()

    # Cheating
    #res_mu = np.abs(res_mu)
    
    # Calculate upper and lower bounds of error bounds given by 2xsigma
    # Taking average residual over all posterior samples to give observation error. Then also adding observation error sd
    # NOT 100% ON THIS (ALSO - NEED TO SQUARE STANDARD DEVIATION FIRST)
    gp_mu = gp_mu.to_numpy().reshape([-1])
    gp_sigma = gp_sigma.to_numpy().reshape([-1])

    # These calcuations were clearly bollocks...
    #ub = gp_mu + 2*np.sqrt(res_sd**2 + gp_sigma**2)
    #lb = gp_mu - 2*np.sqrt(res_sd**2 + gp_sigma**2)
    #ub = (gp_mu+res_mu+2*res_sd+2*gp_sigma).reshape([-1,1])
    #lb = (gp_mu-res_mu-2*res_sd-2*gp_sigma).reshape([-1,1])
    #bounds = np.concatenate((lb,ub), axis=1)
    #lb = np.min(bounds,axis=1)
    #ub = np.max(bounds,axis=1)

    # Bounds looking at E(eta + epsilon) + 2*Var(eta+epsilon)
    # Ignore correlation
    #ub2 = (gp_mu + res_mu) + 2*np.sqrt(gp_sigma**2 + res_sd**2)
    #lb2 = (gp_mu + res_mu) - 2*np.sqrt(gp_sigma**2 + res_sd**2)
    # Using experesssion for residual mean and variance
    #ub2 = (gp_mu + res_mu) + 2*np.sqrt(gp_sigma**2 + res_sd**2 +2*int_mu - 2*res_mu*gp_mu)
    #lb2 = (gp_mu + res_mu) - 2*np.sqrt(gp_sigma**2 + res_sd**2 +2*int_mu - 2*res_mu*gp_mu)
    #ub2[np.isnan(ub2)] = (gp_mu[np.isnan(ub2)] + res_mu[np.isnan(ub2)])
    #lb2[np.isnan(lb2)] = (gp_mu[np.isnan(lb2)] + res_mu[np.isnan(lb2)])

    # Bounds loking at E(eta)
    #print(gp_mu)
    # These are identical to the GP_sd I already output
    #ub3 = (disp_exp-res_mu) + 2*res_sd
    #lb3 = (disp_exp-res_mu) - 2*res_sd

    # (I think it should actually be plus and minus the res so centred around experimental data - fix when better dataset)
    #print(pd.DataFrame(np.concatenate((gp_mu.reshape([-1,1]), gp_sigma.reshape([-1,1]), res_mu.reshape([-1, 1]), res_sd.reshape([-1, 1]), lb.reshape([-1, 1]), ub.reshape([-1, 1])),axis=1),columns=["GP_mean","GP_sd","res_mean","res_sd","lb","ub"]))
    
    fig, ax = plt.subplots()
    # ax.plot(-gp_sam.to_numpy(),applied_load)
    #for sam in obs_err_dt:
    #    ax.plot(gp_mu-2*(np.sqrt(gp_sigma**2 + sam**2)),load,'-g',linewidth=0.25)
    #    ax.plot(gp_mu+2*(np.sqrt(gp_sigma**2 + sam**2)),load,'-g',linewidth=0.25)
    ax.plot(gp_mu-2*(np.sqrt(gp_sigma**2 + np.mean(obs_err_dt)**2)), load,'-g',linewidth=1.5)
    ax.plot(gp_mu+2*(np.sqrt(gp_sigma**2 + np.mean(obs_err_dt)**2)), load,'-g',linewidth=1.5)
    ax.plot(gp_mu,load,'-k',linewidth=2.0)
    ax.plot(gp_mu-2*gp_sigma,load,'-c',linewidth=1.5)
    ax.plot(gp_mu+2*gp_sigma,load,'-c',linewidth=1.5)
    #ax.plot(ub,load,'--m',linewidth=1.0)
    #ax.plot(lb,load,'--m',linewidth=1.0)
    #ax.plot(ub2,load,'--g',linewidth=1.0)
    #ax.plot(lb2,load,'--g',linewidth=1.0)
    #ax.plot(ub3,load,'--g',linewidth=1.0)
    #ax.plot(lb3,load,'--g',linewidth=1.0)
    # ax.plot(all_exp_data["Mean DIC Displacement"].to_numpy(),all_exp_data["Load"].to_numpy(),'bx')
    ax.plot(disp_exp, load,'rx',markersize=10.0, markeredgewidth =2.0)
    ax.set_ylabel("Force (kN)")
    ax.set_xlabel("Displacement (mm)")

def residuals_plot_2cam(exp_data_file, selected_points, disp_str, load, cam_str = [], prior_sam_file = None, posterior_sam_file = None, gp_pred_file=None, gp_mu_file=None, gp_sd_file=None, obs_err_dt = None):
    # Load in experimental data
    exp_data = pd.read_csv(exp_data_file)
    # Sort into correct format
    # exp_data.columns.values[2] = "point_ind"
    # exp_data["point_ind"] = exp_data["point_ind"].astype("int")
    

    # Also load GP samples
    #gp_sam = pd.read_csv(gp_pred_file, header=None)
    gp_mu = pd.read_csv(gp_mu_file, header=None)
    gp_sigma = pd.read_csv(gp_sd_file, header=None)

    # Select relevant points and displacement from the experimental data
    exp_disp = {}
    for cam in cam_str:
        exp_data_subset = exp_data[(exp_data["point_ind"].isin(selected_points[cam])) & (exp_data["Camera_pair"] == cam)].sort_values("Increment")
        increments = exp_data_subset["Increment"].unique()
        # Get mean displacement
        disp_exp = np.array([(exp_data_subset[exp_data_subset["Increment"] == increment][disp_str].mean()) for increment in increments])
        exp_disp[cam] = -disp_exp

    print(exp_disp)
    #gp_sam = - gp_sam
    gp_mu = - gp_mu
    gp_mu = gp_mu.to_numpy().reshape([-1])
    gp_sigma = gp_sigma.to_numpy().reshape([-1])
    
    fig, ax = plt.subplots()
    # ax.plot(-gp_sam.to_numpy(),applied_load)
    #for sam in obs_err_dt:
    #    ax.plot(gp_mu-2*(np.sqrt(gp_sigma**2 + sam**2)),load,'-g',linewidth=0.25)
    #    ax.plot(gp_mu+2*(np.sqrt(gp_sigma**2 + sam**2)),load,'-g',linewidth=0.25)
    obs_err_col = ['-r','--b']
    dt_col = ['rx', 'bx']
    for i, cam in enumerate(cam_str):
        ax.plot(gp_mu-2*(np.sqrt(gp_sigma**2 + np.mean(obs_err_dt[cam])**2)), load,obs_err_col[i],linewidth=1.5)
        ax.plot(gp_mu+2*(np.sqrt(gp_sigma**2 + np.mean(obs_err_dt[cam])**2)), load,obs_err_col[i],linewidth=1.5)
    ax.plot(gp_mu,load,'-k',linewidth=2.0)
    ax.plot(gp_mu-2*gp_sigma,load,'-c',linewidth=1.5)
    ax.plot(gp_mu+2*gp_sigma,load,'-c',linewidth=1.5)
    #ax.plot(ub,load,'--m',linewidth=1.0)
    #ax.plot(lb,load,'--m',linewidth=1.0)
    #ax.plot(ub2,load,'--g',linewidth=1.0)
    #ax.plot(lb2,load,'--g',linewidth=1.0)
    #ax.plot(ub3,load,'--g',linewidth=1.0)
    #ax.plot(lb3,load,'--g',linewidth=1.0)
    # ax.plot(all_exp_data["Mean DIC Displacement"].to_numpy(),all_exp_data["Load"].to_numpy(),'bx')
    for i, cam in enumerate(cam_str):
        ax.plot(exp_disp[cam], load, dt_col[i],markersize=10.0, markeredgewidth =2.0)
    ax.set_ylabel("Force (kN)")
    ax.set_xlabel("Displacement (mm)")

    # FINE TO CRACK ON WITH 5E-2. CHECK DIC INSPIRED PRIORS OUT OF CURIOSITY BUT LIKELY WORSE
    # THEN DO RESIDUALS PLOT
    # REMEMBER I RERAN NOMINAL PROPERTIES TO GET ABAQUS OUTPUT FOR THIS - USE INTERPOLATE DATA 
    # CODE TO GET PRIOR PLOTS

def prior_post_residuals_2cam(exp_data_file, selected_points, disp_str, load, cam_str = [], prior_sam = None, posterior_sam = None, posterior_mu = None, posterior_sd = None, obs_err_dt = None):
    # Load in experimental data

    exp_data = pd.read_csv(exp_data_file)
    # Sort into correct format
    # exp_data.columns.values[2] = "point_ind"
    # exp_data["point_ind"] = exp_data["point_ind"].astype("int")
    if posterior_mu is None:
        posterior_mu = posterior_sam.mean(axis=1).values
        posterior_sd = posterior_sam.std(axis=1).values
        posterior_lb = posterior_sam.quantile(0.05,axis=1).values
        posterior_ub = posterior_sam.quantile(0.95,axis=1).values
        posterior_sam = posterior_sam.values
    
    prior_mu = prior_sam.mean(axis=1).values
    prior_sd = prior_sam.std(axis=1).values
    prior_ub = prior_sam.quantile(0.95, axis=1).values
    prior_lb = prior_sam.quantile(0.05, axis=1).values
    prior_sam = prior_sam.values
    
    # Select relevant points and displacement from the experimental data
    exp_disp = {}
    for cam in cam_str:
        exp_data_subset = exp_data[(exp_data["point_ind"].isin(selected_points[cam])) & (exp_data["Camera_pair"] == cam)].sort_values("Increment")
        increments = exp_data_subset["Increment"].unique()
        # Get mean displacement
        disp_exp = np.array([(exp_data_subset[exp_data_subset["Increment"] == increment][disp_str].mean()) for increment in increments])
        exp_disp[cam] = -disp_exp


    #gp_sam = - gp_sam
    # gp_mu = - gp_mu
    # gp_mu = gp_mu.to_numpy().reshape([-1])
    # gp_sigma = gp_sigma.to_numpy().reshape([-1])
    
    fig, ax = plt.subplots()
    # ax.plot(-gp_sam.to_numpy(),applied_load)
    #for sam in obs_err_dt:
    #    ax.plot(gp_mu-2*(np.sqrt(gp_sigma**2 + sam**2)),load,'-g',linewidth=0.25)
    #    ax.plot(gp_mu+2*(np.sqrt(gp_sigma**2 + sam**2)),load,'-g',linewidth=0.25)
    obs_err_col = ['g','m']
    dt_col = ['gx', 'mx']
   # for i, cam in enumerate(cam_str):
   # #    ax.plot(gp_mu-2*(np.sqrt(gp_sigma**2 + np.mean(obs_err_dt[cam])**2)), load,obs_err_col[i],linewidth=1.5)
   #     ax.plot(gp_mu+2*(np.sqrt(gp_sigma**2 + np.mean(obs_err_dt[cam])**2)), load,obs_err_col[i],linewidth=1.5)
   # plt.show()
    ax.plot(prior_mu,load,'-', color = "red", linewidth = 1.5)
   # ax.plot(prior_mu+2*prior_sd, prior_force, '-r', linewidth = 2.0)
    #ax.plot(prior_mu-2*prior_sd, prior_force, '-r', linewidth = 2.0)
    ax.plot(prior_lb, load, '-r', linewidth = 0.25)
    ax.plot(prior_ub, load, '-r', linewidth = 0.25)
    ax.fill_betweenx(load, prior_lb, prior_ub, color=(1,0.75,0.75,1))
    for i, cam in enumerate(cam_str):
        ax.errorbar(exp_disp[cam], load, xerr=2*np.mean(obs_err_dt[cam]), fmt="none", capsize = 5.0, markersize=10.0, elinewidth=1.0,color=obs_err_col[i])#,markeredgewidth=2.0)
        #ax.errorbar(posterior_mu.reshape(-1), load, xerr=2*np.mean(obs_err_dt[cam]), fmt=dt_col[i], capsize = 5.0, markersize=10.0, elinewidth=1.0)
    ax.plot(posterior_mu,load,'-',color='b',linewidth=2.0)
    ax.plot(posterior_mu-2*posterior_sd,load,'-b',linewidth=0.5)
    ax.plot(posterior_mu+2*posterior_sd,load,'-b',linewidth=0.5)
    ax.fill_betweenx(load, posterior_mu-2*posterior_sd, posterior_mu+2*posterior_sd, color=(0.75,0.75,1,1))
    for i, cam in enumerate(cam_str):
        ax.plot(exp_disp[cam], load, dt_col[i],markersize=10.0, markeredgewidth =2.0)
    ax.set_ylabel("Force (kN)")
    ax.set_xlabel("Displacement (mm)")
    ax.set_ylim((0, 150))
    ax.set_yticks([0,25,50,75,100,125,150])
    ax.set_xlim((-0.03,1.4)) # for w
    #ax.set_xlim((-0.3,7.5)) # for u

if __name__ == "__main__":
    set_plot_params()


    # STICK WITH 1E-2 FOR NOW - EXPERIMENTAL DIC WAS WORSE - THOUGH W AND U MORE COMPARABLE.
    # MAYBE CONSIDER ON RP_CALIBRATION PLOT CODE AS JUST LOADING MC SAMPLES AND ADDING ERROR
    # BARS
    # Also prior MC samples for error bars on prior plot? ABAQUS or emulator predictions
    # Calculates sample residuals force-displacement
    gp_pred_file = "E:\\Calibration_outputs_for_paper\\Predictions_obserr5e-2\\LHSDesign100x8_max_eta_sam_w_410.csv"
    gp_pred_file_u = "E:\\Calibration_outputs_for_paper\\Predictions_obserr5e-2\\LHSDesign100x8_max_eta_sam_u.csv"
    exp_data_file = "..\\inputs\\Interpolated_DIC_multistep_150kN_combined.csv"
    gp_mu_file = "E:\\Calibration_outputs_for_paper\\Predictions_obserr5e-2\\LHSDesign100x8_max_eta_sam_mu_w_DIC.csv"   
    gp_mu_file_u = "E:\\Calibration_outputs_for_paper\\Predictions_obserr5e-2\\LHSDesign100x8_max_eta_sam_mu_u.csv"   
    gp_sd_file = "E:\\Calibration_outputs_for_paper\\Predictions_obserr5e-2\\LHSDesign100x8_max_eta_sam_sigma_w_DIC.csv"
    gp_sd_file_u = "E:\\Calibration_outputs_for_paper\\Predictions_obserr5e-2\\LHSDesign100x8_max_eta_sam_sigma_u.csv"   
    all_exp_file = "..\\inputs\\CS02P_mean_w_tip_"
    all_exp_file_u = "..\\inputs\\CS02P_mean_u_mid_"
    prior_sam_file = "E:\\Calibration_outputs_for_paper\\LHSDesign250x8_1_interp_w_tip.csv"
    prior_sam_file_u = "E:\\Calibration_outputs_for_paper\\LHSDesign250x8_1_u_max_displacements.csv"

    obs_err_file = "E:\\Calibration_outputs_for_paper\\Predictions_obserr5e-2\\LHSDesign100x8_observation_error_samples.csv"
    # Load observation error sample
    obs_err = pd.read_csv(obs_err_file).to_numpy()

    #selected_points = np.array([26805, 2483, 34982, 8697, 32996])
    # DOUBLE CHECK THESE
    point_subset_w = {"LC" : [26805, 34982, 32996], "RC" : [12734 ,33496, 18927]}
    point_subset_u = {"LC" : [22841, 23197, 23275, 23294, 23717], "RC" : [8513, 8529, 8566, 15030, 28731]}
    obs_err_col_w = {"LC" : obs_err[:,1], "RC" :  obs_err[:,3]}
    obs_err_col_u = {"LC" : obs_err[:,0], "RC" :  obs_err[:,2]}

    # applied_load = np.linspace(0.0, 200.0, 17)
    applied_load = np.concatenate((np.array([5.39]), np.linspace(10.0, 150.0, 15)))
    max_force = 150.0
    #for cam_pair, points in point_subset_w.items():
    #    residuals_plot(exp_data_file, all_exp_file+cam_pair+".csv", points, "w_rot", applied_load, cam_str=cam_pair, gp_pred_file=gp_pred_file, gp_mu_file=gp_mu_file, gp_sd_file=gp_sd_file, obs_err_dt=obs_err_col_w[cam_pair])

    #for cam_pair, points in point_subset_u.items():
    #    residuals_plot(exp_data_file, all_exp_file_u+cam_pair+".csv", points, "u_rot", applied_load, cam_str=cam_pair, gp_pred_file=gp_pred_file_u, gp_mu_file=gp_mu_file_u, gp_sd_file=gp_sd_file_u, obs_err_dt=obs_err_col_u[cam_pair])
    
    prior_sam = pd.read_csv(prior_sam_file)
    prior_sam_u = pd.read_csv(prior_sam_file_u,header=None)
    prior_sam_u = -prior_sam_u.T

    posterior_sam = -pd.read_csv(gp_pred_file,header=None)

    posterior_mu = -pd.read_csv(gp_mu_file, header = None).values.reshape([-1])
    posterior_sd = pd.read_csv(gp_sd_file, header = None).values.reshape([-1])
    posterior_mu_u = -pd.read_csv(gp_mu_file_u, header=None).values.reshape([-1])
    posterior_sd_u = pd.read_csv(gp_sd_file_u, header = None).values.reshape([-1])
    prior_post_residuals_2cam(exp_data_file, point_subset_w, "w_rot", applied_load, cam_str = [key for key in point_subset_w.keys()], prior_sam = -prior_sam[["disp_"+str(i) for i in range(prior_sam.shape[1]//3)]], posterior_mu = posterior_mu, posterior_sd=posterior_sd, obs_err_dt = obs_err_col_w)
    prior_post_residuals_2cam(exp_data_file, point_subset_u, "u_rot", applied_load, cam_str = [key for key in point_subset_u.keys()], prior_sam = prior_sam_u, posterior_mu = posterior_mu_u, posterior_sd=posterior_sd_u, obs_err_dt = obs_err_col_u)
    plt.show()
    ssdfsdf
    print(posterior_sam)
    posterior_mu = posterior_sam.mean(axis=1).values    
    posterior_sd = posterior_sam.std(axis=1).values
    posterior_lb = posterior_sam.quantile(0.05,axis=1).values
    posterior_ub = posterior_sam.quantile(0.95,axis=1).values
    posterior_sam = posterior_sam.values
    prior_force = prior_sam["force_0"].values
    prior_sam = -prior_sam[["disp_"+str(i) for i in range(prior_sam.shape[1]//3)]]
    prior_mu = prior_sam.mean(axis=1).values
    prior_sd = prior_sam.std(axis=1).values
    prior_ub = prior_sam.quantile(0.95, axis=1).values
    prior_lb = prior_sam.quantile(0.05, axis=1).values
    prior_sam = prior_sam.values
    fig, ax = plt.subplots()
    print(prior_sd)
    #ax.plot(prior_sam, prior_force, linewidth = 0.25)
    ax.plot(prior_mu, prior_force, '-b', linewidth = 2.0)
   # ax.plot(prior_mu+2*prior_sd, prior_force, '-r', linewidth = 2.0)
    #ax.plot(prior_mu-2*prior_sd, prior_force, '-r', linewidth = 2.0)
    ax.plot(prior_lb, prior_force, '-r', linewidth = 2.0)
    ax.fill_betweenx(prior_force, prior_lb, prior_ub, color=(1,0.75,0.75,1))
    ax.plot(prior_ub, prior_force, '-r', linewidth = 2.0)
    ax.plot(posterior_sam, prior_force, linewidth = 0.25)
    ax.plot(posterior_mu+2*posterior_sd, prior_force,'-g', linewidth = 2.0)
    ax.plot(posterior_mu-2*posterior_sd, prior_force,'-g', linewidth = 2.0)
    ax.fill_betweenx(prior_force, posterior_lb, posterior_ub, color=(0.75,0.75,1,1))
    ax.plot(posterior_lb, prior_force, '-c', linewidth = 2.0)
    ax.plot(posterior_ub, prior_force, '-c', linewidth = 2.0)
    ax.plot(posterior_mu, prior_force, '-m', linewidth=2.0)
    ax.set_ylim((0.0,150.0))
    
    
    plt.show()
    asdsad

    residuals_plot_2cam(exp_data_file, point_subset_w, "w_rot", applied_load, cam_str= [key for key in point_subset_w.keys()], gp_pred_file=gp_pred_file, gp_mu_file=gp_mu_file, gp_sd_file=gp_sd_file, obs_err_dt = obs_err_col_w)
    residuals_plot_2cam(exp_data_file, point_subset_u, "u_rot", applied_load, cam_str= [key for key in point_subset_u.keys()], gp_pred_file=gp_pred_file_u, gp_mu_file=gp_mu_file_u, gp_sd_file=gp_sd_file_u, obs_err_dt = obs_err_col_u)
    plt.show()