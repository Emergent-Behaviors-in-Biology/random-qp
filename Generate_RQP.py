import time
import pandas as pd
import matplotlib
from matplotlib import cm
from matplotlib import colors
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import numpy as np
import itertools
from Eco_function.eco_lib import *
from Eco_function.eco_plot import *
from Eco_function.eco_func import *
from Eco_function.Model_cavity import *
import pdb
import os.path
import pickle
from scipy.integrate import odeint

# Output results columns
columns=['mean_phiN','mean_phiR','mean_N', 'q_N','mean_R','q_R', 'survive','mu','sigc','m','sigm','K','sigK','S','M','gamma','N_bar','qN_bar','R_bar','qR_bar', 'survive_bar','phiN_bar','phiR_bar', 'consumed_power','power_bar','opt_f','opt_f_bar']

# Initializa the parameters
parameters = {}

parameters['sample_size']=50;

parameters['S'] =100;
parameters['M']=100;

parameters['K']=1;
parameters['sigma_K']=1.0;

parameters['mu']=1.;
parameters['sigma_c']=1.0; 

parameters['m']=1.;
parameters['sigma_m']=0.1;

# Initializa the parameters for ODE solver
parameters['t0']=0;
parameters['t1']=6000;
parameters['Nt']=10000;
parameters['loop_size']=10;


# Initializa the filename and varaibles range
file_name='RQP_sigc'+str(parameters['sigma_c'])+'.csv'
gamma_range=np.linspace(0.1, 10., num=20)


# pass variables to simulations
def func_simulation(para):
    parameters['sample_size']=para[0];
    parameters['S'] =para[1];
    parameters['M']=para[2];

    parameters['K']=para[3];
    parameters['sigma_K']=para[4];

    parameters['mu']=para[5];
    parameters['sigma_c']=para[6]; 

    parameters['m']=para[7];
    parameters['sigma_m']=para[8];
    parameters['loop_size']=para[9];


    parameters['t0']=para[10];
    parameters['t1']=para[11];
    parameters['Nt']=para[12];
    index=para[13]
    para_df = pd.DataFrame(columns=columns) 
    Model=Cavity_simulation(parameters)
    Model.gamma_flag='S/M'
    if Model.gamma_flag=='M/S':
            gamma=float(parameters['M'])/float(parameters['S'])
    if Model.gamma_flag=='S/M':
            gamma=float(parameters['S'])/float(parameters['M'])
    Model.gamma_flag='S/M'
    Model.initialize_random_variable()
    mean_var=Model.ode_simulation(Dynamics='quadratic', Simulation_type='QP')
    para_df_current = pd.DataFrame([[mean_var['phi_N'],mean_var['phi_R'],mean_var['mean_N'], mean_var['q_N'],mean_var['mean_R'], mean_var['q_R'], mean_var['Survive']/parameters['S'],parameters['mu'],\
                                parameters['sigma_c'],parameters['m'],parameters['sigma_m'],parameters['K'], parameters['sigma_K'], parameters['S'], parameters['M'],gamma,mean_var['mean_N_bar'], mean_var['q_N_bar'],\
                                mean_var['mean_R_bar'],mean_var['q_R_bar'],mean_var['Survive_bar']/parameters['S'],mean_var['phi_N_bar'],mean_var['phi_R_bar'],mean_var['power'],mean_var['power_bar'],mean_var['opti_f'],mean_var['opti_f_bar'] ]],columns=columns)
    para_df =pd.concat([para_df,para_df_current],ignore_index=True) 
    return para_df

 # Save differnet varibles as a job set 
index=0;
jobs=[];
for par in gamma_range:  
    parameters['S']=int(parameters['M']*par)
    jobs.append([parameters['sample_size'],parameters['S'],parameters['M'],parameters['K'],parameters['sigma_K'], parameters['mu'], parameters['sigma_c'],parameters['m'],parameters['sigma_m'],parameters['loop_size'],parameters['t0'],parameters['t1'],parameters['Nt']  ,index])
    index=index+1

# Input varibles set into simulations
results_df=pd.DataFrame(columns=columns) 
for job in jobs:
    df=func_simulation(job)
    results_df =pd.concat([results_df,df],ignore_index=True) 

# Save results
with open(file_name, 'a') as f:
        results_df.to_csv(f, index=False,encoding='utf-8')