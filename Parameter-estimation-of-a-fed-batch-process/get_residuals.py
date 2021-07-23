#code to get residuals and simulated values on t_grid
import numpy as np
import pandas as pd
from lmfit import minimize, Parameters
from model_funcs import sim_single_exp


def get_residuals(p, y0_dict, c_dict, datasets_dict):


    exp_names = y0_dict.keys() # experiment names

    res_and_data= {} # will contain residuals

    for exp in exp_names: # loop over experiments
            y0 = y0_dict[exp]
            c = c_dict[exp]     
            datasets = datasets_dict[exp]

            residuals_list = residuals_single_exp2(p, c, y0, datasets)
            res_and_data[exp] = residuals_list

    return res_and_data



def residuals_single_exp2(p, c, y0, datasets):

    residuals_list = []  

    for dat in datasets: # loop over datasets
        residuals_dict = {}
        t_grid = dat.index.values  # index of "dat" = time grid of measurements = time grid for simulation   
        sim_exp = sim_single_exp(t_grid, y0, p, c) # simulate experiment with this time grid
            
        for var in dat: # loop over all measured variables
            res_var = (sim_exp[var] - dat[var]).values # weighted residuals for this measured variable
            residuals_dict[var] = dat[var].values
            residuals_dict["{0}_fitted".format(var)] = sim_exp[var].values
            residuals_dict["{0}_residuals".format(var)] = res_var
            residuals_dict["{0}_abs_residuals".format(var)] = abs(res_var)
        
        
        residuals_df = pd.DataFrame (residuals_dict, index = t_grid)
        residuals_list.append(residuals_df)

    return residuals_list