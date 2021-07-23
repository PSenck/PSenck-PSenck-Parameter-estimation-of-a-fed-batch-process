# PSenck-PSenck-Parameter-estimation-of-a-fed-batch-process

This is a toolbox to estimate model parameters of a fed-batch process. The toolbox works if you have all modules in one folder. 
You only need to adjust the paths to load the attached data in the measurement data folder (Messdaten). 
For this, paths for the metadata, the online data, the off-line data and the CO2 data have to be adapted in control_parameter_estimation_jupyter and synthetic_data. 
The module control_parameter_estimation_jupyter performs the parameter estimation with real measurement data and the module synthetic_data performs a parameter estimation with synthetic data. 
The module model_funcs contains the kinetic model and parest_funcs contains the parameter estimation routines.
The module load_and_process contains functions to load raw data and to pre-process them accordingly.
The function get_residuals calculates the residuals between measurement data and simulated data in order to perform a subsequent statistical analysis.
