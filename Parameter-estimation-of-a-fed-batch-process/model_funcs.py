# implements the kinetic model

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp

def model_rhs(t, y, p, c): 

      """ r.h.s. (right hand side) function of the ODE model 

      INPUT:
      t ... current time [h]
      y ... state vector:
            y[0] ... substrate mass (mS) [g]
            y[1] ... bio dry mass (mX) [g]
            y[2] ... ethanol mass (mE) [g]
            y[3] ... volume of fermentation broth [L]
      p ... structure with parameter values to be estimated, cf. lmfit.Parameters
      c ... dict of control values
            c["feed_on"] ... time point when feed was switched on [h]
            c["feed_rate"] ... feed rate [L/h]
            c["csf"] ... substrate concentration in feed [g/L]
            c["M_base"] ... base concentration [mol/L]
            c["gas_flow"] ... aeration rate [L/h]
            c["T"] ... temperature in reactor [°C]
            c["wf_on"] ... weighting factor for online data [dimensionless]
            c["wf_CO2"] ... weighting factor for CO2 data [dimensionless]
            c["pressure"] ... mean pressure in CO2 sensor [bar]
            
      OUTPUT:
      dy_dt ... time derivative of state vector
      """
      #fix parameters 
      dV_gas_dt = c["gas_flow"]     # air flow in each experiment L/h
      R = 0.08314 #bar*l/mol*K
      T = c["T"] + 273.15     #Kelvin  
      pressure = c["pressure"] #bar
      M_CO2 = 44.01     #g/mol

      # (potential) fit parameters
      qsmax = p["qsmax"].value
      qemax = p["qemax"].value
      Ks    = p["Ks"].value
      Ke = p["Ke"].value
      Ki = p["Ki"].value
      

      qO2max = p["qO2max"].value
      qm_max = p["qm_max"].value
      g_e = p["g_e"].value
      #YCO2s_m = 1.4657061593696088
      YCO2s_m = p["YCO2s_m"].value
      
      Yxs_ox = p["Yxs_ox"].value
      Yxs_red = p["Yxs_red"].value
      Yxe_et = p["Yxe_et"].value
      Yxg_glyc = p["Yxg_glyc"].value   

      Yes_red = p["Yes_red"].value
      Ygs_red = p["Ygs_red"].value
      
      #could be used in other models
      # mumaxE = p['mumaxE'].value
      # YCO2x_ox = p["YCO2x_ox"].value
      # YCO2x_red = p["YCO2x_red"].value
      # YCO2x_et = p["YCO2x_et"].value
      # YCO2x_glyc = p["YCO2x_glyc"].value

      YCO2s_ox = p["YCO2s_ox"].value
      YCO2s_red = p["YCO2s_red"].value
      YCO2e_et = p["YCO2e_et"].value
      YCO2g_glyc = p["YCO2g_glyc"].value

      YO2s_ox = p["YO2s_ox"].value
      YO2e_et = p["YO2e_et"].value
      YO2g_glyc = p["YO2g_glyc"].value
       
      # controls
      feed_on = c["feed_on"] # time point when feed was switched on [h]
      feed_rate = c["feed_rate"] # feed rate [L/h]
      Fin = feed_rate * (t > feed_on) # becomes 0 if t < feed_on
      csf = c["csf"] # substrate concentration in feed [g/L]

      # masses and concentrations
      mS, mX, mE, V = y
      cS, cX, cE = [mS, mX, mE] / V    


      #kinetics
      
      qs = qsmax * cS / (cS + Ks)

      if qs > qm_max:
            qm = qm_max
      else:
            qm = qs

      qs = qsmax * cS / (cS + Ks) - qm

      qO2_s_max = qs*YO2s_ox

      Monod_O2 = 1   #cO/(cO+Ko)
 
      qO2 = qO2max * Monod_O2


      if qO2_s_max <= qO2:
            
            qO2_s = qO2_s_max
            qsOx = qO2_s / YO2s_ox                  
            qsRed = 0
            
            qO2res = qO2 - qO2_s_max                       
            qe_Emax = qemax * cE/(cE + Ke) * Ki / (cS+ Ki)
            qO2_e = qe_Emax * YO2e_et
            qO2_g = g_e * qe_Emax*YO2g_glyc  

            if qO2_e + qO2_g <= qO2res:

                  qe = qe_Emax
            else:
                  qe = qO2res / (YO2e_et + g_e * YO2g_glyc)
                  
            qg = g_e*qe

      else:
            qO2_s = qO2 #before, there was qO2max
            qsOx = qO2_s / YO2s_ox         
            qsRed = qs - qsOx 
            qe = 0.0
            qg = 0.0
       

      muOx = qsOx * Yxs_ox
      muRed = qsRed * Yxs_red
      muE = qe * Yxe_et
      muG = qg * Yxg_glyc

      muTotal = muOx+muRed+muE+muG
      qCO2 = qm*YCO2s_m + qsOx* YCO2s_ox + qsRed* YCO2s_red + qe* YCO2e_et + qg*YCO2g_glyc

      dmS_dt = cX *V*(- qsOx- qsRed)+ csf *Fin    #-qm substraction ?
      dmX_dt = muTotal * cX * V
      dmE_dt = (qsRed * Yes_red - qe)*cX*V
      
      dmCO2_dt = qCO2 * cX * V

      #1
      # dnAbgas_dt = pressure * dV_gas_dt / (R*T)
      # dnCO2_dt = dmCO2_dt/M_CO2
      # CO2_percent = 100 * dnCO2_dt / dnAbgas_dt

      #2 equivalent to #1
      dnCO2_dt = dmCO2_dt/M_CO2
      dvCO2_dt = (dnCO2_dt*R*T)/pressure
      CO2_percent = 100 * dvCO2_dt / dV_gas_dt

      dV_dt  = + Fin

      return dmS_dt, dmX_dt, dmE_dt, dV_dt, CO2_percent

    

def wrapper_model_rhs(t, y, p, c):
      """ wrapper function of model_rhs
      Only there to call up the actual model and to return dmS_dt, dmX_dt, dmE_dt, dV_dt since only these should be integrated. 
      """
      dmS_dt, dmX_dt, dmE_dt, dV_dt, CO2_percent = model_rhs(t, y, p, c)
      return dmS_dt, dmX_dt, dmE_dt, dV_dt


def sim_single_exp(t_grid, y0, p, c):  
      """ simulates single experiment and calculates (measured) quantities
    
      INPUT:
      t_grid ... time grid on which to generate the solution [h]
      y0 ... initial state vector:
            y0[0] ... substrate mass (mS) [g]
            y0[1] ... bio dry mass (mX) [g]
            y0[2] ... ethanol mass (mE) [g]
            y0[3] ... volume of fermentation broth [L]
      p ... structure with parameter values to be estimated, cf. lmfit.Parameters
      c ... dict of control values
            c["feed_on"] ... time point when feed was switched on [h]
            c["feed_rate"] ... feed rate [L/h]
            c["csf"] ... substrate concentration in feed [g/L]
            c["M_base"] ... base concentration [mol/L]
            c["gas_flow"] ... aeration rate [L/h]
            c["T"] ... temperature in reactor [°C]
            c["wf_on"] ... weighting factor for online data [dimensionless]
            c["wf_CO2"] ... weighting factor for CO2 data [dimensionless]
            c["pressure"] ... mean pressure in CO2 sensor [bar]      
      OUTPUT:
      sim_exp ... data frame with simulated cS, cX, V, CO2 and base consumption rate over time (for single experiment)
      """
    
      # run ODE solver to get solution y(t)
      y_t = solve_ivp(wrapper_model_rhs, [np.min(t_grid), np.max(t_grid)], y0, t_eval=t_grid, args = (p, c), method= "Radau", first_step = 0.0000001, max_step= 0.1).y.T

      # unpack solution into vectors
      mS, mX, mE, V = [y_t[:,i] for i in range(4)]
      cS, cX, cE = [mS, mX, mE] / V

      # for base consumption rate: get value of dmX_dt at all times t
      dmX_dt = np.array([model_rhs(t_grid[i], y_t[i,:], p, c) for i in range(len(t_grid))])[:,1]
      base_rate = p['base_coef'].value /c["M_base"] * dmX_dt * 1000
      CO2_percent = np.array([model_rhs(t_grid[i], y_t[i,:], p, c) for i in range(len(t_grid))])[:,4]

      # pack everything neatly together to a pandas df
      sim_exp = pd.DataFrame(
            {'t': t_grid,
            'cS': cS, 'cX': cX,
            'V': V, 'base_rate': base_rate,
            'CO2' : CO2_percent, 'cE' : cE }
            ).set_index('t') # make time column the index column

      return sim_exp
