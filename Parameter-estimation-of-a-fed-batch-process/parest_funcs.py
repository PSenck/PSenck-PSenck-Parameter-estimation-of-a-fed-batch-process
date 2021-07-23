# implements the parameter estimation routines
import numpy as np
import pandas as pd
from lmfit import minimize, Parameters
from model_funcs import sim_single_exp

import sympy as sp
from sympy import Eq
from sympy import symbols




def residuals_single_exp(p, c, y0, datasets):
      """ Calculate residuals for a single experiment, 
      INPUT:
      p ... structure with parameter values to be estimated, cf. lmfit.Parameters
      c ... list of control values, to be passed to model function
            c["feed_on"] ... time point when feed was switched on [h]
            c["feed_rate"] ... feed rate [L/h]
            c["csf"] ... substrate concentration in feed [g/L]
            c["M_base"] ... base concentration [mol/L]
            c["gas_flow"] ... aeration rate [L/h]
            c["T"] ... temperature in reactor [°C]
            c["wf_on"] ... weighting factor for online data [dimensionless]
            c["wf_CO2"] ... weighting factor for CO2 data [dimensionless]
            c["pressure"] ... mean pressure in CO2 sensor [bar]

      y0 ... initial state vector:
            y0[0] ... substrate mass (mS) [g]
            y0[1] ... bio dry mass (mX) [g]
            y0[2] ... ethanol mass (mE) [g]
            y0[3] ... volume of fermentation broth [L]
      datasets ... list of data frames with measurement data(online,offline,CO2)

      OUTPUT:
      res ... long vector with all residuals for this experiment
      """ 

      res_single = np.array([]) # empty array, will contain residuals

      weighting_factor = {"cX": 1.0, "cS": 1.0, "cE" : 1.0, "base_rate": 1.0*c["wf_on"], "CO2" : 1.0*c["wf_CO2"]} # individual weighting factor c[5]len(off/on) c[6]len(off/CO2)

      for dat in datasets: # loop over datasets
            t_grid = dat.index.values  # index of "dat" = time grid of measurements = time grid for simulation   
            sim_exp = sim_single_exp(t_grid, y0, p, c) # simulate experiment with this time grid

            for var in dat: # loop over all measured variables
                  res_var = weighting_factor[var]*(sim_exp[var] - dat[var]).values # weighted residuals for this measured variable
                  res_single = np.append(res_single, res_var) # append to long residual vector

      return res_single

def residuals_all_exp(p, y0_dict, c_dict, datasets_dict):

      """ Calcule different Yield coefficients, determined by chemical reactions solved by linear equation systems.
      Calculate residuals for all experiment
      INPUT:
      p ... structure with parameter values to be estimated, cf. lmfit.Parameters
      y0_dict ... dict: keys: experiment names, values: initial state vector y0
      c_dict ... dict: keys: experiment names, values: dict with control variables
      datasets_dict ... dictionary: keys: experiment names, values: list of data frame with measurement data   

      OUTPUT:
      Yield coefficients required in the model are added to the lmfit.paramater object 
      res ... super long vector with all residuals for all experiment
      """

      #Code for chemical balancing 

      #Order: C,H,O,N
      gluc = np.array([6.0,12.0,6.0,0.0])
      O2 = np.array([0.0, 0.0, 2.0, 0.0])
      NH3 = np.array([0.0,3.0,0.0,1.0])
      biomass = np.array([1.0,p["HX"].value, p["OX"].value, p["NX"].value])
      CO2 = np.array([1.0,0.0,2.0,0.0])
      H2O = np.array([0.0,2.0,1.0,0.0])
      etoh = np.array([2.0,6.0,1.0,0.0])
      glyc = np.array([3.0,8.0,3.0,0.0])

      


      MW_element_dict = {"C": 12.011, "H": 1.0079, "O": 15.999, "N": 14.007}        #molar masses
      molecule = {"gluc": gluc, "O2": O2, "NH3" : NH3, "biomass": biomass, "CO2" : CO2, "H2O":  H2O, "etoh": etoh, "glyc": glyc}

      MW = {}           #creating dict with masses of molecules
      for key, mol in molecule.items():
            molecule_MW_array = ([])
            for vectorvalue, weight in zip (mol, MW_element_dict.values()):
                  vw = vectorvalue*weight
                  molecule_MW_array= np.append(molecule_MW_array, vw)
            MW[key] = sum(molecule_MW_array)

      NX1 = p["NX"].value
      GE = (p["g_e"]/MW["glyc"]) * MW["etoh"]   #from mass ratio(p["g_e"]) to molar ratio. Glycerol per ethanol

      #1. oxidative glucose consumption: gluc+ a*O2 + b*NX*NH3 = b*biomass + c*CO2 + d*H2O 
      a,b,c,d, NX = symbols("a b c d NX")
      Yxs_ox = p["Yxs_ox"].value
      b1 = Yxs_ox* MW["gluc"]/MW["biomass"]     #calculate stoichiometric coefficient

      eqOx_list = []
      for num in range(3):
            eqOx = Eq(gluc[num]+ a*O2[num]+ b*NX*NH3[num], b*biomass[num]+ c*CO2[num]+ d*H2O[num])
            eqOx = eqOx.subs({b: b1, NX: NX1})
            eqOx_list.append(eqOx)
      
      solution_Ox = sp.solve(eqOx_list, (a, c, d), dict= True)
      a1, c1, d1 = np.float(solution_Ox[0][a]), np.float(solution_Ox[0][c]), np.float(solution_Ox[0][d])
      YCO2x_ox = c1/b1 * MW["CO2"]/MW["biomass"]
      YCO2s_ox = c1/1 * MW["CO2"]/MW["gluc"]
      YO2s_ox = a1/1 * MW["O2"]/MW["gluc"]

      #Yield coefficients showing up later in fit results. Not actual changeable parameters. Vary have to be False
      p.add("YCO2x_ox", value=YCO2x_ox, vary=False)
      p.add("YCO2s_ox", value=YCO2s_ox, vary=False)
      p.add("YO2s_ox", value=YO2s_ox, vary=False)
      
      # stop

      #2. reductive glucose consumption:  gluc+ g*NX*NH3 = g*biomass + h*CO2 + i*H2O + j*etOH + GE*j*glyc
      g,h,i,j, NX = symbols("g h i j NX")
      Yxs_red = p["Yxs_red"].value
      g1 = Yxs_red* MW["gluc"]/MW["biomass"]

      eqRed_list = []
      for num in range(3):
            eqRed = Eq(gluc[num]+ g*NX*NH3[num] , g*biomass[num]+ h*CO2[num]+ i*H2O[num]+ j*etoh[num]+ GE*j*glyc[num] )
            eqRed = eqRed.subs({g: g1, NX: NX1})
            eqRed_list.append(eqRed)
      
      solution_Red = sp.solve(eqRed_list, (h, i, j), dict= True)
      h1,i1,j1 = np.float(solution_Red[0][h]), np.float(solution_Red[0][i]), np.float(solution_Red[0][j])

      Yes_red = j1/1 * MW["etoh"]/MW["gluc"]
      Ygs_red = GE*j1/1 * MW["glyc"]/MW["gluc"]   
      YCO2x_red = h1/g1 * MW["CO2"]/MW["biomass"]
      YCO2s_red = h1/1 * MW["CO2"]/MW["gluc"]

      p.add("Yes_red", value=Yes_red, vary=False)
      p.add("Ygs_red", value=Ygs_red, vary=False)
      p.add("YCO2x_red", value=YCO2x_red, vary=False)
      p.add("YCO2s_red", value=YCO2s_red, vary=False)


      #3. oxidative ethanol consumption: etoh + k*O2 + l*NX*NH3 = l*biomass + m*CO2 + n*H2O
      k,l,m,n, NX = symbols("k l m n NX")
      Yxe_et = p["Yxe_et"].value      
      l1 = Yxe_et* MW["etoh"]/MW["biomass"]

      eqEt_list = []
      for num in range(3):
            eqEt = Eq(etoh[num]+ k*O2[num]+ l*NX*NH3[num], + l*biomass[num]+ m*CO2[num]+ n*H2O[num])
            eqEt = eqEt.subs({l: l1, NX: NX1})
            eqEt_list.append(eqEt)
      
      solution_Et = sp.solve(eqEt_list, (k, m, n), dict= True)
      k1, m1, n1 = np.float(solution_Et[0][k]), np.float(solution_Et[0][m]), np.float(solution_Et[0][n])

      YCO2x_et = m1/l1 * MW["CO2"]/MW["biomass"]
      YCO2e_et = m1/1 * MW["CO2"]/MW["etoh"]
      YO2e_et = k1 * MW["O2"]/MW["etoh"]

      p.add("YCO2x_et", value=YCO2x_et, vary=False)
      p.add("YCO2e_et", value=YCO2e_et, vary=False)
      p.add("YO2e_et", value=YO2e_et, vary=False)



      #4. oxidative glycerol consumption: glyc+ u*O2 + v*NX*NH3 = v*biomass + w*CO2 + x*H2O
      u,v,w,x,NX = symbols("u v w x NX")
      Yxg_glyc = p["Yxg_glyc"].value
      v1 = Yxg_glyc * MW["glyc"]/MW["biomass"]

      eqGlyc_list = []
      for num in range(3):
            eqGlyc = Eq(glyc[num]+ u*O2[num]+ v*NX*NH3[num] , v*biomass[num]+ w*CO2[num]+ x*H2O[num])
            eqGlyc = eqGlyc.subs({v: v1, NX: NX1})
            eqGlyc_list.append(eqGlyc)
      
      solution_glyc = sp.solve(eqGlyc_list, (u, w, x), dict= True)
      u1, w1, x1 = np.float(solution_glyc[0][u]), np.float(solution_glyc[0][w]), np.float(solution_glyc[0][x])

      YCO2x_glyc = w1/v1 * MW["CO2"]/MW["biomass"]
      YCO2g_glyc = w1/1 * MW["CO2"]/MW["glyc"]
      YO2g_glyc = u1/1 * MW["O2"]/MW["etoh"]
      
      p.add("YCO2x_glyc", value=YCO2x_glyc, vary=False)
      p.add("YCO2g_glyc", value=YCO2g_glyc, vary=False)
      p.add("YO2g_glyc", value=YO2g_glyc, vary=False)

      #5. maintenance metabolism : gluc + 6*O2 = 6*CO2 + 6*H2O

      YCO2s_m = 6 * MW["CO2"]/MW["gluc"]
      p.add("YCO2s_m", value=YCO2s_m, vary=False)
      
      exp_names = y0_dict.keys() # experiment names

      res_all= [] # empty (list which will be an array), will contain residuals

      for exp in exp_names: # loop over experiments
            y0 = y0_dict[exp]
            c = c_dict[exp]     
            datasets = datasets_dict[exp]

            res_this_exp = residuals_single_exp(p, c, y0, datasets)
            res_all = np.append(res_all, res_this_exp)

      return res_all


