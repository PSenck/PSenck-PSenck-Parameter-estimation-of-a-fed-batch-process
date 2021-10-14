import numpy as np
import pandas as pd
import glob
import os.path
from copy import deepcopy

def load_online(path, exp_numbers, cols):

    """load online rawdata as pd.DataFrame, measured by the Sartorius device. The CSV data must be in a folder with file names: online_number
    Example: online_4 for online measurement in experiment number 4.
   
    INPUT:
    path ... path string of folder where data is present 
    exp_numbers ...  list of numbers considered out of all possible experiments in the folder
    cols ... list of strings for columns of interest

    OUTPUT:
    online_raw: dictionary with key: fermentation number, value: pd.DataFrame

    """

    online_files = {"ferm{0}".format(i) : os.path.join(path,"online_{0}.CSV".format(i)) for i in exp_numbers}
    online_raw = {}

    for ferm, f in online_files.items():

        df = pd.read_csv(f ,sep=";",encoding= "unicode_escape",decimal=",", skiprows=[1,2] , skipfooter=1, usecols = cols, engine="python")
        online_raw[ferm] = df
    
    return online_raw


def process_online (online_raw, start_dict, end_dict):
    """ process online rawdata as pd.DataFrame. 
    Timestamps are created, Dataframes gets filtered according to start and end timepoints. 
    Fermentation time as decimal number and first time derivative of the base consumption are created.

    INPUT:
    online_raw ... dict with Key: fermentation number, value: pd.Dataframes
    start_dict ... dict with start time points
    end_dict ... dict with end time points

    OUTPUT:
    online_dict: dictionary with key: fermentation number, value: processed pd.DataFrame
    """
    online_dict = {}
    for ferm, df  in online_raw.items():
        df = deepcopy(df)   #to avoid pandas annoying SettingWithCopyWarning
        #another way to avoid this would be: pd.options.mode.chained_assignment = None  # default='warn'
        df["PDatTime"] = pd.to_datetime( df["PDatTime"], format = "%d.%m.%Y  %H:%M:%S" )
        df = df[(df["PDatTime"] >= start_dict[ferm] ) &  (df["PDatTime"] <= end_dict[ferm])]  
        df["t"] = (df["PDatTime"] - start_dict[ferm]) / pd.Timedelta(1,"h")        #time as decimal number   
        if "BASET" in df:
            df["BASET"] = pd.to_numeric(df["BASET"] , downcast="float" , errors="coerce") # some values were recognized as string
            df["base_rate"] = df["BASET"].diff() / df["t"].diff()       
        df.set_index("t", inplace = True, drop = True) 
        online_dict[ferm] = df

    return online_dict

def load_offline(path, exp_numbers, cols):

    """load offline rawdata as pd.DataFrame, measured by bio dry mass determination and HPLC. The CSV data must be in a folder with file names: offline_number
    Example: offline_4 for offline measurement in experiment number 4.
   
    INPUT:
    path ... path string of folder where data is present 
    exp_numbers ...  list of numbers considered out of all possible experiments in the folder
    cols ... list of strings for columns of interest

    OUTPUT:
    offline_raw: dictionary with key: fermentation number, value: pd.DataFrame

    """

    offline_files = {"ferm{0}".format(i) : os.path.join(path,"offline_{0}.CSV".format(i)) for i in exp_numbers}
    #CSV should be converted to csv on Linux/Mac
    #offline_files = {"ferm{0}".format(i) : os.path.join(path,"offline_{0}.csv".format(i)) for i in exp_numbers}
    offline_raw = {}

    for ferm, f in offline_files.items():

        df = pd.read_csv(f,sep=";", encoding= 'unicode_escape', header = 0, usecols = cols)
        offline_raw[ferm] = df

    return offline_raw


def process_offline(offline_raw, start_dict, end_dict):

    """ process offline rawdata as pd.DataFrame. 
    Timestamps are created, Dataframes gets filtered according to start and end timepoints. 
    Fermentation time as decimal number is created.

    INPUT:
    offline_raw ... dict with Key: fermentation number, value: pd.Dataframes
    start_dict ... dict with start time points
    end_dict ... dict with end time points

    OUTPUT:
    offline_dict: dictionary with key: fermentation number, value: processed pd.DataFrame
    """
    
    offline_dict = {}
    for ferm, df  in offline_raw.items():
        df = deepcopy(df) #to avoid pandas annoying SettingWithCopyWarning
        df["ts"] = pd.to_datetime( df["ts"], format = "%d.%m.%Y %H:%M") 
        df = df[(df["ts"] >= start_dict[ferm] ) &  (df["ts"] <= end_dict[ferm])]
        df["t"] = (df["ts"] - start_dict[ferm]) / pd.Timedelta(1,"h")
        df.set_index("t", inplace = True, drop = True) 
        offline_dict[ferm] = df

    return offline_dict


def load_CO2(path, exp_numbers):

    """load CO2 rawdata as pd.DataFrame, measured by the BlueSens CO2 sensor. The  .dat data must be in a folder with file names: CO2_number.
    Example: CO2_4 for CO2 measurement in experiment number 4.
    
    INPUT:
    path ... path string of folder where data is present 
    exp_numbers ...  list of numbers considered out of all experiments in the folder

    OUTPUT:
    CO2_raw: dictionary with key: fermentation number, value: pd.DataFrame

    """

    CO2_files = {"ferm{0}".format(i) : os.path.join(path,"CO2_{0}.dat".format(i)) for i in exp_numbers}           
    CO2_raw = {}

    for ferm, f in CO2_files.items():

        df = pd.read_csv(f, sep=";", encoding= "unicode_escape", header = 0, skiprows=[0], usecols=[0,2,4], names =["ts","CO2", "p"])
        CO2_raw[ferm] = df

    return CO2_raw


def process_CO2(CO2_raw, start_dict, end_dict):

    """ 
    Process CO2 rawdata as pd.DataFrame. 
    Timestamps are created, Dataframes gets filtered according to start and end timepoints. 
    Fermentation time as decimal number is created.

    INPUT:
    CO2_raw ... dict with Key: fermentation number, value: pd.Dataframes
    start_dict ... dict with start time points
    end_dict ... dict with end time points

    OUTPUT:
    CO2_dict: dictionary with key: fermentation number, value: processed pd.DataFrame
    """

    CO2_dict = {}

    for ferm, df  in CO2_raw.items():  

        df = deepcopy(df)   #to avoid pandas annoying SettingWithCopyWarning
        try:
            df["ts"] = pd.to_datetime( df["ts"], format = "%d.%m.%Y %H:%M:%S", exact= False)        #sometimes this format is required sometimes not, depending on single specific rows
        except:
            df["ts"] = pd.to_datetime( df["ts"] ) #, exact= False

        #df["ts"] = pd.to_datetime( df["ts"], format = "%d.%m.%Y %H:%M", exact= False, errors= "coerce") # this could maybe be an alternative solution instead of try/except

        df = df[(df["ts"] >= start_dict[ferm] ) &  (df["ts"] <= end_dict[ferm])]
        df["t"] = (df["ts"] - start_dict[ferm]) / pd.Timedelta(1,'h')
        df.set_index("t", inplace = True, drop = True)
        CO2_dict[ferm] = df

    return CO2_dict



def process_online_optional_filter(online_raw, start_dict = None, end_dict = None):
    "same function as process_online, except start and/or end_dicts can be None"

    online_dict_optional = {}

    for ferm, df  in online_raw.items():
        df = deepcopy(df)   #to avoid pandas annoying SettingWithCopyWarning
        #another way to avoid this would be: pd.options.mode.chained_assignment = None  # default='warn'
        df["PDatTime"] = pd.to_datetime( df["PDatTime"], format = "%d.%m.%Y  %H:%M:%S" )
        
        if start_dict == None and end_dict == None:
            df["t"] = (df["PDatTime"] - df["PDatTime"][0]) / pd.Timedelta(1,"h")
        
        elif end_dict == None and start_dict != None:
            df = df[(df["PDatTime"] >= start_dict[ferm] )]  
            df["t"] = (df["PDatTime"] - start_dict[ferm]) / pd.Timedelta(1,"h") 

        elif start_dict == None and end_dict != None:
            df = df[(df["PDatTime"] <= end_dict[ferm])]  
            df["t"] = (df["PDatTime"] - df["PDatTime"][0]) / pd.Timedelta(1,"h")
        else:
            df = df[(df["PDatTime"] >= start_dict[ferm] ) &  (df["PDatTime"] <= end_dict[ferm])]  
            df["t"] = (df["PDatTime"] - start_dict[ferm]) / pd.Timedelta(1,"h")        #time as decimal number   
            
        if "BASET" in df:
            df["BASET"] = pd.to_numeric(df["BASET"] , downcast="float" , errors="coerce") # some values were recognized as string
            df["base_rate"] = df["BASET"].diff() / df["t"].diff()  

        df.set_index("t", inplace = True, drop = True) 
        online_dict_optional[ferm] = df

    return online_dict_optional



def process_offline_optional_filter(offline_raw, start_dict = None, end_dict = None):

    "same function as process_offline, except start and/or end_dicts can be None"
    
    offline_dict_optional = {}

    for ferm, df  in offline_raw.items():
        df = deepcopy(df) #to avoid pandas annoying SettingWithCopyWarning
        df["ts"] = pd.to_datetime( df["ts"], format = "%d.%m.%Y %H:%M") 

        if start_dict == None and end_dict == None:
            df["t"] = (df["ts"] - df["ts"][0]) / pd.Timedelta(1,"h")
        
        elif end_dict == None and start_dict != None:
            df = df[(df["ts"] >= start_dict[ferm] )]  
            df["t"] = (df["ts"] - start_dict[ferm]) / pd.Timedelta(1,"h") 

        elif start_dict == None and end_dict != None:
            df = df[(df["ts"] <= end_dict[ferm])]  
            df["t"] = (df["ts"] - df["ts"][0]) / pd.Timedelta(1,"h")
        else:
            df = df[(df["ts"] >= start_dict[ferm] ) &  (df["ts"] <= end_dict[ferm])]  
            df["t"] = (df["ts"] - start_dict[ferm]) / pd.Timedelta(1,"h")        #time as decimal number 

        df.set_index("t", inplace = True, drop = True) 
        offline_dict_optional[ferm] = df

    return offline_dict_optional



def process_CO2_optional_filter(CO2_raw, start_dict = None, end_dict = None):

    "same function as process_CO2, except start and/or end_dicts can be None"

    CO2_dict_optional = {}

    for ferm, df  in CO2_raw.items():  

        df = deepcopy(df)   #to avoid pandas annoying SettingWithCopyWarning
        try:
            df["ts"] = pd.to_datetime( df["ts"], format = "%d.%m.%Y %H:%M:%S", exact= False)        #sometimes this format is required sometimes not, depending on single specific rows
        except:
            df["ts"] = pd.to_datetime( df["ts"] ) #, exact= False

        #df["ts"] = pd.to_datetime( df["ts"], format = "%d.%m.%Y %H:%M", exact= False, errors= "coerce") # this could maybe be an alternative solution instead of try/except

        if start_dict == None and end_dict == None:
            df["t"] = (df["ts"] - df["ts"][0]) / pd.Timedelta(1,"h")
        
        elif end_dict == None and start_dict != None:
            df = df[(df["ts"] >= start_dict[ferm] )]  
            df["t"] = (df["ts"] - start_dict[ferm]) / pd.Timedelta(1,"h") 

        elif start_dict == None and end_dict != None:
            df = df[(df["ts"] <= end_dict[ferm])]  
            df["t"] = (df["ts"] - df["ts"][0]) / pd.Timedelta(1,"h")
        else:
            df = df[(df["ts"] >= start_dict[ferm] ) &  (df["ts"] <= end_dict[ferm])]  
            df["t"] = (df["ts"] - start_dict[ferm]) / pd.Timedelta(1,"h")        #time as decimal number 

        df.set_index("t", inplace = True, drop = True)
        CO2_dict_optional[ferm] = df

    return CO2_dict_optional
