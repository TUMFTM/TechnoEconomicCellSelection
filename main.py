"""
Created on Fri Apr 30 10:42:58 2021

@author: olaf_
"""

from pandas import read_excel
from modules.cellPostprocessing import postprocess
from modules.method import Method
from modules.myplots import Myplots

#%% Load & postprocess cell data
cells_raw = read_excel("inputs/CellDatabase_v5.xlsx") #Load cell database
cells = postprocess(cells_raw) #Post process cell data

#%% Evaluate default scenario (350kW charging)
method_350kW = Method(350) #Initialize method
results_350kW = method_350kW.eval_cells(cells) #Execute for cells in database
topcell_350kW = results_350kW.sort_values("cpar", ascending=False).iloc[0,:] #Select best performing cell

#%% Evaluate 1MW charging scenario
method_1MW = Method(1000) #Initialize method
results_1MW = method_1MW.eval_cells(cells) #Execute for cells in database
topcell_1MW = results_1MW.sort_values("cpar", ascending=False).iloc[0,:] #Select best performing cell

#%% write to result struct
scenarios = {"350kW": {"method": method_350kW, 
                       "results": results_350kW, 
                       "topcell": topcell_350kW}, 
             "1MW": {"method": method_1MW, 
                     "results": results_1MW, 
                     "topcell": topcell_1MW}}

#%% Generate plots
myplots = Myplots() #Initialize plotting class
myplots.massimpact(method_1MW) #Results vehicle simulation
myplots.batterysizing(scenarios["1MW"]) #Results battery sizing
myplots.waterfall(scenarios["1MW"]) #Results cost model
myplots.cell_assessment(scenarios) #Cost parity price over payload & battery volume
myplots.sensitivity(scenarios, cells, ["ncycle", "mspec", "rho_bat", "Ccharge"]) # Parameter sensitivity of cell parameters
myplots.sensitivity(scenarios, cells, ["c_cell2system", "m_cell2system", "vol_cell2system"]) # Parameter sensitivity of pack parameters
myplots.sensitivity(scenarios, cells, ["annual_mileage", "servicelife", "c_diesel", "c_elec", "consumption"]) # Parameter sensitivity of system parameters
