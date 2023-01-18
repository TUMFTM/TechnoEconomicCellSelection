"""
Created on Sat May  1 15:09:42 2021

@author: olaf_
"""
import numpy as np
import pandas as pd
def postprocess(cells_raw):

    # Set index and convert undefined values to NaN
    cells_raw.set_index("Index #", inplace = True)        
    cells_raw.replace("not defined", np.nan, inplace=True)   
    cells_raw.replace('not applicable', np.nan, inplace=True)
    
    # Calculate gravimetric and volumetric energy densities
    for i, cell in cells_raw.iterrows():
        
        E_cell = cell["Max Capacity (AH)"]*cell["Nominal Voltage (V)"] # Cell energy in Wh
        
        if pd.notna(cell["vol. Energy Density (Wh/l)"]):
            cells_raw.at[i,"rho_bat"] = cell["vol. Energy Density (Wh/l)"] # volumetric density in Wh/l
        else:
            if cell["Format"] == "Cyl":
                v_cell = (0.25*np.pi*cell["Diameter (mm)"]**2*cell["Height (mm)"]*1e-6) #cell volume in liter
            else:            
                v_cell = (cell["Length (mm)"]*cell["Height (mm)"]*cell["Width (mm)"]*1e-6) #cell volume in liter
            cells_raw.at[i,"rho_bat"] = E_cell/v_cell # volumetric density in Wh/l
        
        if pd.notna(cell["grav. Energy Density (Wh/kg)"]):
            cells_raw.at[i,"mspec"] = cell["grav. Energy Density (Wh/kg)"] # gravimetric density in Wh/kg
        else: 
            cells_raw.at[i,"mspec"] = E_cell/cell["Weight (gr)"]*1000 # gravimetric density in Wh/kg
    
    # Determine cell chemistry based on nominal voltage
    cells_raw.loc[cells_raw["Nominal Voltage (V)"] > 3.4, "Chemistry"] = "NMC/NCA"
    cells_raw.loc[cells_raw["Nominal Voltage (V)"] <= 3, "Chemistry"] = "LTO"
    cells_raw.loc[[3 < u <= 3.4 for u in cells_raw["Nominal Voltage (V)"]], "Chemistry"] = "LFP"

    # Limit to relevant columns
    relevant_columns = ["Company Name","Part #",
                        "Chemistry",
                        "Nominal Voltage (V)",
                        "Format",
                        'Max Capacity (AH)',
                        "Cycle Life",
                        "Last Cycle % OF Initial Capacity",
                        "mspec",
                        "rho_bat",
                        "Crate_charge_cyclelife"]
    cells = cells_raw[relevant_columns] 
    
    # Update names
    updated_names = ["Company","Name","Chemistry","Unom","Format","Capacity",
                     "ncycle","EOL","mspec","rho_bat","Ccharge"]
    cells.columns = updated_names
    
    # Remove incomplete data
    cells = cells.dropna()
    
    # Convert EOL from percent to per unit
    cells.loc[:, "EOL"] = cells["EOL"]/100
    
    # Add display name
    cells["displayname"] = [f"{comp} {name}" for comp,name in zip(cells.Company, cells.Name)]
    cells.sort_values("displayname", inplace=True)
    
    # Add columns for results 
    cells["Ebat"] = ""
    cells["Econs"] = ""
    cells["t_bat"] = ""
    cells["Vbat"] = ""
    cells["mbat"] = ""
    cells["maxpayload"] = ""
    cells["cpar"] = ""

    return cells