"""
Created on Thu Aug 26 08:49:02 2021

@author: ga78tas
"""

from modules.vehicle_simulation import Vehicle_simulation
from modules.sizing import Sizing
from modules.costmodel import Costmodel

class Method():
    
    def __init__(self, p_fc):

        # Initialize submodels
        self.vehicle = Vehicle_simulation() #Initialize vehicle simulation
        self.sizing = Sizing(p_fc) #Initialize battery sizing algorithm
        self.costmodel = Costmodel(self.vehicle.motor_power, self.sizing.dt_con, p_fc) #Initialize cost model
        
        # Generate energy consumption lookup
        self.f_con, self.f_v_avg = self.vehicle.simulate_consumption(
            self.sizing.m_nobat, self.sizing.m_bet_max)
        
    def eval_cells(self, cells_input): 
        print("Sizing battery and performing cost parity analysis")
        
        # Create copy of dataframe
        cells = cells_input.copy()
        for i, cell in cells.iterrows():    
            
            # Run method for single cell
            cbat_par, payload_max, vol_bat, con, con_max, e_bat, bat_life_km, m_bat, \
                c_pt, c_tax, c_toll, c_maint, c_ene, c_bat = self.eval_cell(cell)
                    
            # Write results to dataframe
            cells.at[i,"Ebat"] = e_bat
            cells.at[i,"Econs"] = con
            cells.at[i,"bat_life_km"] = bat_life_km
            cells.at[i,"Vbat"] = vol_bat 
            cells.at[i,"mbat"] = m_bat
            cells.at[i,"maxpayload"] = payload_max
            cells.at[i,"cpar"] = cbat_par
            cells.at[i,"c_pt"] = c_pt
            cells.at[i,"c_tax"] = c_tax
            cells.at[i,"c_toll"] = c_toll
            cells.at[i,"c_maint"] = c_maint
            cells.at[i,"c_ene"] = c_ene            
            cells.at[i,"c_bat"] = c_bat
            
        cells.dropna(inplace=True)
        return cells
    
    def eval_cell(self, cell):
        
        # Battery sizing algorithm
        e_bat, con, con_max, bat_life_km, vol_bat, m_bat, payload_max = \
            self.sizing.size_battery(self.f_con, self.f_v_avg, cell)
            
        # If the reference load can be transported
        if payload_max > self.sizing.ref_load:
            
            #determine cost parity price (break-even price)
            cbat_par = self.costmodel.costparityanalysis(e_bat, con, bat_life_km)
                        
            #and calculate other cost components 
            _, c_pt, c_tax, c_toll, c_maint, c_ene, c_bat = \
                self.costmodel.calculate_cost("bet", con, cbat_par, e_bat, bat_life_km)
        else: 
            cbat_par = None
            c_pt = None 
            c_tax = None
            c_toll = None 
            c_maint = None 
            c_ene = None 
            c_bat = None
            
        return cbat_par, payload_max, vol_bat, con, con_max, e_bat, bat_life_km, m_bat, \
                c_pt, c_tax, c_toll, c_maint, c_ene, c_bat