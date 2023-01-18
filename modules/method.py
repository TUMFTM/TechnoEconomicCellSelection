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
        
        # Parameters relevant to multiple functios
        self.s_annual = 116_000 # Annual mileage in km [VECTO 5-LH]
        self.p_fc = p_fc #Installed charging power in kW
        
    def eval_cells(self, cells_input): 
        print("Sizing battery and performing cost parity analysis")
        
        # Create copy of dataframe
        cells = cells_input.copy()
        for i, cell in cells.iterrows():    
            
            # Run method for single cell
            cbat_par, payload_max, vol_bat, con, con_max, e_bat, t_bat, m_bat, \
            c_pt, c_tax, c_toll, c_maint, c_ene, c_bat_init, c_bat_replacement, \
                c_bat_residual, c_imputed_interest= self.eval_cell(cell)
            
            # Write results to dataframe
            cells.at[i,"Ebat"] = e_bat
            cells.at[i,"Econs"] = con
            cells.at[i,"t_bat"] = t_bat
            cells.at[i,"Vbat"] = vol_bat 
            cells.at[i,"mbat"] = m_bat
            cells.at[i,"maxpayload"] = payload_max
            cells.at[i,"cpar"] = cbat_par
            cells.at[i,"c_pt"] = c_pt
            cells.at[i,"c_tax"] = c_tax
            cells.at[i,"c_toll"] = c_toll
            cells.at[i,"c_maint"] = c_maint
            cells.at[i,"c_ene"] = c_ene            
            cells.at[i,"c_bat_init"] = c_bat_init
            cells.at[i,"c_bat_repl"] = c_bat_replacement
            cells.at[i,"c_bat_intr"] = c_imputed_interest
            cells.at[i,"c_bat_resi"] = c_bat_residual
            
        return cells
    
    def eval_cell(self, cell):
        
        # Battery sizing algorithm
        e_bat, con, con_max, t_bat, vol_bat, m_bat, payload_max = \
            self.sizing.size_battery(self.s_annual, self.f_con, self.f_v_avg, cell)
            
        # If the reference load can be transported
        if payload_max > self.sizing.ref_load:
            
            #determine cost parity price (break-even price)
            cbat_par = self.costmodel.costparityanalysis(self.s_annual, e_bat, con, t_bat)
                        
            #and calculate other cost components 
            _, c_pt, c_tax, c_toll, c_maint, c_ene, c_bat_init, c_bat_replacement, c_bat_residual, c_imputed_interest = \
                self.costmodel.calculate_cost("bet", self.s_annual, con, cbat_par, e_bat, t_bat)
        else: 
            cbat_par = None
            c_pt = None 
            c_tax = None
            c_toll = None 
            c_maint = None 
            c_ene = None 
            c_bat_init = None
            c_bat_replacement = None
            c_bat_residual = None
            c_imputed_interest = None
            
        return cbat_par, payload_max, vol_bat, con, con_max, e_bat, t_bat, m_bat, \
                c_pt, c_tax, c_toll, c_maint, c_ene, c_bat_init, c_bat_replacement, \
                c_bat_residual, c_imputed_interest