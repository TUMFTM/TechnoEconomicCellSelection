# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 13:08:03 2022

@author: ga78tas
"""

import numpy as np
from scipy.optimize import bisect

class Costmodel():
    
    def __init__(self, motor_power, dt_con, p_fc):
        
        # Mobility parameters
        self.r = 1.095 # Implicit discount rate (derived from observed investment behavior, not rational economics)  [EU Commission impact assessment on CO2 emission performance standards for new heavy duty vehicles]
        self.servicelife = 5 # Service life in years [Wietschel et al. 2019]
        
        #Powertrain
        dmc_ice = 72.0 # Engine direct manufacturing costs in €/kW [Link et al. 2021]
        dmc_tank = 2.0 # Tank direct manufacturing costs in €/L [Link et al. 2021]
        dmc_eat = 19.8 # Exhaust after treatmenst direct manufacturing cost in €/kW [Link et al. 2021]
        vol_tank = 800 # Tank volume in L [Assumption]
        dmc_motor = 32 # Motor direct manufacturing cost in €/kW [Link et al. 2021]
        dmc_pe_hv = 35 # Power electronics direct manufacturing cost in €/kW [Link et al. 2021]
        imc = 1.45 # indirect manufacturing costs.
        self.c_pt = {
            "dt": imc*(motor_power/1000*(dmc_ice+dmc_eat)+vol_tank*dmc_tank), # Cost of DT powertrain components
            "bet": imc*motor_power/1000*(dmc_motor+dmc_pe_hv) # Cost of BET powertrain excl. battery
            }
        
        # Toll costs
        c_toll_road = 0.169 # €/km Mautsatz für Infrastruktur für LKW>18t & ab 4 Achsen [BFStrMG]
        c_toll_air = 0.011 # €/km Mautsatz für Luftverschmutzung EURO6 Diesel [BFStrMG]
        c_toll_noise = 0.002 # €/km Mautsatz für Lärmbelästigung EURO6 Diesel [BFStrMG]
        self.z_toll = 0.92 # Share of toll roads in total mileage % [Hülsmann]
        self.c_toll = {
            "dt": c_toll_road + c_toll_air + c_toll_noise,
            "bet": 0 # Electric vehicles are exempt from toll [BFStrMG]
                }
        
        # Taxes
        self.c_tax = {
            "dt": 556, # €/a [Kraftstg §9 Ab. (1) No. 4a]
            "bet": 556/2 # Electric vehicles get a 50% tax deduction [Kraftstg §9 Ab. (2)]
            }
        
        # Maintenance
        self.c_maintenance = {
            "dt": 0.147, # €/km [Kleiner et al.] For a 40t tractor-trailer in long haul operation
            "bet": 0.098 # €/km [Kleiner et al.] (excluding replacement of the battery system, electric machine and the power electronic)
        }
        
        # Energy consumption
        self.c_elec_sc = 0.2515 # Slow charging electricity cost in EUR/kWh excl. VAT [ISI - LFS III 2021]
        self.c_elec_fc = {350: 0.35/1.19, # 350kW charging electricity cost in EUR/kWh excl. VAT [Ionity price]
                     1000: 0.44/1.19} # 1MW charging electricity cost in EUR/kWh excl. VAT [Assumed to be upper limit announced by Scheuer]
        self.z_fc = 0.2 #Share of fast charged energy [Basma et al.]
        self.dt_con = dt_con #DT consumption in liter/km
        self.c_ene = {
            "dt": 1.609/1.19, # Cost of diesel fuel in €/liter}
            "bet": self.z_fc*self.c_elec_fc[p_fc] + (1-self.z_fc)*self.c_elec_sc # Cost of charging in €/kWh
            }
        
        # Battery
        self.c_cell2system = 1.3 # Cost scaling factor cell to system [EUCAR State of the art 2019] 
        self.z_scr = 0.15 # Residual value of battery at end of life [Burke and Fulton 2019]        
        self.c_par_avg = 105.77 #Average BEV cell price 2022 [BNEF 2022: https://about.bnef.com/blog/lithium-ion-battery-pack-prices-rise-for-first-time-to-an-average-of-151-kwh/]         

    def residualValue (self, total_mileage): 
        # Empiric relation found by Friedrich und Kleiner (2017)
        result = 0.951 * np.exp(-0.002 * total_mileage / 1000) # Share or original value
        return result
    
    def calculate_cost(self, veh, s_annual, con, cbat_spec, e_bat, t_bat):
        
        # Powertrain cost
        c_pt_residual = (self.c_pt[veh]*self.residualValue(s_annual*self.servicelife)
                         *self.r**-self.servicelife) #residual value
        c_pt = self.c_pt[veh] - c_pt_residual #total powertrain cost
        
        # Taxes (due once a year)        
        c_tax = sum([self.c_tax[veh]*self.r**-t 
                      for t in range(1, self.servicelife+1)]) 
        
        # Maintenance costs (due twice a year)      
        c_maint = sum([self.c_maintenance[veh]*s_annual/2*self.r**(-t/2)
                       for t in range(0,2*self.servicelife)])
        
        # Toll costs (due weekly)      
        c_toll = sum([self.z_toll*self.c_toll[veh]*s_annual/52*self.r**(-t/52)
                     for t in range(0,52*self.servicelife)])
        
        # Energy consumption costs (due weekly)
        c_ene = sum([self.c_ene[veh]*con*s_annual/52*self.r**(-t/52) 
                     for t in range(0, 52*self.servicelife)])
        
        ## Battery Costs
        if veh=="bet": # only for BET
            
            # Initial battery costs
            c_bat_init = cbat_spec*self.c_cell2system*e_bat
            
            # Battery replacement costs
            n_replacements = int(self.servicelife/t_bat) # number of replacements
            t_replacement = [t_bat*n for n in range(1,n_replacements+1)] # time of replacments in years            
            c_bat_replacements = [cbat_spec*self.c_cell2system*e_bat*self.r**-t for t in t_replacement]
            c_bat_replacement = sum(c_bat_replacements)
            
            # Residual value
            c_bat_scrappage = [self.z_scr*cbat_spec*self.c_cell2system*e_bat*self.r**-t for t in t_replacement] #scrappage values
            bat_eol_soh = (1-self.servicelife/t_bat%1) #remaining capacity at end of lifetime
            c_bat_eol = ((self.z_scr + bat_eol_soh*(1-self.z_scr))
                         *cbat_spec*self.c_cell2system*e_bat* self.r**-self.servicelife) # Add resale value of remaining capacity at EOL
            c_bat_residuals = c_bat_scrappage + [c_bat_eol]
            c_bat_residual = -sum(c_bat_residuals)
            
            # Imputed interest
            t_operation = [t_bat]*n_replacements + [self.servicelife - t_bat*n_replacements]#Vector of investment durations
            c_bat_inv = [c_bat_init] + c_bat_replacements
            avg_bound_investment = [(c_purchase+c_res)/2 for (c_purchase,c_res) 
                                    in zip(c_bat_inv, c_bat_residuals)] #investment till scrappage
            c_bat_imputed_interest = sum([invest*(self.r**t-1) for (invest,t)
                                       in zip(avg_bound_investment, t_operation)])
        else:
            c_bat_init = 0      
            c_bat_replacement = 0
            c_bat_residual = 0
            c_bat_imputed_interest = 0
        
        #Imputed interest
        c_pt_imputed_interest = (self.c_pt[veh]+c_pt_residual)/2*(self.r**self.servicelife-1) #imputed interest
        c_imputed_interest = c_pt_imputed_interest + c_bat_imputed_interest 

        #Total costs
        c_tot = sum([c_pt, c_tax, c_toll, c_maint, c_ene, c_bat_init, c_bat_replacement, c_bat_residual, c_imputed_interest])
        return c_tot, c_pt, c_tax, c_toll, c_maint, c_ene, c_bat_init, c_bat_replacement, c_bat_residual, c_imputed_interest
        
    def costparityanalysis(self, s_annual, e_bat, bet_con, t_bat):
        
        #Determine DT cost
        dt_costs = self.calculate_cost("dt", s_annual, self.dt_con, 0, 0, 0)[0]
        
        #Find cost parity price using bisection
        c_bat_min = 0 #minimum considered specific battery cost in EUR/kWh
        c_bat_max = 500 #maximum considered specific battery cost in EUR/kWh   
        c_bat_par = bisect(lambda c_bat: dt_costs - self.calculate_cost(
            "bet", s_annual, bet_con, c_bat, e_bat, t_bat)[0], c_bat_min, c_bat_max)
        
        return c_bat_par