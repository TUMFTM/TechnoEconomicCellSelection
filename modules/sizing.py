# -*- coding: utf-8 -*-
"""
Created on Fri Jan 28 13:58:38 2022

@author: ga78tas
"""

class Sizing():
    
    def __init__(self, p_fc):
        
        # Diesel truck payload capability and powertrain volume
        m_dt_szm = 7753.136555182494 # Mass of semi-truck tractor in kg [average of all DT registered between the 1st January 2019 and the 30 June 2020]
        m_trailer = 7500 # trailer mass in kg [Norris and Escher 2017]
        m_dt_max = 40e3 # maximum gross vehicle weight in kg [§ 34 StVZO]   
        self.dt_payload_max = m_dt_max - m_dt_szm - m_trailer #Maximum payload diesel truck
        self.dt_vol_pt = 3250 # DT powertrain volume [Master thesis M. Bstieler]

        # Diesel truck consumption
        dt_con_LL = 24.72226283/100 # Energy consumption in L/km [average of all Class 5-LH 4x2 DT registered between the 1st January 2019 and the 30 June 2020]
        dt_con_RL = 32.35729932/100 # Energy consumption in L/km [average of all Class 5-LH 4x2 DT registered between the 1st January 2019 and the 30 June 2020]
        self.low_load = 2.6e3 # Low load mass in kg [VECTO]
        self.ref_load = 19.3e3 # Reference load mass in kg [VECTO]
        self.share_low_load = 0.27/0.9 # Share of low load operation for Category 5-LH [VECTO]
        self.share_ref_load = 0.63/0.9 # Share of reference load operation for Category 5-LH [VECTO]
        self.dt_con = self.share_low_load*dt_con_LL+self.share_ref_load*dt_con_RL
        
        # BET
        m_empty_chassis = 0.75*m_dt_szm # chassis mass in kg [Phadke 2021]
        m_bet_drivetrain = 450 # drivetrain mass in kg [Mareev et al. 2017]  
        self.m_bet_max = 42e3 # maximum gross vehicle weight in kg [§ 34 StVZO, Punkt 6a, Satz 2]    
        self.m_nobat = m_empty_chassis + m_bet_drivetrain + m_trailer # Mass of BET excluding battery
        self.z_usable = 0.93 # Share of usable energy [Wassiliadis 2022]
        self.cp_share = 0.8 # share of battery that can be charged with max. C-rate [Nykvist 2019]
        self.t_cal = {"NMC/NCA": 13, "LFP": 15, "LTO": 20} #Generic calendar life until 80% of initial capacity for NMC and LFP cells [Hesse 2017] 

        #Volumetric packaging efficiency: Ratio of energy density from cell [Wh/l] to system level [kWh/l]     
        vol_c2s_id3 = 0.4145985401459854 #Volumetric packaging efficienciy of VW ID.3 cell [Wassiliadis et al.]
        vol_c2s_lfp = 0.553 #Tesla Model 3 0.64 #Volumetric packaging efficienciy of LFP cell [Frith et al.]
        vol_c2s_round_to_prism = 0.295/0.353 # Reduced volumetric packaging efficiency of cyclindrical cell compared to pouch/prismatic [Löbberding]
        self.vol_cell2system = {
            "NMC/NCA": {
                "Pouch": vol_c2s_id3,
                "Prismatic": vol_c2s_id3,
                "Cylindrical": vol_c2s_id3 * vol_c2s_round_to_prism
                },
            "LFP": {
                "Pouch": vol_c2s_lfp,
                "Prismatic": vol_c2s_lfp,
                "Cylindrical": vol_c2s_lfp * vol_c2s_round_to_prism
                },
            "LTO": {
                "Pouch": vol_c2s_lfp,
                "Prismatic": vol_c2s_lfp,
                "Cylindrical": vol_c2s_lfp * vol_c2s_round_to_prism
                }
            }
        
        #Gravimetric packaging efficiency: Ratio of specific energy from cell [Wh/kg] to system level [kWh/kg]
        m_c2s_id3 = 0.6336996336996337 #Gravimetric packaging efficienciy of VW ID.3 cell [Wassiliadis et al.]
        m_c2s_lfp = 0.716# Tesla model 3 0.84 #Gravimetric packaging efficienciy of LFP cell [Frith et al.]
        m_c2s_round_to_prism = 0.552/0.575 # Reduced gravimetric packaging efficiency of cyclindrical cell compared to pouch/prismatic [Löbberding]
        self.m_cell2system = {
            "NMC/NCA": {
                "Pouch": m_c2s_id3,
                "Prismatic": m_c2s_id3,
                "Cylindrical": m_c2s_id3 * m_c2s_round_to_prism
                },
            "LFP": {
                "Pouch": m_c2s_lfp,
                "Prismatic": m_c2s_lfp,
                "Cylindrical": m_c2s_lfp * m_c2s_round_to_prism
                },
            "LTO": {
                "Pouch": m_c2s_lfp,
                "Prismatic": m_c2s_lfp,
                "Cylindrical": m_c2s_lfp * m_c2s_round_to_prism
                }
            }
        
        # Operating parameters
        self.p_fc = p_fc # Available charging power during rest period in kW
        self.t_trip_max = 4.5 # Maximum trip driving time in h [StVZO Artikel 7]
        self.t_drive_max = 9 # Longest driving time with single rest period in h [StVZO Artikel 7]
        t_break = 0.75 # Mandatory driver break in h [StVZO Artikel 7]
        t_plug_unplug = 5/60 # Time required to connect and disconnect to the charger [Assumption]
        self.t_fc = t_break-t_plug_unplug # Time available for charging during the rest period in h
        
    def size_battery(self, s_annual, f_con, f_v_avg, cell):
        
        #Determine maximum energy consumption
        con_max = f_con(self.m_bet_max).item()
        v_avg_max = f_v_avg(self.m_bet_max).item()
        e_max = con_max*v_avg_max*self.t_drive_max #Daily max. energy demand in kWh

        #Size battery
        ene_req_day = e_max/(1 + min(cell["Ccharge"]*self.t_fc/cell["EOL"]/self.z_usable, self.cp_share)) #Required energy based on ability to recharge in kWh
        ene_req_charger = (e_max - self.p_fc*self.t_fc) #Required energy based on available charger power in kWh
        e_bat = max(ene_req_day, ene_req_charger)/cell["EOL"]/self.z_usable #Required battery size in kWh
        
        # Determine reulsting battery & vehicle properties
        m_bat = e_bat/cell["mspec"]/self.m_cell2system[cell["Chemistry"]][cell["Format"]]*1000 #Battery mass in kg
        m_tot = m_bat + self.m_nobat #Vehicle mass in kg
        payload_max = self.m_bet_max-m_tot #Maximum payload in kg
        vol_bat = e_bat/cell["rho_bat"]/self.vol_cell2system[cell["Chemistry"]][cell["Format"]]*1000 #Battery volume in l
        
        #If RL payload is feasible, calculate characteristic energy consumption and battery lifetime in km
        if payload_max >= self.ref_load:
            # Determine energy consumption at usecase loads
            con = (self.share_low_load*f_con(m_tot+self.low_load) 
                   + self.share_ref_load*f_con(m_tot+self.ref_load))
            
            # Determine battery life            
            t_cal = self.t_cal[cell["Chemistry"]]*(1-cell["EOL"])/0.2
            n_annual = s_annual*con/e_bat #Annual full equivalent cycles
            t_bat = 1/(1/t_cal+n_annual/cell["ncycle"])
        else:
            con = None
            t_bat = None
        
        return e_bat, con, con_max, t_bat, vol_bat, m_bat, payload_max