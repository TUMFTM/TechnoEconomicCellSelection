"""
Created on Tue Jun 15 14:49:16 2021

@author: ga78tas
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import rc
from matplotlib.ticker import StrMethodFormatter
from copy import deepcopy

class Myplots:
    
    def __init__(self):
        self.width_wide = 12.39 #full ppt slide
        self.width = 5.33 #half ppt slide
        self.height = 5.33
        self.font = {'weight': 'normal', 'size': 14}      
        rc('font', **self.font)
        
    def massimpact(self, method):

        # Generate figure
        fig, (ax0, ax1) = plt.subplots(1,2, 
                                        tight_layout = True, 
                                        figsize=(self.width_wide,self.height))
        
        #impact on energy consumption
        ax0.plot(method.f_con.x/1000, method.f_con.y, label="LongHaul")
        ax0.set_xlabel("Battery weight in t")
        ax0.set_ylabel("Consumption in kWh/km")
        ax0.set_ylim(bottom=0)
        ax0.set_xlim(left=0)
        ax0.grid()
        
        #impact on average speed
        ax1.plot(method.f_v_avg.x/1000, method.f_v_avg.y, label="LongHaul")
        ax1.set_xlabel("Battery weight in t")
        ax1.set_ylabel("Average speed in km/h")
        ax1.set_ylim(0, 90)
        ax1.set_xlim(left=0)
        ax1.grid()
             
    def batterysizing(self, scenario):
        
        method = scenario["method"]
        topcell = scenario["topcell"]
        
        #Define Crate and charging power ranges
        tch = 24-method.sizing.t_drive_max
        Pmin = method.f_con(method.sizing.m_bet_max)*method.f_v_avg(method.sizing.m_bet_max)*method.sizing.t_drive_max/tch
        Cmin = 1/tch
        Crates = np.linspace(Cmin,1.5,110)
        Pchs = np.linspace(Pmin, 1000, 100)
        
        #Generate method and cell objects
        methods = []
        for p in Pchs:
            method_var = deepcopy(method)
            method_var.sizing.p_fc = p
            methods.append(method_var)
        
        cells = [{"EOL": topcell.EOL,
                  "mspec": 1000 , #random irrelevant number
                  "rho_bat": 1000, #random irrelevant number
                  "ncycle": 1000, #random irrelevant number
                  "Format": "Pouch", #random irrelevant number
                  "Chemistry": "NMC/NCA", #random cell chemistry
                  "Ccharge": C} for C in Crates]
        
        #Find required battery size
        Ebats = [[method.sizing.size_battery(method.s_annual, method.f_con, method.f_v_avg, c)[0] for method in methods] for c in cells]
        
        #Generate figure
        fig, ax = plt.subplots(figsize=(self.width, 0.3*self.width_wide))
        CS = ax.contourf(Pchs, Crates, Ebats, 7)
        cbar = fig.colorbar(CS)
        cbar.ax.set_ylabel("Required battery size in kWh")
        plt.xlabel("Charger power in kW")
        plt.ylabel("Charging C-rate in 1/h")
        plt.xlim(left=0)
        plt.ylim(bottom=0)
        fig.tight_layout()
        
    def waterfall(self, scenario): 
        
        method = scenario["method"]
        cell = scenario["topcell"]
        
        #Write to dataframe
        cost_components = ["powertrain", "taxes", "toll", "maintenance", "energy", "battery"]
        colors = ["red", "orange", "yellow", "green", "blue", "indigo"]
        dt_costs = method.costmodel.calculate_cost("dt", method.s_annual, method.costmodel.dt_con, 0, 0, 0)[1:7]        
        bet_costs = [cell.c_pt, cell.c_tax, cell.c_toll, cell.c_maint, cell.c_ene, 
                     cell.c_bat_init+cell.c_bat_repl+cell.c_bat_resi+cell.c_bat_intr]
        data = pd.DataFrame([dt_costs, bet_costs, colors],
                      index = ["dt_costs","bet_costs","color"], 
                      columns = cost_components)
        
        #Calculate differences
        data.loc["delta"] = data.loc["dt_costs"]-data.loc["bet_costs"]
        data = data.loc[:,data.loc["delta"]!=0] #only non-zero differences are considered
        data.loc["bottom"] = sum(data.loc["dt_costs"]) - data.loc["delta"].cumsum() #bottomlines of waterfall bars
        ctot = sum(data.loc["dt_costs"]) #total cost
        
        #Calculate variables for stacked barcharts
        data.loc["dt_cumsum"] = data.loc["dt_costs"].cumsum().shift(1).fillna(0)
        data.loc["bet_cumsum"] = data.loc["bet_costs"].cumsum().shift(1).fillna(0)
        
        #generate plot
        textoffset = 1e4
        fig, ax = plt.subplots(figsize=(12.4,5.33)) 
        component_labels = {
            "powertrain": "Powertrain", 
            "taxes": "Taxes", 
            "toll": "Toll", 
            "maintenance": "Maintenance", 
            "energy": "Energy",
            "battery": "Battery", 
            "subsidy": "Subsidy"
            }            
        for i, param in enumerate(data):
            # Diesel truck bar plot
            ax.bar(0, data[param]["dt_costs"], 
                    bottom = data[param]["dt_cumsum"],
                    color = data[param]["color"]) 
            
            # Cost component differences
            ax.bar(i+1, data[param]["delta"],
                    bottom = data[param]["bottom"], 
                    color = data[param]["color"], 
                    label=component_labels[param]) 
            
            # Battery electric truck bar plot
            ax.bar(len(data.columns)+1, data[param]["bet_costs"], 
                    bottom = data[param]["bet_cumsum"],
                    color = data[param]["color"]) 
                
            # Add annotations
            if data[param]["delta"] >= 0: 
                ax.annotate("-€{:,.0f}".format(data[param]["delta"]),
                        (i+1, data[param]["bottom"]-textoffset),
                        ha = "center", 
                        va = "top")
            else:
                ax.annotate("+€{:,.0f}".format(abs(data[param]["delta"])), 
                        (i+1, data[param]["bottom"]+textoffset),
                        ha = "center", 
                        va = "bottom")
        
        # Draw connecting horizontal lines
        connectinglines = [ctot] + data.loc["bottom"].to_list()
        [ax.hlines(connectinglines[i], i-0.4, i+1.4, 
                    color = "k", 
                    linewidth = 0.5) for i in range(0, len(connectinglines))]
        
        # Finetuning
        ax.set_xticks([0, len(data.columns)+1])
        ax.set_xticklabels(["DT costs","BET costs"])
        ax.set_ylim(0, 1.2*max(data.loc["bet_cumsum"][-1],
                                data.loc["dt_cumsum"][-1]))
        ax.yaxis.set_major_formatter(StrMethodFormatter('€{x:,.0f}')) 
        ax.legend(ncol = len(data.columns), 
                  loc='upper center',
                  bbox_to_anchor=(0.5, -0.1))  
        ax.set_axisbelow(True)
        ax.grid(axis='y')
        fig.tight_layout()
        fig.subplots_adjust(bottom=0.2)
        
    def cell_assessment(self, scenarios):
                
        #Define marker shapes and colors for all plots
        markershapes = {"Cylindrical": "o", #Round
                        "Pouch": "d", #Thin diamond
                        "Prismatic": "s"} #Square
        markercolors = {"NMC/NCA": "blue",
                  "LFP": "gray", 
                  "LTO": "red"}
        shape = {"Cylindrical": "Cylindrical",
                "Pouch": "Pouch",
                "Prismatic": "Prismatic"}
        chems = {"NMC/NCA": "NMC/NCA", 
                "LFP": "LFP", 
                "LTO": "LTO"}
        
        # Generate plot
        fig, axs = plt.subplots(len(scenarios), 2, 
                                  figsize=(self.width_wide,0.6*self.width_wide),
                                  sharex = "col",
                                  sharey = "all")        
        
        for ax_row, scenario in zip(axs, scenarios.values()):
            result = scenario["results"]
            method = scenario["method"]
            topcell = scenario["topcell"]
            
            #Draw scatters of different cell formats and chemistries
            result.dropna(inplace=True)
            result.sort_values("Chemistry", ascending=False, inplace=True)
            groups = result.groupby(["Chemistry", "Format"], sort = False)
            for name, group in groups:                                
                ax_row[0].scatter(group.maxpayload, 
                            group.cpar, 
                            marker = markershapes[name[1]],
                            color = markercolors[name[0]],
                            s = 80,
                            label=f"{chems[name[0]]} {shape[name[1]]}")
                
                ax_row[1].scatter(group.Vbat, 
                            group.cpar, 
                            marker = markershapes[name[1]],
                            color = markercolors[name[0]],
                            s = 80,
                            label=f"{chems[name[0]]} {shape[name[1]]}")
            
            #Draw reference cell
            ax_row[0].scatter(topcell.maxpayload, topcell.cpar, 
                              marker = markershapes[topcell.Format], 
                              color = "orange", s=80, label="reference cell")
            ax_row[1].scatter(topcell.Vbat, topcell.cpar, 
                              marker = markershapes[topcell.Format], 
                              color = "orange", s=80, label="reference cell")
            
            #Draw lines for reference load, DT payload and DT volume
            ax_row[0].set_xlim(0, 28e3)
            ax_row[0].set_ylim(0, 360)
            ax_row[0].axvline(method.sizing.ref_load, color = "red", linestyle = "dashed")
            ax_row[0].axvline(method.sizing.dt_payload_max, color = "red", linestyle = "dashed")
            ax_row[0].text(0.95*method.sizing.ref_load, ax_row[0].get_ylim()[1]*0.5, "Reference load",
                            color = "red", rotation=90,
                            ha = "center", va = "center")
            ax_row[0].text(1.075*method.sizing.dt_payload_max, ax_row[0].get_ylim()[1]*0.5, "Max. payload\n diesel truck",
                            color = "red", rotation=90,
                            ha = "center", va = "center")            
            ax_row[0].set_xlim(0, method.sizing.dt_payload_max*1.13)
            
            ax_row[1].set_xlim(0, 8750)
            ax_row[1].set_ylim(0, 360)
            ax_row[1].axvline(method.sizing.dt_vol_pt, color = "red", linestyle = "dashed")
            ax_row[1].text(0.75*method.sizing.dt_vol_pt,  ax_row[0].get_ylim()[1]*0.5, "Powertrain volume\n Diesel truck",
                            color = "red", rotation=90,
                            ha = "center", va = "center")
            ax_row[1].set_xlim(left = 0)
        
        # Finetuning
        for ax in axs.flat: ax.grid()
        handles, labels = plt.gca().get_legend_handles_labels()        
        fig.legend(handles, labels, loc='lower center', ncol = len(handles))
        axs[-1][0].set_xlabel("Maximum payload in kg")
        axs[-1][1].set_xlabel("Battery volume in liter")
        axs[0][0].set_ylabel("Cost parity price in €/kWh")
        axs[1][0].set_ylabel("Cost parity price in €/kWh")
        
        # insert annotation of usecase
        for ax, s in zip(axs[:,0], scenarios):
            ax.annotate(s, xy=(-0.25, 0.5),
                        xycoords='axes fraction',
                        ha='right', va='center')
        fig.tight_layout()
        
    def sensitivity(self, scenarios, cells, selection):
        
        fig, axs = plt.subplots(len(scenarios),2,
                                sharey = "all",
                                sharex = "col",
                                figsize=(self.width_wide, 0.5*self.width_wide))
        
        #Define markers
        markers = ["o", "v", "D", "s", "*"]        
        for ax_row, (scenario_name, scenario) in zip(axs, scenarios.items()):           
            topcell = scenario["topcell"]
            method = scenario["method"]

            #Set parameter variations
            parameter_variation = self.define_parameters(topcell, method, cells, selection)
            
            #Calculate cost parity price
            c_pars, payloads, volumes, con_maxs, con_refs = self.param_sensivitiy(method, topcell, parameter_variation)
            
            #Plot results
            for param, m in zip(selection, markers):
                if param in ["m_cell2system", "vol_cell2system"]:
                    paramlabel = f"{param}: {parameter_variation[param][0][topcell.Chemistry][topcell.Format]*1000:.2f}-{parameter_variation[param][-1][topcell.Chemistry][topcell.Format]*1000:.2f}"
                elif param == "consumption":
                    paramlabel = f"{param}: {con_maxs[0]:.2f}//{con_refs[0]:.2f}-{con_maxs[-1]:.2f}//{con_refs[-1]:.2f}"
                else:
                    paramlabel = f"{param}: {parameter_variation[param][0]:.2f}-{parameter_variation[param][-1]:.2f}"
                ax_row[0].plot(payloads[param], c_pars[param], 
                               marker = m, markevery = [0,-1], 
                               label = paramlabel)
                ax_row[1].plot(volumes[param], c_pars[param], 
                               marker = m, markevery = [0,-1])

            #Draw lines for reference load, DT payload and DT volume
            ax_row[0].legend(numpoints=2)
            ax_row[0].set_xlim(0, 28e3)
            ax_row[0].set_ylim(0, 360)
            ax_row[0].axvline(method.sizing.ref_load, color = "red", linestyle = "dashed")
            ax_row[0].axvline(method.sizing.dt_payload_max, color = "red", linestyle = "dashed")
            ax_row[0].text(0.95*method.sizing.ref_load, ax_row[0].get_ylim()[1]*0.5, "Reference load",
                            color = "red", rotation=90,
                            ha = "center", va = "center")
            ax_row[0].text(1.075*method.sizing.dt_payload_max, ax_row[0].get_ylim()[1]*0.5, "Max. payload\n diesel truck",
                            color = "red", rotation=90,
                            ha = "center", va = "center")            
            ax_row[0].set_xlim(0, method.sizing.dt_payload_max*1.13)
            
            ax_row[1].set_xlim(0, 8750)
            ax_row[1].set_ylim(0, 360)
            ax_row[1].axvline(method.sizing.dt_vol_pt, color = "red", linestyle = "dashed")
            ax_row[1].text(0.75*method.sizing.dt_vol_pt,  ax_row[0].get_ylim()[1]*0.5, "Powertrain volume\n Diesel truck",
                            color = "red", rotation=90,
                            ha = "center", va = "center")
        
        # Finetuning
        for ax in axs.flat: ax.grid()
        axs[-1][0].set_xlabel("Maximum payload in kg")
        axs[-1][1].set_xlabel("Battery volume in liter")
        axs[0][0].set_ylabel("Cost parity price in €/kWh")
        axs[1][0].set_ylabel("Cost parity price in €/kWh")
        
        # insert annotation of usecase
        for ax, s in zip(axs[:,0], scenarios):
            ax.annotate(s, xy=(-0.25, 0.5),
                        xycoords='axes fraction',
                        ha='right', va='center')
        fig.tight_layout()
        
    def define_parameters(self, topcell, method, cells, selection):
        
        #Define parameter variations
        n_var = 50 # Number of variations
        n_var_cons = 10 # Fewer variations for energy consumption to limit computation time
        
        #Calculate minimum and maximum values
        m_bat_max = method.sizing.m_bet_max - method.sizing.m_nobat-method.sizing.ref_load #Maximum allowed battery weight to transport reference load in kg
        m_cell2system_min = topcell.Ebat/m_bat_max/topcell.mspec  + 1e-9 #Minimum gravimetric packaging efficiency
        mspec_min = 1000*topcell.Ebat/m_bat_max/method.sizing.m_cell2system[topcell.Chemistry][topcell.Format] + 1e-9 #Minimum specific energy in Wh/kg
        Emax = topcell.mspec*method.sizing.m_cell2system[topcell.Chemistry][topcell.Format]*m_bat_max/1000 #Maximum battery size that can transport the reference load
        Eday = method.f_con(method.sizing.m_bet_max)*method.f_v_avg(method.sizing.m_bet_max)*method.sizing.t_drive_max #Daily energy consumption
        Crate_min = (Eday-Emax*method.sizing.z_usable*topcell.EOL)/Emax/method.sizing.t_fc + 1e-9 #Minimum Crate
        Crate_max = method.sizing.z_usable*topcell.EOL*min(
                method.sizing.cp_share/method.sizing.t_fc, 
                method.sizing.p_fc/(Eday-method.sizing.p_fc*method.sizing.t_fc)) #C-rate limit: higher C-rates do not result in a smaller required battery size

        #Define parameter variations
        variations = {"ncycle": np.linspace(min(cells["ncycle"]), max(cells["ncycle"]), n_var), 
                      "mspec": np.linspace(mspec_min, max(cells["mspec"]), n_var), 
                      "rho_bat": np.linspace(min(cells["rho_bat"]), max(cells["rho_bat"]), n_var), 
                      "Ccharge": np.linspace(max(Crate_min, min(cells["Ccharge"])), Crate_max, n_var),
                      "vol_cell2system": [
                          {"NMC/NCA": {"Cylindrical":cy, "Pouch":po, "Prismatic":pr},
                           "LFP": {"Cylindrical":cy, "Pouch":po, "Prismatic":pr},
                           "LTO": {"Cylindrical":cy, "Pouch":po, "Prismatic":pr}}
                       for cy,po,pr in zip(
                          np.linspace(0.247, 0.356, n_var), #Min & max values for cylindrical cells [Löbberding]
                          np.linspace(0.168, 0.528, n_var), #Min & max values for pouch cells [Löbberding]
                          np.linspace(0.168, 0.528, n_var) #Min & max values for prismatic cells [Löbberding]
                          )
                          ],               
                      "m_cell2system": [
                          {"NMC/NCA": {"Cylindrical":cy, "Pouch":po, "Prismatic":pr},
                           "LFP": {"Cylindrical":cy, "Pouch":po, "Prismatic":pr},
                           "LTO": {"Cylindrical":cy, "Pouch":po, "Prismatic":pr}}
                       for cy,po,pr in zip(
                          np.linspace(max(m_cell2system_min,0.484), 0.672, n_var), #Min & max values for cylindrical cells [Löbberding]
                          np.linspace(max(m_cell2system_min,0.495), 0.742, n_var), #Min & max values for pouch cells [Löbberding]
                          np.linspace(max(m_cell2system_min,0.495), 0.742, n_var), #Min & max values for prismatic cells [Löbberding]
                          )
                          ],
                      "c_cell2system": np.linspace(1.2, 2.21, n_var), #Min & max values for 2020 [König] 
                      "annual_mileage": np.linspace(96_850, 180_000, n_var), # maximum annual mileage
                      "servicelife": np.arange(2,11), #Payback period considered by large fleets [IEA 2017] until average truck life in Germany [Wolff 2021] 
                      "c_diesel": np.linspace(1.233/1.19, 2.16/1.19, n_var), #Lowest and highest Diesel cost in Germany between Jan 2021 and March 2022
                      "c_elec": np.linspace(0.05, 0.44/1.19, n_var), #Value reported by ... et al. and maximum electricity price 
                      "consumption": [{"crr": crr, "cd_a": cd_a} for crr, cd_a in 
                                      zip(np.linspace(0.00385, method.vehicle.fr, n_var_cons), #Minimum to average reported coefficient of rolling resistance
                                          np.linspace(4.325, method.vehicle.cd_a, n_var_cons))] #Minimum to average reported drag area
                       }
        
        selection_variations = {s: variations[s] for s in selection}
        
        return selection_variations
    
    def param_sensivitiy(self, method, cell, parameters):
        
        c_pars = {}
        payloads = {}
        volumes = {}
        con_maxs = []
        con_refs = []
        for param, variations in parameters.items():
            c_pars[param] = []
            payloads[param] = []
            volumes[param] = []
            if param in ["ncycle", "mspec", "Ccharge", "rho_bat"]: #Cell parameters
                for var in variations:
                    varcell = cell.copy()
                    varcell[param] = var
                    c_par, payload, volume = method.eval_cell(varcell)[0:3]
                    c_pars[param].append(c_par)
                    payloads[param].append(payload)
                    volumes[param].append(volume)   
            elif param in ["vol_cell2system", "m_cell2system", "p_fc"]: #Sizing parameters
                for var in variations:
                    var_method = deepcopy(method)
                    setattr(var_method.sizing, param, var)
                    c_par, payload, volume = var_method.eval_cell(cell)[0:3]
                    c_pars[param].append(c_par)
                    payloads[param].append(payload)
                    volumes[param].append(volume)
            elif param in ["annual_mileage", "servicelife","c_cell2system"]: #Costmodel parameters
                for var in variations:
                    var_method = deepcopy(method)
                    setattr(var_method.costmodel, param, var)
                    c_par, payload, volume = var_method.eval_cell(cell)[0:3]
                    c_pars[param].append(c_par)
                    payloads[param].append(payload)
                    volumes[param].append(volume)   
            elif param == "c_elec": #Special case electricity costs
                for var in variations:
                    var_method = deepcopy(method)
                    var_method.costmodel.c_ene["bet"] = var
                    c_par, payload, volume = var_method.eval_cell(cell)[0:3]
                    c_pars[param].append(c_par)
                    payloads[param].append(payload)
                    volumes[param].append(volume)  
            elif param == "c_diesel": #Special case diesel costs
                for var in variations:
                    var_method = deepcopy(method)
                    var_method.costmodel.c_ene["dt"] = var
                    c_par, payload, volume = var_method.eval_cell(cell)[0:3]
                    c_pars[param].append(c_par)
                    payloads[param].append(payload)
                    volumes[param].append(volume)   
            elif param == "consumption": #Special case energy consumption
                for var in variations:
                    var_method = deepcopy(method)
                    var_method.vehicle.fr = var["crr"]
                    var_method.vehicle.cd_a = var["cd_a"]
                    var_method.f_con, var_method.f_v_avg  = var_method.vehicle.simulate_consumption(
                        var_method.sizing.m_nobat, var_method.sizing.m_bet_max)
                    c_par, payload, volume, con_ref, con_max = var_method.eval_cell(cell)[0:5]
                    c_pars[param].append(c_par)
                    payloads[param].append(payload)
                    volumes[param].append(volume) 
                    con_refs.append(con_ref)
                    con_maxs.append(con_max)
            else:
                print(f"Please update the parameter sensitivity function to account for the parameter variation of {param}")                        

        return c_pars, payloads, volumes, con_maxs, con_refs