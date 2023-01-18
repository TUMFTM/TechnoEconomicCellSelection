# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 11:20:57 2022

@author: ga78tas
"""
import numpy as np
from pandas import read_csv
from scipy.interpolate import interp1d
from numba import njit

class Vehicle_simulation():

    def __init__(self):
        
        # Vehicle parameters
        self.motor_power = 352.3675316051114e3 #Max power in W [average of all DT registered between the 1st January 2019 and the 30 June 2020]
        self.fr = 0.00385 # Rolling friction coefficient [Lowest of all DT registered between the 1st January 2019 and the 30 June 2020]
        self.cd_a = 4.325 # Drag coefficient x Front surface area in m^2 [Lowest of all DT registered between the 1st January 2019 and the 30 June 2020]
        self.p_aux = 2.3e3 # auxiliary power consumption in W [Zhao 2013]
        self.eta = 0.85 # overall powertrain efficiency [Earl 2018]
        
        # Driver parameters
        self.decmax = 1 # Maximum deceleration in m/s/s [VECTO]
        self.accmax = 1 # Maximum acceleration in m/s/s [VECTO]
        
        # Environmental constants
        self.g = 9.81 # Gravitational acceleration in m/s/s [VECTO]
        self.rho = 1.188 # Air density in kg/m/m/m [VECTO]
        
        # Setup simulation
        drivingcycle = read_csv("inputs/LongHaul_drivingcycle.vdri") #Load VECTO Longhaul driving cycle
        self.preprocessing(drivingcycle) #Generate target speed profile with deceleration phases
        
    def preprocessing(self, missionProfile):
        
        # Unpack missionprofile        
        s = missionProfile["<s>"].values # Distance in m
        v = missionProfile["<v>"].values/3.6 # Speed in m/s
        grad = missionProfile["<grad>"].values # Gradient in %
        stop = missionProfile["<stop>"].values # Stop duration in s
        
        # Calculate distance step along route
        ds = np.diff(s)
        
        # Generate array with dec phases
        i2 = np.where(np.diff(v)<0)[0]+1 #ends of deceleration phases
        i1 = np.append([0], i2[:-1]) #start of deceleration phases
        v_target = np.zeros(len(v))
        for i1, i2 in zip(i1,i2):
            v_target[i1:i2] = np.minimum(
                v[i1:i2], # drive cycle speed
                np.sqrt(v[i2]**2+2*self.decmax*(s[i2]-s[i1:i2]))) #deceleration speed
                
        self.s = s
        self.v = v
        self.grad = grad
        self.stop = stop
        self.vtarget = v_target
        self.ds = ds

    def simulate_consumption(self, m_min, m_max):
        
        # Generate lookup 
        n_points = 4
        gvws = np.linspace(m_min, m_max, n_points)
        cons, v_avgs = zip(*[self.run_sim(gvw) for gvw in gvws])        
        f_con = interp1d(gvws,cons)
        f_v_avg = interp1d(gvws, v_avgs)
        
        return f_con, f_v_avg
        
    def run_sim(self, gvw):
        
        return self.run_sim_fast(
            self.stop, self.ds, self.vtarget, self.grad, self.g, self.rho, self.accmax, 
            self.fr, self.cd_a, self.motor_power, self.eta, self.p_aux, gvw)
    
    @staticmethod
    @njit
    def run_sim_fast(stop_arr, ds_arr, vtarget_arr, grad_arr, g, rho, accmax, fr, cd_a, motor_power, eta, p_aux, gvw):
        
        #Initial simulation states
        s = [0]
        t = [0]
        v = [0]
        p = []
        
        for (stop, ds, vtarget, grad) in zip(stop_arr, ds_arr, vtarget_arr[1:], grad_arr):
            
            # If the vehicle is stopped, add extra entry to list
            if stop>0:
                s.append(s[-1])
                t.append(t[-1]+stop)
                v.append(0)
                p.append(0)   
            
            # Determine target acceleration
            atarget = (vtarget**2-v[-1]**2)/2/ds #target acceleration
            
            # Determine power limited maximum acceleration
            f_roll = gvw*g*fr*np.cos(np.arctan(grad/100)) #Rolling resistance in N
            f_drag = 0.5*rho*cd_a*v[-1]**2 #Drag resistance in N
            f_incl = gvw*g*np.sin(np.arctan(grad/100)) #Inclination resistance in N
            f_max = motor_power/v[-1] if v[-1]>0 else 1e9  #Max driving force in N
            apower = (f_max*eta-f_roll-f_drag-f_incl)/gvw #Max acceleration in m/s/s
            
            # Determine applied acceleration and new states
            a = min(atarget, apower, accmax) #Applied acceleration in m/s/s
            v_new = np.sqrt(v[-1]**2+2*a*ds) #New vehicle speed in m/s
            t_new = t[-1] + 2*ds/(v[-1]+v_new) #New time since start in s
            f_res = gvw*a+f_roll+f_drag+f_incl
            p_new = f_res*v[-1]*eta**np.sign(-f_res) #Applied power in W
            
            # Append new states to lists
            s.append(s[-1]+ds)
            t.append(t_new)
            v.append(v_new)
            p.append(p_new)
            
        p_bat = [max(p_mot, -motor_power) for p_mot in p] # Limit regeneration to motor power (remainder uses mechanical breaking)
        t_diff = [t2-t1 for t1, t2 in zip(t[:-1], t[1:])] # Calculate step durations
        e_con = [p*dt for p,dt in zip(p_bat, t_diff)] # Calculate energy consumption per timestep
        e_tot = sum(e_con)+p_aux*t[-1] # Calculate total energy consumption
        con = e_tot/3600/s[-1] # Convert energy consumption to kWh/km
        v_avg = s[-1]/t[-1]*3.6 # Average speed in km/h
        
        return con, v_avg
    