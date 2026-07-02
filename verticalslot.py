# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 15:37:14 2017

@author: tewing
"""

import numpy as np

gravity = 9.81
water_density = 997
import time

class FishwayDefinition:
    
    discharge_coefficient = 0.7
    discharge_coefficient_escalation = 0.0
    discharge_exponent = 0.5
    non_slot_loss = 0.0
    
    flow_rate = 0
    
    slot_inverts = []
    slot_discharge_coefficient = []
    slot_differential = []
    slot_flow_depth = []
    
    pool_inverts = []
    pool_areas = []
    pool_levels = []
    pool_diss_rate = []
    pool_supp_flows = []
    pool_supp_velocities = []
    
    vertical_slot_heights = []
    vertical_slot_areas = []
    vertical_slot_type = []
    
    upstream_water_level = 0
    downstream_water_level = 0
    
    has_been_run = False
    has_passed_sanity_check = False
    upstream_water_level_error = 1.0
    time_to_solve = 0
    
    database = []
    
    minimum_upstream_height_to_top_slot_invert_threshold = 0.25
    
    def __init__(self,\
                slot_inverts,\
                pool_inverts,\
                pool_areas,\
                upstream_water_level,\
                downstream_water_level):
                    
        self.slot_inverts = np.array(slot_inverts)
        self.pool_inverts = np.array(pool_inverts)
        self.pool_areas = np.array(pool_areas)
        self.upstream_water_level = upstream_water_level
        self.downstream_water_level = downstream_water_level
        self.vertical_slot_type = np.zeros(self.slot_inverts.size, dtype=int)
        self.pool_levels = np.zeros(self.slot_inverts.size-1)
        self.pool_supp_flows = np.zeros(self.slot_inverts.size-1)
        self.pool_supp_velocities = np.zeros(self.slot_inverts.size-1)
        self.slot_discharge_coefficient = self.discharge_coefficient*np.ones(self.slot_inverts.size)
    
        self.slot_differential = []
        self.slot_flow_depth = []
        self.pool_diss_rate = []
        self.vertical_slot_heights = []
        self.vertical_slot_areas = []
        self.database = []
    
    def add_slot_definition(self, height_array, area_array):
        self.vertical_slot_heights.append(np.array(height_array))
        self.vertical_slot_areas.append(np.array(area_array))
        return len(self.vertical_slot_heights)-1
    
    def set_slot_type(self, slot_index, slot_type):
        self.vertical_slot_type[slot_index] = slot_type
    
    def set_supplemental_flow(self, pool_index, flow, velocity):
        self.pool_supp_flows[pool_index] = flow
        self.pool_supp_velocities[pool_index] = velocity
    
    def get_area_for_height(self, slot_index, height_over_invert):
        return np.interp(height_over_invert,\
            self.vertical_slot_heights[self.vertical_slot_type[slot_index]],\
            self.vertical_slot_areas[self.vertical_slot_type[slot_index]])
        
    def get_area_for_elevation(self, slot_index, elevation):
        return self.get_area_for_height(slot_index,\
            elevation-self.slot_inverts[slot_index]);
    
    def slot_flow_rate(self, slot_index, us_elevation, ds_elevation):
        
        slot_invert_elevation = self.slot_inverts[slot_index]
        
        if(us_elevation<slot_invert_elevation or us_elevation<ds_elevation):
            # No flow - water levels below slot invert
            return 0.0;
        else:
            # Flow is downstream (positive)
            deltaH = us_elevation-max([ds_elevation,slot_invert_elevation])
            spec_discharge_coefficient = self.slot_discharge_coefficient[slot_index] + \
                (us_elevation-slot_invert_elevation)*self.discharge_coefficient_escalation
            #orificeVelocity = spec_discharge_coefficient*np.power(2.0*9.81*(deltaH-self.non_slot_loss), self.discharge_exponent)
            orificeVelocity = self.slot_discharge_coefficient[slot_index]*np.sqrt(2.0*9.81*deltaH)
            return orificeVelocity*self.get_area_for_elevation(slot_index,\
                (max([ds_elevation,slot_invert_elevation])+us_elevation)/2.0)

    def calc_upstream_water_level(self, slot_index, ds_elevation, flow_rate):
    	 
        if(flow_rate==0): return [ds_elevation, 0.0]
        
        dtype = [('flow',float),('elevation',float)]
        
        # Add zero flow point
        
        database = np.array([(0.0, ds_elevation)],dtype=dtype)
        
        # Add previous results to table
        
        if(slot_index<self.pool_levels.size):
            guessElevation = self.pool_levels[slot_index]
        else:
            guessElevation = self.upstream_water_level
            
        if(guessElevation>ds_elevation):
            guessFlow = self.slot_flow_rate(slot_index,\
                guessElevation, ds_elevation)
            database = np.append(database,\
                np.array([(guessFlow,guessElevation)], dtype=dtype))
        
        counter=0
        
        while(counter<200):
            
            # If within the existing data, interpolate, otherwise 
            # add another point
            
            if(flow_rate>=np.max(database['flow'])):
                guessElevation = np.max(database['elevation'])+1.0;
            else:
                guessElevation = np.interp(pow(flow_rate,2),\
                    np.square(database['flow']), database['elevation'])
            
            guessFlow = self.slot_flow_rate(slot_index,\
                guessElevation, ds_elevation)
                
            flowError = guessFlow-flow_rate
            
            if((abs(flowError)/flow_rate<0.0001) and (flow_rate>0)):
                return [guessElevation, flowError];
            
            database = np.append(database,\
                np.array([(guessFlow,guessElevation)], dtype=dtype))
                
            database = np.sort(database, order='elevation')
            counter = counter + 1;
            
        return [guessElevation, flowError]
    
    def update_flow(self, flow_rate):
        
        self.flow_rate = flow_rate;
        n_pools = self.pool_levels.size
        n_slots = self.slot_inverts.size
        supp_flow = np.sum(self.pool_supp_flows);
        
        # Calculate elevation of first pool, using the downstream
        # water level as a starting point
        
        [self.pool_levels[0], dummy] = self.calc_upstream_water_level(
            0,
            self.downstream_water_level,
            flow_rate+supp_flow)
        
        # Calculate elevations in intermediate pools
        
        for pool_index in range(1,n_pools):
            
            # Calculate the cumulative supplemental flow            
            
            supp_flow = np.sum(self.pool_supp_flows[pool_index:n_pools])
            
            # Calculate the water of this pool based on the last one            
            
            [self.pool_levels[pool_index], dummy] =\
                self.calc_upstream_water_level(
                    pool_index,
                    self.pool_levels[pool_index-1],
                    flow_rate+supp_flow)
        
        # Calculate estimate of upstream level
        
        return self.calc_upstream_water_level(
            n_slots-1,
            self.pool_levels[n_pools-1],
            flow_rate)
        
    def solve_hydraulics(self):
        
        start_time = time.time()
    
        # Mark fishway as run, even if it doesn't run with success
        
        self.has_been_run = True
        
        # Sanity check fishway configuration
        
        highest_slot_invert = np.min(self.slot_inverts)
        
        for search_slot in range(len(self.slot_inverts)):
            this_slot_type = self.vertical_slot_type[search_slot]
            this_slot_height = self.vertical_slot_heights[this_slot_type]
            this_slot_area = self.vertical_slot_areas[this_slot_type]
            this_slot_invert = this_slot_height[np.argmax(this_slot_area>0)-1] + self.slot_inverts[search_slot]
            if(this_slot_invert>highest_slot_invert):
                highest_slot_invert = this_slot_invert
                        
        if(self.upstream_water_level<=self.downstream_water_level):
            # Reverse/no flow
            self.has_passed_sanity_check = False;
            return;
        elif(self.upstream_water_level<highest_slot_invert+self.minimum_upstream_height_to_top_slot_invert_threshold):
            # Upstream level below maximum slot invert
            self.has_passed_sanity_check = False;
            return;
        elif(self.upstream_water_level<np.max(self.slot_inverts)):
            # Upstream level below maximum slot invert level
            self.has_passed_sanity_check = False;
            return;
        else:
            # Flow is going to go through the fishway
            self.has_passed_sanity_check = True;
            
        # Solve hydraulics, init with even water surface
        
        average_drop = (self.upstream_water_level -\
            self.downstream_water_level)/self.slot_inverts.size
            
        for pool_index in range(0,self.pool_levels.size):
            self.pool_levels[pool_index] = self.downstream_water_level +\
                (pool_index + 1)*average_drop
        
        # Calculate a guesstimate flow, based on the even water surface
        
        init_flows = np.zeros(self.slot_inverts.size)
        init_flows[0] = self.slot_flow_rate(0, self.pool_levels[0],\
            self.downstream_water_level)
            
        init_flows[self.slot_inverts.size-1] = self.slot_flow_rate(\
            self.slot_inverts.size-1,\
            self.upstream_water_level,\
            self.pool_levels[-1])
        
        for slot_index in range(1,self.slot_inverts.size-2):
            init_flows[slot_index] = self.slot_flow_rate(\
                slot_index, self.pool_levels[slot_index],\
                self.pool_levels[slot_index-1])
        
        flow_guess = np.mean(init_flows)
        if(flow_guess <= 0.25): flow_guess = 0.25
    
        [guessed_upstream_level, guessed_upstream_level_error] =\
            self.update_flow(flow_guess)
        
        dtype = [('flow',float),('upstreamLevel',float)]
        
        self.database = np.array([(flow_guess, guessed_upstream_level)],dtype=dtype)
        
        self.database = np.append(self.database,\
            np.array([(0.0, np.max(self.slot_inverts))], dtype=dtype))
        
        height_error = 1.0;
        counter = 0;
    
        while(counter<500 and height_error>0.0001):
            
            self.database = np.sort(self.database, order='flow');
            
            if(self.upstream_water_level>=np.max(self.database['upstreamLevel'])):
                flow_guess = np.max(self.database['flow'])*1.4
            elif(self.upstream_water_level<=np.min(self.database['upstreamLevel'])):
                flow_guess = np.min(self.database['flow'])*0.7
            else:
                flow_guess = np.interp(self.upstream_water_level,\
                    self.database['upstreamLevel'], self.database['flow'])
            
            [guessed_upstream_level, guessed_upstream_level_error] =\
                self.update_flow(flow_guess)
                
            self.database = np.append(self.database,\
                np.array([(flow_guess, guessed_upstream_level)], dtype=dtype))
            
            height_error = abs(guessed_upstream_level -\
                self.upstream_water_level)
                
            counter = counter + 1
            
        self.upstream_water_level_error = height_error
        
        # Calculate volumes and other derived arrays

        n_slots = self.slot_inverts.size
        n_pools = self.pool_levels.size
        
        self.pool_volumes = (self.pool_levels - self.pool_inverts) *\
            self.pool_areas        
        
        self.slot_differential = np.zeros(n_slots);
        
        self.slot_differential[0] = self.pool_levels[0] -\
            self.downstream_water_level
            
        self.slot_differential[n_slots-1] = self.upstream_water_level -\
            self.pool_levels[n_pools-1]
        
        self.slot_flow_depth = np.zeros(n_slots);        
                
        self.slot_flow_depth[0] = ((self.pool_levels[0] +\
            self.downstream_water_level)/2.0) - self.slot_inverts[0]
            
        self.slot_flow_depth[n_slots-1] = ((self.upstream_water_level +\
            self.pool_levels[n_pools-1])/2.0) - self.slot_inverts[n_slots-1]
    
        for slot_index in range(1,n_slots-1):
            self.slot_differential[slot_index] =\
                self.pool_levels[slot_index] - self.pool_levels[slot_index-1]
                
            self.slot_flow_depth[slot_index] =\
                ((self.pool_levels[slot_index] +\
                self.pool_levels[slot_index-1])/2.0) -\
                self.slot_inverts[slot_index]
        
        self.pool_diss_rate = np.zeros(n_pools);
        
        for pool_index in range(0,n_pools):
            
            total_us_slot_flow = self.get_slot_flow_rate(pool_index+1)
            
            if(self.pool_volumes[pool_index]>0):
            
                self.pool_diss_rate[pool_index] = ((water_density * gravity *\
                    total_us_slot_flow * self.slot_differential[pool_index+1]) +\
                    (water_density * 0.5 *\
                    np.power(self.pool_supp_velocities[pool_index],2) *\
                    self.pool_supp_flows[pool_index])) /\
                    self.pool_volumes[pool_index]
                    
            else:
                self.pool_diss_rate[pool_index] = np.nan
       
        self.time_to_solve = time.time() - start_time
    
    def get_slot_flow_rate(self, slot_index):
        n_pools = self.pool_levels.size
        return np.sum(self.pool_supp_flows[slot_index:n_pools]) +\
            self.flow_rate    
    
    def criteria_run_complete(self):
        if(self.upstream_water_level_error<=0.001 and\
                self.has_been_run and self.has_passed_sanity_check):
            return True
        else:
            return False
        
    def criteria_dissipation(self, diss_limit):
        if(self.criteria_run_complete() and\
                np.max(self.pool_diss_rate)<diss_limit):
            return True
        else:
            return False
        
    def criteria_min_entry_differential(self, diff_limit):
        if(self.criteria_run_complete() and\
                self.slot_differential[0]>diff_limit):
            return True
        else:
            return False

    def criteria_max_differential(self, diff_limit):
        if(self.criteria_run_complete() and\
                np.max(self.slot_differential)<diff_limit):
            return True
        else:
            return False
    
    def criteria_min_slot_level(self, level_limit):
        if(self.criteria_run_complete() and\
                np.min(self.slot_flow_depth)>level_limit):
            return True
        else:
            return False
        
    # Convenience method for calculating a combined slot height vs. area array from components
    @staticmethod
    def combineSubSlots(height_array_list, area_array_list, max_height):

        # Compile combined height array, starting with zero and max height in case it's not in any of the source arrays

        combined_heights = np.array([0.0, max_height])

        for i in range(len(area_array_list)):
            if len(height_array_list[i])>0:
                combined_heights = np.append(combined_heights, height_array_list[i])

        unique_heights = np.unique(combined_heights)

        total_area = np.zeros(len(unique_heights))

        for i in range(len(area_array_list)):
            if len(height_array_list[i])>0:
                total_area = total_area + np.interp(unique_heights, height_array_list[i], area_array_list[i])

        return unique_heights, total_area