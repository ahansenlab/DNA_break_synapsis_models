# When state changes from "apart" to "contact", the change is only register if a minimum duration of first_protein_arrival_time has elapsed

from polychrom.hdf5_format import list_URIs, load_URI
import pandas as pd
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
import pickle
import os

with open(snakemake.input.units_file,'r') as infile:
    units_file = json.load(infile)
        
dt_min = units_file['seconds_per_sim_dt']['internal TAD']/60

# parse file name for the number of chromosomes
folder_name = snakemake.input.simulation_folder

with open(snakemake.input.monomer_displacement_per_timestep,'rb') as infile:
    displacement = pickle.load(infile)   
    ends_displacement = displacement['ends_displacement']
    
URI_start = snakemake.config['DSB_simulation_parameters']["initial_URIs_skipped"]
last_URIs = snakemake.config['DSB_simulation_parameters']["last_URIs_used_for_synapsis"]
if "DSB_sim" in snakemake.input.units_file:
    URI_start = 1000 # restartBondUpdaterEveryBlocks * BondUpdaterCountBeforeDSB

R = snakemake.config['DSB_simulation_parameters']["contact_radii_for_first_passage_times"]
R += [1+np.around(np.mean(ends_displacement),decimals=2)] #use the one monomer size + step size the last contact radius to plot

first_protein_arrival_time = 1/dt_min # use ku70/80 arrival time,~1 minute to estimate when the two DSB ends become sticky (in simulation time unit)
print(first_protein_arrival_time)
colors = ['lime','cyan','magenta','orangered','royalblue']
# Get difference of two lists
def Diff(li1, li2):
    return list(set(li1) - set(li2)) + list(set(li2) - set(li1))

inner_encountered_separated_count = {r:0 for r in R} # number of encountered and separated events for the inner-chromosome end pairs
between_encountered_separated_count = {r:0 for r in R}  # number of encountered and separated events for the between-chromosome end pairs
inner_unrecaptured_time = {r:[] for r in R} # dict of lists of time window between exit and end of observation window for inner-chromosome end pairs
between_unrecaptured_time = {r:[] for r in R} # dict of list of time window between exit and end of observation window for between-chromosome end pairs

# get list of all pairwise positions of monomer ends 
files = list_URIs(snakemake.input.simulation_folder)

        
if 'multichain' in snakemake.input.simulation_folder:
        
    # look for simulations extended from the original run, if there is any
    num_chroms = int(folder_name.split('_C')[1].split('_collisionrate')[0])
    
    file_list = []
    end_point_dict={3:39999999,7:19999999,10:9999999,14:9999999}

    date = folder_name.split('sim_')[1].split('_dens')[0]
    parameter_string = folder_name.split(date)[1].split('_v5')[0]
    date = date.replace("_","")
    for path in os.listdir("results/"):
        if "ContinuedFrom"+date in path and parameter_string in path:
            file_list.extend(files[:end_point_dict[num_chroms]+1])
            additional_folder = path
            additional_files = list_URIs("results/"+additional_folder)
            file_list.extend(additional_files[:end_point_dict[num_chroms]+1])

            additional_date = additional_folder.split('sim_')[1].split('_dens')[0]
            additional_date = additional_date.replace("_","")
            for second_path in os.listdir("results/"):
                if "ContinuedFrom"+additional_date in second_path and parameter_string in second_path:
                    additional_folder = second_path
                    additional_files = list_URIs("results/"+additional_folder)
                    file_list.extend(additional_files)
        else:
            file_list.extend(files)

    initial_state = {r:None for r in R} # initialize
    old_contacts_dict = {r:None for r in R} # initialize
    
    all_index = np.arange((num_chroms)*(num_chroms*2-1)).tolist()
    inner_chrom_index = [num_chroms*i+(i+num_chroms)-((i+2)*(i+1))//2 for i in range(num_chroms)]
    between_chrom_index = Diff(inner_chrom_index,all_index)
        
    times_changed = {r:{idx:[0] for idx in range(num_chroms)} for r in R}
    constrained = {r:0 for r in R}
    unconstrained = {r:0 for r in R}
    if num_chroms == 3:
        last_URIs = last_URIs*5
    elif num_chroms == 5:
        last_URIs = last_URIs*4
    observation_window = len(files[-last_URIs:]) # the time window of observation
    
    fig1 = plt.figure(figsize=(6,15))
    ax_set1 = fig1.subplots(3, 1)
    
    fig2 = plt.figure(figsize=(6,15))
    ax_set2 = fig2.subplots(3, 1)
    
    for fi,f in enumerate(files[-last_URIs:]): 
        u = load_URI(f)  
        x0_ends = []
        y0_ends = []
        z0_ends = []
        x1_ends = []
        y1_ends = []
        z1_ends = []
        for c in range(0,num_chroms):
            x1_ends.append(u['pos'][num_chroms+c][0]) # get the chromosome end positions
            y1_ends.append(u['pos'][num_chroms+c][1]) # get the chromosome end positions
            z1_ends.append(u['pos'][num_chroms+c][2]) # get the chromosome end positions
            x0_ends.append(u['pos'][c][0]) # get the chromosome end positions
            y0_ends.append(u['pos'][c][1]) # get the chromosome end positions
            z0_ends.append(u['pos'][c][2]) # get the chromosome end positions     
            
            if fi<=600 and c<=3:
                ax_set1[0].scatter(fi,u['pos'][num_chroms+c][0],s=10,c=colors[c],marker='o')
                ax_set1[0].scatter(fi,u['pos'][c][0],s=10,c=colors[c],marker='s', alpha=0.5)
                ax_set1[1].scatter(fi,u['pos'][num_chroms+c][1],s=10,c=colors[c],marker='o')
                ax_set1[1].scatter(fi,u['pos'][c][1],s=10,c=colors[c],marker='s', alpha=0.5)
                ax_set1[2].scatter(fi,u['pos'][num_chroms+c][2],s=10,c=colors[c],marker='o')
                ax_set1[2].scatter(fi,u['pos'][c][2],s=10,c=colors[c],marker='s', alpha=0.5)

        # compute distances between all monomer ends    
        dists= np.sqrt((np.asarray(x0_ends)-np.asarray(x1_ends))**2+(np.asarray(y0_ends)-np.asarray(y1_ends))**2+(np.asarray(z0_ends)-np.asarray(z1_ends))**2)

        if fi == 0:
            # compute initial state
            # 1 if in contact, 0 if no contact
            for r in R:
                new_contacts = dists<r
                initial_state[r] =  new_contacts
                old_contacts_dict[r] = new_contacts

        for r in R:    
            new_contacts = (dists<r)
            old_contacts = old_contacts_dict[r]

            has_change = [ix for ix, x in enumerate(new_contacts^old_contacts) if x==True and ((old_contacts[ix]==False and (fi-times_changed[r][ix][-1])>=first_protein_arrival_time) or (old_contacts[ix]==True))]
            for idx in has_change:
                times_changed[r][idx].append(fi)
                ## only update the contact that has "truly" changed:
                old_contacts_dict[r][idx] = new_contacts[idx]


    # compute the time in each state
    no_contact_intervals = {r: {idx:[0] for idx in range(num_chroms)} for r in R}
    with_contact_intervals = {r: {idx:[0] for idx in range(num_chroms)} for r in R}
    for r in R:
        for idx in range(num_chroms):
            # compute interval lengths
            times_changed[r][idx].pop(0) # remove the 0 time point
            if len(times_changed[r][idx])>0:
                time_changed_array = np.asarray(times_changed[r][idx])
                last_inner_encountered_separated_index = np.argmax(time_changed_array>(observation_window-2*60/dt_min))
                intervals = np.diff(times_changed[r][idx][:last_inner_encountered_separated_index+1]).astype(int)
            else:
                intervals = np.asarray([])

            # split intervals into contact and no contact
            if initial_state[r][idx] == True: # has contact at time zero
                with_contact_intervals[r][idx] = intervals[1::2] # first interval does not contain time zero
                no_contact_intervals[r][idx] = intervals[0::2] # first interval does not contain time zero
                if len(times_changed[r][idx])>0:
                    inner_encountered_separated_count[r] += len(times_changed[r][idx][:last_inner_encountered_separated_index][0::2])
                    if np.mod(len(times_changed[r][idx]),2)==1:
                        inner_unrecaptured_time[r].append(observation_window-times_changed[r][idx][-1])

            else:
                with_contact_intervals[r][idx] = intervals[0::2]            
                no_contact_intervals[r][idx] = intervals[1::2]
                if len(intervals)>0:
                    inner_encountered_separated_count[r] += len(times_changed[r][idx][:last_inner_encountered_separated_index][1::2])
                    if np.mod(len(times_changed[r][idx]),2)==0:
                        inner_unrecaptured_time[r].append(observation_window-times_changed[r][idx][-1])
                
elif 'brokenloop' in snakemake.input.simulation_folder:
    
    num_repeats = int(folder_name.split('_Repeat')[1].split('_DSB')[0])
    num_breaks = int(int(folder_name.split('_DSB')[1].split('_LoopSize')[0])/num_repeats) # number of DSBs per repeat
    
    idx_length = num_breaks*num_repeats
    times_changed = {r:{idx:[0] for idx in range(idx_length)} for r in R}
    
    initial_state = {r:None for r in R} # initialize
    old_contacts_dict = {r:None for r in R}# initialize
    observation_window = len(files[-last_URIs:]) # the time window of observation
    for fi,f in enumerate(files[-last_URIs:]): 
        u = load_URI(f)  
        x0_ends = []
        y0_ends = []
        z0_ends = []
        x1_ends = []
        y1_ends = []
        z1_ends = []
        for c in range(0,num_repeats):
            for i in range(0,num_breaks):
                x0_ends.append(u['pos'][c*2*num_breaks+i][0]) # get the broken DSB end positions
                y0_ends.append(u['pos'][c*2*num_breaks+i][1]) # get the broken DSB end positions
                z0_ends.append(u['pos'][c*2*num_breaks+i][2]) # get the broken DSB end positions
                x1_ends.append(u['pos'][c*2*num_breaks+num_breaks+i][0]) # get the broken DSB end positions
                y1_ends.append(u['pos'][c*2*num_breaks+num_breaks+i][1]) # get the broken DSB end positions
                z1_ends.append(u['pos'][c*2*num_breaks+num_breaks+i][2]) # get the broken DSB end positions        
         

        # compute distances between all monomer ends
        dists = np.sqrt((np.asarray(x0_ends)-np.asarray(x1_ends))**2+(np.asarray(y0_ends)-np.asarray(y1_ends))**2+(np.asarray(z0_ends)-np.asarray(z1_ends))**2)
        
        if fi == 0:
            # compute initial state
            # 1 if in contact, 0 if no contact
            for r in R:
                new_contacts = dists<r
                initial_state[r] = new_contacts
                old_contacts_dict[r] = new_contacts

        for r in R:    
            new_contacts = (dists<r)
            old_contacts = old_contacts_dict[r]

            has_change = [ix for ix, x in enumerate(new_contacts^old_contacts) if x==True and ((old_contacts[ix]==False and (fi-times_changed[r][ix][-1])>=first_protein_arrival_time) or (old_contacts[ix]==True))]
            for idx in has_change:
                times_changed[r][idx].append(fi)
                ## only update the contact that has "truly" changed:
                old_contacts_dict[r][idx] = new_contacts[idx]


    # compute the time in each state
    no_contact_intervals = {r:{idx:[0] for idx in range(idx_length)} for r in R}
    with_contact_intervals = {r:{idx:[0] for idx in range(idx_length)} for r in R}
    for r in R:
        for idx in range(idx_length):
            # compute interval lengths
            times_changed[r][idx].pop(0) # remove the 0 time point
            intervals = np.diff(times_changed[r][idx]).astype(int)

            # split intervals into contact and no contact
            if initial_state[r][idx] == True: # has contact at time zero
                with_contact_intervals[r][idx] = intervals[1::2] # first interval does not contain time zero
                no_contact_intervals[r][idx] = intervals[0::2] # first interval does not contain time zero
                if len(times_changed[r][idx])>0:
                    inner_encountered_separated_count[r] += len(times_changed[r][idx][0::2])
                    if np.mod(len(times_changed[r][idx]),2)==1:
                        inner_unrecaptured_time[r].append(observation_window-times_changed[r][idx][-1])
            else:
                with_contact_intervals[r][idx] = intervals[0::2]            
                no_contact_intervals[r][idx] = intervals[1::2]
                if len(intervals)>0:
                    inner_encountered_separated_count[r] += len(times_changed[r][idx][1::2])
                    if np.mod(len(times_changed[r][idx]),2)==0:
                        inner_unrecaptured_time[r].append(observation_window-times_changed[r][idx][-1])
                        
elif 'DSB_sim' in snakemake.input.simulation_folder:
    
    with open(snakemake.input.simulation_folder+'/DSB_boundary_coordinates.npy',"rb") as f:
        DSB_coordinates = np.load(f)
        boundary_coordinates = np.load(f)
    initial_state = {r:None for r in R} # initialize
    old_contacts_dict = {r:None for r in R} # initialize
    initial_state = {r:None for r in R} # initialize
    num_DSBs = int(int(folder_name.split('_DSB')[1].split('_collisionrate')[0]))
    
    if num_DSBs==0:
        num_DSBs = 23
       
    times_changed = {r:{idx:[0] for idx in range(num_DSBs)} for r in R}
    constrained = {r:0 for r in R}
    unconstrained = {r:0 for r in R}
    synapsed = {r:{idx:False for idx in range(num_DSBs)} for r in R}
    observation_window = len(files[URI_start:]) # the time window of observation

    fig1 = plt.figure(figsize=(6,15))
    ax_set1 = fig1.subplots(3, 1)
    
    fig2 = plt.figure(figsize=(6,15))
    ax_set2 = fig2.subplots(3, 1)
    
    for fi,f in enumerate(files[URI_start:]): 
#         if np.mod(fi-56,96)!=0:
        u = load_URI(f)  
        x0_ends = []
        y0_ends = []
        z0_ends = []
        x1_ends = []
        y1_ends = []
        z1_ends = []
        if num_DSBs>0:
            for c in range(0,num_DSBs):
                x0_ends.append(u['pos'][2*c][0]) # get the broken DSB end positions
                y0_ends.append(u['pos'][2*c][1]) # get the broken DSB end positions
                z0_ends.append(u['pos'][2*c][2]) # get the broken DSB end positions
                x1_ends.append(u['pos'][2*c+1][0]) # get the broken DSB end positions
                y1_ends.append(u['pos'][2*c+1][1]) # get the broken DSB end positions
                z1_ends.append(u['pos'][2*c+1][2]) # get the broken DSB end positions    

                if fi<=600 and c<=3:
                    ax_set1[0].scatter(fi,u['pos'][2*c][0],s=10,c=colors[c],marker='o')
                    ax_set1[0].scatter(fi,u['pos'][2*c+1][0],s=10,c=colors[c],marker='s', alpha=0.5)
                    ax_set1[1].scatter(fi,u['pos'][2*c][1],s=10,c=colors[c],marker='o')
                    ax_set1[1].scatter(fi,u['pos'][2*c+1][1],s=10,c=colors[c],marker='s', alpha=0.5)
                    ax_set1[2].scatter(fi,u['pos'][2*c][2],s=10,c=colors[c],marker='o')
                    ax_set1[2].scatter(fi,u['pos'][2*c+1][2],s=10,c=colors[c],marker='s', alpha=0.5)
        else:
             for c in range(0,23):
                x0_ends.append(u['pos'][2*c][0]) # get the broken DSB end positions
                y0_ends.append(u['pos'][2*c][1]) # get the broken DSB end positions
                z0_ends.append(u['pos'][2*c][2]) # get the broken DSB end positions
                x1_ends.append(u['pos'][2*c+1][0]) # get the broken DSB end positions
                y1_ends.append(u['pos'][2*c+1][1]) # get the broken DSB end positions
                z1_ends.append(u['pos'][2*c+1][2]) # get the broken DSB end positions  


        # compute distances between all monomer ends
        dists= np.sqrt((np.asarray(x0_ends)-np.asarray(x1_ends))**2+(np.asarray(y0_ends)-np.asarray(y1_ends))**2+(np.asarray(z0_ends)-np.asarray(z1_ends))**2)

        for r in R:    
            new_contacts = (dists<r)
            has_change = [ix for ix, x in enumerate(new_contacts) if x==True and synapsed[r][ix]==False and (fi-times_changed[r][ix][-1])>=first_protein_arrival_time]
            for idx in has_change:
                times_changed[r][idx].append(fi)
                synapsed[r][idx]= True
                left_extruders = u['SMCs'][:,0]
                right_extruders = u['SMCs'][:,1]
                break_site = DSB_coordinates[idx*2]
                test_left = left_extruders - break_site
                test_right = right_extruders - (break_site+1)
                test = np.multiply(test_left,test_right) 
                constraining_LEF_index = np.argwhere(test<0)
                if len(constraining_LEF_index)>0:
                    constrained[r] += 1
                else:
                    unconstrained[r] += 1
                    
    # check the DSB end movement 
    for fi,f in enumerate(files[URI_start-100:URI_start]): 
        u = load_URI(f)
        for c in range(0,num_DSBs):
            if fi<=300 and c<=3:
                ax_set2[0].scatter(fi,u['pos'][2*c][0],s=10,c=colors[c],marker='o')
                ax_set2[0].scatter(fi,u['pos'][2*c+1][0],s=10,c=colors[c],marker='s', alpha=0.5)
                ax_set2[1].scatter(fi,u['pos'][2*c][1],s=10,c=colors[c],marker='o')
                ax_set2[1].scatter(fi,u['pos'][2*c+1][1],s=10,c=colors[c],marker='s', alpha=0.5)
                ax_set2[2].scatter(fi,u['pos'][2*c][2],s=10,c=colors[c],marker='o')
                ax_set2[2].scatter(fi,u['pos'][2*c+1][2],s=10,c=colors[c],marker='s', alpha=0.5)

    # compute the time in each state
    no_contact_intervals = {r: {idx:[0] for idx in range(num_DSBs)} for r in R}
    with_contact_intervals = {r: {idx:[0] for idx in range(num_DSBs)} for r in R}
    for r in R:
        for idx in range(num_DSBs):
            # compute interval lengths
            intervals = np.diff(times_changed[r][idx]).astype(int)
            no_contact_intervals[r][idx] = intervals[0::2] 
            if len(times_changed[r][idx])==1:
                inner_unrecaptured_time[r].append(observation_window-times_changed[r][idx][-1])
                    
with open(snakemake.output.first_passage_times_v2,'wb') as outfile:
    pickle.dump({'initial_state':initial_state, 
                 'times_changed': times_changed, 
                 'no_contact_intervals':no_contact_intervals,
                 'with_contact_intervals':with_contact_intervals,
                'inner_encountered_separated_count':inner_encountered_separated_count,
                'between_encountered_separated_count':between_encountered_separated_count,
                'inner_unrecaptured_time':inner_unrecaptured_time,
                'between_unrecaptured_time':between_unrecaptured_time,
                'constrained':constrained,
                'unconstrained':unconstrained},
                outfile)               
fig1.savefig(snakemake.output.DSB_ends_traces,bbox_inches='tight')
fig2.savefig(snakemake.output.DSB_ends_traces_preDSB,bbox_inches='tight')
   