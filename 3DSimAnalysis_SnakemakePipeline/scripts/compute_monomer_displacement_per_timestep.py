from polychrom.hdf5_format import list_URIs, load_URI
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
import pickle

# parse file name for the number of chromosomes
folder_name = snakemake.input.simulation_folder


URI_start = snakemake.config['DSB_simulation_parameters']["initial_URIs_skipped"]
if "DSB_sim" or "frozenloop_stochastic" in snakemake.input.simulation_folder:
    URI_start = 1000
R = snakemake.config['DSB_simulation_parameters']["contact_radii_for_first_passage_times"]

ends_displacement = []  # initialize
internal_displacement = [] # initialize

# get list of all pairwise positions of monomer ends 
files = list_URIs(snakemake.input.simulation_folder)

if 'multichain' in snakemake.input.simulation_folder:
    num_chroms = int(folder_name.split('_C')[1].split('_collisionrate')[0])
    
    if num_chroms>50:
        for fi,f in enumerate(files[URI_start:]): 
            u = load_URI(f)  
            x0_ends = []
            y0_ends = []
            z0_ends = []
            x0_internal = []
            y0_internal = []
            z0_internal = []
            for c in range(0,num_chroms):
                x0_ends.append(u['pos'][num_chroms+c][0]) # get the chromosome end positions
                y0_ends.append(u['pos'][num_chroms+c][1]) # get the chromosome end positions
                z0_ends.append(u['pos'][num_chroms+c][2]) # get the chromosome end positions
                x0_ends.append(u['pos'][c][0]) # get the chromosome end positions
                y0_ends.append(u['pos'][c][1]) # get the chromosome end positions
                z0_ends.append(u['pos'][c][2]) # get the chromosome end positions 
                x0_internal.append(u['pos'][num_chroms+c][0]) # get the chromosome internal positions
                y0_internal.append(u['pos'][num_chroms+c][1]) # get the chromosome internal positions
                z0_internal.append(u['pos'][num_chroms+c][2]) # get the chromosome internal positions
                x0_internal.append(u['pos'][c][0]) # get the chromosome internal positions
                y0_internal.append(u['pos'][c][1]) # get the chromosome internal positions
                z0_internal.append(u['pos'][c][2]) # get the chromosome internal positions        


            if fi != 0:
                ends_displacement.append(((old_x0_ends-x0_ends)**2+(old_y0_ends-y0_ends)**2+(old_z0_ends-z0_ends)**2)**0.5)
                internal_displacement.append(((old_x0_internal-x0_internal)**2+(old_y0_internal-y0_internal)**2+(old_z0_internal-z0_internal)**2)**0.5)

            old_x0_ends = np.copy(x0_ends)
            old_y0_ends = np.copy(y0_ends)
            old_z0_ends = np.copy(z0_ends)
            old_x0_internal = np.copy(x0_internal)
            old_y0_internal = np.copy(y0_internal)
            old_z0_internal = np.copy(z0_internal)
    else:
        for fi,f in enumerate(files[-10000:]): 
            u = load_URI(f)  
            x0_ends = []
            y0_ends = []
            z0_ends = []
            x0_internal = []
            y0_internal = []
            z0_internal = []
            for c in range(0,num_chroms):
                x0_ends.append(u['pos'][num_chroms+c][0]) # get the chromosome end positions
                y0_ends.append(u['pos'][num_chroms+c][1]) # get the chromosome end positions
                z0_ends.append(u['pos'][num_chroms+c][2]) # get the chromosome end positions
                x0_ends.append(u['pos'][c][0]) # get the chromosome end positions
                y0_ends.append(u['pos'][c][1]) # get the chromosome end positions
                z0_ends.append(u['pos'][c][2]) # get the chromosome end positions 
                x0_internal.append(u['pos'][num_chroms*2+c][0]) # get the chromosome internal positions
                y0_internal.append(u['pos'][num_chroms*2+c][1]) # get the chromosome internal positions
                z0_internal.append(u['pos'][num_chroms*2+c][2]) # get the chromosome internal positions
                x0_internal.append(u['pos'][num_chroms*3+c][0]) # get the chromosome internal positions
                y0_internal.append(u['pos'][num_chroms*3+c][1]) # get the chromosome internal positions
                z0_internal.append(u['pos'][num_chroms*3+c][2]) # get the chromosome internal positions        


            if fi != 0:
                ends_displacement.append(((old_x0_ends-x0_ends)**2+(old_y0_ends-y0_ends)**2+(old_z0_ends-z0_ends)**2)**0.5)
                internal_displacement.append(((old_x0_internal-x0_internal)**2+(old_y0_internal-y0_internal)**2+(old_z0_internal-z0_internal)**2)**0.5)

            old_x0_ends = np.copy(x0_ends)
            old_y0_ends = np.copy(y0_ends)
            old_z0_ends = np.copy(z0_ends)
            old_x0_internal = np.copy(x0_internal)
            old_y0_internal = np.copy(y0_internal)
            old_z0_internal = np.copy(z0_internal)

elif 'DSB_sim' or 'frozenloop' in snakemake.input.simulation_folder:
    with open(snakemake.input.simulation_folder+'/DSB_boundary_coordinates.npy',"rb") as f:
        DSB_coordinates = np.load(f)
    
    if 'frozenloop_sim' in snakemake.input.simulation_folder:
        URI_start = 3380000
    num_DSBs = int(len(DSB_coordinates)/2)
    for fi,f in enumerate(files[URI_start:]): 
        u = load_URI(f)  
        x0_ends = []
        y0_ends = []
        z0_ends = []
        x0_internal = []
        y0_internal = []
        z0_internal = []
        if num_DSBs>0:
            for c in range(0,num_DSBs):
                x0_ends.append(u['pos'][2*c][0]) # get the chromosome end positions
                y0_ends.append(u['pos'][2*c][1]) # get the chromosome end positions
                z0_ends.append(u['pos'][2*c][2]) # get the chromosome end positions
                x0_ends.append(u['pos'][2*c+1][0]) # get the chromosome end positions
                y0_ends.append(u['pos'][2*c+1][1]) # get the chromosome end positions
                z0_ends.append(u['pos'][2*c+1][2]) # get the chromosome end positions 
                x0_internal.append(u['pos'][2*c][0]) # get the chromosome internal positions
                y0_internal.append(u['pos'][2*c][1]) # get the chromosome internal positions
                z0_internal.append(u['pos'][2*c][2]) # get the chromosome internal positions
                x0_internal.append(u['pos'][2*c+1][0]) # get the chromosome internal positions
                y0_internal.append(u['pos'][2*c+1][1]) # get the chromosome internal positions
                z0_internal.append(u['pos'][2*c+1][2]) # get the chromosome internal positions  
        else:
            num_fbn2_tad = 35 # number of fbn2 tad tracked
            u_test = load_URI(files[URI_start]) 
            max_index = min(num_fbn2_tad,int(len(u_test['pos'])/2))
            for c in range(0,max_index):
                x0_ends.append(u['pos'][-2*c-1][0]) # get the chromosome end positions
                y0_ends.append(u['pos'][-2*c-1][1]) # get the chromosome end positions
                z0_ends.append(u['pos'][-2*c-1][2]) # get the chromosome end positions
                x0_ends.append(u['pos'][-2*c-2][0]) # get the chromosome end positions
                y0_ends.append(u['pos'][-2*c-2][1]) # get the chromosome end positions
                z0_ends.append(u['pos'][-2*c-2][2]) # get the chromosome end positions                
                x0_internal.append(u['pos'][-2*c-1][0]) # get the chromosome end positions
                y0_internal.append(u['pos'][-2*c-1][1]) # get the chromosome end positions
                z0_internal.append(u['pos'][-2*c-1][2]) # get the chromosome end positions
                x0_internal.append(u['pos'][-2*c-2][0]) # get the chromosome end positions
                y0_internal.append(u['pos'][-2*c-2][1]) # get the chromosome end positions
                z0_internal.append(u['pos'][-2*c-2][2]) # get the chromosome end positions  



        if fi != 0:
            ends_displacement.append(((old_x0_ends-x0_ends)**2+(old_y0_ends-y0_ends)**2+(old_z0_ends-z0_ends)**2)**0.5)
            internal_displacement.append(((old_x0_internal-x0_internal)**2+(old_y0_internal-y0_internal)**2+(old_z0_internal-z0_internal)**2)**0.5)

        old_x0_ends = np.copy(x0_ends)
        old_y0_ends = np.copy(y0_ends)
        old_z0_ends = np.copy(z0_ends)
        old_x0_internal = np.copy(x0_internal)
        old_y0_internal = np.copy(y0_internal)
        old_z0_internal = np.copy(z0_internal)
        
elif 'brokenloop' in snakemake.input.simulation_folder:
    
    num_breaks = int(folder_name.split('_DSB')[1].split('_LoopSize')[0]) # number of DSBs
    URI_start = 3380000
    for fi,f in enumerate(files[URI_start:]):  
        u = load_URI(f)  
        x0_ends = []
        y0_ends = []
        z0_ends = []
        x1_ends = []
        y1_ends = []
        z1_ends = []
        for c in range(0,num_breaks):
            x0_ends.append(u['pos'][2*c][0]) # get the broken DSB end positions
            y0_ends.append(u['pos'][2*c][1]) # get the broken DSB end positions
            z0_ends.append(u['pos'][2*c][2]) # get the broken DSB end positions
            x1_ends.append(u['pos'][2*c+1][0]) # get the broken DSB end positions
            y1_ends.append(u['pos'][2*c+1][1]) # get the broken DSB end positions
            z1_ends.append(u['pos'][2*c+1][2]) # get the broken DSB end positions           


        if fi != 0:
            ends_displacement.append(((old_x0_ends-x0_ends)**2+(old_y0_ends-y0_ends)**2+(old_z0_ends-z0_ends)**2)**0.5)
            ends_displacement.append(((old_x1_ends-x1_ends)**2+(old_y1_ends-y1_ends)**2+(old_z1_ends-z1_ends)**2)**0.5)
            

        old_x0_ends = np.copy(x0_ends)
        old_y0_ends = np.copy(y0_ends)
        old_z0_ends = np.copy(z0_ends)
        old_x1_ends = np.copy(x1_ends)
        old_y1_ends = np.copy(y1_ends)
        old_z1_ends = np.copy(z1_ends)
    
with open(snakemake.output.monomer_displacement_per_timestep,'wb') as outfile:
    pickle.dump({'ends_displacement':ends_displacement, 
                 'internal_displacement': internal_displacement},
                outfile)               
        
   