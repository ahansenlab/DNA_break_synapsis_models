# When state changes from "apart" to "contact", the change is only register if a minimum duration of first_protein_arrival_time has elapsed

from polychrom.hdf5_format import list_URIs, load_URI
import pandas as pd
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
import pickle
import os


# parse file name for the number of chromosomes
folder_name = snakemake.input.simulation_folder

# get list of all pairwise positions of monomer ends 
files = list_URIs(snakemake.input.simulation_folder)

# list for constraining LEFs
constraining_LEF_size = []
if 'DSB_sim' in snakemake.input.simulation_folder:
    URI_start = 1000 # restartBondUpdaterEveryBlocks * BondUpdaterCountBeforeDSB
    with open(snakemake.input.simulation_folder+'/DSB_boundary_coordinates.npy',"rb") as f:
        DSB_coordinates = np.load(f)
        boundary_coordinates = np.load(f)

    num_DSBs = int(int(folder_name.split('_DSB')[1].split('_collisionrate')[0]))
       
    f = files[URI_start]
    u = load_URI(f)
    left_extruders = u['SMCs'][:,0]
    right_extruders = u['SMCs'][:,1]
    for i in range(num_DSBs):
        break_site = DSB_coordinates[i*2]
        test_left = left_extruders - break_site
        test_right = right_extruders - (break_site+1)
        test = np.multiply(test_left,test_right) 
        constraining_LEF_index = np.argwhere(test<0)
        # identify the inner most constraining LEF
        loop_size = np.abs(right_extruders[constraining_LEF_index] - left_extruders[constraining_LEF_index])
        if len(constraining_LEF_index)>0:
            innermost_index = np.argmin(loop_size)
            constraining_left = left_extruders[constraining_LEF_index[innermost_index]] 
            constraining_right = right_extruders[constraining_LEF_index[innermost_index]] 
            constraining_LEF_size.append(constraining_right-constraining_left)
    
file = snakemake.output.constraining_loop_size
with open(file, 'wb') as g:
    np.save(g, np.asarray(constraining_LEF_size))#capture radius 4
