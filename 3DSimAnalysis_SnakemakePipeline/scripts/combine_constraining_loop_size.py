# When state changes from "apart" to "contact", the change is only register if a minimum duration of first_protein_arrival_time has elapsed

from polychrom.hdf5_format import list_URIs, load_URI
import pandas as pd
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist
import pickle
import os

# list for combining constraining LEFs size from independent runs 
all_constraining_LEF_size = []
for i, constraining_LEF_size_file in enumerate(snakemake.input.constraining_loop_size):
    file = constraining_LEF_size_file
    with open(file, 'rb') as infile:
        constraining_LEF_size = np.load(infile)
    constraining_LEF_size_list = constraining_LEF_size.tolist()
    all_constraining_LEF_size.extend(constraining_LEF_size_list)

file = snakemake.output.all_constraining_loop_size
with open(file, 'wb') as g:
    np.save(g, np.asarray(all_constraining_LEF_size))#capture radius 4
