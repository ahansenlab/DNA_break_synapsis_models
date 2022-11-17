import json
import numpy as np
from scipy import special, stats
from matplotlib import pyplot as plt

import tracklib as tl
from tracklib.analysis import bild, msdfit
from tracklib.analysis.msdfit.lib import TwoLocusRouseFit 


from MSD_fit_shared import get_named_params, run_fit, plot_fit_results

import pandas as pd

import io

tag_name = 'chromosome ends'
in_file = snakemake.input[0]
out_figure_files = [f for f in snakemake.output[0]]
out_file = snakemake.output[0]

single_locus_data = tl.TaggedSet()
single_locus_df = pd.read_csv(in_file,usecols=['id', 't', 'x0', 'y0', 'z0'])
single_locus_df['x1']=0
single_locus_df['y1']=0
single_locus_df['z1']=0
csv_stream = io.StringIO(single_locus_df.to_csv())
this_data = tl.io.load.csv(csv_stream,['None','id', 't', 'x', 'y', 'z', 'x2', 'y2', 'z2'],delimiter=',', skip_header=1,tags=[tag_name])
single_locus_data |= this_data

# combine the simulation data into a single tagged set
two_loci_data = tl.TaggedSet()
this_data = tl.io.load.csv(in_file, 
                        ['id', 't', 'x', 'y', 'z', 'x2', 'y2', 'z2'], 
                        skip_header=1,
                        delimiter=',',
                        tags=[tag_name])
two_loci_data |= this_data 

# compute the displacements in time   
single_locus_data = single_locus_data.process(lambda traj: traj.relative())   
two_loci_data = two_loci_data.process(lambda traj: traj.relative())    

fig1 = plt.figure(figsize=(6,4))

single_locus_msd = tl.analysis.MSD(single_locus_data)
two_loci_msd = tl.analysis.MSD(two_loci_data)
single_locus_dt = np.arange(1, len(single_locus_msd))
two_loci_dt = np.arange(1, len(two_loci_msd))


plt.plot(single_locus_dt, single_locus_msd[1:], color='b',
            linewidth=2, 
            label='single-locus',
           )

plt.plot(two_loci_dt, two_loci_msd[1:], color='g',
             linewidth=2,
            label='two-loci',
           )

plt.legend()
plt.xscale('log')
plt.xlabel('time')
plt.yscale('log')
plt.ylabel('MSD')

fig1.savefig(out_file,bbox_inches='tight')