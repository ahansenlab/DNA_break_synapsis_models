import json
import numpy as np
from scipy import special, stats
from matplotlib import pyplot as plt

import tracklib as tl
from tracklib.analysis import bild, msdfit
from tracklib.analysis.msdfit.lib import TwoLocusRouseFit 

from MSD_fit_shared import get_named_params, run_fit, plot_fit_results

tag_names = 'internal TAD'

in_files = snakemake.input[0]
out_figure_files = snakemake.output[0]
out_fit_file = snakemake.output[-1]

if 'DSB_sim' not in in_files[0] and 'C25' not in in_files[0]:
    # combine the simulation data into a single tagged set
    data = tl.TaggedSet()
    this_data = tl.io.load.csv(in_files, 
                            ['id', 't', 'x', 'y', 'z', 'x2', 'y2', 'z2'], 
                            skip_header=1,
                            delimiter=',',
                               tags=[tag_names])
    data |= this_data 

    # compute the displacements in time
    data = data.process(lambda traj: traj.relative())    

    out = {}
    data.makeSelection()
    data.makeSelection(tags=tag_names)
    print("Running fit for", tag_names)

    # fit the MSD
    mci, profiler = run_fit(data, fit_with_localization_error=False)  
    params = profiler.point_estimate['params']
    out[tag_names] = get_named_params(params,fit_with_localization_error=False)

    # plot the fit results
    plot_fit_results(data, mci, profiler, out_figure_files, title = tag_names, fit_with_localization_error=False)
else:
    out={'msd_fit':'unit conversion already done in calibration simulation'}
    fig, axs = plt.subplots(1, 2, figsize=[13, 4])
    plt.savefig(out_figure_files,bbox_inches='tight')
    plt.close()

with open(out_fit_file, 'w') as outfile:
    json.dump(out, outfile)         



