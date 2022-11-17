import json
import numpy as np
from scipy import special, stats
from matplotlib import pyplot as plt

import tracklib as tl
from tracklib.analysis import bild, msdfit
from tracklib.analysis.msdfit.lib import TwoLocusRouseFit 

from MSD_fit_shared import get_named_params, run_fit, plot_fit_results

if 'DSB_sim' not in snakemake.input[0] and 'C25'not in snakemake.input[0]:
    # load in the experimental trajectory   
    tag = snakemake.wildcards.group
    data = tl.TaggedSet()

    # combine all the data from the sample group into a single tagged set
    for tag_file in snakemake.input.tagged_set:
        this_data = tl.io.load.csv(tag_file, 
                                ['id', 't', 'x', 'y', 'z', 'x2', 'y2', 'z2'], 
                                skip_header=1,delimiter='\t',tags=[tag])

        data |= this_data

    # compute the displacements in time
    data = data.process(lambda traj: traj.relative())    

    # fit the MSD
    print("Running fit for", tag)
    mci, profiler = run_fit(data)

    # retrieve the best fit parameters
    params = profiler.point_estimate['params']

    # save the fit parameters as a text readable json file
    out = get_named_params(params)
    #        'log_confidence_intervals':mci,
    #       'logL':profiler.point_estimate['logL']}

    # plot the fit results
    plot_fit_results(data, mci, profiler, snakemake.output.figure, title = tag)
    fig, axs = plt.subplots(1, 2, figsize=[13, 4])
    plt.savefig(snakemake.output.place_holder_figure,bbox_inches='tight')
    plt.close()
    
else:
    out={'msd_fit':'unit conversion already done in calibration simulation'}
    fig, axs = plt.subplots(1, 2, figsize=[13, 4])
    plt.savefig(snakemake.output.figure,bbox_inches='tight')
    plt.close()
    fig, axs = plt.subplots(1, 2, figsize=[13, 4])
    plt.savefig(snakemake.output.place_holder_figure,bbox_inches='tight')
    plt.close()
with open(snakemake.output.fit_result, 'w') as outfile:
    json.dump(out, outfile)         
       


