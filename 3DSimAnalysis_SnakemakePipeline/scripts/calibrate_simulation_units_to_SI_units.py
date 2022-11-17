import os
import json
import numpy as np

if 'multichain_sim' in snakemake.input[0] and '_C1_'in snakemake.input[0]:
    expt_file = snakemake.input.expt_fit
    sim_file = snakemake.input.sim_fit

    MSD_fit_expt =  json.load(open(expt_file,'r'))
    MSD_fit_sim = json.load(open(sim_file,'r'))

    dt_expt = snakemake.config['DSB_simulation_parameters']['experiment_time_per_sample_in_sec'] # seconds
    dx_expt = snakemake.config['DSB_simulation_parameters']['experiment_spatial_units_in_nanometers']

    expt_steady_state = np.mean(MSD_fit_expt['Js_1d']) # in microns squared
    expt_diffusivity = np.sum(MSD_fit_expt['Gs_1d']) # in microns squared


    dx = {}
    dt = {}

    tag = 'internal TAD'
    sim_steady_state = np.mean(MSD_fit_sim[tag]['Js_1d']) # in monomers squared
    sim_diffusivity = np.sum(MSD_fit_sim[tag]['Gs_1d']) # in monomers squared
    dx_sim = np.sqrt(expt_steady_state/sim_steady_state)*dx_expt # nm per monomer
    dt_sim = ( (sim_diffusivity*dx_sim**2)/(expt_diffusivity*dx_expt**2) )**2 * dt_expt # secs per sim MSD sample
    dt_sim = dt_sim/snakemake.config['DSB_simulation_parameters']['skip_every_X_URIs'] # secs per sim step

    dt[tag] = dt_sim
    dx[tag] = dx_sim
else:
    calibration_units_file_dict ={70000:'multichain_sim_dens0.2_N70000_C1_collisionrate_1_trunc50_stepsperblock30.unit_conversion.json'}
    if 'multichain_sim' in snakemake.input[0]:
        num_monomers = int(snakemake.input[0].split('_N')[1].split('_C')[0])
    elif 'DSB_sim' or 'frozenloop_stochasticConstrainingLEF' in snakemake.input[0]:
        num_monomers = int(snakemake.input[0].split('_N')[1].split('_DSB')[0])
    elif 'brokenloop' in snakemake.input[0]:
        num_monomers = int(snakemake.input[0].split('_N')[1].split('_Repeat')[0])
    elif 'frozenloop_sim' in snakemake.input[0]:
        num_monomers = int(snakemake.input[0].split('_N')[1].split('_bondsFrom')[0])
    calibration_units_file_name = '../results/ΔRAD21/'+calibration_units_file_dict[num_monomers]
    
    with open(calibration_units_file_name,'r') as infile:
        units_file = json.load(infile)
    if 'DSB_sim' in snakemake.input[0] or 'frozenloop_stochasticConstrainingLEF' in snakemake.input[0]:
        steps_per_time_step = int(snakemake.input[0].split('_stepsPerBlock')[1].split('_trunc')[0])
    elif 'multichain_sim' in snakemake.input[0]:
        steps_per_time_step = int(snakemake.input[0].split('_stepsperblock')[1].split('_v5')[0])
    elif 'brokenloop' in snakemake.input[0]:
        steps_per_time_step = 30
    elif 'frozenloop_sim' in snakemake.input[0]:
        steps_per_time_step = int(snakemake.input[0].split('_stepsperblock')[1].split('_v5')[0])
    dx = {}
    dt = {}
    dt['internal TAD'] = units_file['seconds_per_sim_dt']['internal TAD']*steps_per_time_step/30 # in calibration simulation the step size is 30
    dx['internal TAD'] =  units_file['nanometers_per_monomer']['internal TAD']

with open(snakemake.output.units_file,'w') as outfile:
    json.dump({'nanometers_per_monomer':dx, 'seconds_per_sim_dt':dt},outfile)