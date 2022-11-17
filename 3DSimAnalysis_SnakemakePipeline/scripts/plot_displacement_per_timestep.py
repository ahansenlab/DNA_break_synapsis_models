import json
import pickle
import matplotlib.pyplot as plt
import numpy as np

with open(snakemake.input.units_file,'r') as infile:
    units_file = json.load(infile)

with open(snakemake.input.monomer_displacement_per_timestep,'rb') as infile:
    displacement = pickle.load(infile)   
    ends_displacement = displacement['ends_displacement']
    internal_displacement = displacement['internal_displacement']

R = snakemake.config['DSB_simulation_parameters']['contact_radii_for_first_passage_times']

dt_min = units_file['seconds_per_sim_dt']['internal TAD']/60
dx_nm =  units_file['nanometers_per_monomer']['internal TAD']

plt.figure(figsize=(10,10))

# now, create histograms    
plt.subplot(2,2,1)
plt.hist(np.asarray([i for i in ends_displacement]).flatten(),bins='auto',
         cumulative=True,density=True,histtype='step',label=f'mean = '+str(np.mean(ends_displacement))) 
# plt.xscale('log')
plt.xlabel('displacement/simulation_step at polymer ends (sim units)')
plt.ylabel('Cumulative probability')
plt.ylim([0,1])
plt.xlim([0,np.max(ends_displacement)*0.99])
plt.legend(loc='lower right')

if len(internal_displacement)>0:
    plt.subplot(2,2,2)
    plt.hist(np.asarray([i*dx_nm for i in ends_displacement]).flatten(),bins='auto',
             cumulative=True,density=True,histtype='step',label=f'mean = '+str(dx_nm*np.mean(ends_displacement))+' nm') 
    # plt.xscale('log')
    plt.xlabel('displacement/simulation_step at polymer ends (nm)')
    plt.ylabel('Cumulative probability')
    plt.ylim([0,1])
    plt.xlim([0,(np.max(internal_displacement)*0.99)*dx_nm])
    plt.legend(loc='lower right')

    plt.subplot(2,2,3)
    plt.hist(np.asarray([i for i in internal_displacement]).flatten(),bins='auto',
             cumulative=True,density=True,histtype='step', label=f'mean = '+str(np.mean(internal_displacement)))  
    # plt.xscale('log')
    plt.xlabel('displacement/simulation_step at internal TAD (sim units)')
    plt.ylim([0,1])
    plt.xlim([0,np.max(internal_displacement)*0.99])
    plt.legend(loc='lower right')


    plt.hist(np.asarray([i*dx_nm for i in internal_displacement]).flatten(),bins='auto',
             cumulative=True,density=True,histtype='step',label=f'mean = '+str(dx_nm*np.mean(internal_displacement))+' nm') 
    # plt.xscale('log')
    plt.xlabel('displacement/simulation_step at internal TAD (nm)')
    plt.ylim([0,1])
    plt.xlim([0,(np.max(internal_displacement)*0.99)*dx_nm])
    plt.legend(loc='lower right')

plt.savefig(snakemake.output.displacement_per_timestep_distribution,bbox_inches='tight')

