import json
import pickle
import matplotlib.pyplot as plt
import numpy as np

with open(snakemake.input.units_file,'r') as infile:
    units_file = json.load(infile)

with open(snakemake.input.first_passage_times_v2,'rb') as infile:
    passage_times =  pickle.load(infile)   
    no_contact_intervals = passage_times['no_contact_intervals']
    with_contact_intervals = passage_times['with_contact_intervals']

color_list =['dodgerblue', 'orangered','limegreen','magenta']
    
R = snakemake.config['DSB_simulation_parameters']['contact_radii_for_first_passage_times']
loop_size = ['all'] # loop size in which DSB occurs, in the unit of # monomers

dt_min = units_file['seconds_per_sim_dt']['internal TAD']/60
dx_nm =  units_file['nanometers_per_monomer']['internal TAD']

first_protein_arrival_time = 1/dt_min # use ku70/80 arrival time to estimate when the two DSB ends become sticky (in simulation time unit)

bin_edge_0 = np.arange(1,10)
bin_edge= np.hstack((bin_edge_0,bin_edge_0*10,bin_edge_0*100,bin_edge_0*1000,bin_edge_0*10000,bin_edge_0*100000))
bin_edge=bin_edge.tolist()
try:
    x = no_contact_intervals[R[0]][loop_size[0]][0]
except KeyError:
    plt.figure(figsize=(8,5*len(R)))
    for ri, r in enumerate(R):                                                      
        # first need to consolidate all the times
        intervals = []
        for idx in no_contact_intervals[r].keys():
            intervals.extend(no_contact_intervals[r][idx])   

        # now, create histograms    
        plt.subplot(len(R)+1,2,ri*2+1)
        min_len = 0
        plt.hist([i for i in intervals if i>min_len],bins=bin_edge,
                 cumulative=True,density=True,histtype='step',label=f'Radius {r} (monomers)') 
        plt.xscale('log')
        plt.xlabel('Time to recapture (sim units)')
        plt.ylabel('Cumulative probability')
        plt.ylim([0,1])
        if len(intervals)>0:
            plt.xlim([first_protein_arrival_time,np.max(intervals)+1])
        plt.legend(loc='lower right')


        plt.subplot(len(R)+1,2,ri*2+2)
        plt.hist([i*dt_min for i in intervals if i>min_len],bins=bin_edge,
                 cumulative=True,density=True,histtype='step', label=f'Radius {int(r*dx_nm)} nm') 
        plt.xscale('log')
        plt.xlabel('Time to recapture (minutes)')
        plt.ylim([0,1])
        if len(intervals)>0:
            plt.xlim([1,(np.max(intervals)+1)*dt_min])
        plt.legend(loc='lower right')

    plt.savefig(snakemake.output.first_recapture_v2,bbox_inches='tight')


    plt.figure(figsize=(8,5*len(R)))
    for ri, r in enumerate(R):                                                      
        # first need to consolidate all the times
        intervals = []
        for idx in with_contact_intervals[r].keys():
            intervals.extend(with_contact_intervals[r][idx])   

        # now, create histograms    
        plt.subplot(len(R)+1,2,ri*2+1)
        min_len = 0
        plt.hist([i for i in intervals if i>min_len],bins=bin_edge,
                 cumulative=True,density=True,histtype='step',label=f'Radius {r} (monomers)') 
        plt.xscale('log')
        plt.xlabel('Time to exit (sim units)')
        plt.ylabel('Cumulative probability')
        plt.ylim([0,1])
        if len(intervals)>0:
            plt.xlim([first_protein_arrival_time,np.max(intervals)+1])
        plt.legend(loc='lower right')


        plt.subplot(len(R)+1,2,ri*2+2)
        plt.hist([i*dt_min for i in intervals if i>min_len],bins=bin_edge,
                 cumulative=True,density=True,histtype='step', label=f'Radius {int(r*dx_nm)} nm') 
        plt.xscale('log')
        plt.xlabel('Time to exit (minutes)')
        plt.ylim([0,1])
        if len(intervals)>0:
            plt.xlim([1,(np.max(intervals)+1)*dt_min])
        plt.legend(loc='lower right')
        
    plt.savefig(snakemake.output.first_exit_v2,bbox_inches='tight')
    
else:
    plt.figure(figsize=(8*len(loop_size),6.8*len(R)))
    for ri, r in enumerate(R): 
        for li,l in enumerate(loop_size): 
            # first need to consolidate all the times
            intervals = []
            for idx in no_contact_intervals[r][l].keys():
                intervals.extend(no_contact_intervals[r][l][idx])   

            # now, create histograms    
            plt.subplot(len(R)+1,2*len(loop_size),li*2+ri*2*len(loop_size)+1)
            min_len = 0
            plt.hist([i for i in intervals if i>min_len],bins=bin_edge,
                     cumulative=True,density=True,histtype='step',label=f'Radius {r} (monomers)',color=color_list[li]) 
            plt.xscale('log')
            plt.xlabel('Time to recapture (sim units)')
            plt.ylabel('Cumulative probability')
            plt.ylim([0,1])
            plt.xlim([first_protein_arrival_time,np.max(intervals)+1])
            plt.legend(loc='lower right')
            plt.title('loop size: '+ str(l)+ ' monomers')


            plt.subplot(len(R)+1,2*len(loop_size),li*2+ri*2*len(loop_size)+2)
            plt.hist([i*dt_min for i in intervals if i>min_len],bins=bin_edge,
                     cumulative=True,density=True,histtype='step', label=f'Radius {int(r*dx_nm)} nm',color=color_list[li]) 
            plt.xscale('log')
            plt.xlabel('Time to recapture (minutes)')
            plt.ylim([0,1])
            plt.xlim([1,(np.max(intervals)+1)*dt_min])
            plt.legend(loc='lower right')
            plt.title('loop size: '+ str(l)+ ' monomers')

    plt.savefig(snakemake.output.first_recapture_v2,bbox_inches='tight')


    plt.figure(figsize=(8*len(loop_size),6.8*len(R)))
    for ri, r in enumerate(R): 
        for li,l in enumerate(loop_size): 
            # first need to consolidate all the times
            intervals = []
            for idx in with_contact_intervals[r][l].keys():
                intervals.extend(with_contact_intervals[r][l][idx])   

            # now, create histograms    
            plt.subplot(len(R)+1,2*len(loop_size),li*2+ri*2*len(loop_size)+1)
            min_len = 0
            plt.hist([i for i in intervals if i>min_len],bins=bin_edge,
                     cumulative=True,density=True,histtype='step',label=f'Radius {r} (monomers)',color=color_list[li]) 
            plt.xscale('log')
            plt.xlabel('Time to exit (sim units)')
            plt.ylabel('Cumulative probability')
            plt.ylim([0,1])
            plt.xlim([first_protein_arrival_time,np.max(intervals)+1])
            plt.legend(loc='lower right')
            plt.title('loop size: '+ str(l)+ ' monomers')


            plt.subplot(len(R)+1,2*len(loop_size),li*2+ri*2*len(loop_size)+2)
            plt.hist([i*dt_min for i in intervals if i>min_len],bins=bin_edge,
                     cumulative=True,density=True,histtype='step', label=f'Radius {int(r*dx_nm)} nm',color=color_list[li]) 
            plt.xscale('log')
            plt.xlabel('Time to exit (minutes)')
            plt.ylim([0,1])
            plt.xlim([1,(np.max(intervals)+1)*dt_min])
            plt.legend(loc='lower right')
            plt.title('loop size: '+ str(l)+ ' monomers')

    
    plt.savefig(snakemake.output.first_exit_v2,bbox_inches='tight')