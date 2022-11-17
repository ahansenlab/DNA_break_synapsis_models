# This version plots the simulation unit figures in the first column and the SI unit figures in the second column
import json
import pickle
import matplotlib.pyplot as plt
import numpy as np
import bootstrapped.bootstrap as bs
import bootstrapped.stats_functions as bs_stats
import tracklib as tl
from tracklib.util.stats import KM_survival
import tracklib as tl
from tracklib.util.stats import MLE_censored_exponential
num_chrom = []
inner_chrom_num_recapture =[]
between_chrom_num_recapture =[]
fraction_inner_chrom_num_recapture =[]
fraction_inner_chrom_with_recapture=[]
fraction_between_chrom_num_recapture =[]
inner_chrom_num_exit =[]
between_chrom_num_exit =[]
fraction_inner_chrom_num_exit =[]
fraction_between_chrom_num_exit =[]


# Get difference of two lists
def Diff(li1, li2):
    return list(set(li1) - set(li2)) + list(set(li2) - set(li1))

# recapture_times_1D = np.asarray([ 1.03333333,  1.41666667,  1.65      ,  1.65      ,  1.65      ,
#         1.8       ,  1.85      ,  1.88333333,  2.06666667,  2.48333333,
#         2.63333333,  2.66666667,  3.38333333,  3.46666667,  3.53333333,
#         3.63333333,  3.68333333,  3.7       ,  3.75      ,  4.01666667,
#         4.11666667,  4.16666667,  4.2       ,  4.28333333,  4.35      ,
#         4.41666667,  4.63333333,  4.8       ,  5.05      ,  5.31666667,
#         5.33333333,  5.33333333,  5.43333333,  5.5       ,  5.56666667,
#         5.58333333,  5.6       ,  5.66666667,  5.7       ,  5.8       ,
#         5.98333333,  6.05      ,  6.33333333,  6.93333333,  7.05      ,
#         7.13333333,  7.26666667,  7.28333333,  7.41666667,  7.98333333,
#         8.03333333,  8.18333333,  8.3       ,  8.38333333,  8.43333333,
#         8.63333333,  9.36666667,  9.4       ,  9.51666667,  9.61666667,
#         9.65      ,  9.68333333, 10.16666667, 10.2       , 10.23333333,
#        10.28333333, 10.31666667, 10.4       , 10.41666667, 10.63333333,
#        10.76666667, 10.8       , 10.98333333, 11.13333333, 11.26666667,
#        11.43333333, 12.11666667, 12.15      , 12.46666667, 12.76666667,
#        12.86666667, 13.2       , 13.51666667, 13.6       , 13.71666667,
#        14.13333333, 15.46666667, 15.88333333, 17.33333333, 17.58333333,
#        18.21666667, 19.06666667, 23.03333333, 23.3       , 23.31666667,
#        25.53333333, 27.58333333, 29.03333333, 29.41666667, 45.9       ,
#        47.18333333,  1.16666667,  1.48333333,  1.71666667,  1.8       ,
#         1.88333333,  1.88333333,  1.93333333,  2.01666667,  2.11666667,
#         2.18333333,  2.21666667,  2.3       ,  2.33333333,  2.85      ,
#         2.86666667,  2.95      ,  3.7       ,  3.73333333,  3.81666667,
#         4.05      ,  4.1       ,  4.11666667,  4.15      ,  4.23333333,
#         4.28333333,  4.28333333,  4.4       ,  4.41666667,  4.5       ,
#         4.53333333,  4.75      ,  4.78333333,  4.9       ,  4.95      ,
#         5.        ,  5.05      ,  5.08333333,  5.11666667,  5.2       ,
#         5.2       ,  5.48333333,  5.78333333,  5.93333333,  5.95      ,
#         5.98333333,  6.15      ,  6.15      ,  6.31666667,  6.33333333,
#         6.38333333,  6.51666667,  6.63333333,  6.85      ,  6.93333333,
#         7.13333333,  7.25      ,  7.33333333,  7.55      ,  7.6       ,
#         7.6       ,  7.96666667,  8.11666667,  8.43333333,  8.61666667,
#         8.93333333,  9.51666667,  9.53333333,  9.65      ,  9.95      ,
#        10.23333333, 10.8       , 10.91666667, 10.96666667, 11.01666667,
#        11.1       , 11.43333333, 11.46666667, 11.8       , 11.81666667,
#        11.86666667, 12.03333333, 12.1       , 12.23333333, 13.43333333,
#        14.03333333, 14.28333333, 14.45      , 15.2       , 16.68333333,
#        17.81666667, 18.06666667, 19.2       , 19.61666667, 31.71666667,
#        31.73333333, 35.05      ,  1.15      ,  1.33333333,  1.43333333,
#         1.61666667,  1.78333333,  1.95      ,  1.98333333,  2.16666667,
#         2.26666667,  2.51666667,  2.65      ,  2.7       ,  2.73333333,
#         3.1       ,  3.2       ,  3.25      ,  3.28333333,  3.38333333,
#         3.53333333,  3.56666667,  3.6       ,  3.63333333,  3.83333333,
#         3.88333333,  3.98333333,  4.        ,  4.1       ,  4.11666667,
#         4.16666667,  4.2       ,  4.28333333,  4.31666667,  4.35      ,
#         4.51666667,  4.63333333,  4.68333333,  4.73333333,  5.21666667,
#         5.25      ,  5.61666667,  6.03333333,  6.31666667,  6.38333333,
#         6.51666667,  6.61666667,  6.65      ,  6.73333333,  7.        ,
#         7.33333333,  7.5       ,  7.61666667,  7.66666667,  7.7       ,
#         7.73333333,  7.75      ,  7.8       ,  8.01666667,  8.06666667,
#         8.15      ,  8.16666667,  8.28333333,  8.31666667,  8.46666667,
#         8.48333333,  8.51666667,  8.63333333,  8.8       ,  9.75      ,
#         9.88333333,  9.91666667, 10.45      , 11.01666667, 11.36666667,
#        11.36666667, 11.78333333, 12.23333333, 12.61666667, 12.66666667,
#        12.98333333, 12.98333333, 13.35      , 13.65      , 13.83333333,
#        13.85      , 14.48333333, 15.51666667, 15.86666667, 16.71666667,
#        18.76666667, 19.65      , 20.35      , 20.46666667, 21.38333333,
#        21.5       , 24.48333333, 32.8    ])
recapture_times_1D = np.asarray([])
# success_rate = 0.4677
success_rate = 0

recapture_times_3D =[[],[],[]] 
uncapture_times_3D =[[],[],[]] 
inner_chrom_intervals=[[],[],[]] 
inner_recapture_count_global =np.zeros(3)
DSB_num_global =np.zeros(3)
constrained_global=np.zeros(3)
unconstrained_global=np.zeros(3)
for i, (units_file_name, first_passage_times_file, monomer_displacement_file) in enumerate(zip(snakemake.input.units_file,snakemake.input.first_passage_times_v2,snakemake.input.monomer_displacement_per_timestep)):
    inner_chrom_num_recapture.append([])
    inner_chrom_num_exit.append([])
    fraction_inner_chrom_num_recapture.append([])
    fraction_inner_chrom_with_recapture.append([])
    fraction_inner_chrom_num_exit.append([])
    if "multichain" in first_passage_times_file:
        between_chrom_num_exit.append([])
        between_chrom_num_recapture.append([])
        fraction_between_chrom_num_recapture.append([])
        fraction_between_chrom_num_exit.append([])

    

    with open(units_file_name,'r') as infile:
        units_file = json.load(infile)
    dt_min = units_file['seconds_per_sim_dt']['internal TAD']/60
    dx_nm =  units_file['nanometers_per_monomer']['internal TAD']


    with open( first_passage_times_file,'rb') as infile:
        passage_times =  pickle.load(infile)   
        no_contact_intervals = passage_times['no_contact_intervals']
        with_contact_intervals = passage_times['with_contact_intervals']
        inner_encountered_separated_count = passage_times['inner_encountered_separated_count']
        between_encountered_separated_count = passage_times['between_encountered_separated_count']
        inner_unrecaptured_time = passage_times['inner_unrecaptured_time']
        between_unrecaptured_time = passage_times['between_unrecaptured_time']
        constrained = passage_times['constrained']
        unconstrained = passage_times['unconstrained']
    cmap = plt.get_cmap("tab10")
    
    if "multichain" in first_passage_times_file:
        current_num_chrom = units_file_name.split('_C')[1].split('_collisionrate')[0]
        num_chrom.append("num_chrom = "+ current_num_chrom)

        current_num_chrom = int(current_num_chrom)
        all_index = np.arange((current_num_chrom)*(current_num_chrom*2-1)).tolist()
        inner_chrom_index = [current_num_chrom*i+(i+current_num_chrom)-((i+2)*(i+1))//2 for i in range(current_num_chrom)]
        between_chrom_index = Diff(inner_chrom_index,all_index)
    elif "brokenloop" in first_passage_times_file:
        current_num_loops = units_file_name.split('_DSB')[1].split('_LoopSize')[0]
    elif "DSB_sim" in first_passage_times_file:
        current_num_DSBs = units_file_name.split('_DSB')[1].split('_collisionrate')[0]

    with open( monomer_displacement_file,'rb') as infile:
        displacement = pickle.load(infile)   
        ends_displacement = displacement['ends_displacement']
    
    R = snakemake.config['DSB_simulation_parameters']['contact_radii_for_first_passage_times']
    R += [1+np.around(np.mean(ends_displacement),decimals=2)] #use the one monomer size + step size the last contact radius to plot
    loop_size = ['all'] # loop size in which DSB occurs, in the unit of # monomers
    
    if "multichain" in first_passage_times_file:
        N_monomers = int(units_file_name.split('_N')[1].split('_C')[0])
    elif "brokenloop" in first_passage_times_file:
        N_monomers = int(units_file_name.split('_N')[1].split('_Repeat')[0])
    elif "DSB_sim" in first_passage_times_file:
        N_monomers = int(units_file_name.split('_N')[1].split('_DSB')[0])
    volume_density = float(units_file_name.split('_dens')[1].split('_N')[0])
    print(dx_nm)
    radius = 0.5*(N_monomers / volume_density) ** (1/3) * dx_nm # in nm

    first_protein_arrival_time = 1/dt_min # use ku70/80 arrival time to estimate when the two DSB ends become sticky (in minutes)

    if i == 0:
        fig1 = plt.figure(figsize=(12,5*len(R)))
        ax_set1 =  fig1.subplots(len(R), 2)
        fig2 = plt.figure(figsize=(12,5*len(R)))
        ax_set2 =  fig2.subplots(len(R), 2) 
        fig3 = plt.figure(figsize=(12,5*len(R)))
        ax_set3 =  fig3.subplots(len(R), 2)
        fig4 = plt.figure(figsize=(12,5*len(R)))
        ax_set4 =  fig4.subplots(len(R), 2)

    for ri, r in enumerate(no_contact_intervals.keys()):                                                      
        # first need to consolidate all the times

        between_chrom_intervals = []
        inner_index_counter = []
        between_index_counter = []
        inner_recapture_count = 0
        if "multichain" in first_passage_times_file:
            for idx in no_contact_intervals[r].keys():

                inner_chrom_intervals.extend(no_contact_intervals[r][idx]) 
                recapture_times_3D[ri].extend( no_contact_intervals[r][idx])
                inner_recapture_count+=len(no_contact_intervals[r][idx])
                if len(no_contact_intervals[r][idx])>0:
                    inner_index_counter.append(idx)

        if "brokenloop" in first_passage_times_file:
            for idx in no_contact_intervals[r].keys():
                inner_chrom_intervals.extend(no_contact_intervals[r][idx]) 
                inner_recapture_count+=len(no_contact_intervals[r][idx])
                if len(no_contact_intervals[r][idx])>0:
                    inner_index_counter.append(idx)
                    
        if "DSB_sim" in first_passage_times_file:
            for idx in no_contact_intervals[r].keys():
                inner_chrom_intervals[ri].extend(no_contact_intervals[r][idx])
                recapture_times_3D[ri].extend( no_contact_intervals[r][idx])
                inner_recapture_count+=len(no_contact_intervals[r][idx])
                if len(no_contact_intervals[r][idx])>0:
                    inner_index_counter.append(idx)


        inner_chrom_num_recapture[-1].append(len(inner_chrom_intervals[ri]))
        if inner_encountered_separated_count[r]>0:
            fraction_inner_chrom_num_recapture[-1].append(np.around(inner_recapture_count/inner_encountered_separated_count[r]*100))
        elif "DSB_sim" in units_file_name:
            fraction_inner_chrom_num_recapture[-1].append(np.around(inner_recapture_count/int(current_num_DSBs)*100))
            inner_recapture_count_global[ri]+=inner_recapture_count
            DSB_num_global[ri] += int(current_num_DSBs)
            constrained_global[ri] += constrained[r]
            unconstrained_global[ri] += unconstrained[r]
        else:
            fraction_inner_chrom_num_recapture[-1].append(-1)
        if "brokenloop" in units_file_name:
            fraction_inner_chrom_with_recapture[-1].append(int(len(set(inner_index_counter))/int(current_num_loops)*100))
        elif "multichain" in units_file_name:
            fraction_inner_chrom_with_recapture[-1].append(int(len(set(inner_index_counter))/int(current_num_chrom)*100))
            inner_recapture_count_global[ri]+=inner_recapture_count
        elif "DSB_sim" in units_file_name:
            if int(current_num_DSBs)>0:
                fraction_inner_chrom_with_recapture[-1].append(int(len(set(inner_index_counter))/int(current_num_DSBs)*100))
            else:
                fraction_inner_chrom_with_recapture[-1].append(-1)
            
for ri, r in enumerate(with_contact_intervals.keys()): 
    min_len = 0
    print(len(recapture_times_3D[ri]))
    list_to_plot = np.asarray([j for j in recapture_times_3D[ri] if j>min_len])
    print(list_to_plot)
    bins=np.arange(np.amin(list_to_plot),np.amax(list_to_plot)+1)
    if "DSB_sim" in first_passage_times_file:
        ax_set1[ri,0].hist(list_to_plot,
                 cumulative=True,density=False,weights=np.ones_like(list_to_plot)/len(list_to_plot)*inner_recapture_count_global[ri]/DSB_num_global[ri],histtype='step', label=f'Radius {int(r*dx_nm)} nm',bins=bins) 
    else:
        ax_set1[ri,0].hist(list_to_plot,
                 cumulative=True,density=False,weights=np.ones_like(list_to_plot)/len(list_to_plot)*fraction_inner_chrom_num_recapture[i][ri]/100,histtype='step', label=f'Radius {int(r*dx_nm)} nm',bins=bins) 

    ax_set1[ri,0].set_xscale('log')
    ax_set1[ri,0].set_xlabel('Time to recapture (sim time step)')
    ax_set1[ri,0].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')
    ax_set1[ri,0].set_ylim([0,1.1])

    all_intervals = np.asarray(recapture_times_3D[ri]+inner_unrecaptured_time[r])

    censored = np.concatenate((np.zeros(len(recapture_times_3D[ri])),np.ones(len(inner_unrecaptured_time[r]))),axis=None)
    
    KM_survival_curve = KM_survival(all_intervals,censored,conf=0.95, Tmax=np.inf, S1at=1/dt_min)#calculating KM survival curve
    if 0.5 in KM_survival_curve[:,1]:
        median_idx = np.argwhere(KM_survival_curve[:,1]==0.5)[0][0]
        median = KM_survival_curve[median_idx,0]
    else:
        a= np.argwhere(KM_survival_curve[:,1]-0.5<0)
        print(KM_survival_curve[:,1])
        print(censored)
        median_idx = np.argwhere(KM_survival_curve[:,1]-0.5<0)[0][0]
        median = (KM_survival_curve[median_idx,0]+KM_survival_curve[median_idx-1,0])/2

    if 0.5 in KM_survival_curve[:,2]:
        median_lbnd_idx = np.argwhere(KM_survival_curve[:,2]==0.5)[0][0]
        median_lbnd = KM_survival_curve[median_lbnd_idx,0]
    else:
        median_lbnd_idx = np.argwhere(KM_survival_curve[:,2]-0.5<0)[0][0]
        median_lbnd = (KM_survival_curve[median_lbnd_idx,0]+KM_survival_curve[median_lbnd_idx-1,0])/2

    if 0.5 in KM_survival_curve[:,3]:
        median_ubnd_idx = np.argwhere(KM_survival_curve[:,3]==0.5)[0][0]
        median_ubnd = KM_survival_curve[median_ubnd_idx,0]
    else:
        median_ubnd_idx = np.argwhere(KM_survival_curve[:,3]-0.5<0)[0][0]
        median_ubnd = (KM_survival_curve[median_ubnd_idx,0]+KM_survival_curve[median_ubnd_idx-1,0])/2

    mean,mean_lbnd,mean_ubnd = MLE_censored_exponential(all_intervals,censored,conf=0.95)
#     mean, mean_lbnd, mean_ubnd = (0,0,0)


    ax_set3[ri,0].errorbar(radius,mean,yerr=[np.asarray([np.abs(mean_lbnd-mean)]),np.asarray([np.abs(mean_ubnd-mean)])],fmt='o', capthick=2, color = cmap(i))
    ax_set3[ri,0].set_xlabel('Confinement radius (nm)')
    ax_set3[ri,0].set_ylabel('Mean time to recapture (sim time step)')
    ax_set3[ri,0].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')

    
    list_to_plot = [j*dt_min for j in recapture_times_3D[ri] if j>min_len]
    bins=np.arange(np.amin(list_to_plot),np.amax(list_to_plot)+1,step=0.01)
    if "DSB_sim" in first_passage_times_file:
        ax_set1[ri,1].hist(list_to_plot,
                 cumulative=True,density=False,weights=np.ones_like(list_to_plot)/len(list_to_plot)*inner_recapture_count_global[ri]/DSB_num_global[ri],histtype='step', label=f'Radius {int(r*dx_nm)} nm',bins=bins) 
    else:
        ax_set1[ri,1].hist(list_to_plot,
                 cumulative=True,density=False,weights=np.ones_like(list_to_plot)/len(list_to_plot)*fraction_inner_chrom_num_recapture[i][ri]/100,histtype='step', label=f'Radius {int(r*dx_nm)} nm',bins=bins)
    ax_set1[ri,1].hist(recapture_times_1D,
             cumulative=True,density=False,weights=np.ones_like(recapture_times_1D)/len(recapture_times_1D)*success_rate,histtype='step',color = 'grey') 
    ax_set1[ri,1].set_xscale('log')
    ax_set1[ri,1].set_xlabel('Time to recapture (minutes)')
    ax_set1[ri,1].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')
    ax_set1[ri,1].set_ylim([0,1.1])

    
    ax_set3[ri,1].errorbar(radius,mean*dt_min,yerr=[np.asarray([np.abs(mean_lbnd-mean)*dt_min]),np.asarray([np.abs(mean_ubnd-mean)*dt_min])],fmt='o', capthick=2, color = cmap(i))
    ax_set3[ri,1].set_xlabel('Confinement radius (nm)')
    ax_set3[ri,1].set_ylabel('Mean time to recapture (minutes)')
    ax_set3[ri,1].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')

    if i == len(snakemake.input.units_file)-1:

        for ri, r in enumerate(no_contact_intervals.keys()):  
            if "DSB_sim" in first_passage_times_file:
                legend_text = ['n = '+str(inner_recapture_count_global[ri]) +', %synapsis = '+str(np.around(inner_recapture_count_global[ri]/DSB_num_global[ri]*100,2))+', %constrained|synapsed = '+str(np.around(constrained_global[ri]/(constrained_global[ri]+unconstrained_global[ri])*100,2))+', # of DSB = '+str(DSB_num_global[ri]) for j,a in enumerate(snakemake.input.units_file)]
            else:
                legend_text = ['n = '+str(inner_recapture_count_global[ri]) +', %synapsis = '+str(np.around(fraction_inner_chrom_num_recapture[i][ri],2))+', # of DSB = '+str(inner_encountered_separated_count[r]) for j,a in enumerate(snakemake.input.units_file)]
            ax_set1[ri,0].legend(legend_text,loc='lower right',framealpha=0.2)
            ax_set3[ri,0].legend(legend_text,loc='lower right',framealpha=0.2)

    min_len = 0
    list_to_plot=[j for j in inner_chrom_intervals[ri] if j>min_len]
    ax_set2[ri,0].hist(list_to_plot,
             cumulative=True,density=False,weights=np.ones_like(list_to_plot)/len(list_to_plot)*fraction_inner_chrom_num_recapture[i][ri]/100,histtype='step', label=f'Radius {int(r*dx_nm)} nm') 
    ax_set2[ri,0].set_xscale('log')
    ax_set2[ri,0].set_xlabel('Time to exit (sim time step)')
    ax_set2[ri,0].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')
    ax_set2[ri,0].set_ylim([0,1.1])

    ax_set4[ri,0].scatter(radius,np.mean(inner_chrom_intervals[ri]),
             s=25) 
    ax_set4[ri,0].set_xlabel('Confinement radius (nm)')
    ax_set4[ri,0].set_ylabel('Mean time to exit (sim time step)')
    ax_set4[ri,0].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')
    
    list_to_plot = [j*dt_min for j in inner_chrom_intervals[ri] if j>min_len]
    ax_set2[ri,1].hist(list_to_plot,
             cumulative=True,density=False,weights=np.ones_like(list_to_plot)/len(list_to_plot)*fraction_inner_chrom_num_recapture[i][ri]/100,histtype='step', label=f'Radius {int(r*dx_nm)} nm') 
    ax_set2[ri,1].set_xscale('log')
    ax_set2[ri,1].set_xlabel('Time to exit (minutes)')
    ax_set2[ri,1].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')
    ax_set2[ri,1].set_ylim([0,1.1])

    ax_set4[ri,1].scatter(radius,np.mean(inner_chrom_intervals[ri])*dt_min ,
             s=25) 
    ax_set4[ri,1].set_xlabel('Confinement radius (nm)')
    ax_set4[ri,1].set_ylabel('Mean time to exit (minutes)')
    ax_set4[ri,1].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')

    if i == len(snakemake.input.units_file)-1:

        for ri, r in enumerate(no_contact_intervals.keys()):  
            legend_text = ['n = '+str(inner_chrom_num_recapture[j][ri]) +', %A = '+str(fraction_inner_chrom_num_recapture[j][ri])+', %B = '+str(fraction_inner_chrom_with_recapture[j][ri]) for j,a in enumerate(snakemake.input.units_file)]
            ax_set2[ri,0].legend(legend_text,loc='lower right',framealpha=0.2)
            ax_set4[ri,0].legend(legend_text,loc='lower right',framealpha=0.2)
            ax_set2[ri,1].legend(legend_text,loc='lower right',framealpha=0.2)
            ax_set4[ri,1].legend(legend_text,loc='lower right',framealpha=0.2)

fig1.savefig(snakemake.output.first_recapture_combined_v2,bbox_inches='tight')
fig2.savefig(snakemake.output.first_exit_combined_v2,bbox_inches='tight')
fig3.savefig(snakemake.output.first_recapture_vs_confinementR_v2,bbox_inches='tight')
fig4.savefig(snakemake.output.first_exit_vs_confinementR_v2,bbox_inches='tight')