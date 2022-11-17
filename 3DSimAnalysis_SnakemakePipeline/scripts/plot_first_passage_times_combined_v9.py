# This version plots the simulation unit figures in the first column and the SI unit figures in the second column
import json
import pickle
import matplotlib.pyplot as plt
import numpy as np
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

# DKW inequality, a non-parametric upper and lower confidence interval
def edf(data, alpha=.05, scaling_factor=1,x0=None, x1=None ):
    x0 = data.min() if x0 is None else x0
    x1 = data.max() if x1 is None else x1
    x = np.linspace(x0, x1, 1000)
    N = data.size
    y = np.zeros_like(x)
    l = np.zeros_like(x)
    u = np.zeros_like(x)
    e = np.sqrt(1.0/(2*N) * np.log(2./alpha))
    for i, xx in enumerate(x):
        y[i] = np.sum(data <= xx)/N*scaling_factor
        l[i] = np.maximum( y[i] - e, 0 )
        u[i] = np.minimum( y[i] + e, 1 )
    return x, y, l, u

# Get difference of two lists
def Diff(li1, li2):
    return list(set(li1) - set(li2)) + list(set(li2) - set(li1))

file = snakemake.input.recapture_times_1D
with open(file, 'rb') as infile:
    recapture_times_1D = np.load(infile)
    success_rate = np.load(infile)
# convert to minutes
recapture_times_1D = recapture_times_1D/30


recapture_times_3D =[[],[],[]] 
exit_times_3D =[[],[],[]] 
recapture_times_3D_constrained =[[],[],[]] 
recapture_times_3D_unconstrained =[[],[],[]] 
exit_times_3D_constrained =[[],[],[]] 
exit_times_3D_unconstrained =[[],[],[]] 
uncapture_times_3D =[[],[],[]] 
uncapture_times_3D_constrained =[[],[],[]] 
uncapture_times_3D_unconstrained =[[],[],[]] 
inner_chrom_intervals=[[],[],[]] 
inner_recapture_count_global =np.zeros(3)
inner_recapture_count_global_constrained =np.zeros(3)
inner_recapture_count_global_unconstrained =np.zeros(3)
DSB_num_global =np.zeros(3)
DSB_num_global_constrained =np.zeros(3)
DSB_num_global_unconstrained =np.zeros(3)
constrained_global=np.zeros(3)
unconstrained_global=np.zeros(3)
for i, (units_file_name, first_passage_times_file) in enumerate(zip(snakemake.input.units_file,snakemake.input.first_passage_times_v2)):
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
    first_protein_arrival_time = 1/dt_min # use ku70/80 arrival time,~1 minute to estimate when the two DSB ends become sticky (in simulation time unit)
    dx_nm =  units_file['nanometers_per_monomer']['internal TAD']


    with open( first_passage_times_file,'rb') as infile:
        passage_times =  pickle.load(infile)   
        initial_state = passage_times['initial_state']
        print(initial_state)
        no_contact_intervals = passage_times['no_contact_intervals']
        with_contact_intervals = passage_times['with_contact_intervals']
        inner_encountered_separated_count = passage_times['inner_encountered_separated_count']
        between_encountered_separated_count = passage_times['between_encountered_separated_count']
        inner_unrecaptured_time = passage_times['inner_unrecaptured_time']
        between_unrecaptured_time = passage_times['between_unrecaptured_time']
        print(inner_unrecaptured_time)
        constrained = passage_times['constrained']
        unconstrained = passage_times['unconstrained']
        constrained_time0 = passage_times['constrained_time0']
        if "frozenloop" in first_passage_times_file:
            inner_encountered_separated_count_constrained = passage_times['inner_encountered_separated_count_constrained']
            inner_encountered_separated_count_unconstrained = passage_times['inner_encountered_separated_count_unconstrained']
        
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
    elif "DSB_sim" in first_passage_times_file or "frozenloop_stochastic" in first_passage_times_file:
        current_num_DSBs = units_file_name.split('_DSB')[1].split('_collisionrate')[0]

    R = snakemake.config['DSB_simulation_parameters']['contact_radii_for_first_passage_times']
    loop_size = ['all'] # loop size in which DSB occurs, in the unit of # monomers
    
    if "multichain" in first_passage_times_file:
        N_monomers = int(units_file_name.split('_N')[1].split('_C')[0])
    elif "brokenloop" in first_passage_times_file:
        N_monomers = int(units_file_name.split('_N')[1].split('_Repeat')[0])
    elif "DSB_sim" in first_passage_times_file or "frozenloop_stochastic" in first_passage_times_file:
        N_monomers = int(units_file_name.split('_N')[1].split('_DSB')[0])
    elif "frozenloop_sim" in first_passage_times_file:
        N_monomers = int(units_file_name.split('_N')[1].split('_bondsFrom')[0])
        
    volume_density = float(units_file_name.split('_dens')[1].split('_N')[0])

    radius = 0.5*(N_monomers / volume_density) ** (1/3) * dx_nm # in nm

    if i == 0:
        fig1 = plt.figure(figsize=(12,5*len(R)))
        ax_set1 =  fig1.subplots(len(R), 2)
        fig2 = plt.figure(figsize=(12,5*len(R)))
        ax_set2 =  fig2.subplots(len(R), 2) 
        fig3 = plt.figure(figsize=(12,5*len(R)))
        ax_set3 =  fig3.subplots(len(R), 2)
        fig4 = plt.figure(figsize=(12,5*len(R)))
        ax_set4 =  fig4.subplots(len(R), 2)
        fig5 = plt.figure(figsize=(12,5*len(R)))
        ax_set5 =  fig5.subplots(len(R), 2)
        fig6 = plt.figure(figsize=(12,5*len(R)))
        ax_set6 =  fig6.subplots(len(R), 2)

    for ri, r in enumerate(no_contact_intervals.keys()):                                                      
        # first need to consolidate all the times

        between_chrom_intervals = []
        inner_index_counter = []
        between_index_counter = []
        inner_recapture_count = 0
        inner_recapture_count_constrained = 0
        inner_recapture_count_unconstrained = 0
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
                recapture_times_3D[ri].extend(no_contact_intervals[r][idx])
                inner_recapture_count+=len(no_contact_intervals[r][idx])
                if len(no_contact_intervals[r][idx])>0:
                    inner_index_counter.append(idx)
                    
        if "DSB_sim" in first_passage_times_file or "frozenloop_stochastic" in first_passage_times_file:
            for idx in no_contact_intervals[r].keys():
                inner_chrom_intervals[ri].extend(no_contact_intervals[r][idx])
                recapture_times_3D[ri].extend( no_contact_intervals[r][idx])
                uncapture_times_3D[ri].extend(inner_unrecaptured_time[r][idx])
                if constrained_time0[idx]==1:
                    recapture_times_3D_constrained[ri].extend( no_contact_intervals[r][idx])
                    inner_recapture_count_constrained+=len(no_contact_intervals[r][idx])
                    uncapture_times_3D_constrained[ri].extend(inner_unrecaptured_time[r][idx])
                else:
                    recapture_times_3D_unconstrained[ri].extend( no_contact_intervals[r][idx])
                    inner_recapture_count_unconstrained+=len(no_contact_intervals[r][idx])
                    uncapture_times_3D_unconstrained[ri].extend(inner_unrecaptured_time[r][idx])
                inner_recapture_count+=len(no_contact_intervals[r][idx])
                if len(no_contact_intervals[r][idx])>0:
                    inner_index_counter.append(idx)
        print(constrained_time0)
        if "frozenloop_sim" in first_passage_times_file:
            for idx in no_contact_intervals[r].keys():
                inner_chrom_intervals[ri].extend(no_contact_intervals[r][idx])
                recapture_times_3D[ri].extend( no_contact_intervals[r][idx])
                exit_times_3D[ri].extend( with_contact_intervals[r][idx])
                uncapture_times_3D[ri].extend(inner_unrecaptured_time[r][idx])
                if constrained_time0[idx]==1:
                    recapture_times_3D_constrained[ri].extend( no_contact_intervals[r][idx])
                    exit_times_3D_constrained[ri].extend( with_contact_intervals[r][idx])
                    inner_recapture_count_constrained+=len(no_contact_intervals[r][idx])
                    uncapture_times_3D_constrained[ri].extend(inner_unrecaptured_time[r][idx])
                else:
                    print(len(no_contact_intervals[r][idx]))
                    recapture_times_3D_unconstrained[ri].extend( no_contact_intervals[r][idx])
                    exit_times_3D_unconstrained[ri].extend( with_contact_intervals[r][idx])
                    inner_recapture_count_unconstrained+=len(no_contact_intervals[r][idx])
                    uncapture_times_3D_unconstrained[ri].extend(inner_unrecaptured_time[r][idx])
                inner_recapture_count+=len(no_contact_intervals[r][idx])
                if len(no_contact_intervals[r][idx])>0:
                    inner_index_counter.append(idx)


        inner_chrom_num_recapture[-1].append(len(inner_chrom_intervals[ri]))
        if inner_encountered_separated_count[r]>0:
            fraction_inner_chrom_num_recapture[-1].append(np.around(inner_recapture_count/inner_encountered_separated_count[r]*100))
            if "frozenloop_sim" in units_file_name:
                inner_recapture_count_global[ri]+=inner_recapture_count
                inner_recapture_count_global_constrained[ri]+=inner_recapture_count_constrained
                inner_recapture_count_global_unconstrained[ri]+=inner_recapture_count_unconstrained
                DSB_num_global[ri] += inner_encountered_separated_count[r]
                DSB_num_global_constrained[ri] += inner_encountered_separated_count_constrained[r]
                DSB_num_global_unconstrained[ri] += inner_encountered_separated_count_unconstrained[r]
                constrained_global[ri] += inner_encountered_separated_count_constrained[r]
                unconstrained_global[ri] += inner_encountered_separated_count_unconstrained[r]
        elif "DSB_sim" in units_file_name or "frozenloop_stochastic" in units_file_name:
            fraction_inner_chrom_num_recapture[-1].append(np.around(inner_recapture_count/int(current_num_DSBs)*100))
            inner_recapture_count_global[ri]+=inner_recapture_count
            inner_recapture_count_global_constrained[ri]+=inner_recapture_count_constrained
            inner_recapture_count_global_unconstrained[ri]+=inner_recapture_count_unconstrained
            DSB_num_global[ri] += int(current_num_DSBs)
            DSB_num_global_constrained[ri] += np.count_nonzero(constrained_time0)
            DSB_num_global_unconstrained[ri] += int(current_num_DSBs)-np.count_nonzero(constrained_time0)
            constrained_global[ri] += np.count_nonzero(constrained_time0)
            unconstrained_global[ri] += int(current_num_DSBs)-np.count_nonzero(constrained_time0)
        else:
            fraction_inner_chrom_num_recapture[-1].append(-1)
        if "brokenloop" in units_file_name:
            fraction_inner_chrom_with_recapture[-1].append(int(len(set(inner_index_counter))/int(current_num_loops)*100))
            inner_recapture_count_global[ri]+=inner_recapture_count
        elif "multichain" in units_file_name:
            fraction_inner_chrom_with_recapture[-1].append(int(len(set(inner_index_counter))/int(current_num_chrom)*100))
            inner_recapture_count_global[ri]+=inner_recapture_count
        elif "DSB_sim" in units_file_name or "frozenloop_stochastic" in units_file_name:
            if int(current_num_DSBs)>0:
                fraction_inner_chrom_with_recapture[-1].append(int(len(set(inner_index_counter))/int(current_num_DSBs)*100))
            else:
                fraction_inner_chrom_with_recapture[-1].append(-1)
        elif "frozenloop_sim" in units_file_name:
            fraction_inner_chrom_with_recapture[-1].append(int(len(set(inner_index_counter))/len(constrained_time0)*100))
            
for ri, r in enumerate(with_contact_intervals.keys()): 
    min_len = 0
    list_to_plot = np.asarray([j for j in recapture_times_3D[ri] if j>min_len])
    if "DSB_sim" in first_passage_times_file or "frozenloop" in first_passage_times_file:
        scaling_factor = inner_recapture_count_global[ri]/DSB_num_global[ri]
    else: 
        scaling_factor = fraction_inner_chrom_num_recapture[i][ri]/100
        
    x, y, l, u = edf(np.asarray(list_to_plot), alpha=0.05, scaling_factor = scaling_factor)
    ax_set1[ri,0].fill_between(x, l, u, facecolor = 'royalblue',alpha=0.5)
    ax_set1[ri,0].plot(x, y, '-',color = 'royalblue')
    ax_set1[ri,0].set_xscale('log')
    ax_set1[ri,0].set_xlabel('Time to recapture (sim time step)')
    ax_set1[ri,0].set_ylabel('Density')
    ax_set1[ri,0].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')
    ax_set1[ri,0].set_ylim([0,1.05])
    
    # convert to SI units
    list_to_plot = [j*dt_min for j in recapture_times_3D[ri] if j>min_len]
    x, y, l, u = edf(np.asarray(list_to_plot), alpha=0.05, scaling_factor = scaling_factor)
    ax_set1[ri,1].fill_between(x, l, u, facecolor = 'royalblue',alpha=0.5)
    ax_set1[ri,1].plot(x, y, '-',color = 'royalblue')
    ax_set1[ri,1].set_xscale('log')
    ax_set1[ri,1].set_xlabel('Time to recapture (minutes)')
    ax_set1[ri,1].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')
    ax_set1[ri,1].set_ylim([0,1.05])

    # Plots conditioned on whether the ends are constrained
    if "DSB_sim" in first_passage_times_file or "frozenloop" in first_passage_times_file:
        min_len = 0
        list_to_plot = np.asarray([j for j in recapture_times_3D_constrained[ri] if j>min_len])
        
        scaling_factor = inner_recapture_count_global_constrained[ri]/DSB_num_global_constrained[ri]
        x, y, l, u = edf(np.asarray(list_to_plot), alpha=0.05, scaling_factor = scaling_factor)
        ax_set5[ri,0].fill_between(x, l, u, facecolor = 'darkviolet',alpha=0.5)
        ax_set5[ri,0].plot(x, y, '-',color = 'darkviolet')
        if len(recapture_times_3D_unconstrained[ri])>0:
            list2_to_plot = np.asarray([j for j in recapture_times_3D_unconstrained[ri] if j>min_len])
            scaling_factor = inner_recapture_count_global_unconstrained[ri]/DSB_num_global_unconstrained[ri]
            x, y, l, u = edf(np.asarray(list2_to_plot), alpha=0.05, scaling_factor = scaling_factor)
            ax_set5[ri,0].fill_between(x, l, u, facecolor = 'forestgreen',alpha=0.5)
            ax_set5[ri,0].plot(x, y, '-',color = 'forestgreen')
        ax_set5[ri,0].set_xscale('log')
        ax_set5[ri,0].set_xlabel('Time to recapture (sim time step)')
        ax_set5[ri,0].set_ylabel('Density')
        ax_set5[ri,0].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')
        ax_set5[ri,0].set_ylim([0,1.05])

        # convert to SI units
        list_to_plot = np.asarray([j*dt_min for j in recapture_times_3D_constrained[ri] if j>min_len])
        print(np.sum(list_to_plot))
        scaling_factor = inner_recapture_count_global_constrained[ri]/DSB_num_global_constrained[ri]
        x, y, l, u = edf(np.asarray(list_to_plot), alpha=0.05, scaling_factor = scaling_factor)
        ax_set5[ri,1].fill_between(x, l, u, facecolor = 'darkviolet',alpha=0.5)
        ax_set5[ri,1].plot(x, y, '-',color = 'darkviolet')
        if len(recapture_times_3D_unconstrained[ri])>0:
            list2_to_plot = np.asarray([j*dt_min for j in recapture_times_3D_unconstrained[ri] if j>min_len])
            print(np.sum(list2_to_plot))
            scaling_factor = inner_recapture_count_global_unconstrained[ri]/DSB_num_global_unconstrained[ri]
            x, y, l, u = edf(np.asarray(list2_to_plot), alpha=0.05, scaling_factor = scaling_factor)
            ax_set5[ri,1].fill_between(x, l, u, facecolor = 'forestgreen',alpha=0.5)
            ax_set5[ri,1].plot(x, y, '-',color = 'forestgreen')
        ax_set5[ri,1].set_xscale('log')
        ax_set5[ri,1].set_xlabel('Time to recapture (sim time step)')
        ax_set5[ri,1].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')
        ax_set5[ri,1].set_ylim([0,1.05])

    all_intervals = np.asarray(recapture_times_3D[ri]+uncapture_times_3D[ri])

    censored = np.concatenate((np.zeros(len(recapture_times_3D[ri])),np.ones(len(uncapture_times_3D[ri]))),axis=None)
    
    KM_survival_curve = KM_survival(all_intervals,censored,conf=0.95, Tmax=np.inf, S1at=1/dt_min)#calculating KM survival curve
    if 0.5 in KM_survival_curve[:,1]:
        median_idx = np.argwhere(KM_survival_curve[:,1]==0.5)[0][0]
        median = KM_survival_curve[median_idx,0]
    else:
        a= np.argwhere(KM_survival_curve[:,1]-0.5<0)
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


    ax_set3[ri,0].errorbar(radius,mean,yerr=[np.asarray([np.abs(mean_lbnd-mean)]),np.asarray([np.abs(mean_ubnd-mean)])],fmt='o', capthick=2, color = 'royalblue')
    ax_set3[ri,0].set_xlabel('Confinement radius (nm)')
    ax_set3[ri,0].set_ylabel('Mean time to recapture (sim time step)')
    ax_set3[ri,0].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')

    
    ax_set3[ri,1].errorbar(radius,mean*dt_min,yerr=[np.asarray([np.abs(mean_lbnd-mean)*dt_min]),np.asarray([np.abs(mean_ubnd-mean)*dt_min])],fmt='o', capthick=2, color = 'royalblue')
    ax_set3[ri,1].set_xlabel('Confinement radius (nm)')
    ax_set3[ri,1].set_ylabel('Mean time to recapture (minutes)')
    ax_set3[ri,1].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')
    
    # compute average recapture time conditioned on constrained or not
    if "DSB_sim" in first_passage_times_file or "frozenloop" in first_passage_times_file:
        all_intervals = np.asarray(recapture_times_3D_constrained[ri]+uncapture_times_3D_constrained[ri])

        censored = np.concatenate((np.zeros(len(recapture_times_3D_constrained[ri])),np.ones(len(uncapture_times_3D_constrained[ri]))),axis=None)

        KM_survival_curve = KM_survival(all_intervals,censored,conf=0.95, Tmax=np.inf, S1at=1/dt_min)#calculating KM survival curve
        if 0.5 in KM_survival_curve[:,1]:
            median_idx = np.argwhere(KM_survival_curve[:,1]==0.5)[0][0]
            median = KM_survival_curve[median_idx,0]
        else:
            a= np.argwhere(KM_survival_curve[:,1]-0.5<0)
#             print(KM_survival_curve[:,1])
#             print(censored)
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


        ax_set6[ri,0].errorbar(radius,mean,yerr=[np.asarray([np.abs(mean_lbnd-mean)]),np.asarray([np.abs(mean_ubnd-mean)])],fmt='o', capthick=2, color = 'darkviolet')
        ax_set6[ri,0].set_xlabel('Confinement radius (nm)')
        ax_set6[ri,0].set_ylabel('Mean time to recapture (sim time step)')
        ax_set6[ri,0].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')


        ax_set6[ri,1].errorbar(radius,mean*dt_min,yerr=[np.asarray([np.abs(mean_lbnd-mean)*dt_min]),np.asarray([np.abs(mean_ubnd-mean)*dt_min])],fmt='o', capthick=2, color = 'darkviolet')
        ax_set6[ri,1].set_xlabel('Confinement radius (nm)')
        ax_set6[ri,1].set_ylabel('Mean time to recapture (minutes)')
        ax_set6[ri,1].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')
        
        if len(recapture_times_3D_unconstrained[ri])>0:
            all_intervals = np.asarray(recapture_times_3D_unconstrained[ri]+uncapture_times_3D_unconstrained[ri])

            censored = np.concatenate((np.zeros(len(recapture_times_3D_unconstrained[ri])),np.ones(len(uncapture_times_3D_unconstrained[ri]))),axis=None)

            KM_survival_curve = KM_survival(all_intervals,censored,conf=0.95, Tmax=np.inf, S1at=1/dt_min)#calculating KM survival curve
            if 0.5 in KM_survival_curve[:,1]:
                median_idx = np.argwhere(KM_survival_curve[:,1]==0.5)[0][0]
                median = KM_survival_curve[median_idx,0]
            else:
                a= np.argwhere(KM_survival_curve[:,1]-0.5<0)
    #             print(KM_survival_curve[:,1])
    #             print(censored)
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


            ax_set6[ri,0].errorbar(radius,mean,yerr=[np.asarray([np.abs(mean_lbnd-mean)]),np.asarray([np.abs(mean_ubnd-mean)])],fmt='o', capthick=2, color = 'forestgreen')
            ax_set6[ri,0].set_xlabel('Confinement radius (nm)')
            ax_set6[ri,0].set_ylabel('Mean time to recapture (sim time step)')
            ax_set6[ri,0].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')


            ax_set6[ri,1].errorbar(radius,mean*dt_min,yerr=[np.asarray([np.abs(mean_lbnd-mean)*dt_min]),np.asarray([np.abs(mean_ubnd-mean)*dt_min])],fmt='o', capthick=2, color = 'forestgreen')
            ax_set6[ri,1].set_xlabel('Confinement radius (nm)')
            ax_set6[ri,1].set_ylabel('Mean time to recapture (minutes)')
            ax_set6[ri,1].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')
        
    if i == len(snakemake.input.units_file)-1:

        for Ri, R in enumerate(no_contact_intervals.keys()):  
            if "DSB_sim" in first_passage_times_file or "frozenloop" in first_passage_times_file:
                legend_text = ['%synapsis = '+str(np.around(inner_recapture_count_global[Ri]/DSB_num_global[Ri]*100,2))+', %constrained = '+str(np.around(constrained_global[Ri]/DSB_num_global[Ri]*100,2))+', # of DSB = '+str(int(DSB_num_global[Ri]))]
                legend_text2 = ['%synapsis = '+str(np.around(inner_recapture_count_global_constrained[Ri]/DSB_num_global_constrained[Ri]*100,2))+', # of constrained DSB = '+str(int(DSB_num_global_constrained[Ri]))]
                
                legend_text2.append('%synapsis = '+str(np.around(inner_recapture_count_global_unconstrained[Ri]/DSB_num_global_unconstrained[Ri]*100,2))+', # of unconstrained DSB = '+str(int(DSB_num_global_unconstrained[Ri])))
                ax_set5[Ri,1].legend(legend_text2,loc='lower right',framealpha=0.2)
                ax_set6[Ri,1].legend(legend_text2,loc='lower right',framealpha=0.2)
                
            else:
                legend_text = ['n = '+str(inner_recapture_count_global[Ri]) +', %synapsis = '+str(np.around(fraction_inner_chrom_num_recapture[i][Ri],2))+', # of DSB = '+str(int(inner_encountered_separated_count[r]))]
            ax_set1[Ri,1].legend(legend_text,loc='lower right',framealpha=0.2)
            ax_set3[Ri,1].legend(legend_text,loc='lower right',framealpha=0.2)
            
                

    if "frozenloop_sim" in first_passage_times_file:
        min_len = 0
        list_to_plot = np.asarray([j for j in exit_times_3D_constrained[ri] if j>min_len])
        scaling_factor = 1
        x, y, l, u = edf(np.asarray(list_to_plot), alpha=0.05, scaling_factor = scaling_factor)
        ax_set2[ri,0].fill_between(x, l, u, facecolor = 'darkviolet',alpha=0.5)
        ax_set2[ri,0].plot(x, y, '-',color = 'darkviolet')
        if len(exit_times_3D_unconstrained[ri])>0:
            list2_to_plot = np.asarray([j for j in exit_times_3D_unconstrained[ri] if j>min_len])
            scaling_factor = inner_recapture_count_global_unconstrained[ri]/DSB_num_global_unconstrained[ri]
            x, y, l, u = edf(np.asarray(list2_to_plot), alpha=0.05, scaling_factor = scaling_factor)
            ax_set2[ri,0].fill_between(x, l, u, facecolor = 'forestgreen',alpha=0.5)
            ax_set2[ri,0].plot(x, y, '-',color = 'forestgreen')
        ax_set2[ri,0].set_xscale('log')
        ax_set2[ri,0].set_xlabel('Time to exit (sim time step)')
        ax_set2[ri,0].set_ylabel('Density')
        ax_set2[ri,0].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')
        ax_set2[ri,0].set_ylim([0,1.05])

        # convert to SI units
        list_to_plot = np.asarray([j*dt_min for j in exit_times_3D_constrained[ri] if j>min_len])
        scaling_factor = 1
        x, y, l, u = edf(np.asarray(list_to_plot), alpha=0.05, scaling_factor = scaling_factor)
        ax_set2[ri,1].fill_between(x, l, u, facecolor = 'darkviolet',alpha=0.5)
        ax_set2[ri,1].plot(x, y, '-',color = 'darkviolet')
        if len(exit_times_3D_unconstrained[ri])>0:
            list2_to_plot = np.asarray([j*dt_min for j in exit_times_3D_unconstrained[ri] if j>min_len])
            scaling_factor = inner_recapture_count_global_unconstrained[ri]/DSB_num_global_unconstrained[ri]
            x, y, l, u = edf(np.asarray(list2_to_plot), alpha=0.05, scaling_factor = scaling_factor)
            ax_set2[ri,1].fill_between(x, l, u, facecolor = 'forestgreen',alpha=0.5)
            ax_set2[ri,1].plot(x, y, '-',color = 'forestgreen')
        ax_set2[ri,1].set_xscale('log')
        ax_set2[ri,1].set_xlabel('Time to exit (sim time step)')
        ax_set2[ri,1].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')
        ax_set2[ri,1].set_ylim([0,1.05])
        
        all_intervals = np.asarray(exit_times_3D_constrained[ri])

        censored = np.concatenate((np.zeros(len(exit_times_3D_constrained[ri]))),axis=None)

        KM_survival_curve = KM_survival(all_intervals,censored,conf=0.95, Tmax=np.inf, S1at=1/dt_min)#calculating KM survival curve
        if 0.5 in KM_survival_curve[:,1]:
            median_idx = np.argwhere(KM_survival_curve[:,1]==0.5)[0][0]
            median = KM_survival_curve[median_idx,0]
        else:
            a= np.argwhere(KM_survival_curve[:,1]-0.5<0)
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


        ax_set4[ri,0].errorbar(radius,mean,yerr=[np.asarray([np.abs(mean_lbnd-mean)]),np.asarray([np.abs(mean_ubnd-mean)])],fmt='o', capthick=2, color = 'darkviolet')
        ax_set4[ri,0].set_xlabel('Confinement radius (nm)')
        ax_set4[ri,0].set_ylabel('Mean time to exit (sim time step)')
        ax_set4[ri,0].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')


        ax_set4[ri,1].errorbar(radius,mean*dt_min,yerr=[np.asarray([np.abs(mean_lbnd-mean)*dt_min]),np.asarray([np.abs(mean_ubnd-mean)*dt_min])],fmt='o', capthick=2, color = 'darkviolet')
        ax_set4[ri,1].set_xlabel('Confinement radius (nm)')
        ax_set4[ri,1].set_ylabel('Mean time to exit (minutes)')
        ax_set4[ri,1].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')
        
        if len(exit_times_3D_unconstrained[ri])>0:
            all_intervals = np.asarray(exit_times_3D_unconstrained[ri])

            censored = np.concatenate((np.zeros(len(exit_times_3D_unconstrained[ri]))),axis=None)

            KM_survival_curve = KM_survival(all_intervals,censored,conf=0.95, Tmax=np.inf, S1at=1/dt_min)#calculating KM survival curve
            if 0.5 in KM_survival_curve[:,1]:
                median_idx = np.argwhere(KM_survival_curve[:,1]==0.5)[0][0]
                median = KM_survival_curve[median_idx,0]
            else:
                a= np.argwhere(KM_survival_curve[:,1]-0.5<0)
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


            ax_set4[ri,0].errorbar(radius,mean,yerr=[np.asarray([np.abs(mean_lbnd-mean)]),np.asarray([np.abs(mean_ubnd-mean)])],fmt='o', capthick=2, color = 'forestgreen')
            ax_set4[ri,0].set_xlabel('Confinement radius (nm)')
            ax_set4[ri,0].set_ylabel('Mean time to exit (sim time step)')
            ax_set4[ri,0].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')


            ax_set4[ri,1].errorbar(radius,mean*dt_min,yerr=[np.asarray([np.abs(mean_lbnd-mean)*dt_min]),np.asarray([np.abs(mean_ubnd-mean)*dt_min])],fmt='o', capthick=2, color = 'forestgreen')
            ax_set4[ri,1].set_xlabel('Confinement radius (nm)')
            ax_set4[ri,1].set_ylabel('Mean time to exit (minutes)')
            ax_set4[ri,1].set_title('Capture radius = '+ str(np.around(r,decimals=1))+ ' monomer')

    if i == len(snakemake.input.units_file)-1:

        for ri, r in enumerate(with_contact_intervals.keys()):  
            legend_text = ['# of exits = '+ str(len( exit_times_3D_constrained[ri]))+', # of unconstrained exits = '+str(len( exit_times_3D_unconstrained[ri]))]
            ax_set2[ri,0].legend(legend_text,loc='lower right',framealpha=0.2)
            ax_set4[ri,0].legend(legend_text,loc='lower right',framealpha=0.2)
            ax_set2[ri,1].legend(legend_text,loc='lower right',framealpha=0.2)
            ax_set4[ri,1].legend(legend_text,loc='lower right',framealpha=0.2)

recapture_times_3D_radius3 = np.asarray([j*dt_min for j in recapture_times_3D[0]])
uncapture_times_3D_radius3 = np.asarray([j*dt_min for j in uncapture_times_3D[0]])

recapture_times_3D_radius4 = np.asarray([j*dt_min for j in recapture_times_3D[1]])
uncapture_times_3D_radius4 = np.asarray([j*dt_min for j in uncapture_times_3D[1]])

recapture_times_3D_radius5 = np.asarray([j*dt_min for j in recapture_times_3D[2]])
uncapture_times_3D_radius5 = np.asarray([j*dt_min for j in uncapture_times_3D[2]])
if "DSB_sim" in first_passage_times_file or "frozenloop" in first_passage_times_file:
    recapture_times_3D_constrained_radius3 = np.asarray([j*dt_min for j in recapture_times_3D_constrained[0]])
    recapture_times_3D_unconstrained_radius3 = np.asarray([j*dt_min for j in recapture_times_3D_unconstrained[0]])
    uncapture_times_3D_constrained_radius3 = np.asarray([j*dt_min for j in uncapture_times_3D_constrained[0]])
    uncapture_times_3D_unconstrained_radius3 = np.asarray([j*dt_min for j in uncapture_times_3D_unconstrained[0]])
    
    recapture_times_3D_constrained_radius4 = np.asarray([j*dt_min for j in recapture_times_3D_constrained[1]])
    recapture_times_3D_unconstrained_radius4 = np.asarray([j*dt_min for j in recapture_times_3D_unconstrained[1]])
    uncapture_times_3D_constrained_radius4 = np.asarray([j*dt_min for j in uncapture_times_3D_constrained[1]])
    uncapture_times_3D_unconstrained_radius4 = np.asarray([j*dt_min for j in uncapture_times_3D_unconstrained[1]])
    
    recapture_times_3D_constrained_radius5 = np.asarray([j*dt_min for j in recapture_times_3D_constrained[2]])
    recapture_times_3D_unconstrained_radius5 = np.asarray([j*dt_min for j in recapture_times_3D_unconstrained[2]])
    uncapture_times_3D_constrained_radius5 = np.asarray([j*dt_min for j in uncapture_times_3D_constrained[2]])
    uncapture_times_3D_unconstrained_radius5 = np.asarray([j*dt_min for j in uncapture_times_3D_unconstrained[2]])
    
file = snakemake.output.first_recapture_times
with open(file, 'wb') as g:
    np.save(g, recapture_times_3D_radius4)#capture radius 4
    if "DSB_sim" in first_passage_times_file or "frozenloop" in first_passage_times_file:
        np.save(g, recapture_times_3D_constrained_radius4)#capture radius 4
        np.save(g, recapture_times_3D_unconstrained_radius4)#capture radius 4
    np.save(g, uncapture_times_3D_radius4)#capture radius 4
    if "DSB_sim" in first_passage_times_file or "frozenloop" in first_passage_times_file:
        np.save(g, uncapture_times_3D_constrained_radius4)#capture radius 4
        np.save(g, uncapture_times_3D_unconstrained_radius4)#capture radius 4
        
file = snakemake.output.first_recapture_times_radius3
with open(file, 'wb') as g:
    np.save(g, recapture_times_3D_radius3)#capture radius 3
    if "DSB_sim" in first_passage_times_file or "frozenloop" in first_passage_times_file:
        np.save(g, recapture_times_3D_constrained_radius3)#capture radius 3
        np.save(g, recapture_times_3D_unconstrained_radius3)#capture radius 3
    np.save(g, uncapture_times_3D_radius3)#capture radius 3
    if "DSB_sim" in first_passage_times_file or "frozenloop" in first_passage_times_file:
        np.save(g, uncapture_times_3D_constrained_radius3)#capture radius 3
        np.save(g, uncapture_times_3D_unconstrained_radius3)#capture radius 3
        
file = snakemake.output.first_recapture_times_radius5
with open(file, 'wb') as g:
    np.save(g, recapture_times_3D_radius5)#capture radius 5
    if "DSB_sim" in first_passage_times_file or "frozenloop" in first_passage_times_file:
        np.save(g, recapture_times_3D_constrained_radius5)#capture radius 5
        np.save(g, recapture_times_3D_unconstrained_radius5)#capture radius 5
    np.save(g, uncapture_times_3D_radius5)#capture radius 5
    if "DSB_sim" in first_passage_times_file or "frozenloop" in first_passage_times_file:
        np.save(g, uncapture_times_3D_constrained_radius5)#capture radius 5
        np.save(g, uncapture_times_3D_unconstrained_radius5)#capture radius 5


fig1.savefig(snakemake.output.first_recapture_combined_v2,bbox_inches='tight')
fig2.savefig(snakemake.output.first_exit_combined_v2,bbox_inches='tight')
fig3.savefig(snakemake.output.first_recapture_vs_confinementR_v2,bbox_inches='tight')
fig4.savefig(snakemake.output.first_exit_vs_confinementR_v2,bbox_inches='tight')
fig6.savefig(snakemake.output.first_recapture_vs_confinementR_conditioned_v2,bbox_inches='tight')
fig5.savefig(snakemake.output.first_recapture_combined_conditioned_v2,bbox_inches='tight')