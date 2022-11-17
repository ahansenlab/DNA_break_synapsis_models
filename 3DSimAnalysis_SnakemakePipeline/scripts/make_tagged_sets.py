from polychrom.hdf5_format import list_URIs, load_URI
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# parse file name for the number of chromosomes
folder_name = snakemake.input.simulation_folder

URI_start = snakemake.config['DSB_simulation_parameters']["initial_URIs_skipped"]
last_URI_used = snakemake.config['DSB_simulation_parameters']["last_URIs_used"]
skip_every = snakemake.config['DSB_simulation_parameters']["skip_every_X_URIs"]


files = list_URIs(folder_name)

fig0 = plt.figure(figsize=(6,6))
ax_set0 = fig0.subplots(1, 1)

fig = plt.figure(figsize=(6,6))
ax_set = fig.subplots(1, 1)

if 'multichain' in folder_name:
    # look for additional URIs from extended runs
    num_monomers = int(folder_name.split('_N')[1].split('_C')[0])

    end_point_dict={70000:9999999}
    if num_monomers in end_point_dict.keys():
        files = files[:end_point_dict[num_monomers]+1]

    date = folder_name.split('sim_')[1].split('_dens')[0]
    parameter_string = folder_name.split(date)[1].split('_v5')[0]
    date = date.replace("_","")

    for path in os.listdir("/mnt/md0/simulations_ouput/DSB_repair/"):
        if "ContinuedFrom"+date in path and parameter_string in path:
            additional_folder = path
            additional_files = list_URIs("results/"+additional_folder)
            files.extend(additional_files)

    print(len(files))
    
# this grabs the "ends of the chromosomes"
x0_ends = []
y0_ends = []
z0_ends = []
x1_ends = []
y1_ends = []
z1_ends = []
# we want to use the Fbn2_tad_length to calibrate the MSD from simulations
x0_terminalTAD = []
y0_terminalTAD = []
z0_terminalTAD = []
x1_terminalTAD = []
y1_terminalTAD = []
z1_terminalTAD = []
x0_internalTAD = []
y0_internalTAD = []
z0_internalTAD = []
x1_internalTAD = []
y1_internalTAD = []
z1_internalTAD = []

idx = []
frame = []

idx_2 = []
frame_2 = []
if 'multichain' in snakemake.input.simulation_folder:
    
    num_chroms = int(folder_name.split('_C')[1].split('_collisionrate')[0])
    
    if num_chroms>100:
        for fi,f in enumerate(files[URI_start::skip_every]):
            u = load_URI(f)  
            for c in range(0,num_chroms):
                x0_ends.append(u['pos'][num_chroms+c][0]) # get the chromosome end positions
                y0_ends.append(u['pos'][num_chroms+c][1]) # get the chromosome end positions
                z0_ends.append(u['pos'][num_chroms+c][2]) # get the chromosome end positions
                x1_ends.append(u['pos'][c][0]) # get the chromosome end positions
                y1_ends.append(u['pos'][c][1]) # get the chromosome end positions
                z1_ends.append(u['pos'][c][2]) # get the chromosome end positions                


                idx.append(c)
                frame.append(fi)           

        df = pd.DataFrame({'id':idx,
                           't':frame,
                           'x0':x0_ends,
                           'y0':y0_ends,
                           'z0':z0_ends, 
                           'x1':x1_ends,
                           'y1':y1_ends,
                           'z1':z1_ends, })
        df.to_csv(snakemake.output.chrom_ends,index=False) 

        df.to_csv(snakemake.output.terminal_TAD,index=False) 

        df.to_csv(snakemake.output.internal_TAD,index=False)
        
    elif num_chroms == 1:
        for fi,f in enumerate(files[-last_URI_used::skip_every]):

            u = load_URI(f)  
            for c in range(0,num_chroms):
                x0_ends.append(u['pos'][num_chroms+c][0]) # get the chromosome end positions
                y0_ends.append(u['pos'][num_chroms+c][1]) # get the chromosome end positions
                z0_ends.append(u['pos'][num_chroms+c][2]) # get the chromosome end positions
                x1_ends.append(u['pos'][c][0]) # get the chromosome end positions
                y1_ends.append(u['pos'][c][1]) # get the chromosome end positions
                z1_ends.append(u['pos'][c][2]) # get the chromosome end positions                
                x0_terminalTAD.append(u['pos'][2*num_chroms+c][0])
                y0_terminalTAD.append(u['pos'][2*num_chroms+c][1])
                z0_terminalTAD.append(u['pos'][2*num_chroms+c][2])
                x1_terminalTAD.append(u['pos'][3*num_chroms+c][0])
                y1_terminalTAD.append(u['pos'][3*num_chroms+c][1])
                z1_terminalTAD.append(u['pos'][3*num_chroms+c][2])
                x0_internalTAD.append(u['pos'][2*num_chroms+c][0])
                y0_internalTAD.append(u['pos'][2*num_chroms+c][1])
                z0_internalTAD.append(u['pos'][2*num_chroms+c][2])
                x1_internalTAD.append(u['pos'][3*num_chroms+c][0])
                y1_internalTAD.append(u['pos'][3*num_chroms+c][1])
                z1_internalTAD.append(u['pos'][3*num_chroms+c][2])

                idx.append(c)
                frame.append(fi)           

        df = pd.DataFrame({'id':idx,
                           't':frame,
                           'x0':x0_ends,
                           'y0':y0_ends,
                           'z0':z0_ends, 
                           'x1':x1_ends,
                           'y1':y1_ends,
                           'z1':z1_ends, })
        df.to_csv(snakemake.output.chrom_ends,index=False) 

        df = pd.DataFrame({'id':idx,
                           't':frame,
                           'x0':x0_terminalTAD,
                           'y0':y0_terminalTAD,
                           'z0':z0_terminalTAD, 
                           'x1':x1_terminalTAD,
                           'y1':y1_terminalTAD,
                           'z1':z1_terminalTAD, })
        df.to_csv(snakemake.output.terminal_TAD,index=False) 

        df = pd.DataFrame({'id':idx,
                           't':frame,
                           'x0':x0_internalTAD,
                           'y0':y0_internalTAD,
                           'z0':z0_internalTAD, 
                           'x1':x1_internalTAD,
                           'y1':y1_internalTAD,
                           'z1':z1_internalTAD, })
        df.to_csv(snakemake.output.internal_TAD,index=False)         
        
    else:
        cm = plt.cm.get_cmap('Purples')
        fig_counter = 1
        fig_legend = []
        for fi,f in enumerate(files[-last_URI_used::skip_every]):
            u = load_URI(f)  
            for c in range(0,num_chroms):
                x0_ends.append(u['pos'][num_chroms+c][0]) # get the chromosome end positions
                y0_ends.append(u['pos'][num_chroms+c][1]) # get the chromosome end positions
                z0_ends.append(u['pos'][num_chroms+c][2]) # get the chromosome end positions
                x1_ends.append(u['pos'][c][0]) # get the chromosome end positions
                y1_ends.append(u['pos'][c][1]) # get the chromosome end positions
                z1_ends.append(u['pos'][c][2]) # get the chromosome end positions                
                x0_terminalTAD.append(u['pos'][2*num_chroms+c][0])
                y0_terminalTAD.append(u['pos'][2*num_chroms+c][1])
                z0_terminalTAD.append(u['pos'][2*num_chroms+c][2])
                x1_terminalTAD.append(u['pos'][3*num_chroms+c][0])
                y1_terminalTAD.append(u['pos'][3*num_chroms+c][1])
                z1_terminalTAD.append(u['pos'][3*num_chroms+c][2])
                x0_internalTAD.append(u['pos'][2*num_chroms+c][0])
                y0_internalTAD.append(u['pos'][2*num_chroms+c][1])
                z0_internalTAD.append(u['pos'][2*num_chroms+c][2])
                x1_internalTAD.append(u['pos'][3*num_chroms+c][0])
                y1_internalTAD.append(u['pos'][3*num_chroms+c][1])
                z1_internalTAD.append(u['pos'][3*num_chroms+c][2])


                idx.append(c)
                frame.append(fi)      
            if np.mod(fi,10)==0 and fi<=150:
                # compute distances between all monomer ends    
                dists= np.sqrt((np.asarray(x0_ends)-np.asarray(x1_ends))**2+(np.asarray(y0_ends)-np.asarray(y1_ends))**2+(np.asarray(z0_ends)-np.asarray(z1_ends))**2)
                ax_set.hist(dists,cumulative=True,density=True,histtype='step', color = cm(15*fig_counter),linewidth=4) 
                fig_counter += 1
                fig_legend.append('time = '+str(fi*1000)+' time steps, average distance = '+str(np.around(np.mean(dists),2)))
                
                if fig_counter == 1:
                    ax_set1[ri,0].set_xlabel('Distance between DSB ends (monomer)')
                    ax_set1[ri,0].set_title('Distribution of DSB end distance over time')
                    ax_set1[ri,0].set_ylim([0,1.05])
            ax_set.legend(fig_legend,loc='lower right',framealpha=0.2)

        df = pd.DataFrame({'id':idx,
                           't':frame,
                           'x0':x0_ends,
                           'y0':y0_ends,
                           'z0':z0_ends, 
                           'x1':x1_ends,
                           'y1':y1_ends,
                           'z1':z1_ends, })
        df.to_csv(snakemake.output.chrom_ends,index=False) 

        df = pd.DataFrame({'id':idx,
                           't':frame,
                           'x0':x0_terminalTAD,
                           'y0':y0_terminalTAD,
                           'z0':z0_terminalTAD, 
                           'x1':x1_terminalTAD,
                           'y1':y1_terminalTAD,
                           'z1':z1_terminalTAD, })
        df.to_csv(snakemake.output.terminal_TAD,index=False) 

        df = pd.DataFrame({'id':idx,
                           't':frame,
                           'x0':x0_internalTAD,
                           'y0':y0_internalTAD,
                           'z0':z0_internalTAD, 
                           'x1':x1_internalTAD,
                           'y1':y1_internalTAD,
                           'z1':z1_internalTAD, })
        df.to_csv(snakemake.output.internal_TAD,index=False)         



    import tracklib as tl

    data = tl.TaggedSet()

    out_names = [snakemake.output.chrom_ends,snakemake.output.terminal_TAD,snakemake.output.internal_TAD]
    out_tags = ["chromosome_ends","terminal_TAD","internal_TAD"]

elif 'DSB_sim' or 'frozenloop' in snakemake.input.simulation_folder:
    URI_start = 1000 # restartBondUpdaterEveryBlocks * BondUpdaterCountBeforeDSB
    if 'frozenloop_sim' in snakemake.input.simulation_folder:
        URI_start = 3380000
        
    with open(snakemake.input.simulation_folder+'/DSB_boundary_coordinates.npy',"rb") as f:
        DSB_coordinates = np.load(f)

    num_DSBs = int(len(DSB_coordinates)/2)
 
    num_fbn2_tad = 35 # number of fbn2 tad tracked
    u_test = load_URI(files[URI_start]) 
    max_index = min(num_fbn2_tad,int(len(u_test['pos'])/2))
    for fi,f in enumerate(files[URI_start::skip_every]):
        u = load_URI(f)  
        for c in range(0,max_index):
            x0_ends.append(u['pos'][-2*c-1][0]) # get the chromosome end positions
            y0_ends.append(u['pos'][-2*c-1][1]) # get the chromosome end positions
            z0_ends.append(u['pos'][-2*c-1][2]) # get the chromosome end positions
            x1_ends.append(u['pos'][-2*c-2][0]) # get the chromosome end positions
            y1_ends.append(u['pos'][-2*c-2][1]) # get the chromosome end positions
            z1_ends.append(u['pos'][-2*c-2][2]) # get the chromosome end positions                
            x0_terminalTAD.append(u['pos'][-2*c-1][0]) 
            y0_terminalTAD.append(u['pos'][-2*c-1][1])
            z0_terminalTAD.append(u['pos'][-2*c-1][2])
            x1_terminalTAD.append(u['pos'][-2*c-2][0])
            y1_terminalTAD.append(u['pos'][-2*c-2][1])
            z1_terminalTAD.append(u['pos'][-2*c-2][2])
            x0_internalTAD.append(u['pos'][-2*c-1][0]) 
            y0_internalTAD.append(u['pos'][-2*c-1][1])
            z0_internalTAD.append(u['pos'][-2*c-1][2])
            x1_internalTAD.append(u['pos'][-2*c-2][0])
            y1_internalTAD.append(u['pos'][-2*c-2][1])
            z1_internalTAD.append(u['pos'][-2*c-2][2])


            idx.append(c)
            frame.append(fi)           

    df = pd.DataFrame({'id':idx,
                       't':frame,
                       'x0':x0_ends,
                       'y0':y0_ends,
                       'z0':z0_ends, 
                       'x1':x1_ends,
                       'y1':y1_ends,
                       'z1':z1_ends, })
    df.to_csv(snakemake.output.chrom_ends,index=False) 

    df = pd.DataFrame({'id':idx,
                       't':frame,
                       'x0':x0_terminalTAD,
                       'y0':y0_terminalTAD,
                       'z0':z0_terminalTAD, 
                       'x1':x1_terminalTAD,
                       'y1':y1_terminalTAD,
                       'z1':z1_terminalTAD, })
    df.to_csv(snakemake.output.terminal_TAD,index=False) 

    df = pd.DataFrame({'id':idx,
                       't':frame,
                       'x0':x0_internalTAD,
                       'y0':y0_internalTAD,
                       'z0':z0_internalTAD, 
                       'x1':x1_internalTAD,
                       'y1':y1_internalTAD,
                       'z1':z1_internalTAD, })
    df.to_csv(snakemake.output.internal_TAD,index=False)         



    import tracklib as tl

    data = tl.TaggedSet()

    out_names = [snakemake.output.chrom_ends,snakemake.output.terminal_TAD,snakemake.output.internal_TAD]
    out_tags = ["chromosome_ends","terminal_TAD","internal_TAD"]

    
elif 'brokenloop' in snakemake.input.simulation_folder:
    
    num_breaks = int(folder_name.split('_DSB')[1].split('_LoopSize')[0]) # number of DSBs
    URI_start = 3380000
    for fi,f in enumerate(files[URI_start::skip_every]):
        u = load_URI(f)  
        for c in range(0,num_breaks):
            x0_ends.append(u['pos'][2*c][0]) # get the broken DSB end positions
            y0_ends.append(u['pos'][2*c][1]) # get the broken DSB end positions
            z0_ends.append(u['pos'][2*c][2]) # get the broken DSB end positions
            x1_ends.append(u['pos'][2*c+1][0]) # get the broken DSB end positions
            y1_ends.append(u['pos'][2*c+1][1]) # get the broken DSB end positions
            z1_ends.append(u['pos'][2*c+1][2]) # get the broken DSB end positions 
            
            idx.append(c)
            frame.append(fi)    
            
        for c in range(0,int(num_breaks/2)):
            
            x0_terminalTAD.append(u['pos'][2*num_breaks+2*c][0])
            y0_terminalTAD.append(u['pos'][2*num_breaks+2*c][1])
            z0_terminalTAD.append(u['pos'][2*num_breaks+2*c][2])
            x1_terminalTAD.append(u['pos'][2*num_breaks+2*c+1][0])
            y1_terminalTAD.append(u['pos'][2*num_breaks+2*c+1][1])
            z1_terminalTAD.append(u['pos'][2*num_breaks+2*c+1][2])

            idx_2.append(c)
            frame_2.append(fi)           

    df = pd.DataFrame({'id':idx,
                       't':frame,
                       'x0':x0_ends,
                       'y0':y0_ends,
                       'z0':z0_ends, 
                       'x1':x1_ends,
                       'y1':y1_ends,
                       'z1':z1_ends, })
    df.to_csv(snakemake.output.chrom_ends,index=False) 
    
    df = pd.DataFrame({'id':idx_2,
                       't':frame_2,
                       'x0':x0_terminalTAD,
                       'y0':y0_terminalTAD,
                       'z0':z0_terminalTAD, 
                       'x1':x1_terminalTAD,
                       'y1':y1_terminalTAD,
                       'z1':z1_terminalTAD, })

    df.to_csv(snakemake.output.terminal_TAD,index=False) 
    df.to_csv(snakemake.output.internal_TAD,index=False)  
    
    import tracklib as tl

    data = tl.TaggedSet()

    out_names = [snakemake.output.chrom_ends,snakemake.output.chrom_ends,snakemake.output.chrom_ends]
    out_tags = ["chromosome_ends","terminal_TAD","internal_TAD"]
    
for tag, fname in zip(out_tags,out_names):
    # combine all the data from the sample group into a single tagged set
    # for tag_file in snakemake.input.tagged_set:
    this_data = tl.io.load.csv(fname, 
                            ['id', 't', 'x', 'y', 'z', 'x2', 'y2', 'z2'], 
                            skip_header=1,
                            delimiter=',',
                            tags=[tag])
    data |= this_data
    
# compute the displacements in time
data = data.process(lambda traj: traj.relative())  

for tag in out_tags:
    data.makeSelection()
    data.makeSelection(tags=tag)
    msd = tl.analysis.MSD(data)
    ax_set0.loglog(msd,label=tag)

ax_set0.legend()
ax_set0.set_xlabel('Time (sim units)')
ax_set0.set_ylabel('MSD (monomers^2)')
fig0.savefig(snakemake.output.figure,bbox_inches='tight')
fig.savefig(snakemake.output.figure_2,bbox_inches='tight')
    
