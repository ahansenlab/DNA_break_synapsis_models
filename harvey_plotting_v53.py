import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Arc
from matplotlib.patches import FancyArrowPatch
from matplotlib.patches import FancyBboxPatch
from matplotlib.patches import Circle
from matplotlib.patches import Ellipse
from matplotlib.patches import Rectangle
import matplotlib.transforms as mtransforms

from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42

def plot_arch_diagram(ax,xrange,yrange,xlim_low,xlim_high,break_site,scaling_left,scaling_right,
    l, r, n_lef=0,
    height_factor=1.0,
    max_height = 150,
    height=None,
    y_offset=0,
    plot_text=True,
    color=(223.0/255.0,90/255.0,73/255.0),
    lw=2,
    alpha=1,
    ls = '-',
     zorder=0):
    """Visualize an individual loop with an arch diagram.
    """
    #apply scaling to make sure the motors don't go past the DSB ends
    if l <=0:
        l=(l-(xlim_low-break_site))*scaling_left+(xlim_low-break_site)
    else:
        l=xlim_high-break_site-(xlim_high-break_site-l)*scaling_right
            
    if r <=0:
        r=(r-(xlim_low-break_site))*scaling_left+(xlim_low-break_site)
    else:
        r=xlim_high-break_site-(xlim_high-break_site-r)*scaling_right
        
    arc_center = ((l+r)/2,y_offset)
    
    if height is None:
        arc_height = np.nanmin([max_height, np.max(r-l)/2.0*height_factor])
    else:
        arc_height = height
        
    # LEFs represented by arc
    arc = Arc(xy=arc_center,
            width=r-l,
            height=arc_height,
            theta1=0,
            theta2=180,
            alpha=alpha,
            lw=lw,
            color=color,
            capstyle='round',
            ls = ls,  
            zorder = zorder)
    
    ax.add_patch(arc)
    
    left_motor = Ellipse((l,yrange*0.02*1.5+yrange*0.015), 2*xrange*0.015,2*yrange*0.015, color = color,zorder=zorder+1)
    ax.add_patch(left_motor)
    
    right_motor = Ellipse((r,yrange*0.02*1.5+yrange*0.015), 2*xrange*0.015,2*yrange*0.015, color = color,zorder=zorder+2)
    ax.add_patch(right_motor)

    if n_lef > 1 and arc_center[0] < plt.xlim()[1] and plot_text:
        plt.text(x=arc_center[0],
                 y=arc_center[1]+arc_height/2+30,
                 horizontalalignment='center',
                 verticalalignment='center',
                 s=str(int(n_lef)),
                 fontsize=20,
                 color=color
                )
        
def plot_BEs(
    l, r, xrange, yrange,ax,
    color = 'b'):
    """Visualize BEs.
    """
    # Draw BEs
    left_x_tail = l-xrange*0.08*0.5
    left_y_tail = yrange*0.04
    left_x_head = l+xrange*0.08*0.5
    left_y_head = yrange*0.04
    left_BE1 = FancyArrowPatch((left_x_tail, left_y_tail), (left_x_head, left_y_head),
                                 mutation_scale=9,color=color,zorder=2500)
    ax.add_patch(left_BE1)
    
    left_x_tail = l+xrange*0.08*0.5
    left_y_tail = yrange*0.06
    left_x_head = l-xrange*0.08*0.5
    left_y_head = yrange*0.06
    left_BE2 = FancyArrowPatch((left_x_tail, left_y_tail), (left_x_head, left_y_head),
                                 mutation_scale=9,color=color,zorder=2600)
    ax.add_patch(left_BE2)
    
    bb = mtransforms.Bbox([[l-xrange*0.08*0.5, yrange*0.03], [l+xrange*0.08*0.5,yrange*0.07]])
    left_BEbox = FancyBboxPatch((bb.xmin, bb.ymin),
                             abs(bb.width), abs(bb.height),
                             boxstyle="round,pad=0.1",
                             fc="none",
                             ec='black',
                             lw=2,
                             zorder = 2650,
                             mutation_aspect=yrange/xrange)
    left_BEbox.set_boxstyle("round,pad="+str(xrange*0.08*0.1)+", rounding_size="+str(xrange*0.08*0.3))
    ax.add_patch(left_BEbox)
    
    right_x_tail = r+xrange*0.08*0.5
    right_y_tail = yrange*0.04
    right_x_head = r-xrange*0.08*0.5
    right_y_head = yrange*0.04
    right_BE1 = FancyArrowPatch((right_x_tail, right_y_tail), (right_x_head, right_y_head),
                                 mutation_scale=9,color=color,zorder=2700)
    ax.add_patch(right_BE1)
    
    right_x_tail = r-xrange*0.08*0.5
    right_y_tail = yrange*0.06
    right_x_head = r+xrange*0.08*0.5
    right_y_head = yrange*0.06
    right_BE2 = FancyArrowPatch((right_x_tail, right_y_tail), (right_x_head, right_y_head),
                                 mutation_scale=9,color=color,zorder=2800)
    ax.add_patch(right_BE2)
    
    bb = mtransforms.Bbox([[r-xrange*0.08*0.5, yrange*0.03], [r+xrange*0.08*0.5, yrange*0.07]])
    right_BEbox = FancyBboxPatch((bb.xmin, bb.ymin),
                             abs(bb.width), abs(bb.height),
                             boxstyle="round,pad=0.1",
                             fc="none",
                             ec='black',
                             lw=2,
                             zorder = 2850,
                             mutation_aspect=yrange/xrange)
    right_BEbox.set_boxstyle("round,pad="+str(xrange*0.08*0.1)+", rounding_size="+str(xrange*0.08*0.3))
    ax.add_patch(right_BEbox)
    
def plot_DSBs(xrange,
    yrange,ax,
    color = 'r'):
    """Visualize BEs.
    """
    # Draw DSBs
    left_DSB = Ellipse((-xrange*0.025,yrange*0.02), 2*xrange*0.02,2*yrange*0.02, color = color,zorder=2900)
    ax.add_patch(left_DSB)
    
    right_DSB = Ellipse((xrange*0.025,yrange*0.02), 2*xrange*0.02,2*yrange*0.02, color = color,zorder=3000)
    ax.add_patch(right_DSB)
  
        
# internal helper function to help with plotting
def _gridspec_inches(
    wcols,
    hrows,
    wspace=0,
    hspace=0.35,
    fig_kwargs={}):
    
    fig = plt.figure()
    fig_height_inches = (
        sum(hrows)
        )

    fig_width_inches = (
        sum(wcols)
        )

    fig=plt.figure(
        figsize=(fig_width_inches,fig_height_inches),
        subplotpars=matplotlib.figure.SubplotParams(
        left=0,
        right=1,
        bottom=0,
        top=1,
        wspace =0.0,
        hspace = 0.0),
        **fig_kwargs)
    
    fig.set_size_inches(fig_width_inches,fig_height_inches,forward=True)

    
    gs = matplotlib.gridspec.GridSpec(
        len(hrows),
        len(wcols),
        left=0.1,
        right=1,
        top=1,
        bottom=0,
        wspace=wspace,
        hspace=hspace,
        width_ratios=wcols,
        height_ratios=hrows
        )
    return fig, gs

def add_colorbar(gs,count,cmap,max_time_colour=15,time_step=2,label='Extrusion time (sec)'):
    # add colorbar
    plt.subplot(gs[count])
    norm = matplotlib.colors.Normalize(vmin=0,vmax=max_time_colour)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar = plt.colorbar(sm, ticks=np.arange(0,max_time_colour,time_step))
    cbar.ax.tick_params(labelsize=15) 
    cbar.set_label(label,size=20)
    plt.axis('off')
    
def plot_LEFs(plotting_folder,processivity,left_extruders,right_extruders,left_boundary,right_boundary,restraining_left,restraining_right,restraining_left_o,restraining_right_o,break_site,title,superLEFs_len,gapsize_left,gapsize_right,gapsize_left_o,gapsize_right_o,steps_per_frame,extrusion_velocity,upper_time_limit,unconstrained,synapsed, preDSB,gapsize_xmin):                                            
    plt.figure()
    fig, gs = _gridspec_inches([5],[5,2,2])
    max_size_colour = processivity
    height_factor = 1
    cmap=matplotlib.cm.get_cmap('gnuplot2_r')
    ax=plt.subplot(gs[0])
    # the distance between the leftmost BE/constraining LEF subunit and the rightmost BE/constraining LEF subunit 
    xrange = max(restraining_right_o,right_boundary)- min(restraining_left_o,left_boundary)
    xlim_low = min(restraining_left_o,left_boundary)-xrange*0.2
    xlim_high = max(restraining_right_o,right_boundary) + xrange*0.2
    xlim_low_LEF = xlim_low-(xlim_high-xlim_low)*1.5
    xlim_high_LEF = xlim_high+(xlim_high-xlim_low)*1.5
    arch_height_limit = 100
    #find LEFs between the xlim_low and xlim_high
    test_left = np.where((left_extruders-xlim_low_LEF)>=0,1,0) 
    test_right = np.where((right_extruders-xlim_high_LEF)<=0,1,0)
    test_binary = np.multiply(test_left,test_right)
    LEF_index = np.nonzero(test_binary)
    left_LEF = left_extruders[LEF_index]
    right_LEF = right_extruders[LEF_index]
    left_SuperLEF = left_extruders[:superLEFs_len]
    yrange = arch_height_limit/0.6
    
    # small scaling make sure the LEF motors don't move past DSB end visually
    delta = np.sqrt((yrange*(0.015+0.02))**2-(yrange*0.02*1.5+yrange*0.015-yrange*0.02)**2)
    scaling_left = (-xrange*0.025-(xlim_low-break_site)-delta)/(0-(xlim_low-break_site))
    scaling_right = ((xlim_high-break_site)-xrange*0.025-delta)/(xlim_high-break_site-1)
    
    if type(xrange)==np.ndarray or type(xrange)==list:
        xrange=xrange[0]
    if type(yrange)==np.ndarray or type(yrange)==list:
        yrange=yrange[0]
    if type(xlim_low)==np.ndarray or type(xlim_low)==list:
        xlim_low=xlim_low[0]
    if type(xlim_high)==np.ndarray or type(xlim_high)==list:
        xlim_high=xlim_high[0]
    if type(left_boundary)==np.ndarray or type(left_boundary)==list:
        left_boundary=left_boundary[0]
    if type(right_boundary)==np.ndarray or type(right_boundary)==list:
        right_boundary=right_boundary[0]
    if type(break_site)==np.ndarray or type(break_site)==list:
        break_site=break_site[0]
    if type(scaling_left)==np.ndarray or type(scaling_left)==list:
        scaling_left =scaling_left[0]
    if type(scaling_right)==np.ndarray or type(scaling_right)==list:
        scaling_right =scaling_right[0]
    
    for l, r in zip(left_LEF,right_LEF):
        if l != r:
            line_style = '-'
            line_width=2
            if l == restraining_left and r == restraining_right:
                if l in left_SuperLEF:
                    line_style = ':'
                if preDSB:
                    plot_arch_diagram(ax,xrange,yrange,xlim_low,xlim_high,break_site,scaling_left,scaling_right,l-break_site,r-break_site,n_lef=1,y_offset=yrange*0.02*0.75,height_factor=height_factor,\
                                      color='orange',lw=line_width,max_height=arch_height_limit,ls = line_style,zorder=999)
                else:
                    plot_arch_diagram(ax,xrange,yrange,xlim_low,xlim_high,break_site,scaling_left,scaling_right,l-break_site,r-break_site,n_lef=1,y_offset=yrange*0.02*0.75,height_factor=height_factor,\
                                      color='green',lw=line_width,max_height=arch_height_limit,ls = line_style,zorder=999)
            else:
                if l in left_SuperLEF:
                    line_style = ':'
                plot_arch_diagram(ax,xrange,yrange,xlim_low,xlim_high,break_site,scaling_left,scaling_right,l-break_site,r-break_site,n_lef=1,y_offset=yrange*0.02*0.75,height_factor=height_factor,\
                  color='orange',lw=line_width, max_height=arch_height_limit,ls = line_style)  

    stepsize = np.maximum(int(np.ceil((xlim_high-xlim_low)/60)),1)*10
    positive_ticks = np.arange(0,xlim_high-break_site,stepsize)
    negative_ticks = np.arange(0,xlim_low-break_site,-stepsize)
    x_ticks=np.concatenate([np.flip(negative_ticks)[:-1],positive_ticks])
    xlims = [xlim_low-break_site,xlim_high-break_site]
    ylims = [0,yrange]
    
    #Draw DNA
    if synapsed or preDSB:
        DNA = Rectangle((xlim_low-break_site,yrange*0.02*0.5),xlim_high-xlim_low,yrange*0.02,color='dodgerblue',zorder=2000)
        ax.add_patch(DNA)
    else:
        DNAleft = Rectangle((xlim_low-break_site,yrange*0.02*0.5),-xrange*0.025-(xlim_low-break_site)-xrange*0.001,yrange*0.02,color='dodgerblue',zorder=2000)
        ax.add_patch(DNAleft)
        DNAright = Rectangle((xrange*0.026,yrange*0.02*0.5),xlim_high-break_site-xrange*0.026,yrange*0.02,color='dodgerblue',zorder=2050)
        ax.add_patch(DNAright)
        
        plot_DSBs(xrange,yrange,ax,color='purple') 
    
    BE1 = left_boundary-break_site
    BE2 = right_boundary-break_site
    #apply scaling to BEs

    BE1=(BE1-(xlim_low-break_site))*scaling_left+(xlim_low-break_site)
    BE2=xlim_high-break_site-(xlim_high-break_site-BE2)*scaling_right
        
    plot_BEs(BE1,BE2,xrange,yrange,ax,color='red')

    #Make legends
    if superLEFs_len==0:
        DNAlegend_y = yrange*(0.81-0.02)
        DSBlegend_y = yrange*(0.74-0.02)
        BElegend_y = yrange*(0.67-0.02)
    else:
        DNAlegend_y = yrange*(0.81-0.02)-yrange*0.07*2
        DSBlegend_y = yrange*(0.74-0.02)-yrange*0.07*2
        BElegend_y = yrange*(0.67-0.02)-yrange*0.07*2
        
        # long-lived constrainingLEF legend
        constrainingLEFlegend = Arc(xy=(xlim_low-break_site+xrange*0.08,yrange*0.78),
            width=xrange*0.1,
            height=yrange*0.1,
            theta1=0,
            theta2=180,
            alpha=1,
            lw=2,
            color='green',
            capstyle='round',
            ls = ':')
        ax.add_patch(constrainingLEFlegend)
        left_motor = Ellipse((xlim_low-break_site+xrange*0.03,yrange*0.78), 2*xrange*0.015,2*yrange*0.015, color = 'green')
        ax.add_patch(left_motor)
        right_motor = Ellipse((xlim_low-break_site+xrange*0.13,yrange*0.78), 2*xrange*0.015,2*yrange*0.015, color = 'green')
        ax.add_patch(right_motor)
        ax.text(xlim_low-break_site+xrange*0.17,yrange*0.79,'Long-lived constraining LEF',fontsize='x-large')

        #gapbridgingLEF legend
        gapbridgingLEFlegend = Arc(xy=(xlim_low-break_site+xrange*0.08,yrange*0.71),
                width=xrange*0.1,
                height=yrange*0.1,
                theta1=0,
                theta2=180,
                alpha=1,
                lw=2,
                color='orange',
                capstyle='round',
                ls = ':')
        ax.add_patch(gapbridgingLEFlegend)
        ax.add_patch(constrainingLEFlegend)
        left_motor = Ellipse((xlim_low-break_site+xrange*0.03,yrange*0.71), 2*xrange*0.015,2*yrange*0.015, color = 'orange')
        ax.add_patch(left_motor)
        right_motor = Ellipse((xlim_low-break_site+xrange*0.13,yrange*0.71), 2*xrange*0.015,2*yrange*0.015, color = 'orange')
        ax.add_patch(right_motor)
        ax.text(xlim_low-break_site+xrange*0.17,yrange*0.72,'Long-lived gap-bridging and other LEFs',fontsize='x-large')
        
    #DNA legend
    DNAlegend = Rectangle((xlim_low-break_site+xrange*0.03,DNAlegend_y+yrange*0.01),xrange*0.1,yrange*0.02,color='dodgerblue')
    ax.add_patch(DNAlegend)
    ax.text(xlim_low-break_site+xrange*0.17,DNAlegend_y,'Double-stranded DNA',fontsize='x-large')
    #DSB legend
    DSBlegend = Ellipse((xlim_low-break_site+xrange*0.11,DSBlegend_y+yrange*0.02), 2*xrange*0.02,2*yrange*0.02, color = 'purple')
    ax.add_patch(DSBlegend)
    ax.text(xlim_low-break_site+xrange*0.17,DSBlegend_y,'DSB end',fontsize='x-large')
    #BE legend
    left_x_tail = xlim_low-break_site+xrange*0.13
    left_y_tail = BElegend_y+yrange*0.03
    left_x_head = xlim_low-break_site+xrange*0.13-xrange*0.08
    left_y_head = BElegend_y+yrange*0.03

    BElegend1 = FancyArrowPatch((left_x_tail, left_y_tail), (left_x_head, left_y_head),
                                 mutation_scale=9,color='red')
    ax.add_patch(BElegend1)
    
    left_x_tail = xlim_low-break_site+xrange*0.13-xrange*0.08
    left_y_tail = BElegend_y+yrange*0.01
    left_x_head = xlim_low-break_site+xrange*0.13
    left_y_head = BElegend_y+yrange*0.01

    BElegend2 = FancyArrowPatch((left_x_tail, left_y_tail), (left_x_head, left_y_head),
                                 mutation_scale=9,color='red')
    ax.add_patch(BElegend2)
    bb = mtransforms.Bbox([[xlim_low-break_site+xrange*0.13-xrange*0.08, BElegend_y], [xlim_low-break_site+xrange*0.13,BElegend_y+yrange*0.04]])
    BEbox = FancyBboxPatch((bb.xmin, bb.ymin),
                             abs(bb.width), abs(bb.height),
                             boxstyle="round,pad=0.1",
                             fc="none",
                             ec='black',
                             lw=2,
                             mutation_aspect=yrange/xrange)
    BEbox.set_boxstyle("round,pad="+str(xrange*0.08*0.1)+", rounding_size="+str(xrange*0.08*0.3))
    ax.add_patch(BEbox)
    ax.text(xlim_low-break_site+xrange*0.17,BElegend_y,'BE',fontsize='x-large')
    
    #constrainingLEF legend
    constrainingLEFlegend = Arc(xy=(xlim_low-break_site+xrange*0.08,yrange*0.92),
            width=xrange*0.1,
            height=yrange*0.1,
            theta1=0,
            theta2=180,
            alpha=1,
            lw=2,
            color='green',
            capstyle='round',
            ls = '-')
    ax.add_patch(constrainingLEFlegend)
    left_motor = Ellipse((xlim_low-break_site+xrange*0.03,yrange*0.92), 2*xrange*0.015,2*yrange*0.015, color = 'green')
    ax.add_patch(left_motor)
    right_motor = Ellipse((xlim_low-break_site+xrange*0.13,yrange*0.92), 2*xrange*0.015,2*yrange*0.015, color = 'green')
    ax.add_patch(right_motor)
    ax.text(xlim_low-break_site+xrange*0.17,yrange*(0.95-0.02),'Constraining LEF',fontsize='x-large')
    
    #gapbridgingLEF legend
    gapbridgingLEFlegend = Arc(xy=(xlim_low-break_site+xrange*0.08,yrange*0.85),
            width=xrange*0.1,
            height=yrange*0.1,
            theta1=0,
            theta2=180,
            alpha=1,
            lw=2,
            color='orange',
            capstyle='round',
            ls = '-')
    ax.add_patch(gapbridgingLEFlegend)
    ax.add_patch(constrainingLEFlegend)
    left_motor = Ellipse((xlim_low-break_site+xrange*0.03,yrange*0.85), 2*xrange*0.015,2*yrange*0.015, color = 'orange')
    ax.add_patch(left_motor)
    right_motor = Ellipse((xlim_low-break_site+xrange*0.13,yrange*0.85), 2*xrange*0.015,2*yrange*0.015, color = 'orange')
    ax.add_patch(right_motor)
    ax.text(xlim_low-break_site+xrange*0.17,yrange*(0.88-0.02),'Gap-bridging and other LEFs',fontsize='x-large')

    
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    plt.xticks(x_ticks,fontsize = 15)
    axis = plt.gca()
    axis.axes.yaxis.set_visible(False)
    ax.set_title(title,fontsize='xx-large')
    ax.set_xlabel('Distance from DSB (kb)',fontsize='xx-large')
    
    #plot gap size
    ax2=plt.subplot(gs[1])
    x = np.arange(0,len(gapsize_left))*steps_per_frame*2/extrusion_velocity/60+gapsize_xmin

    if not preDSB:
        if unconstrained:
            ax2.text(upper_time_limit*0.1+gapsize_xmin,np.maximum(gapsize_left_o*2,1)*2/3,'No constraining LEF remains, \n two DSB ends may diffuse apart',fontsize='x-large')
        else:
            ax2.scatter(x,gapsize_left, s=15, c='b',marker='d')

    if synapsed:
        ax2.text(upper_time_limit*0.3+gapsize_xmin,np.maximum(gapsize_left_o*2,1)*2/3,'Synapsis achieved',fontsize='x-large')
        
    ax2.set_xlim([gapsize_xmin,gapsize_xmin+upper_time_limit])
    ax2.set_ylim([0,np.maximum(gapsize_left_o*2,1)])
    ax2.set_xlabel('Time elapsed after DSB (min)',fontsize='xx-large')
    ax2.set_ylabel('Left gap size (kb)',fontsize='xx-large')
    ax2.tick_params(axis='both', which='major', labelsize=15)
    
    ax3=plt.subplot(gs[2])
    x = np.arange(0,len(gapsize_right))*steps_per_frame*2/extrusion_velocity/60+gapsize_xmin

    if not preDSB:
        if unconstrained:
            ax3.text(upper_time_limit*0.1+gapsize_xmin,np.maximum(gapsize_right_o*2,1)*2/3,'No constraining LEF remains, \n two DSB ends may diffuse apart',fontsize='x-large')
        else:
            ax3.scatter(x,gapsize_right, s=15, c='b',marker='d')
                
    if synapsed:
        ax3.text(upper_time_limit*0.3+gapsize_xmin,np.maximum(gapsize_right_o*2,1)*2/3,'Synapsis achieved',fontsize='x-large')
    
    ax3.set_xlim([gapsize_xmin,gapsize_xmin+upper_time_limit])
    ax3.set_ylim([0,np.maximum(gapsize_right_o*2,1)])
    ax3.set_xlabel('Time elapsed after DSB (min)',fontsize='xx-large')
    ax3.set_ylabel('Right gap size (kb)',fontsize='xx-large')
    ax3.tick_params(axis='both', which='major', labelsize=15)
    
    plt.savefig(plotting_folder,bbox_inches='tight',format='jpeg')
    plt.close('all')
    
    
def plot_residence(output_folder,residence_list,separation, processivity, boundary_strength, base_stochasticity):                                            
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    plt.figure()
    residence_time = [x*base_stochasticity/60 for x in residence_list]
    x = pd.Series(residence_time, name="residence time(min)")
    residence_array = np.array(residence_time)
    u = np.around(np.mean(residence_array),1)
    s = np.around(np.std(residence_array),1)
    plot_label = 'n = ' + str(len(residence_time)) + ', u = ' + str(u) + ', s = ' + str(s)
    ax = sns.distplot(x,kde_kws={"color": "k", "lw": 3, "label": plot_label})
    ax.set_title('separation = '+ str(separations)+' kb, processivity = ' + str(processivity) + ' kb, stall_prob = '+ str(round(boundary_strength,2)))
    plot_path = output_folder +'/separations '+ str(separations)+'_processivity ' + str(processivity) + '_stall_prob '+ str(round(boundary_strength,2)) + '_residence dist_.png'
    plt.savefig(plot_path,bbox_inches='tight')
    plt.close('all')