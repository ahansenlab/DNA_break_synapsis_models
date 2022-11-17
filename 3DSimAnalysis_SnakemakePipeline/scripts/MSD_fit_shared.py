import json
import numpy as np
from scipy import special, stats
from matplotlib import pyplot as plt

import tracklib as tl
from tracklib.analysis import bild, msdfit
from tracklib.analysis.msdfit.lib import TwoLocusRouseFit 


def run_fit(data,fit_with_localization_error=True):
    fit = TwoLocusRouseFit(data)
    
    if fit_with_localization_error==False:
        fit.fix_values += [(0, 0), (3, 0), (6, 0)] # fixes the "sig_x, sig_y, sig_z" values to zero
    
    profiler = msdfit.Profiler(
        fit,
        verbosity=2,
        max_fit_runs = np.inf,
        profiling = False,
        bracket_strategy = {
            'multiplicative' : False,
            'step' : 0.01,
            'nonidentifiable_cutoffs' : [3, np.inf]
        },
    )
    profiler.restart_on_better_point_estimate = False

    res = profiler.fit.run(show_progress=False,
                           optimization_steps=(dict(method='Nelder-Mead', 
                                                    options={'fatol' : 0.001, 'xatol' : 0.01}),),
                          )
    profiler.point_estimate = res

    mci = profiler.find_MCI(show_progress=False)
    
    return mci, profiler

def plot_fit_results(data, mci, profiler, filename, title=None,fit_with_localization_error=True):
    fig, axs = plt.subplots(1, 2, figsize=[13, 4])

    ax = axs[0]

    msd = tl.analysis.MSD(data)
    dt = np.arange(1, len(msd))
    ax.plot(dt, msd[1:], color='k',
            zorder=5, linewidth=2,
            label='empirical',
           )

    fixer = profiler.fit.get_value_fixer()
    params = profiler.point_estimate['params']
    msdm = profiler.fit.params2msdm(params)
    ax.plot(dt, np.sum([msd(dt) for msd, m in msdm], axis=0),
            label='MLE', color='red',
           )

    ax.legend()
    ax.set_xscale('log')
    ax.set_xlabel('time')
    ax.set_yscale('log')
    ax.set_ylabel('MSD')
    if title is not None:
        ax.set_title(title)

    ax = axs[1]
    if fit_with_localization_error == True:
        for iparam, label in [(0, 'σ_x'), (3, 'σ_y'), (6, 'σ_z'), (1, 'Γ'), (2, 'J'),]:
            logvals = np.array([res['params'][iparam] for res in profiler.ress[iparam]+[profiler.point_estimate]])
            logLs =   np.array([res['logL']           for res in profiler.ress[iparam]+[profiler.point_estimate]])
            ind = np.argsort(logvals)
            ax.plot(logvals[ind] - mci[iparam, 0], logLs[ind],
                    marker='x',
                    label=label)

    else:

        for iparam, label in [(1, 'Γ'), (2, 'J'),]:
            logvals = np.array([res['params'][iparam] for res in profiler.ress[iparam]+[profiler.point_estimate]])
            logLs =   np.array([res['logL']           for res in profiler.ress[iparam]+[profiler.point_estimate]])
            ind = np.argsort(logvals)
            ax.plot(logvals[ind] - mci[iparam, 0], logLs[ind],
                    marker='x',
                    label=label)
        
    ax.legend(loc=(1.02, 0.5))
    ax.set_xlabel('log-parameters, normalized to point estimate')
    ax.set_ylabel('log(L)')
    ax.set_title('Conditional posterior')
    plt.draw()
    plt.savefig(filename,bbox_inches='tight')
    plt.close()

def get_named_params(params,fit_with_localization_error=True):
    out = {}
    out['Gs_1d'] = list(np.exp(params[[1, 4, 7]])) # diffusivity "gamma", per dimension
    out['Js_1d'] = list(np.exp(params[[2, 5, 8]])) # steady state value, per dimension    
    out['G_3d'] = np.sum(out['Gs_1d']) # diffusivity "gamma", total
    out['J_3d'] = np.sum(out['Js_1d']) # steady state value, total    
    
    if fit_with_localization_error==True:
        out['sig2'] = list(np.exp(params[[0, 3, 6]])) # localization error variance, per dimension
        out['sig'] = list(np.sqrt(out['sig2'])) # localization error variance, total
        
    return out