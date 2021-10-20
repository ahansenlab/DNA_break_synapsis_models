import matplotlib.pyplot as plt
import numpy as np
import os

from multiprocessing import Pool
import pickle

from scipy.sparse.csgraph import shortest_path
from scipy.sparse import csr_matrix
#from brandaolib.plotting import zoomArray

def contact_map_generator_from_SMC_list_and_DSBs(smc_pairs,dsb_sites,N=4000,base=10,nsamples=80):
    """
    Parameters
    ----------
    smc_pairs : list of tuples
        List of 2-tuples containing locations of SMCs on the chain of length `N`

    N : int
        Length of the polymer chain.

    base : int
        Division value to convert chain position to index position on `HiC` heatmap

    Returns
    -------

    HiC : ndarray
        Contact proabilities from sampled points of smc_pairs conformation

    HiC_count : ndarray
        Number of sampled points per HiC map position

    """
    row = [x[0] for x in smc_pairs]
    col = [x[1] for x in smc_pairs]
    
    for st, nd in zip(np.r_[0,dsb_sites],np.r_[dsb_sites,N-1]):
        if nd-st <= 1:
            continue
        row += [x for x in np.arange(st,nd) ]
        col += [np.mod(x+1,N) for x in np.arange(st,nd)] 
    
    #row = [x for x in range(N) ] + [x[0] for x in smc_pairs]
    #col = [np.mod(x+1,N) for x in range(N)] + [x[1] for x in smc_pairs]
    deff = shortest_path(csr_matrix((np.ones(len(row)),(row,col)),shape=(N,N)),directed=False)


    # create heatmap
    HiC = np.zeros((N//base,N//base))
    HiC_count = np.zeros((N//base,N//base))

    vals = sorted(np.random.choice(N,nsamples, replace=False))
    for ix in range(len(vals)):
        for iy in range(ix,len(vals)):
            x = vals[ix]
            y = vals[iy]
            if deff[x,y] == 0:
                pc = 1
            else:
                pc = 1/np.sqrt(deff[x,y])**3

            if not np.isnan(pc):
                HiC_count[x//base,y//base] += 1
                HiC_count[y//base,x//base] += 1
                HiC[x//base,y//base] += pc
                HiC[y//base,x//base] += pc
    return HiC, HiC_count



def contact_map_generator_from_SMC_list(smc_pairs,N=4000,base=10,nsamples=80):
    """
    Parameters
    ----------
    smc_pairs : list of tuples
        List of 2-tuples containing locations of SMCs on the chain of length `N`

    N : int
        Length of the polymer chain.

    base : int
        Division value to convert chain position to index position on `HiC` heatmap

    Returns
    -------

    HiC : ndarray
        Contact proabilities from sampled points of smc_pairs conformation

    """

    row = [x for x in range(N) ] + [x[0] for x in smc_pairs]
    col = [np.mod(x+1,N) for x in range(N)] + [x[1] for x in smc_pairs]
    deff = shortest_path(csr_matrix((np.ones(len(row)),(row,col)),shape=(N,N)),directed=False)
    deff = zoomArray(deff,(N//base,N//base))

    HiC = 1/np.sqrt(deff+1)**3

    return HiC

def contact_map_generator_from_SMC_list_slow(smc_pairs,N=4000,base=10,nsamples=80):
    """
    Parameters
    ----------
    smc_pairs : list of tuples
        List of 2-tuples containing locations of SMCs on the chain of length `N`

    N : int
        Length of the polymer chain.

    base : int
        Division value to convert chain position to index position on `HiC` heatmap

    Returns
    -------

    HiC : ndarray
        Contact proabilities from sampled points of smc_pairs conformation

    HiC_count : ndarray
        Number of sampled points per HiC map position

    """

    row = [x for x in range(N) ] + [x[0] for x in smc_pairs]
    col = [np.mod(x+1,N) for x in range(N)] + [x[1] for x in smc_pairs]
    deff = shortest_path(csr_matrix((np.ones(len(row)),(row,col)),shape=(N,N)),directed=False)


    # create heatmap
    HiC = np.zeros((N//base,N//base))
    HiC_count = np.zeros((N//base,N//base))

    vals = sorted(np.random.choice(N,nsamples, replace=False))
    for ix in range(len(vals)):
        for iy in range(ix,len(vals)):
            x = vals[ix]
            y = vals[iy]
            if deff[x,y] == 0:
                pc = 1
            else:
                pc = 1/np.sqrt(deff[x,y])**3

            if not np.isnan(pc):
                HiC_count[x//base,y//base] += 1
                HiC_count[y//base,x//base] += 1
                HiC[x//base,y//base] += pc
                HiC[y//base,x//base] += pc
    return HiC, HiC_count


class heatmap():
    def __init__(self,HiC,HiC_counts,nsamples=80):
        self.HiC = HiC
        self.HiC_count = HiC_counts
        self.numFailed = 0
        self.totsamp = 0
        self.nsamples = nsamples

    def add(self,hmap):
        self.HiC += hmap.HiC
        self.HiC_count += hmap.HiC_count
        self.numFailed += hmap.numFailed
        self.totsamp += hmap.totsamp
