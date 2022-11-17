from __future__ import absolute_import, division, print_function
import numpy as np
import sys
import os
import time
import tempfile
import logging
import warnings

import pickle
import os
import time
import numpy as np
import polychrom

from cooltools.lib.numutils import coarsen, logbins
from functools import reduce

from polychrom import polymerutils
from polychrom import forces
from polychrom import forcekits
from polychrom.simulation import Simulation
from polychrom.starting_conformations import grow_cubic
from polychrom.hdf5_format import HDF5Reporter, list_URIs, load_URI, load_hdf5_file

import simtk.openmm 
import os 
import shutil

import pyximport; pyximport.install()
pyximport.install(setup_args={'include_dirs': np.get_include()})
from DSB_smcTranslocator_v2 import smcTranslocatorDirectional 

import warnings
import h5py 
import glob
import re

from itertools import product
from scipy.ndimage.filters import gaussian_filter
from scipy.sparse import coo_matrix
from scipy.ndimage.filters import gaussian_filter1d
from scipy.stats import expon

from datetime import datetime

from polychrom.hdf5_format import HDF5Reporter, list_URIs, load_URI, load_hdf5_file
import h5py as hp
from matplotlib import pyplot as plt



from collections.abc import Iterable

import simtk.openmm as openmm
import simtk.unit


logging.basicConfig(level=logging.INFO)

class DSB_Simulation(Simulation):
    def do_block(
        self,
        steps=None,
        check_functions=[],
        get_velocities=False,
        save=True,
        save_extras={},
        positions_to_sample=None,
    ):
        """performs one block of simulations, doing steps timesteps,
        or steps_per_block if not specified.

        Parameters
        ----------

        steps : int or None
            Number of timesteps to perform.
        increment : bool, optional
            If true, will not increment self.block and self.steps counters
        """

        if not self.forces_applied:
            if self.verbose:
                logging.info("applying forces")
                sys.stdout.flush()
            self._apply_forces()
            self.forces_applied = True

        a = time.time()
        self.integrator.step(steps)  # integrate!

        self.state = self.context.getState(
            getPositions=True, getVelocities=get_velocities, getEnergy=True
        )

        b = time.time()
        coords = self.state.getPositions(asNumpy=True)
        newcoords = coords / simtk.unit.nanometer
        newcoords = np.array(newcoords, dtype=np.float32)
        if self.kwargs["save_decimals"] is not False:
            newcoords = np.round(newcoords, self.kwargs["save_decimals"])

        self.time = self.state.getTime() / simtk.unit.picosecond

        # calculate energies in KT/particle
        eK = self.state.getKineticEnergy() / self.N / self.kT
        eP = self.state.getPotentialEnergy() / self.N / self.kT
        curtime = self.state.getTime() / simtk.unit.picosecond

        msg = "block %4s " % int(self.block)
        msg += "pos[1]=[%.1lf %.1lf %.1lf] " % tuple(newcoords[0])

        check_fail = False
        for check_function in check_functions:
            if not check_function(newcoords):
                check_fail = True

        if np.isnan(newcoords).any():
            raise IntegrationFailError("Coordinates are NANs")
        if eK > self.eK_critical:
            print(eK)
            print(self.eK_critical)
            raise EKExceedsError("Ek={1} exceeds {0}".format(self.eK_critical, eK))
        if (np.isnan(eK)) or (np.isnan(eP)):
            raise IntegrationFailError("Energy is NAN)")
        if check_fail:
            raise IntegrationFailError("Custom checks failed")

        dif = np.sqrt(np.mean(np.sum((newcoords - self.get_data()) ** 2, axis=1)))
        msg += "dr=%.2lf " % (dif,)
        self.data = coords
        msg += "t=%2.1lfps " % (self.state.getTime() / simtk.unit.picosecond)
        msg += "kin=%.2lf pot=%.2lf " % (eK, eP)
        msg += "Rg=%.3lf " % self.RG()
        msg += "SPS=%.0lf " % (steps / (float(b - a)))

        if (
            self.integrator_type.lower() == "variablelangevin"
            or self.integrator_type.lower() == "variableverlet"
        ):
            dt = self.integrator.getStepSize()
            msg += "dt=%.1lffs " % (dt / simtk.unit.femtosecond)
            mass = self.system.getParticleMass(1)
            dx = simtk.unit.sqrt(2.0 * eK * self.kT / mass) * dt
            msg += "dx=%.2lfpm " % (dx / simtk.unit.nanometer * 1000.0)

        logging.info(msg)

        if positions_to_sample is None:
            result = {
                "pos": newcoords,
                "potentialEnergy": eP,
                "kineticEnergy": eK,
                "time": curtime,
                "block": self.block,
            }
        else: 
            
            result = {
                "pos": [newcoords[pos] for pos in positions_to_sample],
                "potentialEnergy": eP,
                "kineticEnergy": eK,
                "time": curtime,
                "block": self.block,
            }
    
        if get_velocities:
            result["vel"] = self.state.getVelocities() / (
                simtk.unit.nanometer / simtk.unit.picosecond
            )
        result.update(save_extras)
        if save:
            for reporter in self.reporters:
                reporter.report("data", result)

        self.block += 1
        self.step += steps

        return result    
    
def run_simulation(idx, N_monomers, translocator_initialization_steps, \
                  smcStepsPerBlock, \
                  smcBondDist, smcBondWiggleDist,\
                  steps_per_block_beforeDSB, steps_per_block, volume_density, block, \
                  chain, extra_bonds, save_folder, save_every_x_blocks, \
                  total_saved_blocks, restartBondUpdaterEveryBlocks,BondUpdaterCountBeforeDSB, GPU_choice = 0, 
                   overwrite=False,density=0.2,positions_to_sample=None,colrate=0.3,errtol=0.01,trunc=0.001,
                  initial_conformation=None,block_to_save_all = None,save_length=10000):
    
    print(f"Starting simulation. Doing {total_saved_blocks} steps, saving every {save_every_x_blocks} block.")    
    print(f"A total of {(total_saved_blocks-(BondUpdaterCountBeforeDSB*restartBondUpdaterEveryBlocks)) * steps_per_block+(BondUpdaterCountBeforeDSB*restartBondUpdaterEveryBlocks)*steps_per_block_beforeDSB} steps will be performed")
    
    # assertions for easy managing code below 
    assert restartBondUpdaterEveryBlocks % save_every_x_blocks == 0 
    assert (total_saved_blocks * save_every_x_blocks) % restartBondUpdaterEveryBlocks == 0 
    savesPerBondUpdater = restartBondUpdaterEveryBlocks // save_every_x_blocks
    BondUpdaterInitsTotal  = (total_saved_blocks ) * save_every_x_blocks // restartBondUpdaterEveryBlocks
    print("BondUpdater will be initialized {0} times".format(BondUpdaterInitsTotal))    
    
    
    # create SMC translocation object   
    SMCTran = initModel(idx)                        
    SMCTran.steps(translocator_initialization_steps)  # steps to "equilibrate" SMC dynamics

    # now feed bond generators to BondUpdater 
    BondUpdater = simulationBondUpdater(SMCTran)   
       
    # clean up the simulation directory
    folder = save_folder
    
    # set initial configuration of the polymer    
    if initial_conformation is None:
        box = (N_monomers / volume_density) ** 0.33  
        data = grow_cubic(N_monomers, int(box) - 2, method='linear')  # creates a compact conformation, but one can take pre-existing configurations
    else:
        data = initial_conformation
        
    # create the reporter class - this saves polymer configurations and other information you specify
    reporter = HDF5Reporter(folder=folder, max_data_length=save_length,overwrite=True,blocks_only=True)
    
    # Iterate over various BondUpdaterInitializations
    for BondUpdaterCount in range(BondUpdaterInitsTotal):
        
        if BondUpdaterCount >= BondUpdaterCountBeforeDSB:
            set_extra_bonds = None # induce double strand break
        else:
            set_extra_bonds = extra_bonds
        
        # simulation parameters are defined here     
        sim = DSB_Simulation(
                platform="cuda",
                integrator="variableLangevin", 
                error_tol=errtol, 
                GPU = "{}".format(GPU_choice), 
                collision_rate=colrate, 
                N = len(data),
                reporters=[reporter],
                PBCbox=False,
                precision="mixed")

        sim.set_data(data)  # loads a polymer, puts a center of mass at zero
        sim.add_force(forces.spherical_confinement(sim,density=density,k=1))
        # -----------Adding forces ---------------
        sim.add_force(
            forcekits.polymer_chains(
                sim,
                chains=chain, # makes circular chains
                bond_force_func=forces.harmonic_bonds, # adds harmonic bonds
                bond_force_kwargs={
                    'bondLength':1.0,
                    'bondWiggleDistance':0.1, # Bond distance will fluctuate +- 0.1 on average
                    # This is fine because our simulation is "soft" on average (all forces are soft-core)
                    # And this is potentially desirable to make polymer chain softer to avoid jerking 
                    # when a new bond is initialized,
                    'override_checks':True,
                 },

                angle_force_func=forces.angle_force,
                angle_force_kwargs={
                    'k':0.05, # we are making a very flexible polymer, basically not necessary here,
                    'override_checks':True,
                },
                nonbonded_force_func=forces.polynomial_repulsive, # this is the excluded volume potential
                nonbonded_force_kwargs={
                    'trunc':trunc,
                    'radiusMult':1.05, # this is from old code
                },
                except_bonds=True,
                extra_bonds=set_extra_bonds,
                override_checks=True,
            )
        )

        sim.step = block
        kbond = sim.kbondScalingFactor / (smcBondWiggleDist ** 2)
        bondDist = smcBondDist * sim.length_scale
        activeParams = {"length":bondDist,"k":kbond}
        inactiveParams = {"length":bondDist, "k":0}
        BondUpdater.setParams(activeParams, inactiveParams)

        # this step actually puts all bonds in and sets first bonds to be what they should be
        BondUpdater.setup( BondUpdaterCount, BondUpdaterCountBeforeDSB,
                     bondForce=sim.force_dict['harmonic_bonds'],
                    smcStepsPerBlock=smcStepsPerBlock,
                    blocks=restartBondUpdaterEveryBlocks) 
        print("Restarting BondUpdater")
        

        # Minimize energy since we're reconfiguring the random walk into crumples and SMC-dependent loops
        sim.local_energy_minimization() 

        # Iterate over simulation time steps within each BondUpdater 
        for b in range(restartBondUpdaterEveryBlocks):
            # BondUpdater updates bonds at each time step.
            curBonds, pastBonds = BondUpdater.step(sim.context)  
            if (b % save_every_x_blocks == 0):  
                if BondUpdaterCount >= BondUpdaterCountBeforeDSB:
                    # save SMC positions and monomer positions 
                    if b+BondUpdaterCount*restartBondUpdaterEveryBlocks in block_to_save_all:
                        sim.do_block(steps=steps_per_block, save_extras={"SMCs":curBonds,"SMC_step":sim.step},
                                 positions_to_sample=[i for i in range(N_monomers)])
                    else:
                        sim.do_block(steps=steps_per_block, save_extras={"SMCs":curBonds,"SMC_step":sim.step},
                                 positions_to_sample=positions_to_sample)
                else: 
                    sim.do_block(steps=steps_per_block_beforeDSB, save_extras={"SMCs":curBonds,"SMC_step":sim.step},
                             positions_to_sample=positions_to_sample) 
                
            else:
                if BondUpdaterCount >= BondUpdaterCountBeforeDSB:
                    sim.integrator.step(steps_per_block)  # do steps without getting the positions from the GPU (faster)
                else:
                    sim.integrator.step(steps_per_block_beforeDSB)  # do steps without getting the positions from the GPU (faster)

        data = sim.get_data()  # get the polymer positions (i.e. data) 
        block = sim.step       # get the time step
        del sim               # delete the simulation
        time.sleep(0.2)      # wait 200ms to let garbage collector do its magic
    # dump data to  output file!
    reporter.dump_data() 
    done_file = open(os.path.join(folder,'sim_done.txt'),"w+")
    done_file.close()
    
def initModel(idx):

    #unchanging parameters
    BELT_ON=0
    BELT_OFF=1
    switchRate = 0 
    SWITCH_PROB= switchRate # switching rate
    PUSH=0
    PAIRED=0
    SLIDE=1
    SLIDE_PAUSEPROB=0.99 
    loop_prefactor=1.5 
    FULL_LOOP_ENTROPY=1 
    FRACTION_ONESIDED=0

    # Extruder dynamics frequency of sampling parameters 
    numExtruderSteps = 500 # steps taken for each simulation sample
    numInitializationSteps = 10000 # how long we take to equilibrate the simulation

    # Polymer and extruder dynamics parameters #
    processivity = lof[idx][3] # processivity
    separations = lof[idx][4] # separation
    longlived_fraction = lof[idx][8] # fraction of all LEFs that are longlived
    PAUSEPROB= pause_prob_beforeDSB # motor pause probability
    
    normal_sep = separations/(1-longlived_fraction)
    smcNum = int(chrom_size//normal_sep) # number of SMCs loaded
    if longlived_fraction==0:
        longlived_smcNum = 0
    else:
        longlived_sep = separations/longlived_fraction
        longlived_smcNum = int(chrom_size//longlived_sep)
        
    SWITCH =  np.ones(chrom_size,dtype=np.double)*SWITCH_PROB

    birthArray = np.ones(chrom_size)/chrom_size
    deathArray = np.zeros(chrom_size, dtype=np.double) + 1. / (0.5*processivity/(1-PAUSEPROB))  # times 0.5 to account for two-sided extrusion 
    longlived_deathArray = np.zeros(chrom_size, dtype=np.double) + 1. / (0.5*processivity*lof[idx][9]/(1-PAUSEPROB)) 
    deathArray[boundary_coordinates] = 1./ (0.5*processivity/(1-PAUSEPROB))/ lof[idx][7]
    longlived_deathArray[boundary_coordinates] = 1./ (0.5*processivity*lof[idx][9]/(1-PAUSEPROB)) / lof[idx][7]

    boundary_flag = (boundaryStrengthsL>0) + (boundaryStrengthsR>0)


    pauseArray = PAUSEPROB*np.ones(chrom_size, dtype=np.double) 
    
    slidePauseArray = np.zeros(chrom_size, dtype=np.double) + SLIDE_PAUSEPROB
    oneSidedArray = np.zeros(smcNum, dtype=np.int64)
    longlived_oneSidedArray = np.zeros(longlived_smcNum, dtype=np.int64)
    belt_on_array = np.zeros(smcNum, dtype=np.double) + BELT_ON
    belt_off_array = np.zeros(smcNum, dtype=np.double) + BELT_OFF
    spf=slidePauseArray*(1.-(1.-SLIDE_PAUSEPROB)*np.exp(-1.*loop_prefactor))
    spb=slidePauseArray*(1.-(1.-SLIDE_PAUSEPROB)*np.exp(loop_prefactor))
    
    ################### TAD BOUNDARY CONDITION###################
    # stallLeft is from the LEF perspective, boundaryStrength is from the CTCF perspective
    stallLeftArray = boundaryStrengthsR*lof[idx][1] 
    stallRightArray =  boundaryStrengthsL*lof[idx][1]  
    ##################################################################

    transloc = smcTranslocatorDirectional(birthArray, deathArray, longlived_deathArray, stallLeftArray, stallRightArray, pauseArray,
                                         smcNum, longlived_smcNum, oneSidedArray,longlived_oneSidedArray, FRACTION_ONESIDED, slide=SLIDE,
                                         slidepauseForward=spf, slidepauseBackward=spb, switch=SWITCH, pushing=PUSH,
                                        belt_on=belt_on_array, belt_off=belt_off_array,SLIDE_PAUSEPROB=SLIDE_PAUSEPROB) 
    return transloc


class simulationBondUpdater(object):
    """
    This class precomputes simulation bonds for faster dynamic allocation. 
    """

    def __init__(self, smcTransObject):#, plectonemeObject):
        """
        :param smcTransObject: smc translocator object to work with
        """
        self.smcObject = smcTransObject
        self.allBonds = []

    def setParams(self, activeParamDict, inactiveParamDict):
        """
        A method to set parameters for bonds.
        It is a separate method because you may want to have a Simulation object already existing

        :param activeParamDict: a dict (argument:value) of addBond arguments for active bonds
        :param inactiveParamDict:  a dict (argument:value) of addBond arguments for inactive bonds

        """
        self.activeParamDict = activeParamDict
        self.inactiveParamDict = inactiveParamDict


    def setup(self, BondUpdaterCount, BondUpdaterCountBeforeDSB,  bondForce, smcStepsPerBlock, blocks = 100):
        """
        A method that milks smcTranslocator object
        and creates a set of unique bonds, etc.

        :param bondForce: a bondforce object (new after simulation restart!)
        :param blocks: number of blocks to precalculate
        :param smcStepsPerBlock: number of smcTranslocator steps per block
        :return:
        """

        if len(self.allBonds) != 0:
            raise ValueError("Not all bonds were used; {0} sets left".format(len(self.allBonds)))

        self.bondForce = bondForce

        # update smc translocator object
        if BondUpdaterCount == BondUpdaterCountBeforeDSB:
            
            stall_prob_left,stall_prob_right = self.smcObject.get_stall_prob()
            loading_prob = self.smcObject.get_emission()
            # make double strand break site impermeable
            stall_prob_left[DSB_coordinates[1::2]] = 1 # even elements of the DSB coordinates
            stall_prob_right[DSB_coordinates[::2]] = 1 # odd elements of the DSB coordinates
            # boosting by DSB
            processivity = lof[idx][3] 
            PAUSEPROB=pause_prob # motor pause probability
            deathArray = np.zeros(chrom_size, dtype=np.double) + 1. / (0.5*processivity/(1-PAUSEPROB))  # times 0.5 to account for two-sided extrusion 
            longlived_deathArray = np.zeros(chrom_size, dtype=np.double) + 1. / (0.5*processivity*lof[idx][9]/(1-PAUSEPROB)) 
            deathArray[boundary_coordinates] = 1./ (0.5*processivity/(1-PAUSEPROB))/ lof[idx][7]
            longlived_deathArray[boundary_coordinates] = 1./ (0.5*processivity*lof[idx][9]/(1-PAUSEPROB)) / lof[idx][7]
            deathArray[DSB_coordinates] =  1./ (0.5*processivity/(1-PAUSEPROB))/ lof[idx][10] 
            longlived_deathArray[DSB_coordinates] =1./ (0.5*processivity*lof[idx][9]/(1-PAUSEPROB)) / lof[idx][10]

            # DSB flag to ensure correct valency at DSB site: DSB can only stabilize one LEF
            DSB_flag = np.zeros(chrom_size, int)
            DSB_flag[DSB_coordinates] = 1 
            self.smcObject.updateStallprob(stall_prob_left,stall_prob_right)
            self.smcObject.updateDSBFlag(DSB_flag)
            self.smcObject.updateFalloffProb(deathArray,longlived_deathArray)
            # implement targeted loading mechanism
            loading_prob[DSB_coordinates]= 1/chrom_size*lof[idx][11]
            self.smcObject.updateEmissionProb(loading_prob)
            pauseArray = PAUSEPROB*np.ones(chrom_size, dtype=np.double) 
            self.smcObject.updatePauseProb(pauseArray) # increase sampling rate to not miss any synapsis events
            

        #precalculating all bonds
        allBonds = []
        
        left, right = self.smcObject.getSMCs()
        # add SMC bonds
        bonds = [(int(i), int(j)) for i,j in zip(left, right)]

        left, right = self.smcObject.getlonglivedSMCs()
        # add longlivedSMC bonds
        bonds += [(int(i), int(j)) for i,j in zip(left, right)]
        self.curBonds = bonds.copy()
        
        allBonds.append(bonds)
        for dummy in range(blocks):
            self.smcObject.steps(smcStepsPerBlock)
            left, right = self.smcObject.getSMCs()
            # add SMC bonds
            bonds = [(int(i), int(j)) for i,j in zip(left, right)]
            
            left, right = self.smcObject.getlonglivedSMCs()
            # add longlivedSMC bonds
            bonds += [(int(i), int(j)) for i,j in zip(left, right)]

            allBonds.append(bonds)
        
        self.allBonds = allBonds
        self.uniqueBonds = list(set(sum(allBonds, [])))
        
        allBonds.pop(0)
        

        #adding forces and getting bond indices
        self.bondInds = []
        
    
        for bond in self.uniqueBonds:
            paramset = self.activeParamDict if (bond in self.curBonds) else self.inactiveParamDict
            ind = bondForce.addBond(bond[0], bond[1], **paramset) # changed from addBond
            self.bondInds.append(ind)
        self.bondToInd = {i:j for i,j in zip(self.uniqueBonds, self.bondInds)}
        return self.curBonds,[]


    def step(self, context, verbose=False):
        """
        Update the bonds to the next step.
        It sets bonds for you automatically!
        :param context:  context
        :return: (current bonds, previous step bonds); just for reference
        """
        if len(self.allBonds) == 0:
            raise ValueError("No bonds left to run; you should restart simulation and run setup  again")

        pastBonds = self.curBonds
        self.curBonds = self.allBonds.pop(0)  # getting current bonds
        bondsRemove = [i for i in pastBonds if i not in self.curBonds]
        bondsAdd = [i for i in self.curBonds if i not in pastBonds]
        bondsStay = [i for i in pastBonds if i in self.curBonds]
        if verbose:
            print("{0} bonds stay, {1} new bonds, {2} bonds removed".format(len(bondsStay),
                                                                            len(bondsAdd), len(bondsRemove)))
        bondsToChange = bondsAdd + bondsRemove
        bondsIsAdd = [True] * len(bondsAdd) + [False] * len(bondsRemove)
        for bond, isAdd in zip(bondsToChange, bondsIsAdd):
            ind = self.bondToInd[bond]
            paramset = self.activeParamDict if isAdd else self.inactiveParamDict
            self.bondForce.setBondParameters(ind, bond[0], bond[1], **paramset)  # actually updating bonds
        self.bondForce.updateParametersInContext(context)  # now run this to update things in the context
        return self.curBonds, pastBonds
    
## Simulation parameters and model setup
## Define the CTCF motifs
L = [0, 201, 602, 803, 1604, 1805, 2206, 2407] # in the unit of number of monomers (each monomer correspond to 1kb of DNA to emulate 1D simulation)
R = [200, 601, 802, 1603, 1804, 2205, 2406, 3607]
tad_start, tad_end = (0, 3607)
num_segments = tad_end-tad_start+1
stallL = np.zeros(num_segments)
stallR = np.zeros(num_segments)

nsamples = 70000//num_segments 
initialization_steps = 10000
steps_per_sample = 100

## Dummy leftover variables 
num_cores = 1
samps_per_core = 1
####

processivities = [250]
separations = [125]
boundaryPauseProb = [0.5]

ctcf_boost_factors = [16]
longlived_fraction = [0.2]
longlived_boost_factor = [20]
dsb_boost_factor = [2]
targeted_loading_factor = [1000]
llp = [2]
h = [0.4]

for p,s  in zip(R,boundaryPauseProb*len(R)):
    stallR[p] = s
for p,s  in zip(L,boundaryPauseProb*len(L)):
    stallL[p] = s 

pad = int(np.mod(70000,3608)/2) # padding at the two ends of the chromosome to make total number of monomers exactly 70000
boundaryStrengthsL = np.tile(stallL,nsamples)
boundaryStrengthsL = np.hstack((np.zeros(pad),boundaryStrengthsL))
boundaryStrengthsL = np.hstack((boundaryStrengthsL,np.zeros(pad)))
boundaryStrengthsR = np.tile(stallR,nsamples)
boundaryStrengthsR = np.hstack((np.zeros(pad),boundaryStrengthsR))
boundaryStrengthsR = np.hstack((boundaryStrengthsR,np.zeros(pad)))
boundary_coordinates = np.hstack(np.argwhere(boundaryStrengthsL+boundaryStrengthsR>0))
boundary_coordinates = boundary_coordinates.tolist()

chrom_size = len(boundaryStrengthsL)

# randomly determine the coordinate of the DSB in the first TAD
random_pick_0 = int(np.random.random()*(R[1]-L[1]-2)+L[1]+1)+pad
if random_pick_0+1 == R[1]:
    DSB_coordinates = [random_pick_0-1, random_pick_0]
else:
    DSB_coordinates = [random_pick_0, random_pick_0+1]

DSB_TAD_array_size = [R[1]-L[1]]
tad_count = np.asarray([1,0,0,0]) # counter for TAD of each size containing DSB
tad_size = np.asarray([200, 400, 800, 1200]) 
dist_between_break = 3000 # in monomer size(kb), distances between two DSB
within_tile = True
while within_tile:
    # skip dist_between_break 
    next_break = DSB_coordinates[-1] + dist_between_break
    # preventing the next DSB occurring at a boundary element
    if np.count_nonzero(next_break - np.asarray(boundary_coordinates))<len(boundary_coordinates):
        next_break += 2
    # determine the TAD the next DSB occurs in
    test = np.asarray(boundary_coordinates) - next_break
    boundray_index = np.argmax(test>0) # the first BE coordinates to the right of the DSB
    # randomize the coordinate of the DSB within the TAD (so that the distance between DSB and BE is randomized)
    random_pick = int(np.random.random()*(boundary_coordinates[boundray_index]-boundary_coordinates[boundray_index-1]-2)+boundary_coordinates[boundray_index-1]+1)
    if random_pick+1 == boundary_coordinates[boundray_index]:
        DSB_coordinates += [random_pick-1, random_pick]
    else:
        DSB_coordinates += [random_pick, random_pick+1]
    tad_size_index = np.argwhere(tad_size==boundary_coordinates[boundray_index]-boundary_coordinates[boundray_index-1])[0][0]
    tad_count[tad_size_index]+=1
    DSB_TAD_array_size.append(tad_size[tad_size_index])
    if DSB_coordinates[-1] + dist_between_break> chrom_size-1-pad*2:
        within_tile = False

# create chain with breaks
# chain segment lengths 

break_sites_seg = [0] + DSB_coordinates[1::2]+ [chrom_size]
chain = [(x,y,0) for (x,y) in zip(break_sites_seg[0:-1],break_sites_seg[1:])]
extra_bonds = [(x,y) for (x,y) in zip(DSB_coordinates[0::2],DSB_coordinates[1::2])] # add "links" between the DSB ends prior to DSBs

positions_to_sample = DSB_coordinates.copy()

experiment_TAD_size = 515 # fbn2 TAD

internal_to_sample =  [2*x*1000 + int(1000 - int(experiment_TAD_size/2)) for x in range(int(chain[-1][1]/2/1000))]+[2*x*1000+int(1000 + int(experiment_TAD_size/2)+1) for x in range(int(chain[-1][1]/2/1000))] # For calibration with experimental data
positions_to_sample += internal_to_sample

from copy import deepcopy

# get initial conformation
conformation_folder = f'Data/3D_PolymerSimulationEquilibrationrRun/blocks_490000-499999.h5::499999'

from copy import deepcopy
u = load_URI(conformation_folder)  
data = u['pos']

# saving for polymer simulation
GPU_choice = 1

steps_per_block_beforeDSB = 17369 # number of polymer simulation steps per block (17369)
steps_per_block = 434 # number of polymer simulation steps per block (~1/40 of 17369)
total_saved_blocks =  10000000 # total number of blocks saved
save_every_x_blocks = 1 # number of blocks before a 3D configuration is saved
total_blocks_done = total_saved_blocks * save_every_x_blocks
BondUpdaterCountBeforeDSB = int(10/(steps_per_block_beforeDSB/17369))
restartBondUpdaterEveryBlocks = 100 


# block number to store all monomer positions
block_to_save_all = [BondUpdaterCountBeforeDSB*restartBondUpdaterEveryBlocks+1]
block_to_save_all = [i+1 for i in block_to_save_all]

# parameters for SMC translocator
# simulation parameters for smc bonds 
smcBondWiggleDist = 0.2
smcBondDist = 0.5
smcStepsPerBlock = 1#int(1) # this scales the number of SMC substeps per simulation time step

stiff = 2
volume_density = 0.2
block = 0  # starting block 
N_monomers = chain[-1][1]
N_breaks = int(len(DSB_coordinates)/2)
collision_rate = 1
truncated_potentials = 50 # truncation potentials of 10 will prevent most of the polymer crossing

pause_prob_beforeDSB = 1-steps_per_block_beforeDSB/17369 # probability motor pauses
pause_prob = 1-steps_per_block/17369 # probability motor pauses
# parameters for SMC translocator
translocator_initialization_steps = 10000/(1-pause_prob_beforeDSB) # for SMC translocator

date = datetime.today().strftime('%Y_%m_%d')
save_folder = f'Data/DSB_sim_{date}_dens{volume_density}_N{N_monomers}_DSB{N_breaks}_collisionrate_{collision_rate}_stepsPerBlock{steps_per_block}_trunc{truncated_potentials}_restartBondUpdaterEveryBlocks{restartBondUpdaterEveryBlocks}_proc{processivities[0]}_sep{separations[0]}_ctcf{ctcf_boost_factors[0]}_longlivedfraction{longlived_fraction[0]}_longlivedfactor{longlived_boost_factor[0]}_dsb{dsb_boost_factor[0]}_superloading{targeted_loading_factor[0]}_pauseProb{pause_prob}_v5_run98'
name_format_string = save_folder+'/Sim_tests_proc{}_sep{}_lp{}_H{}_boundaryPauseProb{}_ctcf{}_longlivedfraction{}_longlivedfactor{}_dsb{}_superloading{}_nsamples{}_nreps{}'

lof_nums = list(product(llp,boundaryPauseProb,processivities,separations,[''],h,ctcf_boost_factors,longlived_fraction,longlived_boost_factor,dsb_boost_factor,targeted_loading_factor))
lof = [(x[0],x[1],name_format_string.format(x[2],x[3],x[0],x[5],x[1],x[6],x[7],x[8],x[9],x[10],nsamples,int(samps_per_core*num_cores) ),
        x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10]) for x in lof_nums]

if not os.path.exists(save_folder):
    os.makedirs(save_folder)

LogFileName = save_folder+"/log.txt"

o = open(LogFileName, "w") 
o.close()

with open(LogFileName, "a") as f:
    f.write('Tiling factor: \n')
    f.write(str(nsamples)+'\n')
    for i in range(len(tad_size)):
        f.write('number of breaks within a '+str(tad_size[i])+'kb TADs: \n')
        f.write(str(tad_count[i])+'\n')

file = save_folder+'/DSB_boundary_coordinates.npy' 

# saving DSB coordinates and boundary coordinates to robustness 
if not os.path.exists(file):
    with open(file, 'wb') as g:
        np.save(g, np.asarray(DSB_coordinates))
        np.save(g, boundary_coordinates)
        
idx = 0 # dummy for now
run_simulation(idx, N_monomers, translocator_initialization_steps, \
                  smcStepsPerBlock, \
                  smcBondDist, smcBondWiggleDist,\
                  steps_per_block_beforeDSB, steps_per_block, volume_density, block, \
                  chain, extra_bonds, save_folder, save_every_x_blocks, \
                  total_saved_blocks, restartBondUpdaterEveryBlocks,BondUpdaterCountBeforeDSB, GPU_choice = GPU_choice, 
                  overwrite=False,density=0.2,positions_to_sample=positions_to_sample,colrate=collision_rate,errtol=0.01,
                  trunc=truncated_potentials,initial_conformation=data,block_to_save_all = block_to_save_all,save_length=1000)