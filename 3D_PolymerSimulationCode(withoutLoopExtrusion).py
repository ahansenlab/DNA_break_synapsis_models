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

import warnings
import h5py 
import glob

from itertools import product
import re

from datetime import datetime

from polychrom.hdf5_format import HDF5Reporter, list_URIs, load_URI, load_hdf5_file
import h5py as hp
from matplotlib import pyplot as plt



from collections.abc import Iterable

import simtk.openmm as openmm
import simtk.unit

from polychrom import forces

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
    
def run_simulation(idx, N_monomers, \
                  steps_per_block, volume_density, block, \
                  chain, extra_bonds, save_folder, save_every_x_blocks, \
                  skip_x_blocks_at_start, total_saved_blocks, GPU_choice = 0, 
                   overwrite=False,positions_to_sample=None,colrate=0.3,errtol=0.01,trunc=0.001,
                  initial_conformation=None,block_to_save_all = None, save_length=10000):
    
    print(f"Starting simulation. Doing {total_saved_blocks} steps, saving every {save_every_x_blocks} block.")    
    print(f"A total of {total_saved_blocks * steps_per_block} steps will be performed")
       
    # clean up the simulation directory
    folder = save_folder
    if os.path.exists(folder):
        if overwrite == False:
            print(f"Folder {folder} already exists")
            return 0
        else:
            shutil.rmtree(folder)
    
    # set initial configuration of the polymer    
    if initial_conformation is None:
        box = (N_monomers / volume_density) ** 0.33  
        data = grow_cubic(N_monomers, int(box) - 2)  # creates a compact conformation, but one can take pre-existing configurations
    else:
        data = initial_conformation
        
    # create the reporter class - this saves polymer configurations and other information you specify
    reporter = HDF5Reporter(folder=folder, max_data_length=save_length,overwrite=True)

    # simulation parameters are defined here     
    sim = DSB_Simulation(
            platform="cuda",
            integrator="variableLangevin", 
            error_tol=errtol, 
            GPU = "{}".format(GPU_choice), 
            collision_rate=colrate, 
            N = len(data),
            reporters=[reporter],
            PBCbox=False,#[box, box, box],
            precision="mixed")  # timestep not necessary for variableLangevin

    sim.set_data(data)  # loads a polymer, puts a center of mass at zero
    sim.add_force(forces.spherical_confinement(sim,density=volume_density,k=1))
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
                # when a new bond is initialized
             },

            angle_force_func=forces.angle_force,
            angle_force_kwargs={
                'k':0.05 # we are making a very flexible polymer, basically not necessary here
            },
            nonbonded_force_func=forces.polynomial_repulsive, # this is the excluded volume potential
            nonbonded_force_kwargs={
                'trunc':trunc,
                'radiusMult':1.05, # this is from old code
            },
            except_bonds=True,
            extra_bonds=extra_bonds,
            override_checks=True,
        )
    )

    sim.step = block

    # Minimize energy since we're reconfiguring the random walk into crumples and SMC-dependent loops
    sim.local_energy_minimization() 

    # Iterate over simulation time steps within each BondUpdater 
    for b in range(total_blocks_done):
        if (b % save_every_x_blocks == 0) and (b > skip_x_blocks_at_start):  
            # save SMC positions and monomer positions 
            if b in block_to_save_all:
                sim.do_block(steps=steps_per_block, save_extras={"SMC_step":sim.step},
                         positions_to_sample=[i for i in range(N_monomers)]) 
            else:
                sim.do_block(steps=steps_per_block, save_extras={"SMC_step":sim.step},
                         positions_to_sample=positions_to_sample) 
        else:
            sim.integrator.step(steps_per_block)  # do steps without getting the positions from the GPU (faster)

    data = sim.get_data()  # get the polymer positions (i.e. data) 
    block = sim.step       # get the time step
    del sim               # delete the simulation
    time.sleep(0.2)      # wait 200ms to let garbage collector do its magic
    reporter.dump_data() # dump data to  output file!
    done_file = open(os.path.join(folder,'sim_done.txt'),"w+")
    done_file.close()
    
    
# create multiple chains
chain_len = 2500
num_chains = 28
extra_bonds = False
chain = [(x*chain_len,(x+1)*chain_len,0) for x in range(num_chains)]


experiment_TAD_size = 515
ends_to_sample = [x*chain_len for x in range(num_chains)]+[x*chain_len-1 for x in np.arange(1,num_chains+1)]
internal_to_sample = [x*chain_len +int(chain_len/2 - int(experiment_TAD_size/2)) for x in range(num_chains)]+[x*chain_len +int(chain_len/2 + int(experiment_TAD_size/2)) for x in range(num_chains)]
positions_to_sample = ends_to_sample + internal_to_sample

from copy import deepcopy

# get initial conformation
conformation_folder = f'Data/3D_PolymerSimulationEquilibrationrRun/blocks_490000-499999.h5::499999'

from copy import deepcopy
u = load_URI(conformation_folder)  
data = u['pos']
    
# saving for polymer simulation
GPU_choice = 0

steps_per_block = 30 # number of polymer simulation steps per block.
total_saved_blocks = 1000000000 # total number of blocks saved
save_every_x_blocks = 1 # number of blocks before a 3D configuration is saved
skip_x_blocks_at_start = 2 # does not save N initial blocks after restarting the bond positions

total_blocks_done = total_saved_blocks * save_every_x_blocks

# block number to store all monomer positions
block_to_save_all = [1999999,3999999,6999999,9999999,19999999,29999999,39999999,49999999]
block_to_save_all = [i+skip_x_blocks_at_start+1 for i in block_to_save_all]

# simulation parameters for smc bonds  
smcBondWiggleDist = 0.2
smcBondDist = 0.5

stiff = 2
volume_density = 0.2
block = 0  # starting block 
N_monomers = chain[-1][1]
collision_rate = 1
truncated_potentials = 50 # truncation potentials of 10 will prevent most of the polymer crossing

date = datetime.today().strftime('%Y_%m_%d')
save_folder = f'Data/multichain_sim_{date}_dens{volume_density}_N{N_monomers}_C{num_chains}_collisionrate_{collision_rate}_trunc{truncated_potentials}_stepsperblock{steps_per_block}_v5'

if not os.path.exists(save_folder):
    os.makedirs(save_folder)

idx = 0 # dummy for now
run_simulation(idx, N_monomers, \
                  steps_per_block, volume_density, block, \
                  chain, extra_bonds, save_folder, save_every_x_blocks, \
                  skip_x_blocks_at_start, total_saved_blocks, GPU_choice = GPU_choice,
                  overwrite=True, positions_to_sample=positions_to_sample,                             colrate=collision_rate,
                  errtol=0.01,trunc=truncated_potentials,initial_conformation=data,block_to_save_all = block_to_save_all)