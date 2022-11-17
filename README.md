# DNA_break_synapsis_models
This repository contains source code for the article titled _DNA double-strand break end synapsis by DNA loop extrusion_, investigating the role of loop extrusion in DNA double-strand break end synapsis. We will describe how to run 1D simulation with different parameter combinations and generate in silico ChIP-seq metaplots and Hi-C contact maps. 

The jupyter notebooks titled with `_analysis` suffix are for analyzing simulation data.

The `Mathematica Notebooks` folder contains our analytical model to predict synapsis efficiency, the `Data` folder contains the simulation results plotted in the aritcle, and the Figures folder contains output images from the analysis code.

## System requirements

Operating systems should allow running Python. Windows, macOS, and Linux are all supported by Python.

## Installation

Download source code (should only take a minute):

`$ git clone https://github.com/ahansenlab/DNA_break_synapsis_models`

## Installation requirements

### Python dependencies
Python 2.7/3.3+

The following dependencies should be installed correctly before running the code this package.
 * numpy (1.6+)
 * scipy
 * matplotlib
 * cython
 * pandas
 * seaborn
 * multiprocess (if using parallel processing)

We recommend using the conda package manager if you have trouble building these Python dependencies:

`$ conda install numpy scipy matplotlib cython ...`

## Demo

The remaining sections provide instructions to run 1D simulations and analyze the simulation results. Default parameters have been defined for running all the simulation and analysis notebooks so that all notebooks could be run directly (only need to provide user-defined folder names).

Sample output data are provided in the `Data` folder, which were used to generate figures in the manuscript, and can be analyzed using the analysis notebooks for reproducing the figures shown in the manuscript.

Most simulations should finish within 24 hours on a "normal" desktop computer.

## Running 1D simulations

The jupyter notebook `1D_SimulationCode.ipynb` is used for performing 1D simulations of the DNA double-strand break end synapsis by DNA loop extrusion.

**Step 1**. Set up the chromosome

In the second cell, define output folder on line 2, where all simulation results will be stored. Define tiling factor (the number of times the TAD array will be repeated) on line 6. Setting up the TAD array by defining boundary elements on line 9.

**Step 2**. Define simulation parameters

In the third cell, define simulation parameters including spearation, normal LEF processivity, boundary strength (the probability of LEF stalled by an boundary element), the ratio of long-lived LEF processivity over normal LEF processivity, the fold stabilization of LEF at BE, the fold stabilization of LEF at DSB ends, fold increase of loading probability at the DSB ends, the threshold distance for calling synapsis, the probability LEF taking a step, maximum number of simulation time steps post DSB to stop simulation, and whether to take snapshots of the synapsis process.

**Step 3**. Define parallel processing parameters (optional)

In the last cell, define the number of cores the simulation can utilize to parallelize simulations with different parameter combinations.

## Analyzing 1D simulation results

See jupyter notebooks with `_analysis` suffix for examples of analysis of 1D simulation results.

## Running 3D polymer simulations

The python code `3D_PolymerSimulationCode(withLoopExtrusion).py` is used for performing 3D polymer simulations of the DNA double-strand break end synapsis by DNA loop extrusion; the python code `3D_PolymerSimulationCode(withoutLoopExtrusion).py` is used for performing 3D polymer simulations of the DNA double-strand break end synapsis by passive diffusion.

## Analyzing 3D polymer simulation results

Snakemake pipeline in the `3DSimAnalysis_SnakemakePipeline` folder is used for analysis of 3D polymer simulation results.

## Generating in silico ChIP-seq metaplots of LEFs

**Step 1**. Generate the regions file containing DSB coordinates

In the jupyter notebook `WriteDSBcoordinatesToBed.ipynb`, modify the output_folder to be the folder containing 1D simulation results of interest in the second cell, and run the entire notebook to obtain the regions file `DSBcoordinates.bed` containing all the DSB coordinates in the 1D simulations, which will be stored in the same output_folder.

**Step 2**. Generate the score file containing the number of LEFs around the DSB

In the jupyter notebook `WriteLEFcountToBed.ipynb`, modify the output_folder to be the same output folder defined in step 1, and run the entire notebook to obtain the score files containing the LEF count in each bin (default bin size = 5 kb, can be modified on line 3 of the last cell).

**Step 3**. Convert the score file from step 2 into BigWig format

This step requires the bedGraphToBigWig package from UCSC Genome Browser(http://hgdownload.cse.ucsc.edu/admin/exe/), which can be installed by running the following line from comand line:

```bash
conda install -c bioconda ucsc-bedgraphtobigwig
```

1. Create a DSBchrom.sizes text file in the output_folder defined in step 1, with "chr1"	as the first column, and the total length of your chromosome in bp as the second column.
2. Run the following line from command line, with the output_folder as the argument
```bash 
python BedToBigwig.py output_folder
```

**Step 4**. Plot ChIP-seq metaplots

This step requires the deepTools package (https://github.com/deeptools/deepTools), which can be installed by running the following line from comand line:

```bash
conda install -c bioconda deeptools
```

Run the following line from command line, with the output_folder as the argument
```bash 
python heatmap.py output_folder
```

## Generating Hi-C contact maps from LEF coordinates

The jupyter notebook `GeneratingContactMaps.ipynb` is used for generating Hi-C contact maps based on LEF coordiates output from 1D simulations.

Modify the base_folder_string variable to be the common string of the folder names of 1D simulation output folders to extract LEF coordinates from.

Define the index of the time point to pull LEF coordinates from by modifying the time_point_list variable. Index 0-6 correspond to 4min, 10min, 20min, 60min, 90min, 2hr post DSB, respectively.

Run the entire notebook to generate the contact map.