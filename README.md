# Tardanico and Hovestadt 2023 Simulation Model

A spatially explicit individual-based metacommunity simulation model simulating communities of annual, asexual organisms living in continuous fractal landscapes.

## Requirements

The simulation program is written in Julia 1.1.1 and requires this version of Julia to run. In addition, the program requires the following Julia packages: Distributions, ArgParse, DelimitedFiles, Plots, Random, StatsBase, DataFrames, CSV, Printf. As Julia is a just-in-time compiled language and is compiled at runtime, the code itself does not require any particular OS.

## Set-up

Julia 1.1.1 can be downloaded here: https://julialang.org/downloads/oldreleases/
To install the required packages, run the following code, either in the Julia REPL or with a Julia script via the console command line: 

```
using Pkg

Pkg.add("Distributions")
Pkg.add("ArgParse")
Pkg.add("DelimitedFiles")
Pkg.add("Random")
Pkg.add("StatsBase")
Pkg.add("DataFrames")
Pkg.add("CSV")
Pkg.add("Plots")

Pkg.build("Rmath")
ENV["GRDIR"]=""
Pkg.build("GR")
```
## Simulation set-up

The program requires two landscape files and a parameter config file as input. Landscape files are .txt or .csv files containing values for patch environmental attributes. Both landscape files must have matching dimensions. Parameter config files consist of a Julia dictionary. Please refer to the following example:

```
# General Simulation parameters
"burnin" => false, # Burn-in period? Boolean
"bmax" => 50, # Number of timesteps in burn-in period
"tmax" => 10000, # Number of timesteps in main Simulation Run
"rmax" => 1, # Number of replicates. Will repeat simulation with the same landscape. 
             # Keep at 1 and do separate simulation runs.
"cmax" => 200, # Number of timesteps in the climate change period
"c_change" => 0.02, # Rate of change in T
"T_ref" => 12.5, # This is unnecessary for the actual model. Remove this.
"grad" => 0.1, # Gradient steepness
"%*α*)" => 3, # Trade-off strength
"scen" => 1, # Climate scenario. DEFUNCT PARAMETER. KEEP AT 1
#"initial_pop" => 10, # Initial number of individuals in a patch.
"carry_capacity" => 150, # Carring capacity of a patch
"immi" => true, # Is there immigration from outside the landscape? True or false.
"mutate" => false, # Do organisms undergo mutation? True or false
"p_mut" => 0.001, # Probability of mutation.
"mut_sd" => 0.05, # Standard deviation of trait change from mutation.
"mut_decay" => 0, # Change in mutation sd over time. Set to 0 for constant mutation sd.
"p_immi" => 0.005, # Chance of immigration to a patch.
"e_immi" => 2.5, # Expected immigrants in a patch
#"seed" => 1234, # RNG seed
"output_time" => 10000, # Specify when the program should output data
#"output_interval" => 50, # Timestep interval for output of data for individual organisms.
                          # Program will only output once if this is not defined.
#::::::::::::::::::-
# Initial trait distribution parameters
"topt_%*µ*)" => 0, # %*µ*) of T_opt initialization distribution
"topt_%*σ*)" => 1, # %*σ*) of T_opt initialization distribution
"hopt_%*µ*)" => 0, # %*µ*) of H_opt initialization distribution
"hopt_%*σ*)" => 1, # %*σ*) of H_opt initialization distribution
"tsd_%*µ*)" => 0, # %*µ*) of T_sd initialization distribution
"tsd_%*σ*)" => 1, # %*σ*) of T_sd initialization distribution
"hsd_%*µ*)" => 0, # %*µ*) of H_sd initialization distribution
"hsd_%*σ*)" => 1, # %*σ*) of H_sd initialization distribution
#:::::::::::::::::::
# Filename and output directory parameters
"region" => "", # Name of region. Leave as "" if not using real world landscape data. (DEFUNCT; LEAVE AS IS)
"temp_type" => "autocorrelated", # Type of landscape for T attribute. Uniform, clustered, or autocorrelated
"env_type" => "autocorrelated", # Type of landscape for H attribute. Uniform, clustered, or autocorrelated
"autocor_temp" => 0, # Autocorrelation degree (Hurst index) or cluster size of temperature dimension.
"autocor_env" => 0, # Autocorrelation degree (Hurst index) or cluster size of habitat dimension.
"dir" => "/my_output_directory/", # Output directory
)
```

This simulation program is intended to run in the console via a shell script. The script should look like the following example:

```
    # Example shell script

cd /my_install_directory/julia-1.1.1/bin/  # Your install directory for Julia 1.1.1

echo "My scenario, starting script"
source=/my_code_files_directory/LandscapesModel1.jl # The main program

# Your parameter config file for the scenario you want to run
para=/my_parameter_config_directory/my_parameters.jl 

# Optional code for merging scenario output into combined files
merge=/my_code_files_directory/MergeFiles.jl 

# Specify your landscape files to be used with the simulation.
T1=/my_landscape_files/landscape_file.txt 
T2=/my_landscape_files/landscape_file2.txt

H1=/my_landscape_files/landscape_file3.txt
H2=/my_landscape_files/landscape_file4.txt


echo "source = $source"
echo "para = $para"

# Run the program. Each line is a replicate.
./julia -p 1 $source -n $para -t $T1  -e $H1 -r 1 
./julia -p 1 $source -n $para -t $T2  -e $H2 -r 2

echo "Merging output files"

# Run this to merge output files together (Optional)
./julia -p 1 $merge -n $para

echo "Script complete."
```

Note that -r is an optional argument for specifying a seed for the random number generator (RNG). This may be desirable if there is a need to precisely replicate a simulation run. An RNG seed can also be specified in the parameter config file. If a seed is specified in the shell, the program will always use that seed. If no RNG seed is specified, the seed will be based on the system clock.

## Output

The program produces two output files, one file containing data over time for various landscape-wide statistics, and one file containing data on individual organisms in the landscape at a specified time step. The program addtionally creates a third data file containing patch aggregated data calculated from individual data. 

### Data variables for individuals:

|Variable | Explanation |
|---------|-------------|
|Replicate | Replicate number|
|Timestep | The time step of data output|
|H_t | Hurst index for T attribute|
|H_h | Hurst index for H attribute|
alpha| Niche breadth trade-off parameter
clim_scen | Climate scenario (DEFUNCT)
gradient | Strength of compositional heterogeneity (G)
x,y | The x and y indices of the patch an organism is in
T_opt | T niche optimum
H_opt | H niche optimum
T_sd | T niche breadth
H_sd | H niche breadth
disp_l | Dispersal probability (Pdisp)
disp_g | Probability of dispersing via global dispersal (Pglobal)
fert_max | Intrinsic reproduction rate (R0), maximum expected offspring
LineageID | Identifier for organism lineage/species
Origin_x, Origin_y | Lineage patch of first appearance
Origin_time | Lineage time of first appearance
temp_t | Patch T attribute
habitat | Patch H attribute
trend | Fluctuation in T at a given timestep
mean_trend | Mean of T fluctuation,0 if no climate shift
precip_t | Defunct patch attribute from an early iteration of the program

### Patch aggregated data variables:

|Variable | Explanation|
|---------|------------|
x,y | The x and y indices of a patch
Replicate | Replicate number
Timestep | The time step of data output
H_t Hurst | index for T attribute
H_h Hurst | index for H attribute
alpha | Niche breadth trade-off parameter
clim_scen | Climate scenario (DEFUNCT)
gradient | Strength of compositional heterogeneity (G)
pop | patch total population
richness | Patch lineage/species richness
shannon | Patch Shannon-Weiner diversity
simpson | Patch Simpson diverisity
tdiff | Mean T_opt-patch T difference at current timestep
mt_tdiff | Mean T_opt-patch T difference on average over time
tft | Mean fitness for patch T at current timestep[^1]
mft | Mean fitness for patch average T over time[^1]
mfh | Mean fitness for patch H[^2]
time_fit | mean overall fitness at current timestep[^3]
mean_fit | mean overall fitness under patch average conditions over time[^3]
T_opt | patch mean T niche optimum
H_opt | patch mean H niche optimum
T_sd | patch mean T niche breadth
H_sd | patch mean H niche breadth
disp_l | Patch mean dispersal probability (Pdisp)
disp_g | Probability of dispersing via global dispersal (Pglobal)
fert | Patch mean intrinsic reproduction rate (R0), maximum expected offspring
temperature | Patch T attribute
habitat | Patch H attribute
trend | Fluctuation in T at a given timestep

[^1]: Fitness calculated for the T component of an organism’s niche.
[^2]: Fitness calculated for the H component of an organisms’s niche.
[^3]: Combined product of fitness for T and H niche components.

### Data variables for landscape time series:

|Variable | Explanation|
|---------|------------|
Replicate | Replicate number
Timestep | The time step of data output
H_t Hurst | index for T attribute
H_h Hurst | index for H attribute
alpha | Niche breadth trade-off parameter
clim_scen | Climate scenario (DEFUNCT)
grad | Strength of compositional heterogeneity (G)
T_opt | Landscape mean T niche optimum
H_opt | Landscape mean H niche optimum
T_sd | Landscape mean T niche breadth
H_sd | Landscape H niche breadth
disp | Landscape dispersal probability (Pdisp)
disp_g | Landscape mean global dispersal probability (Pglobal)
VT_opt | Landscape T_opt variance
VH_opt | Landscape H_opt variance
VT_sd | Landscape T_sd variance
VH_sd | Landscape H_sd variance
Vdisp | Landscape disp variance
Vdisp_g | Landscape disp_g variance
fert_max | Intrinsic reproduction rate (R0), maximum expected offspring
trend | Fluctuation in T at a given timestep
mean_trend | Mean of T fluctuation, 0 if no climate shift
