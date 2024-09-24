#=
--- VCell Parameter Array Generator ---
This script generates parameter distributions for use with VCell's batch simulation feature
using external, user-defined parameter arrays. It is a Julia implementation of the MATLAB
(live) script of the same name.
=#
## Activate project environment:
# Store current working directories where this script and the Project.toml and Manifest.toml
# files for this project are located:
cwd0 = pwd()                             # path with project environment
cwd = pwd() * "/VCell_batch_results"     # path with VCell data
cd(cwd)

# Activate the project's environment:
using Pkg
Pkg.activate(cwd0)
# Pkg.instantiate()     # installs all (missing) package dependencies for the current project

## Load packages and set random seed:
using Distributions, QuasiMonteCarlo, Random
using CSV, DataFrames
Random.seed!(123)   # set random seed

## Load baseline parameter values:
fn = "fullEGFR_fitted_params_5-12-21.csv"   # file containing base parameter values to perturb
p_in = CSV.File(fn) |> DataFrame     # load parameter names and their baseline values
Base.names(p_in) |> print     # column names in parameter DataFrame


## Select parameters to modify:
use_all_params = true  # choose whether to perturb all parameters or a specific subset of them
if use_all_params
    pnames = p_in.name
else
    pnames = ["kEf","kEr","kdEf","kdEr","kcatE","kdp"]   # names of specific parameters
end
param_inds = in(pnames).(p_in.name) # get indices of parameters we want to change


## Generate parameter distributions:
base_vals = p_in.value[param_inds]  # baseline values of parameters we want to change
lb = base_vals .* 0.1  # lower bound on parameter values
ub = base_vals .* 10.0  # upper bound on parameter values
nspls = 198  # number of parameter samples we want

#=
Available sampling methods:
    → LatinHypercubeSample(), UniformSample(), SobolSample(), LatticeRuleSample(), LowDiscrepancySample(), GridSample()
    → Any distribution defined through `Distribution` or a properly defined custom sampling algorithm
See the `QuasiMonteCarlo.jl` documentation for more information on selecting/defining sampling distributions.
=#
spl_method = LatinHypercubeSample() # select sampling method
params = QuasiMonteCarlo.sample(nspls, lb, ub, spl_method)  # generate parameter distributions


## Format parameter samples for VCell batch input and save file:
pnames_vcell = pnames .* "="  # add "=" to the end of each parameter name
params_vcell = permutedims(pnames_vcell .* string.(params)) # add parameter values as strings after "="
p_out_vcell = DataFrame(params_vcell, :auto)    # put parameter array in correct format for CSV.write
p_out = DataFrame(params', pnames)
fn_out = "param_array2"     # output file name, minus the file identifier
CSV.write(fn_out * ".dat", p_out_vcell, writeheader=false) # generate parameter array file for VCell
CSV.write(fn_out * ".csv", p_out, writeheader=true) # generate parameter array file for loading outside of VCell
