#=
    This scripts contains all necessary code for the model calculations that
    were run on the University of Virginia's Rivanna HPC system. The enclosing
    directory of this script contains the Manifest.toml and Project.toml files
    associated with the package versions that were used for these calculations
    on Rivanna.

    PLEASE NOTE: Julia v1.6.2 was used for these calculations. Attempting to
    run this script with a different version of Julia or different package
    versions may result in unintended behavior and/or other unforeseen errors.
=#

## Load packages:
using Pkg
Pkg.activate(".")
Pkg.instantiate()
Pkg.precompile()

## Launch worker processes and load packages:
using Distributed, ClusterManagers
if nprocs() == 1
    # addprocs_slurm(1000, nodes=25, partition="parallel", time="71:59:59")
    if nprocs() == 1
        addprocs(Threads.nthreads(), exeflags="--project=."); # adds the same number of processes as available CPU cores
    end
end
nprocs()    # check number of available processes
nworkers()  # check number of available workers
@everywhere begin
    using DifferentialEquations, Catalyst, SharedArrays, DataFrames, DataFramesMeta, Statistics
end
using CSV, JLD2, Pipe, QuasiMonteCarlo, StatsBase, Dates
using Random
Random.seed!(123)


## Model reaction network(s):
@everywhere begin
    receptor_rn = @reaction_network begin
        # EGF binding (MEMBRANE):
        kEf, EGF + EGFR → mELₘ + EGF
        kEr, EGFR ← mELₘ

        # EGFR dimerization (MEMBRANE):
        (kdEf, kdEr), 2mELₘ ↔ mELmELₘ

        # EGFR phosphorylation and dephosphorylation (MEMBRANE):
        (kcatE, kdp), mELmELₘ ↔ Eₘ
    end kEf kEr kdEf kdEr kcatE kdp

    adaptor_rn = @reaction_network begin
        # EGFR-GRB2 binding (MEMBRANE):
        (kG2f, kG2r), Eₘ + GRB2 ↔ EG2ₘ

        # EGFR-GRB2:SOS binding (MEMBRANE):
        (kG2SOSf, kG2SOSr), Eₘ + G2SOS ↔ EG2SOSₘ

        # EGFR:GRB2-SOS binding (MEMBRANE):
        (kSOSf, kSOSr), EG2ₘ + SOS ↔ EG2SOSₘ

        # GRB2-SOS binding (CYTOPLASM):
        (kSOSf, kSOSr), GRB2 + SOS ↔ G2SOS
    end kG2f kG2r kG2SOSf kG2SOSr kSOSf kSOSr #kdeg

    mapk_rn = @reaction_network begin
        # Ras guanine nucleotide exchange rxns (MEMBRANE):
        kRgne/(Kmgne + RAS), EG2SOSₘ + RAS → RAStₘ + EG2SOSₘ    # GDP-to-GTP exchange (SOS-catalyzed; Michaelis-Menten kinetics)
        kRhydro, RAStₘ → RAS        # GTP hydrolysis on Ras

        # Ras-RAF1 binding, all (MEMBRANE):
        (kR1f, kR1r), RAStₘ + RAF1 ↔ RAF1ₘ     # Ras-RAF1, forward & reverse
        (kR1r, kR1f), pRAF1ₘ ↔ RAStₘ + pRAF1    # Ras-pRAF1, reverse & forward
        knfpR1r, nfpRAF1ₘ → RAStₘ + nfpRAF1     # Ras-nfpRAF1, unbinding only

        # Ras-BRAF binding, all (MEMBRANE):
        (kBf, kBr), RAStₘ + BRAF ↔ BRAFₘ       # Ras-BRAF, forward & reverse
        knfpBr, nfpBRAFₘ → RAStₘ + nfpBRAF     # Ras-nfpBRAF, unbinding only

        # RAF1 homodimerization (MEMBRANE):
        (kdR1f, kdR1r), RAF1ₘ + pRAF1ₘ ↔ RAF1pRAF1ₘ      # RAF1-pRAF1 , forward & reverse
        (kdR1r, kdR1f), pRAF1dimₘ ↔ 2pRAF1ₘ     # pRAF1-pRAF1, reverse & forward

        # BRAF homodimerization (MEMBRANE):
        (kdR1f, kdR1r), 2BRAFₘ ↔ BRAFdimₘ      # BRAF-BRAF, forward & reverse

        # RAF1-BRAF heterodimerization (MEMBRANE):
        (kdR1f, kdR1r), RAF1ₘ + BRAFₘ ↔ RAF1BRAFₘ    # RAF1-BRAF, forward & reverse
        (kdR1r, kdR1f), pRAF1BRAFₘ ↔ pRAF1ₘ + BRAFₘ      # pRAF1-BRAF dimerization, reverse & forward

        # RAF1 phosphorylation (MEMBRANE):
        kpR1, RAF1BRAFₘ → pRAF1BRAFₘ      # BRAF-catalyzed
        kpR1, RAF1pRAF1ₘ → pRAF1dimₘ      # pRAF1-catalyzed

        # pRAF1 dephosphorylation (MEMBRANE and CYTOPLASM):
        kdpR1, pRAF1ₘ → RAF1ₘ   # pRAF1 dephosphorylation on MEMBRANE
        kdpR1, pRAF1 → RAF1         # pRAF1 dephosphorylation in CYTOPLASM

        # MEK and ERK dephosphorylation (CYTOPLASM):
        kdpMEK, pMEK → MEK  # MEK
        kdpERK, pERK → ERK  # ERK

        # SOS dephosphorylation (CYTOPLASM):
        kdpSOS, nfpSOS → SOS

        # MEK phosphorylation, catalyzed by active RAF species (mass-action) (CYTOPLASM):
        2*kpMEK, BRAFdimₘ + MEK → pMEK + BRAFdimₘ     # BRAFdimₘ-catalzyed
        kpMEK, RAF1BRAFₘ + MEK → pMEK + RAF1BRAFₘ     # RAF1BRAFₘ-catalzyed
        2*kpMEK, pRAF1BRAFₘ + MEK → pMEK + pRAF1BRAFₘ     # pRAF1BRAFₘ-catalzyed
        kpMEK, RAF1pRAF1ₘ + MEK → pMEK + RAF1pRAF1ₘ     # RAF1pRAF1ₘ-catalzyed
        2*kpMEK, pRAF1dimₘ + MEK → pMEK + pRAF1dimₘ     # pRAF1dimₘ-catalzyed
        kpMEK, pRAF1ₘ + MEK → pMEK + pRAF1ₘ     # pRAF1ₘ-catalzyed
        kpMEK, pRAF1 + MEK → pMEK + pRAF1     # pRAF1-catalzyed

        # ERK phosphorylation, catalyed by pMEK (mass-action) (CYTOPLASM):
        kpERK, pMEK + ERK → pERK + pMEK   # pMEK-catalyzed

        # ERK-mediated negative feedback phosphorylation of membrane-bound species (mass-action) (MEMBRANE):
        knfpBR1, pERK + RAF1ₘ → nfpRAF1ₘ + pERK    # RAF1
        knfpBR1, pERK + BRAFₘ → nfpBRAFₘ + pERK    # BRAF
        knfpSOS, pERK + EG2SOSₘ → EG2ₘ + nfpSOS + pERK  # SOS
    end kRgne Kmgne kRhydro kR1f kR1r knfpR1r kBf kBr knfpBr kdR1f kdR1r kpR1 kdpR1 knfpBR1 knfpSOS kdpSOS kdpMEK kdpERK kpMEK kpERK

    rn = merge(receptor_rn, adaptor_rn)
    merge!(rn, mapk_rn)
end


## Variables for getting model outputs of interest:
@everywhere function get_fitting_vars(df)
    rasgtpnames = ["RAStₘ","RAF1ₘ","pRAF1ₘ","nfpRAF1ₘ","BRAFₘ","nfpBRAFₘ"] .* "(t)"
    rasgtpdimnames = ["RAF1pRAF1ₘ","pRAF1dimₘ","BRAFdimₘ","RAF1BRAFₘ","pRAF1BRAFₘ"] .* "(t)"
    raf1m_names = ["RAF1ₘ","pRAF1ₘ","RAF1BRAFₘ","pRAF1BRAFₘ"] .* "(t)";  # membrane species containing just one RAF1 molecule
    raf1mdim_names = ["RAF1pRAF1ₘ","pRAF1dimₘ"] .* "(t)";    # RAF1 dimers at the membrane (2 RAF1 molecules)

    rasgtp = sum.(eachrow(df[:,rasgtpnames])) .+ 2.0.*sum.(eachrow(df[:,rasgtpdimnames]))
    memraf1 = sum.(eachrow(df[:,raf1m_names])) .+ 2.0.*sum.(eachrow(df[:,raf1mdim_names]))
    pmek = df[:,"pMEK(t)"]
    perk = df[:,"pERK(t)"]

    return rasgtp, memraf1, pmek, perk
end

@everywhere function find_names(x, y)
    # Finds the locations of strings in `x` in `y` and returns the indices
    # in `y`, if they exist.

    yinds = Int64[]  # empty array of matches
    for i in eachindex(x)
        append!(yinds, in([x[i]]).(y) |> findall)
    end

    return yinds
end



## ============== Define model parameters ============== ##
# Get species & parameter info from reaction network:
parmap = paramsmap(rn);     # parameter name-index map
smap = speciesmap(rn);    # species name-index map
pnames = string.(rn.ps);  # convert parameter names to strings
@everywhere snames = string.(rn.states);  # convert species names to strings


# Load and set parameter values:
pnamesin = ["kSon","kSoff","kpR1","kdR1f","kR1f","kR1r","kBf","kBr","kdR1r","kdpR1","kiR1f","kiR1r","kpMEK","kpERK","kdpERK",
        "kdpMEK","kiBf","kiBr","knfpR1r","knfpBr", "knfpBR1","kEf","kEr","kdEf","kdEr","kcatE","kdp","kRgne","Kmgne",
        "kRhydro","kG2f","kG2r","kSOSf","kSOSr","kdpSOS","knfpSOS","kG2SOSf","kG2SOSr","knfpiR1r","knfpiBr",
        "kdeg"] # parameter namess
pvals = [2.83E+01,8.72E+00,2.29E+03,9.75E-01,3.31E-05,9.08E+01,5.63E+02,1.73E-02,1.31E+01,2.81E+01,1.48E-06,2.20E-01,1.59E-04,
        1.16E-06,3.77E-01,5.11E-01,1.58E-13,1.56E+00,4.12E-02,5.32E+00,3.81E-13,1.25E+04,2.16E+01,9.79E-01,1.52E+01,6.66E+00,
        5.18E+01,5.75E+00,7.21E+04,2.71E-01,6.52E-01,3.33E+03,5.52E-05,1.91E+00,5.86E-03,7.84E-06,2.59E-01,5.89E+01,1.74E-02,
        1.02E+01,0.0]  # parameter values
paramsin = DataFrame(:names=>names, :values=>pvals)    # load parameter values from external file
pvals = paramsin.values  # parameter values
pdict = Dict(zip(pnamesin, pvals))     # create dictionary of the parameters we just loaded

# Populate parameter vector for solver with parameters in correct order:
p = ones(length(parmap))    # initialize parameter vector with ones
for i in eachindex(pnames)
    p[i] = pdict[pnames[i]] # find parameter names (keys) in the dictionary created from the file we loaded, then get the parameter value
end


## Generate LHS design matrix of parameter sets to run:
n_lhs = 1000
lb = p .* 0.1
ub = p .* 10.0
lhs = QuasiMonteCarlo.sample(n_lhs, lb, ub, LatinHypercubeSample()) |> permutedims
lhs2 = [lhs[i,:] for i in axes(lhs,1)]
currentime = today() |> string
@pipe DataFrame(lhs, pnames) |> CSV.write("LHS_matrix_n=$n_lhs"*"_"*currentime*".csv",_)   # save LHS matrix to file


## Initial conditions:
prot_names = ["EGF","EGFR","GRB2","SOS","RAF1","BRAF","RAS","MEK","ERK"].*"(t)";
prot_abund = [1.70E-03,9.30E+04,6.28E+05,7.70E+03,1.20E+04,1.00E+03,1.35E+05,5.00E+05,6.00E+05]
u0 = zeros(length(snames))    # initial conditions vector

# Find the locations of the initial condition species in the reaction network and set their initial conditions:
ic_locs = find_names(prot_names, snames)
u0[ic_locs] = prot_abund   # set the initial conditions


## ============== Perform stochasticity calculations ============== ##
## Set up JumpProblem:
# Define the stochastic/jump problem:
t0 = 0.0   # initial integration time point
tf = 0.1 # final integration time point
tspan = (t0,tf)    # simulation time bounds

disc_prob = DiscreteProblem(rn, u0, tspan, p)
jump_prob = JumpProblem(rn, disc_prob, SortingDirect(), save_positions=(false,false))

# pmap version of stochasticity function:
@everywhere function compute_stochasticity_pmap(jump_prob, p_in; nreps::Int64=5, tint::Float64=0.1, t0::Float64=0.0, tf::Float64=15.0)
    sols = []   # initialize array to store all solutions for each parameter set

    t = collect(t0:tint:tf)
    prob = remake(jump_prob, p=p_in, tspan=(t0,tf))

    for j in 1:nreps    # run each parameter set n times
        sol_j = solve(prob, SSAStepper(), saveat=t)

        push!(sols, sol_j)
    end

    return sols
end

## Perform stochasticity calculations and save the results:
dt = 0.1   # time step for saving outputs
nreps = 5  # number of simulations to run per parameter set

fn_stochsols = "stochastic_model_Rivanna_solns_LHS-" * string(size(lhs,1)) * "_tf=" * string(tf) * "_dt=" * string(dt) * "_nreps=" * string(nreps) * "_" * currentime * ".jld2"

@time stoch_sols = pmap(p -> compute_stochasticity_pmap(jump_prob, p; nreps=nreps, tint=dt, t0=t0, tf=tf), lhs2);

save(fn_stochsols, "stoch_sols", stoch_sols)
