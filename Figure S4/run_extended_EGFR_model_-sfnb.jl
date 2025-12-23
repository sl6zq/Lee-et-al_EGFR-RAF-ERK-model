## ============== Activate project environment ============== ##
# Store current working directory where this script and the Project.toml and Manifest.toml
# files for this project are located:
cwd = pwd() # get current directory
println("Current working directory: ", cwd)

# Activate and precompile the project's environment:
using Pkg
Pkg.activate(cwd)
Pkg.instantiate()
# Pkg.precompile() # precompiles the packages associated with this project's environment -- may take a few minutes


## Use distributed approach to parallelism (for GSA)?
using Distributed
use_distributed = false
if use_distributed
    # Start workers on local machine using `Distributed` approach to parallelism:
    if nprocs() == 1
        addprocs(Threads.nthreads(), exeflags="--project=."); # adds the same number of processes as available CPU cores
    end
    nprocs()    # check number of available processes
    nworkers()  # check number of available workers
    @everywhere using ProgressMeter, SharedArrays
end

## ============== Load necessary packages/files ============== ##
using DelimitedFiles, CSV, JLD2, DataFrames, DataFramesMeta, Pipe
using Distributions
using AlgebraOfGraphics, CairoMakie
const aog = AlgebraOfGraphics
const cm = CairoMakie

using Catalyst
using SciMLSensitivity
using GlobalSensitivity
using Statistics
using StatsBase
using QuasiMonteCarlo
using Random
@everywhere begin
    using DifferentialEquations, NumericalIntegration

    include("EGFR-Ras-Raf_extended_model_network_-sfnb.jl") # this loads the `rn` variable, which generates the model reaction network
    include("utilities.jl")    # contains miscellaneous functions, including the trapezoidal rule
end

## ============== Load experimental data ============== ##
fn_exp = cwd * "/Exptl_data.csv"
expt = DataFrame(CSV.File(fn_exp))  # load exptl data used to generate model fits
expt_w = unstack(expt)  # exptl data in wide format

# Ras-GTP:
expt_ras = @pipe expt_w |>     # exptl RasGTP data
    filter(row -> row.sorafenib == "N", _) |>
    @select(_, :time, :error, :RasGTP) |> filter(row -> !ismissing(row.RasGTP), _)

# Membrane-bound RAF1:
expt_raf1 = @pipe expt_w |>     # exptl membrane RAF1 data
    filter(row -> row.sorafenib == "N", _) |>
    @select(_, :time, :error, :memRAF1) |> filter(row -> !ismissing(row.memRAF1), _)

# pMEK:
expt_pmek = @pipe expt_w |>     # exptl membrane RAF1 data
    filter(row -> row.sorafenib == "N", _) |>
    @select(_, :time, :error, :pMEK) |> filter(row -> !ismissing(row.pMEK), _)

# pERK:
expt_perk = @pipe expt_w |>     # exptl membrane RAF1 data
    filter(row -> row.sorafenib == "N", _) |>
    @select(_, :time, :error, :pERK) |> filter(row -> !ismissing(row.pERK), _)


## ============== Define model parameters ============== ##
parmap = paramsmap(rn)     # parameter name-index map
smap = speciesmap(rn)    # species name-index map
pnames = [string.(rn.ps)...]  # convert parameter names to strings
snames = [string.(rn.states)...]  # convert species names to strings
fitnames = ["RasGTP","memRAF1","pMEK","pERK"]   # names of model fitting variables

## Load and set parameter values:
# paramsin = readdlm("model_params.txt",'\t');    # load parameter values from external file
paramsin = DataFrame(CSV.File(cwd * "/model_params.txt"))    # load parameter values from external file
pnamesin = paramsin[:,1]        # parameter names
pvals = paramsin[:,2]           # parameter values
pdict = Dict(zip(pnamesin, pvals))     # create dictionary of the parameters we just loaded

# Populate parameter vector for solver with parameters in correct order:
p = ones(length(parmap))    # initialize parameter vector with ones
for i in eachindex(pnames)
    p[i] = pdict[pnames[i]] # find parameter names (keys) in the dictionary created from the file we loaded, then get the parameter value
end
p_df = DataFrame(name = pnames, value = p)

## Species initial conditions:
hela_prots = readdlm(cwd * "/HeLa_prot_abund.txt", '\t')  # HeLa protein abundances (cell⁻¹)
hela_prots = hela_prots[1:end-1,:]    # get rid of sorafenib if running the version without sfnb reactions
prot_names = hela_prots[:,1] .* "(t)"
prot_abund = hela_prots[:,2]
u0 = zeros(length(snames))    # initial conditions vector

# Find the locations of the initial condition species in the reaction network and set their initial conditions:
ic_locs = find_names(prot_names, snames)
u0[ic_locs] = prot_abund   # set the initial conditions

## Define the ODEs and ODE problem:
t0 = 0.0                                # initial integration time point
tf = 60.0                               # final integration time point
tspan = (t0, tf)                        # simulation time bounds
tint = 0.01                             # interval for saving time points
odes = convert(ODESystem, rn)           # convert reaction network to ODESystem
prob = ODEProblem(rn, u0, tspan, p)     # define ODE problem using the model's reaction network


## ============== Compute ODE solution ============== ##
@time sol = solve(prob, Rodas4())
tsol = sol.t
df = DataFrame(sol[:,:]', snames)   # store the solutions in a DataFrame
df.t = tsol # add solution time points to the DataFrame





## ============== Local sensitivity analysis (gradient-based) ============== ##
# Local sensitivity analysis:
lsaprob = ODEForwardSensitivityProblem(odes, u0, tspan, p); # define LSA problem
lsasol = solve(lsaprob, Rodas4(autodiff=false))
t_lsa = lsasol.t    # time vector for LSA solutions


## Extract LSA results and plot:
# Extract sensitivity results:
x1, dp1 = extract_local_sensitivities(lsasol);  # full time series

# Get full LSA timecourse results for model fitting variables (sum local sensitivities across aggregated species):
dp_fitting = []     # allocate array for local sensitivities of fitting variables
dpn_fitting = []    # allocate array for normalized local sensitivities of fitting variables

for i in eachindex(pnames)
    param_i = pnames[i]
    dp1_i = dp1[i]'
    dp1n_i = dp1_i./maximum(dp1_i, dims=1)
    df_i = DataFrame(dp1_i, snames)
    dfn_i = DataFrame(dp1n_i, snames)
    dras_i, dmemraf1_i, dpmek_i, dperk_i = get_fitting_vars(df_i)
    drasn_i, dmemraf1n_i, dpmekn_i, dperkn_i = get_fitting_vars(dfn_i)
    df_out = DataFrame(RasGTP=dras_i, memRAF1=dmemraf1_i, pMEK=dpmek_i, pERK=dperk_i)
    dfn_out = DataFrame(RasGTP=drasn_i, memRAF1=dmemraf1n_i, pMEK=dpmekn_i, pERK=dperkn_i)
    push!(dp_fitting, df_out)
    push!(dpn_fitting, dfn_out)
end

dfit = Dict(zip(pnames, dp_fitting))    # unnormalized local sensitivities
dfitn = Dict(zip(pnames, dpn_fitting))  # normalized local sensitivities


## Plot LSA results for fitting variables:
# # Initialize empty plot array:
# lsahms = []

# # Generate the plots for each parameter using the FULL TIME SERIES:
# for i in eachindex(pnames)
#     param_name = pnames[i]     # name of the parameter you're extracting
#     dpplot = dfitn[param_name] # get sensitivity time series for the parameter
#     @show param_name

#     # Plot a heatmap (per parameter) of the sensitivities:
#     # heatmap(eachindex(lsasol.t), snames, dpplot,
#     hm = @pipe dpplot |>
#         hcat(_, DataFrame(t = t_lsa)) |>
#         stack(_, Not(:t)) |>
#         data(_) * mapping(:t, :variable, :value) * visual(Heatmap, colormap = :bwr) |>
#         draw(_,
#             axis = (title=param_name, xlabel="t", ylabel="Species",
#                     yminorticksvisible = false, yminorgridvisible = false,
#                     xminorticksvisible = false, xminorgridvisible = false,
#                     width = 150, height = 80)
#         )
#     push!(lsahms, hm)
# end

# lsahms_dict = Dict(zip(pnames, lsahms))



## Compute and plot integrated (average) local sensitivities:
# Use trapezoidal rule to integrate the sensitivities over time:
dp1_ave = zeros(length(u0), length(p))

for i in eachindex(dp1)
    dp1_ave[:,i] = integrate(t_lsa, permutedims(dp1[i])) ./ (tspan[end]-tspan[1]);
end

dp1_avenorm = dp1_ave./maximum(abs.(dp1_ave), dims=2)   # normalized by maximum magnitude for each species

# Results for fitting variables only:
dp1ave_fit = @pipe permutedims(dp1_ave) |>
    DataFrame(_, snames) |>
    get_fitting_vars(_)  |>
    hcat(_...)
dp1ave_fit_df = DataFrame(dp1ave_fit, fitnames)
dp1ave_fit_df.params = pnames

dp1aven_fit = dp1ave_fit./maximum(abs.(dp1ave_fit), dims=1)  # normalize by maximum magnitude for each variable
dp1aven_fit_df = DataFrame(dp1aven_fit, fitnames)
dp1aven_fit_df.param = pnames


## Heatmap of normalized sensitivities for fitting variables
hm_norm_lsa = @pipe dp1aven_fit_df |>
    stack(_, Not(:param)) |>
    data(_) * mapping(:param, :variable, :value) * visual(Heatmap, colormap = :bwr) |>
    draw(_,
        axis = (title = "Normalized time-integrated local sensitivities", xlabel = "Parameter", ylabel = "Species",
                xticklabelrotation = pi/2, width = 600, height = 80)
    )
save(cwd*"/plots/LSA-norm_all-species_integrated_heatmap.png", hm_norm_lsa, px_per_unit = 3)
save(cwd*"/plots/LSA-norm_all-species_integrated_heatmap.pdf", hm_norm_lsa, pt_per_unit = 1)












## ============== Global sensitivity analysis ============== ##
## Define GSA problem functions:
# Define allowable parameter ranges:
pbounds = [[p[i].*0.1, p[i].*10.0] for i in eachindex(p)]

# Names of the outputs we're measuring:
gsanames = ["RasGTPₐᵥₑ","memRAF1ₐᵥₑ","pMEKₐᵥₑ","pERKₐᵥₑ","RasGTPₘₐₓ","memRAF1ₘₐₓ","pMEKₘₐₓ","pERKₘₐₓ"]

# GSA function (serial):
gsafunc = function(p_in)
    prob1 = remake(prob; p = p_in)
    sol1 = solve(prob1, Rodas4())
    t1 = sol1.t
    df1 = DataFrame(sol1[:,:]', snames)   # store the solutions in a DataFrame

    # Get solutions of interest:
    rasgtp_1 = sum.(eachrow(df1[:,rasgtpnames])) .+ 2.0.*sum.(eachrow(df1[:,rasgtpdimnames]))
    memraf1_1 = sum.(eachrow(df1[:,raf1m_names])) .+ 2.0.*sum.(eachrow(df1[:,raf1mdim_names]))
    pmek_1 = df1[:, "pMEK(t)"]
    perk_1 = df1[:, "pERK(t)"]
    outputs = [rasgtp_1 memraf1_1 pmek_1 perk_1]

    # Output:
    out1 = integrate(t1, outputs) ./ (t1[end]-t1[1])  # time-averaged species concentrations
    out2 = maximum(outputs', dims=2) # time-averaged species concentrations
    out = [out1; out2]
    return out
end

# Parallelized GSA function:
gsafunc_par = function(p_in)
    @show(size(p_in))
    if use_distributed
        par_prob_func1(prob, i, repeat) = remake(prob; p = p_in[:,i])
        ensemble_prob = EnsembleProblem(prob, prob_func = par_prob_func1)
        sol1 = solve(ensemble_prob, Rodas4(), EnsembleDistributed(); trajectories=size(p_in,2))
        # Now sol1[i] is the solution for the ith set of parameters
        out = zeros(length(gsanames), size(p_in, 2))
        for i in eachindex(sol1)
            soli = sol1[i]
            ti = soli.t
            dfi = DataFrame(soli[:,:]', snames)   # store the solutions in a DataFrame

            # Get solutions of interest:
            rasgtp_i, memraf1_i, pmek_i, perk_i = get_fitting_vars(dfi)
            outputs = [rasgtp_i memraf1_i pmek_i perk_i]

            # Output:
            out1 = integrate(ti, outputs) ./ (ti[end]-ti[1])  # time-averaged species concentrations
            out2 = maximum(outputs', dims=2) # time-averaged species concentrations
            out[:,i] = [out1; out2]
        end
    else
        par_prob_func(prob, i, repeat) = remake(prob; p = p_in[:,i])
        ensemble_prob = EnsembleProblem(prob, prob_func = par_prob_func)
        sol1 = solve(ensemble_prob, Rodas4(), EnsembleThreads(); trajectories=size(p_in,2))

        # Now sol1[i] is the solution for the ith set of parameters
        out = zeros(length(gsanames), size(p_in, 2))
        Threads.@threads for i in eachindex(sol1)
        # for i in eachindex(sol1)
            soli = sol1[i]
            ti = soli.t
            dfi = DataFrame(soli[:,:]', snames)   # store the solutions in a DataFrame

            # Get solutions of interest:
            rasgtp_i, memraf1_i, pmek_i, perk_i = get_fitting_vars(dfi)
            outputs = [rasgtp_i memraf1_i pmek_i perk_i]

            # Output:
            out1 = integrate(ti, outputs) ./ (ti[end]-ti[1])  # time-averaged species concentrations
            out2 = maximum(outputs', dims=2) # time-averaged species concentrations
            out[:,i] = [out1; out2]
        end
    end
    # @show size(out)
    return out
end

gsafunc_par2 = function(p_in)
    @show(size(p_in))
    par_prob_func(prob, i, repeat) = remake(prob; p = p_in[:,i])
    ensemble_prob = EnsembleProblem(prob, prob_func = par_prob_func)
    sol1 = solve(ensemble_prob, Rodas4(), EnsembleThreads(); trajectories=size(p_in,2))
    # Now sol1[i] is the solution for the ith set of parameters
    out = zeros(length(gsanames), size(p_in, 2))
    Threads.@threads for i in eachindex(sol1)
    # for i in eachindex(sol1)
        soli = sol1[i]
        ti = soli.t
        dfi = DataFrame(soli[:,:]', snames)   # store the solutions in a DataFrame

        # Get solutions of interest:
        rasgtp_i, memraf1_i, pmek_i, perk_i = get_fitting_vars(dfi)
        outputs = [rasgtp_i memraf1_i pmek_i perk_i]

        # Output:
        out1 = integrate(ti, outputs) ./ (ti[end]-ti[1])  # time-averaged species concentrations
        out2 = maximum(outputs', dims=2) # time-averaged species concentrations
        out[:,i] = [out1; out2]
    end
    # @show size(out)
    return out[2,:]
end




## ============== eFAST GSA ============== ##
run_efast = true
N_efast = 1000   # number of samples to run per parameter for eFAST → multiply this by # parameters for total number of samples
total_runs_efast = N_efast * length(pbounds)

fn_efast = cwd * "/GSA results/eFAST-GSA-res_" * string(N_efast) * "-spls.jld2"

if run_efast
    @time begin
        efast = GlobalSensitivity.gsa(gsafunc, eFAST(), pbounds, samples=N_efast)
        # efast = GlobalSensitivity.gsa(gsafunc_par, eFAST(), pbounds, n=N_efast, batch=true)
    end
    save(fn_efast, "efast", efast)   # save results
end

if !run_efast && ispath(fn_efast)
    efast = load(fn_efast, "efast") # load eFAST results from disk
end

## Plot eFAST results:
if @isdefined efast
    # # Per-output-normalized eFAST results:
    # efast_S1n = efast.S1 ./ maximum(efast.S1, dims=2)
    # efast_STn = efast.ST ./ maximum(efast.ST, dims=2)

    # eFAST results:
    efast_S1 = efast.S1
    efast_ST = efast.ST

    # First-order indices:
    hm_S1 = @pipe DataFrame(efast_S1', gsanames) |>
        hcat(_, DataFrame(param = pnames)) |>
        stack(_, Not(:param)) |>
        data(_) *
            mapping(:param, :variable => sorter(reverse(gsanames)), :value => "S₁") *
                visual(Heatmap, colormap = :inferno) +
            mapping([4.5]) * visual(HLines, color = :white) |>
        draw(_,
            axis = (title = "First-order eFAST indices", xlabel = "", ylabel = "",
                    xticklabelrotation = pi/2, width = 600, height = 160)
        )
    save(cwd*"/plots/efast_S1_heatmap_"*string(total_runs_efast)*"-samples.png", hm_S1, px_per_unit = 3)
    save(cwd*"/plots/efast_S1_heatmap_"*string(total_runs_efast)*"-samples.pdf", hm_S1, pt_per_unit = 1)

    # Total-order indices:
    hm_ST = @pipe DataFrame(efast_ST', gsanames) |>
        hcat(_, DataFrame(param = pnames)) |>
        stack(_, Not(:param)) |>
        data(_) *
            mapping(:param, :variable => sorter(reverse(gsanames)), :value => "Sₜ") *
                visual(Heatmap, colormap = :inferno) +
            mapping([4.5]) * visual(HLines, color=:white) |>
        draw(_,
            axis = (title = "Total-order eFAST indices", xlabel = "", ylabel = "",
                    xticklabelrotation = pi/2, width = 600, height = 160)
        )

    save(cwd*"/plots/efast_ST_heatmap_"*string(total_runs_efast)*"-samples.png", hm_ST, px_per_unit = 3)
    save(cwd*"/plots/efast_ST_heatmap_"*string(total_runs_efast)*"-samples.pdf", hm_ST, pt_per_unit = 1)
end




## ============== LHS-based GSA ============== ##
# Generate parameter design matrix using Latin Hypercube Sampling (LHS):
Random.seed!(123)
lb = p .* 0.1   # lower bounds on parameter values
ub = p .* 10.0   # upper bounds on parameter values
n_lhs = 3000   # number of parameter combinations to generate (samples)
@time lhs = QuasiMonteCarlo.sample(n_lhs, lb, ub, LatinHypercubeSample());    # generate LHS parameter design matrix
lhs_t = transpose(lhs);    # transpose of LHS parameter matrix
lhs_df = DataFrame(Tables.table(lhs_t, header=pnames));

# ODE solves for LHS parameter design matrix:
@time begin
    lhs_res = gsafunc_par(lhs)
end

# Calculate parameter-output correlations:
corrs = zeros(size(lhs,1), size(lhs_res,1));    # initialize Pearson correlation array
ranks = copy(corrs); # initialize Spearman rank array

Threads.@threads for i in axes(corrs, 2)
    x = lhs_res[i,:]    # result/output vector for which to compute correlations
    y = lhs_t   # parameter matrix
    corrs[:,i] = Statistics.cor(x, y)   # Pearson's R
    ranks[:,i] = StatsBase.corspearman(x, y)    # Spearman's ρ
end


## Plot the LHS GSA results for all outputs:
# Pearson:
hm_pearson = @pipe DataFrame(corrs, gsanames) |>
    hcat(_, DataFrame(param = pnames)) |>
    stack(_, Not(:param)) |>
    data(_) *
        mapping(:param, :variable => sorter(reverse(gsanames)), :value => "Pearson's R") *
            visual(Heatmap, colormap = :bwr) +
        mapping([4.5]) * visual(HLines) |>
    draw(_,
        axis = (title = "LHS-GSA", xlabel = "", ylabel = "",
                xticklabelrotation = pi/2, width = 600, height = 160)
    )
save(cwd*"/plots/LHS-pearson_heatmap_"*string(n_lhs)*"-samples.png", hm_pearson, px_per_unit = 3)
save(cwd*"/plots/LHS-pearson_heatmap_"*string(n_lhs)*"-samples.pdf", hm_pearson, pt_per_unit = 1)

# Spearman:
hm_spearman = @pipe DataFrame(ranks, gsanames) |>
    hcat(_, DataFrame(param = pnames)) |>
    stack(_, Not(:param)) |>
    data(_) *
        mapping(:param, :variable => sorter(reverse(gsanames)), :value => "Spearman's ρ") *
            visual(Heatmap, colormap = :bwr) +
        mapping([4.5]) * visual(HLines) |>
    draw(_,
        axis = (title = "LHS-GSA", xlabel = "", ylabel = "",
                xticklabelrotation = pi/2, width = 600, height = 160)
    )
save(cwd*"/plots/LHS-spearman_heatmap_"*string(n_lhs)*"-samples.png", hm_spearman, px_per_unit = 3)
save(cwd*"/plots/LHS-spearman_heatmap_"*string(n_lhs)*"-samples.pdf", hm_spearman, pt_per_unit = 1)
