#=
--- Analyze VCell Batch Results ---
This script loads in the results from running (deterministic) batch simulations in VCell
using the parameter array input feature and analyzes the results to obtain
estimates of global sensitivities of model outputs with respect to model parameters.
=#

## Activate project environment:
# Store current working directory where this script and the Project.toml and Manifest.toml
# files for this project are located:
cwd0 = pwd()                                       # path with project environment
cwd = pwd() * "/VCell_batch_results"     # path with VCell data

# Activate the project's environment:
using Pkg
Pkg.activate(cwd0)
# Pkg.instantiate()     # installs all (missing) package dependencies for the current project



## Load packages:
using DifferentialEquations # for DiffEq types
using NumericalIntegration
using Pipe
using DataFrames, DataFramesMeta
using CSV, DelimitedFiles
using Statistics, StatsBase
using Dates
using AlgebraOfGraphics, CairoMakie

@isdefined(DataFrames) && return(const dfs = DataFrames)
@isdefined(CairoMakie) && return(const cm = CairoMakie)
@isdefined(AlgebraOfGraphics) && return(const aog = AlgebraOfGraphics)
@isdefined(AlgebraOfGraphics) && set_aog_theme!()

# Set random seed:
using Random
Random.seed!(123)

# Some helpful utility functions:
include("utilities.jl")
include("VCell_batch_funs.jl")

#=
The following file loads the MLJ.jl ecosystem and PLSRegressor object that we need
to use PLSR in our analysis. It also loads some PLSR-related utility functions.
=#
include("utilities_PLSR.jl")



## Load parameter distributions:
fn_params = cwd * "/param_array2.csv"  # file path
params = @pipe CSV.File(fn_params) |> # load parameter values used to generate batch solutions in VCell
    DataFrame(_) |>
    dfs.select(_, Not(r"ki|iB|iR|kSon|kSoff")) |>
    rename(_, :kfpBr => :knfpBr, :kfpR1r => :knfpR1r)
pnames = names(params)



## Load and store batch simulation results:
dir_batchres = cwd * "/raw_batch_results/batchResults_det_-sor"  # directory containing batch results from VCell
batchres, vars, pinds = loadbatch(params, dir_batchres)
vars = replace.(vars, r"_Count"=>"")



## Compute the time-averaged concentration of each species for each simulation and get other outputs of interest:
using NumericalIntegration
cmax = []   # array for maximum concentrations of each species per simulation
cave = []   # array for average concentrations of each species per simulation
cf = []     # array for concentration at final simulation time point

# Calculate maximum and average concentrations and get value at tf:
for i in axes(batchres, 1)
    resi = batchres[i]                  # results from simulation i
    ti = resi[:, 1]                     # time vector
    ci = resi[:, 2:end]                 # matrix of species concentrations
    maxi = maximum(ci, dims = 1)
    avei = Float64[]                    # vector to store average species concentrations
    for j in axes(ci, 2)
        int = integrate(ti, ci[:, j], Trapezoidal()) ./ (ti[end] - ti[1])   # compute average concentration (integrate)
        push!(avei, int)                # store average concentration
    end
    # cave[i,:] = avei
    push!(cf, ci[end, :])
    push!(cmax, maxi)
    push!(cave, avei)
end

# Get matrix of average species concentrations across batch simulations:
cmax = vcat(cmax...)
cave = hcat(cave...) |> permutedims
cf = hcat(cf...) |> permutedims
df_cmax = DataFrame(cmax, vars[2:end])
df_cave = DataFrame(cave, vars[2:end])
df_cf = DataFrame(cf, vars[2:end])

## Get specific outputs of interest:
# Total phosphorylated RAF1:
pRAF1_names = string.(keys(pRAF1_dict))                                             # names of species (columns) containing pRAF1
pRAF1_molec = [get(pRAF1_dict, i, 0) for i in pRAF1_names] |> permutedims           # number of pRAF1 molecules per species
df_cmax.pRAF1_tot = sum.(eachrow(pRAF1_molec.*df_cmax[:, pRAF1_names]));            # max values of total pRAF1 per simulation
df_cave.pRAF1_tot = sum.(eachrow(pRAF1_molec.*df_cave[:, pRAF1_names]));            # average concentration of total pRAF1 per simulation
df_cf.pRAF1_tot = sum.(eachrow(pRAF1_molec.*df_cave[:, pRAF1_names]));              # average concentration of total pRAF1 per simulation

# Membrane-bound RAF1:
memRAF1_names = string.(keys(memRAF1_dict))                                         # names of species (columns) containing membrane RAF1
memRAF1_molec = [get(memRAF1_dict, i, 0) for i in memRAF1_names] |> permutedims     # number of membrane RAF1 molecules per species
df_cmax.memRAF1 = sum.(eachrow(memRAF1_molec.*df_cmax[:, memRAF1_names]));          # max values of total membrane RAF1 per simulation
df_cave.memRAF1 = sum.(eachrow(memRAF1_molec.*df_cave[:, memRAF1_names]));          # average concentration of total membrane RAF1 per simulation
df_cf.memRAF1 = sum.(eachrow(memRAF1_molec.*df_cave[:, memRAF1_names]));            # average concentration of total membrane RAF1 per simulation

# Total Ras-GTP:
RASt_names = string.(keys(RASt_dict))                                               # names of species (columns) containing Ras-GTP
RASt_molec = [get(RASt_dict, i, 0) for i in RASt_names] |> permutedims              # number of Ras-GTP molecules per species
df_cmax.RasGTP_tot = sum.(eachrow(RASt_molec.*df_cmax[:, RASt_names]));             # max values of total Ras-GTP per simulation
df_cave.RasGTP_tot = sum.(eachrow(RASt_molec.*df_cave[:, RASt_names]));             # average concentration of total Ras-GTP per simulation
df_cf.RasGTP_tot = sum.(eachrow(RASt_molec.*df_cave[:, RASt_names]));               # average concentration of total Ras-GTP per simulation

# Total pEGFR:
pEGFR_names = string.(keys(pEGFR_dict))                                             # names of species (columns) containing pEGFR
pEGFR_molec = [get(pEGFR_dict, i, 0) for i in pEGFR_names] |> permutedims           # number of pEGFR molecules per species
df_cmax.pEGFR_tot = sum.(eachrow(pEGFR_molec.*df_cmax[:, pEGFR_names]));            # max values of total pEGFR per simulation
df_cave.pEGFR_tot = sum.(eachrow(pEGFR_molec.*df_cave[:, pEGFR_names]));            # average concentration of total pEGFR per simulation
df_cf.pEGFR_tot = sum.(eachrow(pEGFR_molec.*df_cave[:, pEGFR_names]));              # average concentration of total pEGFR per simulation




## ============== PLSR sensitivity analysis ============== ##
# Get X and Y data:
X = Matrix(params[pinds,:])                                 # parameter values as a matrix

ynames = ["memRAF1", "RasGTP_tot", "pMEK", "pERK"]
ynames_2 = Dict(ynames .=> ["membrane RAF1", "Active RAS", "pMEK", "pERK"])

for (i, yname_i) in enumerate(ynames)
    yname_i2 = ynames_2[yname_i]

    Y = df_cmax[:, yname_i]

    #=
    Standard scale X and Y. MLJ supports standardization of Y but not X in its model
    pipelines, so to make things simpler, we will just scale the data ourselves before
    feeding them to the model.
    =#
    scalerX = StatsBase.fit(ZScoreTransform, X, dims = 1)               # fit standard scaler to X
    scalerY = StatsBase.fit(ZScoreTransform, Y, dims = 1)               # fit standard scaler to Y
    X_plsr = DataFrame(StatsBase.transform(scalerX, X), pnames)         # standardize X
    Y_plsr = StatsBase.transform(scalerY, Y)                            # standardize Y




    ## Train PLSR model:
    # Perform cross-validation on number of PLSR model components and get the best model:
    max_ncomps = 5
    @time df_r2q2 = run_plsr_cv(X_plsr, Y_plsr; ncomps = max_ncomps, cv_type = "LOO")

    # Optimal number of components (where Q²Y is maximal):
    best_ncomp = @subset(df_r2q2, :Q2Y .== maximum(:Q2Y)).ncomp[1]

    # PLSR model with optimal number of components/latent factors:
    best_model = @subset(df_r2q2, :Q2Y .== maximum(:Q2Y)).trained_model[1]


    ## Plot R²Y and Q²Y:
    fig_r2q2 = @pipe df_r2q2 |>
        DataFrames.select(_, Not(:trained_model)) |>
        stack(_, Not(:ncomp)) |>
        @transform(_, :variable = replace(:variable, "R2Y" => "R²Y", "Q2Y" => "Q²Y")) |>
        # Plot variance explained by each component/LV:
        data(_) *
            mapping(:ncomp, :value,
                color = :variable => sorter(["R²Y", "Q²Y"]),
                marker = :variable => sorter(["R²Y", "Q²Y"])) *
                visual(ScatterLines) +
            # Add vline to indicate best/optimal number of LVs:
            mapping([best_ncomp]) * visual(VLines, linestyle = :dash) |>
        draw(_,
            axis = (width = 120, height = 100,
                    xticks = 0:10,
                    title = "LHS-PLSR GSA, $(yname_i2)\n(198 deterministic sims, -sorafenib)",
                    xlabel = "Number of LVs", ylabel = "R²Y or Q²Y",
                    limits = (nothing, nothing, 0, 1))
        )

    display(fig_r2q2)

    save(cwd0 * "/plots/VCell_det_res_-sfnb_PLSR_R2+Q2_$(yname_i).png", fig_r2q2, px_per_unit = 3)
    save(cwd0 * "/plots/VCell_det_res_-sfnb_PLSR_R2+Q2_$(yname_i).pdf", fig_r2q2, pt_per_unit = 1)






    ## ============== Analyze optimal PLSR model ============== ##
    ## Perform bootstrapping on model coefficients and VIP scores:
    @time coefs_boot_df, vip_boot_df = get_bootstrapped_coefs_and_vip_PLS1(best_model; nboot=1980, boot_frac=0.9)

    ## Plot the coefficient boostrapping results:
    fig_coefs_df = @pipe coefs_boot_df |>
        stack(_, Not(:boot_spl)) |>
        groupby(_, :variable) |>
        @combine(_, :μ = mean(:value), :uci = quantile(:value, 0.975), :lci = quantile(:value, 0.025)) |>
        sort(_, :μ, rev=true) |>
        unique(_)

    fig_coefs_boot = @pipe fig_coefs_df |>
        (mapping([0]) * visual(HLines, linestyle = :dash) +     # horizontal line at y = 0
            data(_) * (
                mapping(:variable => sorter(unique(_.variable)), :lci, :uci) * # error bars
                    visual(Rangebars, linewidth = 1.5) +
                mapping(:variable => sorter(unique(_.variable)), :μ) *
                    visual(Scatter, color = aog.wongcolors()[i])
            )
        ) |>
        draw(_,
            axis = (ylabel = "PLSR coefficient",
                    title = "LHS-PLSR SA, $(yname_i2)\n(198 deterministic sims, -sorafenib)",
                    xticklabelrotation = pi/2, height = 200, width = 600))
    display(fig_coefs_boot)
    save(cwd0 * "/plots/VCell_det_res_-sfnb_PLSR_bootstrapped_coefficients_$(yname_i).png", fig_coefs_boot, px_per_unit = 3)
    save(cwd0 * "/plots/VCell_det_res_-sfnb_PLSR_bootstrapped_coefficients_$(yname_i).pdf", fig_coefs_boot, pt_per_unit = 1)



    ## Plot the VIP score boostrapping results:
    fig_vip_df = @pipe vip_boot_df |>
        stack(_, Not(:boot_spl)) |>
        groupby(_, :variable) |>
        @combine(_, :μ = mean(:value), :uci = quantile(:value, 0.975), :lci = quantile(:value, 0.025)) |>
        @transform(_, :vip_thresh = :μ .> 1 .&& :lci .> 1) |>
        sort(_, :μ, rev=true) |>
        unique(_)

    fig_vip_boot = @pipe fig_vip_df |>
        (data(_) *
            (
                mapping(:variable => sorter(unique(_.variable)), :μ, color = :vip_thresh) *
                    visual(BarPlot) +
                mapping(:variable => sorter(unique(_.variable)), :lci, :uci) * # error bars
                    visual(Rangebars, linewidth = 1.5)
            ) +
            mapping([1]) * visual(HLines, linestyle = :dash)
        ) |>
        draw(_,
            axis = (ylabel = "VIP score",
                    title = "LHS-PLSR SA, $(yname_i)\n(198 deterministic sims, -sorafenib)",
                    xticklabelrotation = pi/2, height = 200, width = 600))
    display(fig_vip_boot)
    save(cwd0 * "/plots/VCell_det_res_-sfnb_PLSR_bootstrapped_VIP_scores_$(yname_i).png", fig_vip_boot, px_per_unit = 3)
    save(cwd0 * "/plots/VCell_det_res_-sfnb_PLSR_bootstrapped_VIP_scores_$(yname_i).pdf", fig_vip_boot, pt_per_unit = 1)
end