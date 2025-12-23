#=
--- Analyze stochasticity data from Rivanna ---
This script is for analyzing the stochasticity data that were generated from running the
extended EGFR-RAF model on UVA's Rivanna HPC cluster.
=#

## =============== Activate project environment =============== ##
using Pkg
Pkg.activate(".")
# Pkg.instantiate()   # may need to run this to install necessary packages
# Pkg.precompile()    # run this to speed up subsequent runs of this script


## =============== Load packages and other necessary functions =============== ##
using DifferentialEquations # for DiffEq types
using NumericalIntegration
using Pipe
using DataFrames, DataFramesMeta
using CSV, JLD2, DelimitedFiles
using Statistics, StatsBase
using Dates
using AlgebraOfGraphics, CairoMakie

using StatsModels, GLM
const sm = StatsModels

@isdefined(CairoMakie) && return(const cm = CairoMakie)
@isdefined(AlgebraOfGraphics) && return(const aog = AlgebraOfGraphics)
@isdefined(AlgebraOfGraphics) && set_aog_theme!()

# Set random seed:
using Random
Random.seed!(123)


# Load the reaction network and helpful functions/variables for extracting model outputs:
include("EGFR-Ras-Raf_extended_model_network_-sfnb.jl")

# Some helpful utility functions:
include("utilities.jl")

#=
The following file loads the MLJ.jl ecosystem and PLSRegressor object that we need
to use PLSR in our analysis. It also loads some PLSR-related utility functions.
=#
include("utilities_PLSR.jl")




## =============== Load Rivanna results =============== ##
# File names:
sols_fn = "Rivanna_results/stochastic_model_Rivanna_solns_LHS-1000_tf=15.0_dt=0.1_nreps=5_2021-08-16.jld2"
sols_csv_fn = "Rivanna_results/stochastic_model_Rivanna_solutions.csv"
lhs_fn = "Rivanna_results/LHS_matrix_n=1000_2021-08-16.csv"


## Load parameter data:
lhs_df = CSV.File(lhs_fn) |> DataFrame  # LHS-generated parameter values associated solutions in `sols`
lhs = Matrix(lhs_df[:,:])
pnames = names(lhs_df)  # parameter names
snames = string.(species(rn))  # convert species names to strings from the model's reaction network


## Load stochastic model solutions computed on Rivanna
if (isfile(sols_csv_fn))
    sols_df = CSV.read(sols_csv_fn, DataFrame)
else
    # Get model solutions from JLD2 archive and save to disk as CSV file
    sols = JLD2.load(sols_fn, "stoch_sols");

    sols_df = DataFrame()

    for i in eachindex(sols)
        sols_i = sols[i]
        for j in eachindex(sols_i)
            sols_mat_ij = permutedims(hcat(sols_i[j].u...))
            df_ij = @pipe DataFrame(sols_mat_ij, snames) |>
                @transform(_, :t = sols_i[j].t, :param_idx = i, :sol_idx = j)
            append!(sols_df, df_ij)
        end
    end

    CSV.write(sols_csv_fn, sols_df)
end


## =============== Calculate stochasticity scores =============== ##
# -- Extract membrane Raf1 solutions and calculate stochasticity scores:
stoch_scores = zeros(maximum(sols_df.param_idx))

y = stoch_scores    # alias for stoch_scores

# memraf1 = []
memraf1 = DataFrame()

for i in eachindex(stoch_scores)
    sols_i = @subset(sols_df, :param_idx .== i)

    ## Get integrated and normalized membrane Raf1 solutions:
    memraf1_i = DataFrame()
    for j in unique(sols_i.sol_idx)
        sols_ij = @subset(sols_i, :sol_idx .== j)
        memraf1_ij = get_fitting_vars(sols_ij)[2]   # extract membrane Raf1 values
        t_ij = sols_ij.t

        append!(memraf1_i, DataFrame(t = t_ij, memraf1 = memraf1_ij, param_idx = i, sol_idx = j))
    end

    # Get maximum membrane Raf1 values for each simulation:
    memraf1_i = @pipe memraf1_i |>
        groupby(_, :sol_idx) |>
        @transform(_, :max_memraf1 = maximum(:memraf1))

    # Normalize solutions to maximum:
    @transform!(memraf1_i, :memraf1_norm = :memraf1./:max_memraf1)   # normalize

    if all(unique(memraf1_i.max_memraf1 .== 0.0))
        stoch_score_i = 0.0     # set stochasticity to 0 if membrane Raf1 = 0 at all times for each sim
    else
        # Calculate variance at each time point and integrate over the time course:
        stoch_score_i = @pipe memraf1_i |>
            groupby(_, :t) |>
            @combine(_, :var = var(:memraf1_norm)) |>
            integrate(_.t, _.var, Trapezoidal())[1]
    end

    stoch_scores[i] = stoch_score_i     # store stochasticity score

    append!(memraf1, memraf1_i)         # store membrane Raf1 solutions
end

# Write stochasticity scores to disk:
stoch_scores_fn = "Rivanna_results/stochasticity_scores_Rivanna_n=1000_tf=15.0_dt=0.1_nreps=5_2021-08-06.csv"
if (!isfile(stoch_scores_fn))
    writedlm(stoch_scores_fn, stoch_scores)
end



## =============== Ordinary least squares on stochasticity scores =============== ##
# -- Get data for OLS:
lhsₙ = hcat([zscore(lhs[:,i]) for i in axes(lhs, 2)]...)   # z-scored parameter values
X = lhsₙ    # alias for z-scored LHS values
data1 = @pipe DataFrame([X log10.(y)], [pnames; "y"]) |>   # DataFrame of log-transformed LHS values and stochasticity scores
    filter(row -> !isinf(row.y), _)   # remove Inf values
X_reg = Matrix(data1[!, Not(:y)])   # X data for regression
y_reg = data1.y     # y data for regression


## Fit OLS model:
form1 = sm.term(:y) ~ sum(sm.term.(names(data1, Not(:y))))  # model formula
ols1 = lm(form1, data1)      # intercept is implied in this model
ols1 |> println

# Model results:
coef1_df = @pipe DataFrame(coef = coef(ols1)[2:end], # extract model coefficients
                           lci = confint(ols1)[2:end,1], # lower coefficient CI
                           uci = confint(ols1)[2:end,2], # upper coefficient CI
                           param = pnames) |>    # parameter names
    sort(_, :coef, rev=true)
confint1 = (coef1_df.coef-coef1_df.lci, coef1_df.uci-coef1_df.coef)   # store coefficient confidence intervals for plotting errorbars
R²1 = @pipe r2(ols1) |> round(_, digits=4); println("OLS: R² = ", R²1)
adjR²1 = @pipe adjr2(ols1) |> round(_, digits=4); println("OLS: Adj. R² = ", adjR²1)


## Plot OLS coefficients:
param_order = sort(coef1_df, :coef, rev=false).param

ols_coefs_plot = @pipe coef1_df |>
    sort(_, :coef, rev=false) |>
    @transform(_, @byrow :pval = !(:uci >= 0.0 >= :lci)) |>
    @transform(_, @byrow :uci = :uci - :coef) |>
    @transform(_, @byrow :lci = :coef - :lci) |>
    @transform(_, @byrow :xint = 0.0) |>
    unique(_) |>
    data(_) * (
        mapping(:xint => "", :param => sorter(param_order) => "") * visual(Lines, linestyle = :dash) +
        mapping(:coef => "", :param => sorter(param_order) => "", :lci, :uci, color=:pval => "p < 0.05") * # error bars
            visual(Errorbars, whiskerwidth=0, linewidth=2.5, direction=:x) +
        mapping(:coef => "", :param => sorter(param_order) => "", color=:pval => "p < 0.05") * visual(Scatter) # coefficient estimates
    )

# -- Draw the plot and save to disk:
fig1 = draw(ols_coefs_plot,
            axis = (xlabel="Regression coefficient", ylabel="Model parameter",
                    title="OLS on log-transformed\nmembrane RAF1 stochasticity scores\n(Adj. R²="*string(round(adjR²1, digits=3))*")",
                    width=200, height=550))
display(fig1)
save("Rivanna_results/plots/OLS-model_regression-coefficients_memRAF1-stochasticity.png", fig1, px_per_unit = 3)
save("Rivanna_results/plots/OLS-model_regression-coefficients_memRAF1-stochasticity.pdf", fig1, pt_per_unit = 1)



## =============== STOCHASTICITY PLOTS =============== ##
# -- Interpolate 2D grid of stochasticity values for each desired pair of parameters:
using Dierckx

# Parameters to use in plot of stochasticity scores:
plot_params = [["kR1r","Kmgne"],
               ["kR1r", "kR1f"],
               ["kBr", "kdp"],
               ["kRgne","kR1f"]]

ratio = 3/2
ht = 350
clims = (minimum(data1.y), maximum(data1.y))
clevs = LinRange(clims[1], clims[2], 20)
stoch_fig = Figure(; title="Interpolated stochasticity scores")
for i in eachindex(plot_params)
    par1 = plot_params[i][1]    # param 1 name
    par2 = plot_params[i][2]    # param 2 name

    xplot = data1[!, par1]     # param #1 vals
    yplot = data1[!, par2]     # param #2 vals
    cplot = data1.y     # stochasticity scores
    spl = Spline2D(xplot, yplot, cplot; s = 1e6)  # spline object for interpolation of stochasticity scores

    nitp = Int(200)  # number of grid points along params 1 and 2 for interpolation
    xitp = LinRange(minimum(xplot), maximum(xplot), nitp) |> collect    # grid of param #1 vals to interpolate over
    yitp = LinRange(minimum(yplot), maximum(yplot), nitp) |> collect    # grid of param #2 vals to interpolate over
    zitp = evalgrid(spl, xitp, yitp)    # grid of interpolated stochasticity scores

    # -- Plot filled contour of interpolated stochasticity scores:
    if i <= 2   # set indices for plot layouts
        ni = 1
        mi = i
    elseif i >= 2
        ni = 2
        mi = i-2
    end

    # Filled contour:
    ax, hm = cm.contourf(stoch_fig[ni, mi][1, 1], xitp, yitp, zitp, colormap=:inferno, levels=clevs)
    Colorbar(stoch_fig[ni, mi][1, 2], hm, label = "log₁₀(stochasticity)")
    ax.xlabel = par1 * " z-score"
    ax.ylabel = par2 * " z-score"
end

Label(stoch_fig[1, 1:2, Top()],
      "Interpolated stochasticity scores",         # set title on stochasticity figure
      valign = :bottom,
      font = assetpath("fonts", "NotoSans-Bold.ttf"),
      padding = (0, 0, 5, 0))

display(stoch_fig)              # display stochasticity figure

save("Rivanna_results/plots/Stochasticity_contour_plots.png", stoch_fig, px_per_unit = 3)
save("Rivanna_results/plots/Stochasticity_contour_plots.pdf", stoch_fig, pt_per_unit = 1)



## =============== INDIVIDUAL SOLUTION VISUALIZATIONS =============== ##
# Plot solutions of interest:
n_sols_plot = 30 # 56

sols_plot = @pipe memraf1 |>
    @subset(_, :param_idx .<= n_sols_plot) |>
    unique(_) |>
    data(_) *
        mapping(:t => "t (min)",
                :memraf1 => "membrane Raf1 (cell⁻¹)",
                color = :sol_idx => string,
                layout = :param_idx => x -> (string(x) |> sorter(string.(collect(1:n_sols_plot))))) *
            visual(Stairs) |>
    draw(_,
        axis = (;width = 80, height = 80, limits = (0, nothing, 0, nothing)),
        facet = (linkxaxes = :none, linkyaxes = :none)
    )

display(sols_plot)

save("Rivanna_results/plots/membrane-RAF1_top_stochastic_solutions.png", sols_plot, px_per_unit = 3)
save("Rivanna_results/plots/membrane-RAF1_top_stochastic_solutions.pdf", sols_plot, pt_per_unit = 3)





## =============== LHS-PLSR GSA: All parameter sets =============== ##
# Data for PLSR-based GSA:
x_plsr = lhs
y_plsr = @pipe memraf1  |>
    groupby(_, [:t, :param_idx]) |>
    @combine(_, :mean = mean(:memraf1)) |>
    groupby(_, :param_idx) |>
    @combine(_, :max = maximum(:mean)).max


#=
Standard scale X and Y. MLJ supports standardization of Y but not X in its model
pipelines, so to make things simpler, we will just scale the data ourselves before
feeding them to the model.
=#
scalerX = StatsBase.fit(ZScoreTransform, x_plsr, dims = 1)          # fit standard scaler to X
scalerY = StatsBase.fit(ZScoreTransform, y_plsr, dims = 1)          # fit standard scaler to Y
X_plsr = DataFrame(StatsBase.transform(scalerX, x_plsr), pnames)    # standardize X
Y_plsr = StatsBase.transform(scalerY, y_plsr)                       # standardize Y



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
    unique(_) |>
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
                # title = "LHS-PLSR GSA\n1000 stochastic sims\n(-sorafenib)",
                title = "LHS-PLSR GSA\n(1000 stochastic sims, -sorafenib)",
                xlabel = "Number of LVs", ylabel = "R²Y or Q²Y",
                limits = (nothing, nothing, 0, 1))
    )

display(fig_r2q2)

save("Rivanna_results/plots/PLSR_R2+Q2_all-param-sets.png", fig_r2q2, px_per_unit = 3)
save("Rivanna_results/plots/PLSR_R2+Q2_all-param-sets.pdf", fig_r2q2, pt_per_unit = 1)




## ============== Analyze optimal PLSR model ============== ##
# Perform bootstrapping on model coefficients and VIP scores:
@time coefs_boot_df, vip_boot_df = get_bootstrapped_coefs_and_vip_PLS1(best_model; nboot=10000, boot_frac=0.9)

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
                visual(Scatter, color = aog.wongcolors()[4])
        )
    ) |>
    draw(_,
        axis = (ylabel = "PLSR coefficient",
                title = "LHS-PLSR GSA\n(1000 stochastic sims, -sorafenib)",
                xticklabelrotation = pi/2, height = 200, width = 600))
display(fig_coefs_boot)
save("Rivanna_results/plots/PLSR_bootstrapped_coefficients.png", fig_coefs_boot, px_per_unit = 3)
save("Rivanna_results/plots/PLSR_bootstrapped_coefficients.pdf", fig_coefs_boot, pt_per_unit = 1)



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
                title = "LHS-PLSR GSA\n(1000 stochastic sims, -sorafenib)",
                xticklabelrotation = pi/2, height = 200, width = 600))
display(fig_vip_boot)
save("Rivanna_results/plots/PLSR_bootstrapped_VIP_scores.png", fig_vip_boot, px_per_unit = 3)
save("Rivanna_results/plots/PLSR_bootstrapped_VIP_scores.pdf", fig_vip_boot, pt_per_unit = 1)








## =============== LHS-PLSR GSA: Removing top stochastic parameter sets =============== ##
# Get data for PLSR-based GSA:
q = 0.9     # highest percentile/fraction of stochastic parameter sets to keep
yquant = quantile(y, q)  # calculate desired quantile for filtering out the top stochastic parameter sets
stochi = y .< yquant

# Filter data for the desired samples to analyze:
x_plsr2 = x_plsr[stochi,:]
y_plsr2 = y_plsr[stochi]


# Standard scale X and Y:
scalerX2 = StatsBase.fit(ZScoreTransform, x_plsr2, dims = 1)  # fit standard scaler to X
scalerY2 = StatsBase.fit(ZScoreTransform, y_plsr2, dims = 1)  # fit standard scaler to Y
X2 = DataFrame(StatsBase.transform(scalerX2, x_plsr2), pnames)    # standardize X
Y2 = StatsBase.transform(scalerY2, y_plsr2)                       # standardize Y



## Train PLSR model:
# Perform cross-validation on number of PLSR model components and get the best model:
max_ncomps_filt = 6
@time df_r2q2_filt = run_plsr_cv(X2, Y2; ncomps = max_ncomps_filt, cv_type = "LOO")

# Optimal number of components (where Q²Y is maximal):
best_ncomp_filt = @subset(df_r2q2_filt, :Q2Y .== maximum(:Q2Y)).ncomp[1]

# PLSR model with optimal number of components/latent factors:
best_model_filt = @subset(df_r2q2_filt, :Q2Y .== maximum(:Q2Y)).trained_model[1]


## Plot R²Y and Q²Y:
fig_r2q2_filt = @pipe df_r2q2_filt |>
    DataFrames.select(_, Not(:trained_model)) |>
    stack(_, Not(:ncomp)) |>
    @transform(_, :variable = replace(:variable, "R2Y" => "R²Y", "Q2Y" => "Q²Y")) |>
    # Plot variance explained by each component/LV:
    data(_) *
        mapping(:ncomp, :value,
            color = :variable => sorter(["R²Y", "Q²Y"]),
            marker = :variable => sorter(["R²Y", "Q²Y"])) *
            visual(ScatterLines) +
        # Add vline for best/optimal number of LVs:
        mapping([best_ncomp_filt]) * visual(VLines, linestyle = :dash) |>
    draw(_,
        axis = (width = 120, height = 100,
                xticks = 0:10,
                title = "LHS-PLSR GSA\n(Top "*string(Integer(round((1-q)*100)))*"% stochastic\nsimulations removed)",
                xlabel = "Number of LVs", ylabel = "R²Y or Q²Y",
                limits = (nothing, nothing, 0, 1))
    )
fig_r2q2_filt |> display
save("Rivanna_results/plots/PLSR-filt_R2+Q2_all-param-sets.png", fig_r2q2_filt, px_per_unit = 3)
save("Rivanna_results/plots/PLSR-filt_R2+Q2_all-param-sets.pdf", fig_r2q2_filt, pt_per_unit = 1)




## ============== Analyze optimal PLSR model ============== ##
# Perform bootstrapping on model coefficients and VIP scores:
@time coefs_boot_filt_df, vip_boot_filt_df = get_bootstrapped_coefs_and_vip_PLS1(best_model_filt; nboot=10000, boot_frac=0.9)

## Plot the coefficient boostrapping results:
fig_coefs_filt_df = @pipe coefs_boot_filt_df |>
    stack(_, Not(:boot_spl)) |>
    groupby(_, :variable) |>
    @combine(_, :μ = mean(:value), :uci = quantile(:value, 0.975), :lci = quantile(:value, 0.025)) |>
    sort(_, :μ, rev=true) |>
    unique(_)

fig_coefs_boot_filt = @pipe fig_coefs_filt_df |>
    (mapping([0]) * visual(HLines, linestyle = :dash) +     # horizontal line at y = 0
        data(_) * (
            mapping(:variable => sorter(unique(_.variable)), :lci, :uci) * # error bars
                visual(Rangebars, linewidth = 1.5) +
            mapping(:variable => sorter(unique(_.variable)), :μ) *
                visual(Scatter, color = aog.wongcolors()[4])
        )
    ) |>
    draw(_,
        axis = (ylabel = "PLSR coefficient",
                title = "LHS-PLSR GSA\n(Top "*string(Integer(round((1-q)*100)))*"% stochastic\nsimulations removed)",
                xticklabelrotation = pi/2, height = 200, width = 600))
display(fig_coefs_boot_filt)
save("Rivanna_results/plots/PLSR-filt_bootstrapped_coefficients.png", fig_coefs_boot_filt, px_per_unit = 3)
save("Rivanna_results/plots/PLSR-filt_bootstrapped_coefficients.pdf", fig_coefs_boot_filt, pt_per_unit = 1)



## Plot the VIP score boostrapping results:
fig_vip_filt_df = @pipe vip_boot_filt_df |>
    stack(_, Not(:boot_spl)) |>
    groupby(_, :variable) |>
    @combine(_, :μ = mean(:value), :uci = quantile(:value, 0.975), :lci = quantile(:value, 0.025)) |>
    @transform(_, :vip_thresh = :μ .> 1 .&& :lci .> 1) |>
    sort(_, :μ, rev=true) |>
    unique(_)

fig_vip_boot_filt = @pipe fig_vip_filt_df |>
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
                title = "LHS-PLSR GSA\n(Top "*string(Integer(round((1-q)*100)))*"% stochastic\nsimulations removed)",
                xticklabelrotation = pi/2, height = 200, width = 600))
display(fig_vip_boot_filt)
save("Rivanna_results/plots/PLSR-filt_bootstrapped_VIP_scores.png", fig_vip_boot_filt, px_per_unit = 3)
save("Rivanna_results/plots/PLSR-filt_bootstrapped_VIP_scores.pdf", fig_vip_boot_filt, pt_per_unit = 1)




## ============== Benchmarking other regression models ============== ##
#=
    In this section we test the performance of other regression models for predicting
    membrane RAF1 stochasticity, just to check that none perform better than PLSR.
    The models we'll test include:
        - Ridge regression
        - (kernel) SVR
        - Lasso/elastic net
        - LARS and/or robust regression
        - ϵ-support vector regression (ϵ-SVR) and various kernel functions

    This could easily be coded up in a way to programmatically and automatically
    loop through all of the models and associated hyperparameter tuning we want/need
    to run, but for the sake of getting it done (and since it isn't too hard), we'll
    just manually do training and tuning for each model of interest.
=#

RidgeRegressor = MLJ.@load RidgeRegressor pkg = MLJLinearModels
LassoRegressor = MLJ.@load LassoRegressor pkg = MLJLinearModels
ElasticNetRegressor = MLJ.@load ElasticNetRegressor pkg = MLJLinearModels
HuberRegressor = MLJ.@load HuberRegressor pkg = MLJLinearModels     # robust regressor with Hubero ρ loss that we can tune
EpsilonSVR = MLJ.@load EpsilonSVR pkg = LIBSVM

cv = CV(nfolds = 100, shuffle = true, rng = 123)        # cross-validation to use across all models


## Ridge regression -----------------------------------------------------------
ridgereg = RidgeRegressor()

# Define hyperparams for tuning
r1 = range(ridgereg, :lambda, lower = 1e-3, upper = 1e2, scale = :log10)

# Self-tuning model
self_tuning_rr_model = TunedModel(model = ridgereg,
                                  resampling = cv,
                                  tuning = Grid(resolution = 10),
                                  range = [r1],
                                  measure = [rmse, l2, l2_sum])

# Model machine
self_tuning_ridgereg = machine(self_tuning_rr_model, X_plsr, Y_plsr)

# Fit self-tuning model
fit!(self_tuning_ridgereg, verbosity = 0)

# Model report
display(report(self_tuning_ridgereg))

# Best model after tuning
rr_model = self_tuning_ridgereg.fitresult

# Compute prediction sum of squares
tss = sum((Y_plsr .- mean(Y_plsr)).^2)      # Y total sum of squares
press_ridgereg = 0.0
cvsplit = MLJBase.train_test_pairs(CV(nfolds = length(y)), 1:length(y))

for split in cvsplit
    X_trainᵢ = X_plsr[split[1], :]
    X_testᵢ = X_plsr[split[2], :]
    y_trainᵢ = Y_plsr[split[1]]
    y_testᵢ = Y_plsr[split[2]]
    machᵢ = machine(RidgeRegressor(lambda = rr_model.model.lambda), X_trainᵢ, y_trainᵢ)
    fit!(machᵢ, verbosity = 0)

    pressᵢ = (y_testᵢ .- MLJ.predict(machᵢ, X_testᵢ)).^2
    press_ridgereg += pressᵢ[1]
end

r2y_ridgereg = rsq(MLJ.predict(rr_model, X_plsr), Y_plsr)
q2y_ridgereg = 1 - press_ridgereg/tss




## Lasso -----------------------------------------------------------
lasso = LassoRegressor()

# Self-tuning model
self_tuning_lasso_model = TunedModel(model = lasso,
                                     resampling = cv,
                                     tuning = Grid(resolution = 10),
                                     range = [r1],
                                     measure = [rmse, l2, l2_sum])

# Model machine
self_tuning_lasso = machine(self_tuning_lasso_model, X_plsr, Y_plsr)

# Fit self-tuning model
fit!(self_tuning_lasso, verbosity = 0)

# Model report
display(report(self_tuning_lasso))

# Best model after tuning
lasso_model = self_tuning_lasso.fitresult

# Compute prediction sum of squares
press_lasso = 0.0

for split in cvsplit
    X_trainᵢ = X_plsr[split[1], :]
    X_testᵢ = X_plsr[split[2], :]
    y_trainᵢ = Y_plsr[split[1]]
    y_testᵢ = Y_plsr[split[2]]
    machᵢ = machine(LassoRegressor(lambda = lasso_model.model.lambda), X_trainᵢ, y_trainᵢ)
    fit!(machᵢ, verbosity = 0)

    pressᵢ = (y_testᵢ .- MLJ.predict(machᵢ, X_testᵢ)).^2
    press_lasso += pressᵢ[1]
end

r2y_lasso = rsq(MLJ.predict(lasso_model, X_plsr), Y_plsr)
q2y_lasso = 1 - press_lasso/tss




## Elastic net regression -----------------------------------------------------------
en_reg = ElasticNetRegressor()

# Elastic net hyperparameters
en_r1 = range(en_reg, :lambda, lower = 1e-3, upper = 1e2, scale = :log10)
en_r2 = range(en_reg, :gamma, lower = 1e-3, upper = 1e2, scale = :log10)

# Self-tuning model
self_tuning_en_model = TunedModel(model = en_reg,
                                  resampling = cv,
                                  tuning = Grid(resolution = 10),
                                  range = [en_r1, en_r2],
                                  measure = [rmse, l2, l2_sum])

# Model machine
self_tuning_en = machine(self_tuning_en_model, X_plsr, Y_plsr)

# Fit self-tuning model
fit!(self_tuning_en, verbosity = 0)

# Model report
display(report(self_tuning_en))

# Best model after tuning
en_model = self_tuning_en.fitresult

# Compute prediction sum of squares
press_en = 0.0

for split in cvsplit
    X_trainᵢ = X_plsr[split[1], :]
    X_testᵢ = X_plsr[split[2], :]
    y_trainᵢ = Y_plsr[split[1]]
    y_testᵢ = Y_plsr[split[2]]
    machᵢ = machine(ElasticNetRegressor(lambda = en_model.model.lambda, gamma = en_model.model.gamma), X_trainᵢ, y_trainᵢ)
    fit!(machᵢ, verbosity = 0)

    pressᵢ = (y_testᵢ .- MLJ.predict(machᵢ, X_testᵢ)).^2
    press_en += pressᵢ[1]
end

r2y_en = rsq(MLJ.predict(en_model, X_plsr), Y_plsr)
q2y_en = 1 - press_en/tss




## Robust regression -----------------------------------------------------------
# We'll use robust regression with Huber ρ loss and elastic net penalization.
rob_reg = HuberRegressor(penalty = :en)

# δ hyperparameter
rδ = range(rob_reg, :delta, lower = 1e-3, upper = 1e2, scale = :log10)

# Self-tuning model
self_tuning_rob_model = TunedModel(model = rob_reg,
                                   resampling = cv,
                                   tuning = Grid(resolution = 10),
                                   range = [rδ, en_r1, en_r2],
                                   measure = [rmse, l2, l2_sum])

# Model machine
self_tuning_rob = machine(self_tuning_rob_model, X_plsr, Y_plsr)

# Fit self-tuning model
fit!(self_tuning_rob, verbosity = 0)

# Model report
display(report(self_tuning_rob))

# Best model after tuning
rob_model = self_tuning_rob.fitresult

# Compute prediction sum of squares
press_rob = 0.0

for split in cvsplit
    X_trainᵢ = X_plsr[split[1], :]
    X_testᵢ = X_plsr[split[2], :]
    y_trainᵢ = Y_plsr[split[1]]
    y_testᵢ = Y_plsr[split[2]]
    machᵢ = machine(ElasticNetRegressor(lambda = rob_model.model.lambda, gamma = rob_model.model.gamma), X_trainᵢ, y_trainᵢ)
    fit!(machᵢ, verbosity = 0)

    pressᵢ = (y_testᵢ .- MLJ.predict(machᵢ, X_testᵢ)).^2
    press_rob += pressᵢ[1]
end

r2y_rob = rsq(MLJ.predict(rob_model, X_plsr), Y_plsr)
q2y_rob = 1 - press_rob/tss




## ϵ-SVR -----------------------------------------------------------
# We'll try support vector regression with multiple kernel functions
svr_kernels = [LIBSVM.Kernel.Linear,
               LIBSVM.Kernel.Polynomial,
               LIBSVM.Kernel.RadialBasis,
               LIBSVM.Kernel.Sigmoid]

for kernel in svr_kernels
    svr_reg = EpsilonSVR(kernel = :en)


    # gamma hyperparameter
    if kernel == LIBSVM.Kernel.Polynomial
        rgamma = range(svr_reg, :gamma, lower = 1e-3, upper = 1e2, scale = :log10)
    else
        rgamma = range(svr_reg, :delta, lower = 1e-3, upper = 1e2, scale = :log10)
    end

    # Self-tuning model
    self_tuning_svr_model = TunedModel(model = svr_reg,
    resampling = cv,
    tuning = Grid(resolution = 10),
    range = [rδ, en_r1, en_r2],
    measure = [rmse, l2, l2_sum])

    # Model machine
    self_tuning_svr = machine(self_tuning_svr_model, X_plsr, Y_plsr)

    # Fit self-tuning model
    fit!(self_tuning_svr, verbosity = 0)

    # Model report
    display(report(self_tuning_svr))

    # Best model after tuning
    svr_model = self_tuning_svr.fitresult

    # Compute prediction sum of squares
    press_svr = 0.0

    for split in cvsplit
        X_trainᵢ = X_plsr[split[1], :]
        X_testᵢ = X_plsr[split[2], :]
        y_trainᵢ = Y_plsr[split[1]]
        y_testᵢ = Y_plsr[split[2]]
        machᵢ = machine(ElasticNetRegressor(lambda = svr_model.model.lambda, gamma = svr_model.model.gamma), X_trainᵢ, y_trainᵢ)
        fit!(machᵢ, verbosity = 0)

        pressᵢ = (y_testᵢ .- MLJ.predict(machᵢ, X_testᵢ)).^2
        press_svr += pressᵢ[1]
    end

    r2y_svr = rsq(MLJ.predict(rr_model, X_plsr), Y_plsr)
    q2y_svr = 1 - press_svr/tss
end



